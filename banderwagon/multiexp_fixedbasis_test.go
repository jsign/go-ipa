package banderwagon

import (
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"testing"
	"time"

	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fp"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
)

func TestCorrectness(t *testing.T) {
	t.Parallel()

	points := GenerateRandomPoints(256)

	// Generate random scalars.
	scalars := make([]fr.Element, 256)
	for i := 0; i < len(scalars); i++ {
		scalars[i].SetRandom()
	}

	// Gnark result.
	var gnarkResult bandersnatch.PointProj
	now := time.Now()
	_, err := gnarkResult.MultiExp(points, scalars, bandersnatch.MultiExpConfig{ScalarsMont: true})
	if err != nil {
		t.Fatalf("error in gnark multiexp: %v", err)
	}
	fmt.Printf("gnark took %v\n", time.Since(now))

	// MSM custom result.
	msmEngine, err := bandersnatch.New(points, 8)
	if err != nil {
		t.Fatalf("error in msm engine: %v", err)
	}
	now = time.Now()
	result, err := msmEngine.MSM(scalars)
	if err != nil {
		t.Fatalf("error in msm multiexp: %v", err)
	}
	fmt.Printf("custom msm took %v\n", time.Since(now))

	if !result.Equal(&gnarkResult) {
		t.Fatalf("msm result does not match gnark result")
	}
}

func BenchmarkCustomMSM(b *testing.B) {
	points := GenerateRandomPoints(256)
	// Generate random scalars.
	scalars := make([]fr.Element, 256)
	for i := 0; i < len(scalars); i++ {
		scalars[i].SetRandom()
	}
	msmEngine, _ := bandersnatch.New(points, 8)

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = msmEngine.MSM(scalars)
	}
}

func GenerateRandomPoints(numPoints uint64) []bandersnatch.PointAffine {
	seed := "eth_verkle_oct_2021"

	points := []bandersnatch.PointAffine{}

	var increment uint64 = 0

	for uint64(len(points)) != numPoints {

		digest := sha256.New()
		digest.Write([]byte(seed))

		b := make([]byte, 8)
		binary.BigEndian.PutUint64(b, increment)
		digest.Write(b)

		hash := digest.Sum(nil)

		var x fp.Element
		x.SetBytes(hash)

		increment++

		x_as_bytes := x.Bytes()
		var point_found Element
		err := point_found.SetBytes(x_as_bytes[:])
		if err != nil {
			continue
		}
		var pointAff bandersnatch.PointAffine
		pointAff.FromProj(&point_found.inner)
		points = append(points, pointAff)
	}

	return points
}
