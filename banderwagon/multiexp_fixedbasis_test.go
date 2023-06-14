package banderwagon

import (
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"testing"

	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fp"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
)

func TestCorrectness(t *testing.T) {
	t.Parallel()

	for _, msmLength := range []int{5, 6, 8, 10, 16, 19, 32, 64, 128, 256} {
		points := GenerateRandomPoints(uint64(msmLength))

		windowSize := []int{4, 8}
		for _, w := range windowSize {
			msmEngine, err := bandersnatch.New(points, w)
			if err != nil {
				t.Fatalf("error in msm engine: %v", err)
			}
			for i := 0; i < 10; i++ {
				// Generate random scalars.
				scalars := make([]fr.Element, msmLength)
				for i := 0; i < len(scalars); i++ {
					scalars[i].SetRandom()
				}

				// MSM custom result.
				result, err := msmEngine.MSM(scalars)
				if err != nil {
					t.Fatalf("error in msm multiexp: %v", err)
				}

				// Gnark result.
				var gnarkResult bandersnatch.PointProj
				_, err = gnarkResult.MultiExp(points, scalars, bandersnatch.MultiExpConfig{ScalarsMont: true})
				if err != nil {
					t.Fatalf("error in gnark multiexp: %v", err)
				}

				if !result.Equal(&gnarkResult) {
					t.Fatalf("msm result does not match gnark result")
				}
			}
		}
		fmt.Printf("Correct for msm length %d\n", msmLength)
	}
}

func BenchmarkCustomMSM(b *testing.B) {
	windowSize := []int{4, 8}
	msmLength := []int{1, 2, 4, 8, 16, 32, 64, 128, 256}

	points := GenerateRandomPoints(256)
	// Generate random scalars.
	scalars := make([]fr.Element, len(points))
	for i := 0; i < len(scalars); i++ {
		scalars[i].SetRandom()
	}
	for _, k := range msmLength {
		b.Run(fmt.Sprintf("msm_length=%d", k), func(b *testing.B) {
			for _, w := range windowSize {
				b.Run(fmt.Sprintf("window_size=%d", w), func(b *testing.B) {
					msmEngine, _ := bandersnatch.New(points, w)

					b.ReportAllocs()
					b.ResetTimer()
					for i := 0; i < b.N; i++ {
						_, _ = msmEngine.MSM(scalars[:k])
					}
				})
			}
		})
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
