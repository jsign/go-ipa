package multiproof

import (
	"fmt"
	"testing"

	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
	"github.com/crate-crypto/go-ipa/banderwagon"
	"github.com/crate-crypto/go-ipa/ipa"
)

func BenchmarkCompare(b *testing.B) {
	config := ipa.NewIPASettings()
	pl := banderwagon.NewPrecomputeLagrange(config.SRSPrecompPoints.SRS)
	points := banderwagon.GetAffinePoints(config.SRSPrecompPoints.SRS)

	msmLength := []int{1, 2, 4, 8, 16, 32, 64, 128, 256}
	// Generate random scalars.
	scalars := make([]fr.Element, len(points))
	for i := 0; i < len(scalars); i++ {
		scalars[i].SetRandom()
	}
	for _, k := range msmLength {
		scalarsBench := make([]fr.Element, len(points))
		copy(scalarsBench, scalars[:k])

		b.Run(fmt.Sprintf("msm_length=%d", k), func(b *testing.B) {
			b.Run("custom", func(b *testing.B) {
				msmEngine, _ := bandersnatch.New(points, 8)

				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = msmEngine.MSM(scalars[:k])
				}
			})
			b.Run("gnark", func(b *testing.B) {
				var gnarkResult bandersnatch.PointProj
				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = gnarkResult.MultiExp(points, scalarsBench, bandersnatch.MultiExpConfig{ScalarsMont: true})
				}
			})
			b.Run("precomp", func(b *testing.B) {
				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					pl.Commit(scalarsBench[:k])
				}
			})
		})
	}
}
