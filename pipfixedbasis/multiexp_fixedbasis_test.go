package pipfixedbasis

import (
	"fmt"
	"testing"

	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
	"github.com/crate-crypto/go-ipa/banderwagon"
	"github.com/crate-crypto/go-ipa/ipa"
)

func TestCorrectness(t *testing.T) {
	t.Parallel()

	srs := ipa.GenerateRandomPoints(256)
	points := banderwagon.GetAffinePoints(srs)

	for _, msmLength := range []int{5, 6, 8, 10, 16, 19, 32, 64, 128, 256} {
		points := points[:msmLength]

		windowSize := []int{4, 8}
		for _, w := range windowSize {
			msmEngine, err := New(points, w)
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

	srs := ipa.GenerateRandomPoints(256)
	points := banderwagon.GetAffinePoints(srs)

	// Generate random scalars.
	scalars := make([]fr.Element, len(points))
	for i := 0; i < len(scalars); i++ {
		scalars[i].SetRandom()
	}
	for _, k := range msmLength {
		b.Run(fmt.Sprintf("msm_length=%d", k), func(b *testing.B) {
			for _, w := range windowSize {
				b.Run(fmt.Sprintf("window_size=%d", w), func(b *testing.B) {
					msmEngine, _ := New(points, w)

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
				for _, w := range []int{4, 8} {
					b.Run("window_size="+fmt.Sprintf("%d", w), func(b *testing.B) {
						msmEngine, _ := New(points, w)

						b.ReportAllocs()
						b.ResetTimer()
						for i := 0; i < b.N; i++ {
							_, _ = msmEngine.MSM(scalars[:k])
						}
					})
				}
			})
			b.Run("precomp", func(b *testing.B) {
				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					pl.Commit(scalarsBench[:k])
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

		})
	}
}

func BenchmarkInitialize(b *testing.B) {
	srs := ipa.GenerateRandomPoints(256)
	points := banderwagon.GetAffinePoints(srs)

	for i := 0; i < b.N; i++ {
		New(points, 1)
	}
}
