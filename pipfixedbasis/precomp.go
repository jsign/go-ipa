package pipfixedbasis

import (
	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
)

type PrecompPoint struct {
	numWindows   int
	windowLength int

	windows [][]bandersnatch.PointAffine
}

func NewPrecompPoint(point bandersnatch.PointAffine, windowSize int) PrecompPoint {
	var specialWindow fr.Element
	specialWindow.SetUint64(1 << windowSize)

	res := PrecompPoint{
		numWindows:   256 / precompWindowSize,
		windowLength: 1<<(windowSize-1) + 1,
		windows:      make([][]bandersnatch.PointAffine, 256/precompWindowSize),
	}
	for i := 0; i < len(res.windows); i++ {
		res.windows[i] = make([]bandersnatch.PointAffine, res.windowLength)
		res.windows[i][0].Identity()

		curr := point
		for j := 1; j < len(res.windows[i]); j++ {
			res.windows[i][j] = curr
			curr.Add(&curr, &point)
		}
		point.ScalarMul(&point, &specialWindow)
	}

	return res
}

func (pp *PrecompPoint) ScalarMul(scalar fr.Element, res *bandersnatch.PointProj) {
	scalar.FromMont()
	var carry uint64

	var pNeg bandersnatch.PointAffine
	for l := 0; l < fr.Limbs; l++ {
		const numWindowsInLimb = 64 / precompWindowSize
		for w := 0; w < numWindowsInLimb; w++ {
			windowValue := (scalar[l]>>(precompWindowSize*w))&((1<<precompWindowSize)-1) + carry
			carry = 0
			if windowValue == 0 {
				continue
			}
			if windowValue >= 1<<(precompWindowSize-1) {
				windowValue = (1 << precompWindowSize) - windowValue
				pNeg.Neg(&pp.windows[l*numWindowsInLimb+w][windowValue])
				res.MixedAdd(res, &pNeg)
				carry = 1
			} else {
				res.MixedAdd(res, &pp.windows[l*numWindowsInLimb+w][windowValue])
			}
		}
	}
}
