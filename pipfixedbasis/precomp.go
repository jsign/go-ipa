package pipfixedbasis

import (
	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
)

type PrecompPoint struct {
	windowSize int
	windows    [][]bandersnatch.PointAffine
}

func NewPrecompPoint(point bandersnatch.PointAffine, windowSize int) PrecompPoint {
	var specialWindow fr.Element
	specialWindow.SetUint64(1 << windowSize)

	res := PrecompPoint{
		windowSize: windowSize,
		windows:    make([][]bandersnatch.PointAffine, 256/windowSize),
	}
	for i := 0; i < len(res.windows); i++ {
		res.windows[i] = make([]bandersnatch.PointAffine, 1<<(windowSize-1))
		curr := point
		for j := 0; j < len(res.windows[i]); j++ {
			res.windows[i][j] = curr
			curr.Add(&curr, &point)
		}
		point.ScalarMul(&point, &specialWindow)
	}

	return res
}

func (pp *PrecompPoint) ScalarMul(scalar fr.Element, res *bandersnatch.PointProj) {
	numWindowsInLimb := 64 / pp.windowSize

	scalar.FromMont()
	var carry uint64
	var pNeg bandersnatch.PointAffine
	for l := 0; l < fr.Limbs; l++ {
		for w := 0; w < numWindowsInLimb; w++ {
			windowValue := (scalar[l]>>(pp.windowSize*w))&((1<<pp.windowSize)-1) + carry
			if windowValue == 0 {
				continue
			}
			carry = 0

			if windowValue > 1<<(pp.windowSize-1) {
				windowValue = (1 << pp.windowSize) - windowValue
				if windowValue != 0 {
					pNeg.Neg(&pp.windows[l*numWindowsInLimb+w][windowValue-1])
					res.MixedAdd(res, &pNeg)
				}
				carry = 1
			} else {
				res.MixedAdd(res, &pp.windows[l*numWindowsInLimb+w][windowValue-1])
			}
		}
	}
}
