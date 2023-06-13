package bandersnatch

import (
	"errors"

	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
)

const msmLength = 256
const fieldSizeBits = 256

type MSMFixedBasis struct {
	// windowSize is the number of bits in the window (e.g: 8).
	windowSize int
	// windowMask is a mask to get the c-bits window (e.g: 11111111b).
	windowMask int
	// numWindows is the number of windows (e.g: 256 / 8 = 32).
	numWindows int

	// pointsPowers are the pointsPowers with their powers. pointsPowers[i][j] = pointsPowers[i] * 2^(windowSize * j).
	pointsPowers [msmLength][]PointProj
}

func New(points []PointAffine, windowSize int) (*MSMFixedBasis, error) {
	if len(points) != msmLength {
		return nil, errors.New("points length must be 256")
	}
	if windowSize > fieldSizeBits {
		return nil, errors.New("c must be less than field size")
	}
	numWindows := fieldSizeBits / windowSize

	// Compute the powers of the points.
	var frQ fr.Element
	frQ.SetUint64(1 << windowSize)
	var pointsPowers [msmLength][]PointProj
	for i := range points {
		pointsPowers[i] = make([]PointProj, numWindows)
		pointsPowers[i][0].FromAffine(&points[i])
		for j := 1; j < len(pointsPowers[i]); j++ {
			pointsPowers[i][j].ScalarMul(&pointsPowers[i][j-1], &frQ)
		}
	}

	return &MSMFixedBasis{
		windowSize: windowSize,
		windowMask: (1 << windowSize) - 1,
		numWindows: numWindows,

		pointsPowers: pointsPowers,
	}, nil
}

func (msm *MSMFixedBasis) MSM(scalars []fr.Element) (PointProj, error) {
	// Check that scalar length matches points length.
	if len(scalars) > len(msm.pointsPowers) {
		return PointProj{}, errors.New("more scalars than accepted")
	}

	buckets := make([]PointProj, 1<<msm.windowSize)
	for i := 0; i < len(buckets); i++ {
		buckets[i].Identity()
	}

	for i, scalar := range scalars {
		scalar.FromMont()
		for limbIndex := 0; limbIndex < fr.Limbs; limbIndex++ {
			for w := 0; w < 64/msm.windowSize; w++ {
				windowValue := int(scalar[limbIndex]>>(w*msm.windowSize)) & msm.windowMask
				if windowValue == 0 {
					continue
				}
				buckets[windowValue].Add(&buckets[windowValue], &msm.pointsPowers[i][limbIndex*(64/msm.windowSize)+w])
			}
		}
	}

	var tmp, result PointProj
	tmp.Identity()
	result.Identity()
	for k := len(buckets) - 1; k >= 1; k-- {
		tmp.Add(&tmp, &buckets[k])
		result.Add(&result, &tmp)
	}

	return result, nil
}
