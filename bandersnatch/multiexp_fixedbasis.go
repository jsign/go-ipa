package bandersnatch

import (
	"errors"
	"runtime"
	"sync"

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
	pointsPowers [msmLength][]PointAffine
}

func New(points []PointAffine, windowSize int) (*MSMFixedBasis, error) {
	if len(points) > msmLength {
		return nil, errors.New("max msm length is 256")
	}
	if windowSize > fieldSizeBits {
		return nil, errors.New("c must be less than field size")
	}
	numWindows := fieldSizeBits / windowSize

	// Compute the powers of the points.
	var frQ fr.Element
	frQ.SetUint64(1 << windowSize)
	var pointsPowers [msmLength][]PointAffine
	for i := range points {
		pointsPowers[i] = make([]PointAffine, numWindows)
		pointsPowers[i][0] = points[i]
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

	const workPerRoutine = 8
	count := (len(scalars) + workPerRoutine - 1) / workPerRoutine
	if count > runtime.NumCPU() {
		count = runtime.NumCPU()
	}
	results := make(chan PointProj, count)

	var wg sync.WaitGroup
	doWork := func(idx, start, end int) {
		defer wg.Done()
		buckets := make([]PointProj, 1<<msm.windowSize)
		for i := 0; i < len(buckets); i++ {
			buckets[i].Identity()
		}

		for i, scalar := range scalars[start:end] {
			k := i + start
			scalar.FromMont()
			for limbIndex := 0; limbIndex < fr.Limbs; limbIndex++ {
				for w := 0; w < 64/msm.windowSize; w++ {
					windowValue := int(scalar[limbIndex]>>(w*msm.windowSize)) & msm.windowMask
					if windowValue == 0 {
						continue
					}
					buckets[windowValue].MixedAdd(&buckets[windowValue], &msm.pointsPowers[k][limbIndex*(64/msm.windowSize)+w])
				}
			}
		}

		var tmp, result PointProj
		tmp.Identity()
		result.Identity()
		for k := len(buckets) - 1; k >= 1; k-- {
			if !buckets[k].IsIdentity() {
				tmp.Add(&tmp, &buckets[k])
			}
			result.Add(&result, &tmp)
		}
		results <- result
	}

	wg.Add(count)
	for i := 1; i < count; i++ {
		go doWork(i, i*len(scalars)/count, (i+1)*len(scalars)/count)
	}
	doWork(0, 0, len(scalars)/count)
	result := <-results
	for i := 1; i < count; i++ {
		res := <-results
		result.Add(&result, &res)
	}

	return result, nil
}
