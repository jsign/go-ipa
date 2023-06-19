package pipfixedbasis

import (
	"context"
	"errors"
	"fmt"
	"runtime"

	"github.com/bits-and-blooms/bitset"
	"github.com/crate-crypto/go-ipa/bandersnatch"
	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
	"golang.org/x/sync/errgroup"
)

const (
	precompNumPoints = 256

	pipperMsmLength = 256 - precompNumPoints
	fieldSizeBits   = 256
)

type MSMFixedBasis struct {
	// windowSize is the number of bits in the window (e.g: 8).
	windowSize int
	// windowMask is a mask to get the c-bits window (e.g: 11111111b).
	windowMask int
	// numWindows is the number of windows (e.g: 256 / 8 = 32).
	numWindows int

	// pointsPowers are the pointsPowers with their powers. pointsPowers[i][j] = pointsPowers[i] * 2^(windowSize * j).
	pointsPowers  [pipperMsmLength][]bandersnatch.PointAffine
	precompPoints [precompNumPoints]PrecompPoint
}

func New(points []bandersnatch.PointAffine, windowSize int) (*MSMFixedBasis, error) {
	if len(points) > precompNumPoints+pipperMsmLength {
		return nil, fmt.Errorf("max msm length is %d", precompNumPoints+pipperMsmLength)
	}
	if windowSize > fieldSizeBits {
		return nil, errors.New("c must be less than field size")
	}
	numWindows := fieldSizeBits / windowSize

	// Compute the powers of the points.
	var frQ fr.Element
	frQ.SetUint64(1 << windowSize)
	var pointsPowers [pipperMsmLength][]bandersnatch.PointAffine
	for i := range points[precompNumPoints:] {
		pointsPowers[i] = make([]bandersnatch.PointAffine, numWindows)
		pointsPowers[i][0] = points[i+precompNumPoints]
		for j := 1; j < len(pointsPowers[i]); j++ {
			pointsPowers[i][j].ScalarMul(&pointsPowers[i][j-1], &frQ)
		}
	}

	// Precomputed table.
	group, _ := errgroup.WithContext(context.Background())
	group.SetLimit(runtime.NumCPU())
	var precompPoints [precompNumPoints]PrecompPoint
	for i := 0; i < precompNumPoints; i++ {
		i := i
		group.Go(func() error {
			if i < 5 {
				precompPoints[i] = NewPrecompPoint(points[i], 16)
			} else {
				precompPoints[i] = NewPrecompPoint(points[i], 8)
			}
			return nil
		})
	}
	_ = group.Wait()

	return &MSMFixedBasis{
		windowSize: windowSize,
		windowMask: (1 << windowSize) - 1,
		numWindows: numWindows,

		pointsPowers:  pointsPowers,
		precompPoints: precompPoints,
	}, nil
}

func (msm *MSMFixedBasis) MSM(scalars []fr.Element) (bandersnatch.PointProj, error) {
	// Check that scalar length matches points length.
	if len(scalars) > precompNumPoints+pipperMsmLength {
		return bandersnatch.PointProj{}, errors.New("more scalars than accepted")
	}

	if len(scalars) <= precompNumPoints {
		return msm.msmPrecomp(scalars), nil
	}

	precompRes := msm.msmPrecomp(scalars[:precompNumPoints])
	pipperRes := msm.msmPipper(scalars[precompNumPoints:])
	return *precompRes.Add(&precompRes, &pipperRes), nil
}

func (msm *MSMFixedBasis) msmPrecomp(scalars []fr.Element) bandersnatch.PointProj {
	var res bandersnatch.PointProj
	res.Identity()

	set := bitset.New(uint(len(scalars)))
	for i := range scalars {
		if !scalars[i].IsZero() {
			set.Set(uint(i))
		}
	}

	minScalarsPerRoutine := 4
	if int(set.Count()) <= minScalarsPerRoutine {
		for i := range scalars {
			if set.Test(uint(i)) {
				msm.precompPoints[i].ScalarMul(scalars[i], &res)
			}
		}
		return res
	}

	numBatches := (int(set.Count()) + minScalarsPerRoutine - 1) / minScalarsPerRoutine
	if numBatches > runtime.NumCPU() {
		numBatches = runtime.NumCPU()
	}
	batchSize := (int(set.Count()) + numBatches - 1) / numBatches

	results := make(chan bandersnatch.PointProj, numBatches)
	for start := 0; start < len(scalars); {
		end := start
		size := 0
		for size < batchSize && end < len(scalars) {
			if set.Test(uint(end)) {
				size++
			}
			end++
		}
		go func(start, end int) {
			var res bandersnatch.PointProj
			res.Identity()
			for i := start; i < end; i++ {
				if set.Test(uint(i)) {
					msm.precompPoints[i].ScalarMul(scalars[i], &res)
				}
			}
			results <- res
		}(start, end)
		start = end
	}

	var result bandersnatch.PointProj
	result.Identity()
	for i := 0; i < numBatches; i++ {
		res := <-results
		result.Add(&result, &res)
	}

	return result
}

func (msm *MSMFixedBasis) msmPipper(scalars []fr.Element) bandersnatch.PointProj {
	const minScalarsPerRoutine = 8
	batches := (len(scalars) + minScalarsPerRoutine - 1) / minScalarsPerRoutine
	if batches > runtime.NumCPU() {
		batches = runtime.NumCPU()
	}
	results := make(chan bandersnatch.PointProj, batches)

	for i := 1; i < batches; i++ {
		start := i * len(scalars) / batches
		end := (i + 1) * len(scalars) / batches
		go msm.doWork(msm.pointsPowers[start:end], scalars[start:end], results)
	}
	msm.doWork(msm.pointsPowers[:len(scalars)/batches], scalars[:len(scalars)/batches], results)
	result := <-results
	for i := 1; i < batches; i++ {
		res := <-results
		result.Add(&result, &res)
	}

	return result
}
func (msm *MSMFixedBasis) doWork(points [][]bandersnatch.PointAffine, scalars []fr.Element, results chan<- bandersnatch.PointProj) {
	var bucketz [1<<7 + 1]bandersnatch.PointProj
	buckets := bucketz[:1<<(msm.windowSize-1)+1]

	for i := 0; i < len(buckets); i++ {
		buckets[i].Identity()
	}

	for i, scalar := range scalars {
		if scalar.IsZero() {
			continue
		}

		scalar.FromMont()
		carry := 0

		var pNeg bandersnatch.PointAffine
		for limbIndex := 0; limbIndex < fr.Limbs; limbIndex++ {
			for w := 0; w < 64/msm.windowSize; w++ {
				windowValue := int(scalar[limbIndex]>>(w*msm.windowSize))&msm.windowMask + carry
				if windowValue == 0 {
					continue
				}
				carry = 0

				if windowValue >= (1 << (msm.windowSize - 1)) {
					windowValue = windowValue - (1 << msm.windowSize)
					pNeg.Neg(&points[i][limbIndex*(64/msm.windowSize)+w])
					buckets[-windowValue].MixedAdd(&buckets[-windowValue], &pNeg)
					carry = 1
				} else {
					buckets[windowValue].MixedAdd(&buckets[windowValue], &points[i][limbIndex*(64/msm.windowSize)+w])
				}
			}
		}
	}

	var tmp, result bandersnatch.PointProj
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
