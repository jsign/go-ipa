package bandersnatch

import (
	"context"
	"errors"
	"runtime"

	"github.com/crate-crypto/go-ipa/bandersnatch/fr"
	"golang.org/x/sync/errgroup"
)

const numPrecomputedPoints = 5
const pipperMsmLength = 256 - numPrecomputedPoints
const fieldSizeBits = 256

type MSMFixedBasis struct {
	// windowSize is the number of bits in the window (e.g: 8).
	windowSize int
	// windowMask is a mask to get the c-bits window (e.g: 11111111b).
	windowMask int
	// numWindows is the number of windows (e.g: 256 / 8 = 32).
	numWindows int

	// pointsPowers are the pointsPowers with their powers. pointsPowers[i][j] = pointsPowers[i] * 2^(windowSize * j).
	pointsPowers     [pipperMsmLength][]PointAffine
	firstFivePrecomp [5][16][1 << 16]PointAffine
}

func New(points []PointAffine, windowSize int) (*MSMFixedBasis, error) {
	if len(points) > numPrecomputedPoints+pipperMsmLength {
		return nil, errors.New("max msm length is 256")
	}
	if windowSize > fieldSizeBits {
		return nil, errors.New("c must be less than field size")
	}
	numWindows := fieldSizeBits / windowSize

	// Compute the powers of the points.
	var frQ fr.Element
	frQ.SetUint64(1 << windowSize)
	var pointsPowers [pipperMsmLength][]PointAffine
	for i := range points[5:] {
		pointsPowers[i] = make([]PointAffine, numWindows)
		pointsPowers[i][0] = points[i]
		for j := 1; j < len(pointsPowers[i]); j++ {
			pointsPowers[i][j].ScalarMul(&pointsPowers[i][j-1], &frQ)
		}
	}

	// Precomputed table.
	var specialWindow fr.Element
	specialWindow.SetUint64(1 << 16)
	var firstFivePrecomp [5][16][1 << 16]PointAffine
	group, _ := errgroup.WithContext(context.Background())
	group.SetLimit(runtime.NumCPU())
	for pointIdx := 0; pointIdx < 5; pointIdx++ {
		p := points[pointIdx]
		for windowIdx := 0; windowIdx < 16; windowIdx++ {
			pointIdx := pointIdx
			windowIdx := windowIdx
			base := p
			group.Go(func() error {
				curr := base
				for j := 1; j < 1<<16; j++ {
					firstFivePrecomp[pointIdx][windowIdx][j] = curr
					curr.Add(&curr, &base)
				}
				return nil
			})
			p.ScalarMul(&p, &specialWindow)
		}
	}
	_ = group.Wait()

	return &MSMFixedBasis{
		windowSize: windowSize,
		windowMask: (1 << windowSize) - 1,
		numWindows: numWindows,

		pointsPowers:     pointsPowers,
		firstFivePrecomp: firstFivePrecomp,
	}, nil
}

func (msm *MSMFixedBasis) MSM(scalars []fr.Element) (PointProj, error) {
	// Check that scalar length matches points length.
	if len(scalars) > numPrecomputedPoints+pipperMsmLength {
		return PointProj{}, errors.New("more scalars than accepted")
	}

	if len(scalars) <= numPrecomputedPoints {
		return msm.msmPrecomp(scalars), nil
	}

	precompRes := msm.msmPrecomp(scalars[:numPrecomputedPoints])
	pipperRes := msm.msmPipper(scalars[numPrecomputedPoints:])
	return *precompRes.Add(&precompRes, &pipperRes), nil
}

func (msm *MSMFixedBasis) msmPrecomp(scalars []fr.Element) PointProj {
	var res PointProj
	res.Identity()

	for k, scalar := range scalars {
		scalar.FromMont()
		for l := 0; l < fr.Limbs; l++ {
			for i := 0; i < 4; i++ {
				window := (scalar[l] >> (16 * i)) & 0xFFFF
				if window == 0 {
					continue
				}
				res.MixedAdd(&res, &msm.firstFivePrecomp[k][4*l+i][window])
			}
		}
	}
	return res
}

func (msm *MSMFixedBasis) msmPipper(scalars []fr.Element) PointProj {
	const workPerRoutine = 8
	count := (len(scalars) + workPerRoutine - 1) / workPerRoutine
	if count > runtime.NumCPU() {
		count = runtime.NumCPU()
	}
	results := make(chan PointProj, count)

	for i := 1; i < count; i++ {
		points := msm.pointsPowers[i*len(scalars)/count : (i+1)*len(scalars)/count]
		scalars := scalars[i*len(scalars)/count : (i+1)*len(scalars)/count]
		go msm.doWork(points, scalars, results)
	}
	msm.doWork(msm.pointsPowers[:len(scalars)/count], scalars[:len(scalars)/count], results)
	result := <-results
	for i := 1; i < count; i++ {
		res := <-results
		result.Add(&result, &res)
	}

	return result
}
func (msm *MSMFixedBasis) doWork(points [][]PointAffine, scalars []fr.Element, results chan<- PointProj) {
	var bucketz [1<<7 + 1]PointProj
	buckets := bucketz[:1<<(msm.windowSize-1)+1]

	for i := 0; i < len(buckets); i++ {
		buckets[i].Identity()
	}

	for i, scalar := range scalars {
		scalar.FromMont()
		carry := 0

		var pNeg PointAffine
		for limbIndex := 0; limbIndex < fr.Limbs; limbIndex++ {
			for w := 0; w < 64/msm.windowSize; w++ {
				windowValue := int(scalar[limbIndex]>>(w*msm.windowSize))&msm.windowMask + carry
				carry = 0
				if windowValue == 0 {
					continue
				}

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
