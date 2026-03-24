// Package sampler provides cryptographic sampling primitives for lattice-based proofs.
package sampler

import (
	"encoding/binary"
	"io"
	"math"

	"github.com/KarpelesLab/lnp22/internal/modmath"
	"github.com/KarpelesLab/lnp22/ring"
)

// gaussianCDT stores the cumulative distribution table for a discrete Gaussian
// with a specific standard deviation.
type gaussianCDT struct {
	sigma   float64
	table   []uint64 // CDF values scaled to [0, 2^64)
	tailCut int      // table covers [-tailCut, tailCut]
}

// newGaussianCDT builds a CDT for the centered discrete Gaussian D_σ.
// The tail is cut at ⌈τ*σ⌉ where τ=12 (probability mass in tail is negligible).
func newGaussianCDT(sigma float64) *gaussianCDT {
	tau := 12.0
	tailCut := int(math.Ceil(tau * sigma))

	probs := make([]float64, 2*tailCut+1)
	total := 0.0
	twoSigmaSq := 2.0 * sigma * sigma
	for i := 0; i <= 2*tailCut; i++ {
		x := float64(i - tailCut)
		probs[i] = math.Exp(-x * x / twoSigmaSq)
		total += probs[i]
	}

	table := make([]uint64, 2*tailCut+1)
	cum := 0.0
	for i := 0; i < 2*tailCut; i++ {
		cum += probs[i] / total
		table[i] = uint64(cum * math.Exp2(64))
	}
	table[2*tailCut] = math.MaxUint64

	return &gaussianCDT{sigma: sigma, table: table, tailCut: tailCut}
}

func (g *gaussianCDT) sample(rng io.Reader) int64 {
	var buf [8]byte
	if _, err := io.ReadFull(rng, buf[:]); err != nil {
		panic("sampler: failed to read randomness: " + err.Error())
	}
	r := binary.LittleEndian.Uint64(buf[:])

	lo, hi := 0, len(g.table)-1
	for lo < hi {
		mid := (lo + hi) / 2
		if g.table[mid] < r {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return int64(lo - g.tailCut)
}

// SampleGaussianPoly samples a polynomial with each coefficient drawn from D_σ.
func SampleGaussianPoly(r *ring.Ring, sigma float64, rng io.Reader) ring.Poly {
	cdt := newGaussianCDT(sigma)
	p := r.NewPoly()
	for i := 0; i < r.N; i++ {
		p[i] = modmath.Reduce(cdt.sample(rng), r.Q)
	}
	return p
}

// SampleGaussianVec samples a vector of l polynomials from D_σ.
func SampleGaussianVec(r *ring.Ring, l int, sigma float64, rng io.Reader) ring.PolyVec {
	cdt := newGaussianCDT(sigma)
	v := r.NewPolyVec(l)
	for i := 0; i < l; i++ {
		for j := 0; j < r.N; j++ {
			v[i][j] = modmath.Reduce(cdt.sample(rng), r.Q)
		}
	}
	return v
}
