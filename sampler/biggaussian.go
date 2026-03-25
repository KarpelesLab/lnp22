package sampler

import (
	"io"
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
)

// SampleBigGaussianPoly samples a polynomial with each coefficient from D_σ.
// The Gaussian samples themselves are small (fit in int64), then promoted to *big.Int.
func SampleBigGaussianPoly(r *ring.BigRing, sigma float64, rng io.Reader) ring.BigPoly {
	cdt := newGaussianCDT(sigma)
	p := r.NewPoly()
	for i := 0; i < r.N; i++ {
		p[i] = bigmod.Reduce(big.NewInt(cdt.sample(rng)), r.Q)
	}
	return p
}

// SampleBigGaussianVec samples a vector of l BigPolynomials from D_σ.
func SampleBigGaussianVec(r *ring.BigRing, l int, sigma float64, rng io.Reader) ring.BigPolyVec {
	cdt := newGaussianCDT(sigma)
	v := r.NewPolyVec(l)
	for i := 0; i < l; i++ {
		for j := 0; j < r.N; j++ {
			v[i][j] = bigmod.Reduce(big.NewInt(cdt.sample(rng)), r.Q)
		}
	}
	return v
}
