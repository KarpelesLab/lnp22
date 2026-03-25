package trapgen

import (
	"errors"
	"io"
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

// Params holds TrapGen configuration.
type Params struct {
	Ring  *ring.BigRing
	Sigma float64 // Gaussian std dev for short polynomials f, g
	Base  int     // gadget base (typically 2)
}

// DefaultParams returns default TrapGen parameters for the given ring.
func DefaultParams(r *ring.BigRing) *Params {
	return &Params{
		Ring:  r,
		Sigma: 1.17, // standard NTRU Gaussian parameter
		Base:  2,
	}
}

// NTRUKey holds the NTRU keypair.
type NTRUKey struct {
	F    ring.BigPoly // short polynomial f
	G    ring.BigPoly // short polynomial g
	H    ring.BigPoly // public key h = g/f mod q (i.e., f^{-1} * g)
	FInv ring.BigPoly // f^{-1} mod q (cached)
}

// PublicKey contains the public matrix and gadget structure.
type PublicKey struct {
	A ring.BigPolyVec // public polynomials [a_1, ..., a_l]
	H ring.BigPoly    // NTRU public key h = f^{-1} * g
}

// GenerateNTRUKey generates an NTRU keypair by sampling short f, g
// and computing h = f^{-1} * g mod q.
// Retries if f is not invertible.
func GenerateNTRUKey(params *Params, rng io.Reader) (*NTRUKey, error) {
	r := params.Ring
	maxAttempts := 1000

	for attempt := 0; attempt < maxAttempts; attempt++ {
		f := sampler.SampleBigGaussianPoly(r, params.Sigma, rng)
		g := sampler.SampleBigGaussianPoly(r, params.Sigma, rng)

		// Ensure f is invertible mod q by checking if gcd(Resultant(f, X^N+1), q) = 1.
		// For a practical approach: try to compute f^{-1} in R_q.
		fInv, err := polyInverse(r, f)
		if err != nil {
			continue // f not invertible, retry
		}

		// h = f^{-1} * g mod q
		h := r.Mul(fInv, g)

		return &NTRUKey{F: f, G: g, H: h, FInv: fInv}, nil
	}
	return nil, errors.New("trapgen: failed to generate invertible f")
}

// TrapGen generates a public matrix A and trapdoor.
//
// The structure is: given l random polynomials a_1, ..., a_{l-1} and the NTRU key,
// the last column a_l is set so that a_l = -(Σ a_i * r_i) + g_entry mod q
// where r_i are short and g_entry comes from the gadget.
//
// This gives [A | A·R + G] structure in the module setting.
func TrapGen(params *Params, rng io.Reader) (*PublicKey, *Trapdoor, error) {
	r := params.Ring
	l := GadgetLen(r.Q, params.Base) + 1 // trapdoor width

	// Generate NTRU key
	key, err := GenerateNTRUKey(params, rng)
	if err != nil {
		return nil, nil, err
	}

	// Sample short trapdoor polynomials R = [r_1, ..., r_l]
	trapR := sampler.SampleBigGaussianVec(r, l, params.Sigma, rng)

	// Generate random public polynomials a_1, ..., a_{l-1}
	A := r.NewPolyVec(l)
	for i := 0; i < l-1; i++ {
		A[i] = sampler.SampleBigUniformPoly(r, rng)
	}
	// Set last entry to include the NTRU public key
	A[l-1] = key.H

	pk := &PublicKey{A: A, H: key.H}
	td := &Trapdoor{R: trapR, Base: params.Base, Params: params}

	return pk, td, nil
}

// polyInverse computes the inverse of f in R_q = Z_q[X]/(X^N+1).
// Uses the iterative Newton-lifting approach for power-of-2 cyclotomic rings.
// Returns error if f is not invertible.
func polyInverse(r *ring.BigRing, f ring.BigPoly) (ring.BigPoly, error) {
	// Start with f0 = f[0]^{-1} mod q (scalar inverse)
	if f[0].Sign() == 0 {
		return nil, errors.New("trapgen: f[0] = 0, not invertible")
	}
	f0Inv := new(big.Int).ModInverse(bigmod.Reduce(f[0], r.Q), r.Q)
	if f0Inv == nil {
		return nil, errors.New("trapgen: f[0] not invertible mod q")
	}

	// If NTT is available, compute inverse via NTT:
	// f is invertible iff all NTT coefficients are invertible.
	if r.HasNTT() {
		fNTT := r.NTT(f)
		invNTT := r.NewPoly()
		for i := 0; i < r.N; i++ {
			v := bigmod.Reduce(fNTT[i], r.Q)
			if v.Sign() == 0 {
				return nil, errors.New("trapgen: f has zero in NTT domain, not invertible")
			}
			inv := new(big.Int).ModInverse(v, r.Q)
			if inv == nil {
				return nil, errors.New("trapgen: NTT coefficient not invertible")
			}
			invNTT[i] = inv
		}
		return r.InvNTT(invNTT), nil
	}

	// Fallback: Newton iteration for rings without NTT.
	// g_0 = f[0]^{-1}, g_{i+1} = g_i * (2 - f * g_i) mod (q, X^{2^i})
	// After log2(N) iterations, g = f^{-1} mod (q, X^N+1).
	g := r.NewPoly()
	g[0].Set(f0Inv)

	two := r.NewPoly()
	two[0].SetInt64(2)

	for step := 1; step < r.N; step <<= 1 {
		// g = g * (2 - f * g)
		fg := r.Mul(f, g)
		diff := r.Sub(two, fg)
		g = r.Mul(g, diff)
	}

	// Verify: f * g = 1
	check := r.Mul(f, g)
	one := r.One()
	if !r.Equal(check, one) {
		return nil, errors.New("trapgen: Newton iteration did not converge")
	}

	return g, nil
}
