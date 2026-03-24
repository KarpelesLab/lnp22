// Package commitment implements the BDLOP lattice-based commitment scheme
// (Baum–Damgård–Lyubashevsky–Oechsner–Peikert).
//
// A BDLOP commitment to a vector of messages m = (m_1, ..., m_ℓ) uses
// a public key (B₀, b_1, ..., b_ℓ) and randomness r to produce:
//
//	t₀ = B₀ · r           (binding component, μ polynomials)
//	t_i = ⟨b_i, r⟩ + m_i  (message components, ℓ polynomials)
//
// The scheme is additively homomorphic and supports zero-knowledge opening proofs.
package commitment

import (
	"io"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

// Params holds the commitment scheme parameters.
type Params struct {
	Mu     int     // number of binding equations (rows in B₀)
	Lambda int     // randomness dimension overhead (security parameter component)
	Ell    int     // number of message slots
	Sigma  float64 // Gaussian std dev for randomness
}

// DefaultParams returns default commitment parameters for 128-bit security.
func DefaultParams() *Params {
	return &Params{
		Mu:     3,
		Lambda: 3,
		Ell:    4,
		Sigma:  50,
	}
}

// randDim returns the total dimension of the randomness vector: μ + λ + ℓ.
func (p *Params) randDim() int {
	return p.Mu + p.Lambda + p.Ell
}

// CommitKey is the public commitment key.
type CommitKey struct {
	Params *Params
	B0     ring.PolyMat // μ × (μ+λ+ℓ) matrix
	Bm     ring.PolyMat // ℓ × (μ+λ+ℓ) matrix (row vectors b_1, ..., b_ℓ)
}

// Commitment is the output of the commitment function.
type Commitment struct {
	T0 ring.PolyVec // binding component (μ polynomials)
	Tm ring.PolyVec // message components (ℓ polynomials)
}

// Opening contains the information needed to verify a commitment.
type Opening struct {
	R   ring.PolyVec // randomness vector (μ+λ+ℓ polynomials)
	Msg ring.PolyVec // message vector (ℓ polynomials)
}

// KeyGen generates a random commitment key.
func KeyGen(params *Params, rng io.Reader) *CommitKey {
	d := params.randDim()
	return &CommitKey{
		Params: params,
		B0:     sampler.SampleUniformMat(params.Mu, d, rng),
		Bm:     sampler.SampleUniformMat(params.Ell, d, rng),
	}
}

// Commit creates a commitment to the message vector msg using the given key.
// msg must have length ℓ. Returns the commitment and the opening.
func Commit(ck *CommitKey, msg ring.PolyVec, rng io.Reader) (*Commitment, *Opening) {
	if len(msg) != ck.Params.Ell {
		panic("commitment: message length mismatch")
	}

	// Sample randomness r from D_σ^d
	r := sampler.SampleGaussianVec(ck.Params.randDim(), ck.Params.Sigma, rng)

	// t₀ = B₀ · r
	t0 := ring.MatVecMul(ck.B0, r)

	// t_i = ⟨b_i, r⟩ + m_i for each message slot
	tm := ring.NewPolyVec(ck.Params.Ell)
	for i := 0; i < ck.Params.Ell; i++ {
		ip := ring.InnerProduct(ck.Bm[i], r)
		tm[i] = *ring.Add(ip, &msg[i])
	}

	com := &Commitment{T0: t0, Tm: tm}
	opening := &Opening{R: r, Msg: msg}
	return com, opening
}

// Verify checks that an opening is valid for the given commitment.
func Verify(ck *CommitKey, com *Commitment, opening *Opening) bool {
	// Check t₀ = B₀ · r
	expectedT0 := ring.MatVecMul(ck.B0, opening.R)
	if !ring.VecEqual(com.T0, expectedT0) {
		return false
	}

	// Check t_i = ⟨b_i, r⟩ + m_i
	for i := 0; i < ck.Params.Ell; i++ {
		ip := ring.InnerProduct(ck.Bm[i], opening.R)
		expected := ring.Add(ip, &opening.Msg[i])
		if !ring.Equal(&com.Tm[i], expected) {
			return false
		}
	}

	return true
}

// AddCommitments returns the component-wise sum of two commitments.
// This corresponds to a commitment to (m₁ + m₂) with randomness (r₁ + r₂).
func AddCommitments(a, b *Commitment) *Commitment {
	return &Commitment{
		T0: ring.VecAdd(a.T0, b.T0),
		Tm: ring.VecAdd(a.Tm, b.Tm),
	}
}

// ScalarMulCommitment multiplies a commitment by a scalar polynomial s.
// This corresponds to a commitment to s·m with randomness s·r.
func ScalarMulCommitment(c *Commitment, s *ring.Poly) *Commitment {
	return &Commitment{
		T0: ring.VecScalarMul(c.T0, s),
		Tm: ring.VecScalarMul(c.Tm, s),
	}
}

// AddOpenings returns the component-wise sum of two openings.
func AddOpenings(a, b *Opening) *Opening {
	return &Opening{
		R:   ring.VecAdd(a.R, b.R),
		Msg: ring.VecAdd(a.Msg, b.Msg),
	}
}
