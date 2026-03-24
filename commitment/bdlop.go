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
	Ring   *ring.Ring
	Mu     int     // number of binding equations (rows in B₀)
	Lambda int     // randomness dimension overhead
	Ell    int     // number of message slots
	Sigma  float64 // Gaussian std dev for randomness
}

// DefaultParams returns default commitment parameters for the given ring.
func DefaultParams(r *ring.Ring) *Params {
	return &Params{
		Ring:   r,
		Mu:     3,
		Lambda: 3,
		Ell:    4,
		Sigma:  50,
	}
}

func (p *Params) randDim() int {
	return p.Mu + p.Lambda + p.Ell
}

// CommitKey is the public commitment key.
type CommitKey struct {
	Params *Params
	B0     ring.PolyMat // μ × (μ+λ+ℓ)
	Bm     ring.PolyMat // ℓ × (μ+λ+ℓ)
}

// Commitment is the output of the commitment function.
type Commitment struct {
	T0 ring.PolyVec // binding component (μ polynomials)
	Tm ring.PolyVec // message components (ℓ polynomials)
}

// Opening contains the information needed to verify a commitment.
type Opening struct {
	R   ring.PolyVec // randomness
	Msg ring.PolyVec // message
}

// KeyGen generates a random commitment key.
func KeyGen(params *Params, rng io.Reader) *CommitKey {
	r := params.Ring
	d := params.randDim()
	return &CommitKey{
		Params: params,
		B0:     sampler.SampleUniformMat(r, params.Mu, d, rng),
		Bm:     sampler.SampleUniformMat(r, params.Ell, d, rng),
	}
}

// Commit creates a commitment to the message vector msg.
func Commit(ck *CommitKey, msg ring.PolyVec, rng io.Reader) (*Commitment, *Opening) {
	if len(msg) != ck.Params.Ell {
		panic("commitment: message length mismatch")
	}
	r := ck.Params.Ring

	rv := sampler.SampleGaussianVec(r, ck.Params.randDim(), ck.Params.Sigma, rng)
	t0 := r.MatVecMul(ck.B0, rv)

	tm := r.NewPolyVec(ck.Params.Ell)
	for i := 0; i < ck.Params.Ell; i++ {
		ip := r.InnerProduct(ck.Bm[i], rv)
		tm[i] = r.Add(ip, msg[i])
	}

	return &Commitment{T0: t0, Tm: tm}, &Opening{R: rv, Msg: msg}
}

// Verify checks that an opening is valid for the given commitment.
func Verify(ck *CommitKey, com *Commitment, opening *Opening) bool {
	r := ck.Params.Ring

	expectedT0 := r.MatVecMul(ck.B0, opening.R)
	if !r.VecEqual(com.T0, expectedT0) {
		return false
	}

	for i := 0; i < ck.Params.Ell; i++ {
		ip := r.InnerProduct(ck.Bm[i], opening.R)
		expected := r.Add(ip, opening.Msg[i])
		if !r.Equal(com.Tm[i], expected) {
			return false
		}
	}
	return true
}

// AddCommitments returns the component-wise sum of two commitments.
func AddCommitments(r *ring.Ring, a, b *Commitment) *Commitment {
	return &Commitment{
		T0: r.VecAdd(a.T0, b.T0),
		Tm: r.VecAdd(a.Tm, b.Tm),
	}
}

// ScalarMulCommitment multiplies a commitment by a scalar polynomial.
func ScalarMulCommitment(r *ring.Ring, c *Commitment, s ring.Poly) *Commitment {
	return &Commitment{
		T0: r.VecScalarMul(c.T0, s),
		Tm: r.VecScalarMul(c.Tm, s),
	}
}

// AddOpenings returns the component-wise sum of two openings.
func AddOpenings(r *ring.Ring, a, b *Opening) *Opening {
	return &Opening{
		R:   r.VecAdd(a.R, b.R),
		Msg: r.VecAdd(a.Msg, b.Msg),
	}
}
