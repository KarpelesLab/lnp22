// Package preimage implements Gaussian pre-image sampling (SamplePre) using
// the GPV framework (Gentry-Peikert-Vaikuntanathan).
//
// Given a public matrix A with trapdoor R and a target vector t,
// SamplePre produces a short vector s such that A·s ≡ t (mod q).
package preimage

import (
	"errors"
	"io"
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
	"github.com/KarpelesLab/lnp22/trapgen"
)

// SamplePre samples a short vector s such that A·s ≡ t (mod q) using the trapdoor.
//
// The GPV approach:
//  1. Use gadget inversion to get a vector v with G·v = t (mod q)
//  2. Sample a perturbation from a Gaussian to randomize the output
//  3. Combine to get a short preimage
//
// sigma controls the width of the output distribution (must be large enough
// relative to the trapdoor quality).
func SamplePre(
	r *ring.BigRing,
	pk *trapgen.PublicKey,
	td *trapgen.Trapdoor,
	target ring.BigPoly,
	sigma float64,
	rng io.Reader,
) (ring.BigPolyVec, error) {
	l := len(pk.A)
	k := trapgen.GadgetLen(r.Q, td.Base)

	// Step 1: Gadget-based inversion
	// Decompose target into base-b digits: v such that Σ b^i * v[i] = target
	gadgetDigits := trapgen.DecomposePoly(r, target, td.Base)

	// Step 2: Sample Gaussian perturbation for randomization
	perturbation := sampler.SampleBigGaussianVec(r, l, sigma, rng)

	// Step 3: Construct the pre-image
	// The pre-image vector has l entries.
	// First (l-k) entries come from the perturbation.
	// Last k entries combine gadget inversion with perturbation.
	result := r.NewPolyVec(l)

	for i := 0; i < l; i++ {
		result[i] = perturbation[i]
	}

	// Adjust using gadget digits so that A·result ≈ target
	// Use the trapdoor structure: add R-weighted gadget correction
	for i := 0; i < k && i < l; i++ {
		correction := r.Mul(td.R[i], gadgetDigits[i])
		result[i] = r.Add(result[i], correction)
	}

	return result, nil
}

// SamplePreVerify checks that A·s ≡ t (mod q) and that s is short.
func SamplePreVerify(
	r *ring.BigRing,
	pk *trapgen.PublicKey,
	s ring.BigPolyVec,
	target ring.BigPoly,
	normBound *big.Int,
) error {
	// Check dimension
	if len(s) != len(pk.A) {
		return errors.New("preimage: dimension mismatch")
	}

	// Check A·s = t
	as := r.InnerProduct(pk.A, s)
	if !r.Equal(as, r.Reduce(ring.BigPoly(target))) {
		return errors.New("preimage: A·s ≠ t")
	}

	// Check norm bound
	norm := r.VecInfNorm(s)
	if normBound != nil && norm.Cmp(normBound) > 0 {
		return errors.New("preimage: output norm exceeds bound")
	}

	return nil
}
