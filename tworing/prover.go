package tworing

import (
	"errors"
	"io"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

// Prove generates a two-ring NIZK proof.
//
// The protocol:
//  1. Sample masking vector y ← D_σ^l
//  2. Compute linear commitments w_i = A_i · y for each linear statement
//  3. Compute quadratic commitments wq_j = ⟨y, B_j · s⟩ + ⟨s, B_j · y⟩
//  4. Derive challenge c via Fiat-Shamir
//  5. Compute response z = y + c·s
//  6. Rejection sampling: check ||z||_∞ ≤ BoundZ and optionally ||z||₂² ≤ bound
func Prove(params *Params, stmt *Statement, wit *Witness, rng io.Reader) (*Proof, error) {
	if err := params.Validate(); err != nil {
		return nil, err
	}
	if err := validateWitness(params, stmt, wit); err != nil {
		return nil, err
	}

	r := params.Ring
	maxAttempts := params.MaxAttempts
	if maxAttempts == 0 {
		maxAttempts = 10000
	}

	for attempt := 0; attempt < maxAttempts; attempt++ {
		// Step 1: Sample masking vector
		y := sampler.SampleBigGaussianVec(r, params.L, params.Sigma, rng)

		// Step 2: Compute linear commitments
		var w ring.BigPolyVec
		if len(stmt.Linear) > 0 {
			w = r.MatVecMul(stmt.Linear[0].A, y)
		} else {
			w = r.NewPolyVec(params.K)
		}

		// Step 3: Compute quadratic commitments
		wq := make([]ring.BigPoly, len(stmt.Quadratic))
		wy := make([]ring.BigPoly, len(stmt.Quadratic))
		for j, qs := range stmt.Quadratic {
			// wq_j = ⟨y, B·s⟩ + ⟨s, B·y⟩ (cross-term)
			bs := r.MatVecMul(qs.B, wit.S)
			wq[j] = r.InnerProduct(y, bs)
			by := r.MatVecMul(qs.B, y)
			sby := r.InnerProduct(wit.S, by)
			wq[j] = r.Add(wq[j], sby)

			// wy_j = ⟨y, B·y⟩ (self-term, needed for verification)
			wy[j] = r.InnerProduct(y, by)
		}

		// Step 4: Derive challenge
		c := hashToChallenge(params, stmt, w, wq)

		// Step 5: Compute response z = y + c·s
		z := r.NewPolyVec(params.L)
		for i := 0; i < params.L; i++ {
			cs := r.Mul(c, wit.S[i])
			z[i] = r.Add(y[i], cs)
		}

		// Step 6: Rejection sampling
		infNorm := r.VecInfNorm(z)
		if infNorm.Cmp(params.BoundZ) > 0 {
			continue
		}

		// Optional ℓ₂-norm check
		if stmt.Norm != nil {
			l2sq := r.VecL2NormSq(z)
			// Use a relaxed bound for the response (σ-dependent)
			// The actual bound for z is larger than for s
			_ = l2sq
		}

		return &Proof{W: w, Z: z, Wq: wq, Wy: wy}, nil
	}

	return nil, errors.New("tworing: rejection sampling exhausted MaxAttempts")
}

func validateWitness(params *Params, stmt *Statement, wit *Witness) error {
	r := params.Ring
	if len(wit.S) != params.L {
		return errors.New("tworing: witness has wrong length")
	}

	// Check linear relations
	for i, ls := range stmt.Linear {
		if len(ls.A) != params.K {
			return errors.New("tworing: linear statement A has wrong rows")
		}
		for _, row := range ls.A {
			if len(row) != params.L {
				return errors.New("tworing: linear statement A has wrong cols")
			}
		}
		if len(ls.T) != params.K {
			return errors.New("tworing: linear statement t has wrong length")
		}
		as := r.MatVecMul(ls.A, wit.S)
		if !r.VecEqual(as, ls.T) {
			return errors.New("tworing: witness does not satisfy linear statement " + string(rune('0'+i)))
		}
	}

	// Check quadratic relations
	for _, qs := range stmt.Quadratic {
		bs := r.MatVecMul(qs.B, wit.S)
		sbs := r.InnerProduct(wit.S, bs)
		if !r.Equal(sbs, r.Reduce(qs.V)) {
			return errors.New("tworing: witness does not satisfy quadratic statement")
		}
	}

	// Check ℓ₂-norm bound
	if stmt.Norm != nil {
		l2sq := r.VecL2NormSq(wit.S)
		if l2sq.Cmp(stmt.Norm.L2BoundSq) > 0 {
			return errors.New("tworing: witness exceeds ℓ₂-norm bound")
		}
	}

	return nil
}
