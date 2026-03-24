package nizk

import (
	"errors"
	"io"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

// ProveLinear generates a NIZK proof that the prover knows s such that
// A·s ≡ t (mod q) and ||s||_∞ ≤ β.
//
// The protocol:
//  1. Sample masking vector y ← D_σ^l
//  2. Compute commitment w = A·y (mod q)
//  3. Derive challenge c = H(A, t, w) ∈ C_κ
//  4. Compute response z = y + c·s
//  5. If ||z||_∞ > B_z, reject and restart (rejection sampling)
//  6. Output proof (w, z)
func ProveLinear(params *Params, stmt *Statement, wit *Witness, rng io.Reader) (*LinearProof, error) {
	if err := params.Validate(); err != nil {
		return nil, err
	}
	if err := validateStatement(params, stmt); err != nil {
		return nil, err
	}
	if err := validateWitness(params, stmt, wit); err != nil {
		return nil, err
	}

	maxAttempts := params.MaxAttempts
	if maxAttempts == 0 {
		maxAttempts = 10000
	}

	for attempt := 0; attempt < maxAttempts; attempt++ {
		// Step 1: Sample masking vector y from D_σ^l
		y := sampler.SampleGaussianVec(params.L, params.Sigma, rng)

		// Step 2: Compute commitment w = A·y
		w := ring.MatVecMul(stmt.A, y)

		// Step 3: Derive challenge
		c := HashToChallenge(params, stmt.A, stmt.T, w)

		// Step 4: Compute response z_i = y_i + c·s_i
		z := ring.NewPolyVec(params.L)
		for i := 0; i < params.L; i++ {
			cs := ring.Mul(c, &wit.S[i])
			z[i] = *ring.Add(&y[i], cs)
		}

		// Step 5: Rejection sampling — check norm bound
		if ring.VecInfNorm(z) <= params.BoundZ {
			return &LinearProof{W: w, Z: z}, nil
		}
	}

	return nil, errors.New("nizk: rejection sampling exhausted MaxAttempts")
}

// validateStatement checks that the statement dimensions match the params.
func validateStatement(params *Params, stmt *Statement) error {
	if len(stmt.A) != params.K {
		return errors.New("nizk: matrix A has wrong number of rows")
	}
	for i, row := range stmt.A {
		if len(row) != params.L {
			return errors.New("nizk: matrix A row " + string(rune('0'+i)) + " has wrong number of columns")
		}
	}
	if len(stmt.T) != params.K {
		return errors.New("nizk: target vector t has wrong length")
	}
	return nil
}

// validateWitness checks that the witness is valid: A·s = t and ||s||_∞ ≤ β.
func validateWitness(params *Params, stmt *Statement, wit *Witness) error {
	if len(wit.S) != params.L {
		return errors.New("nizk: witness has wrong length")
	}
	// Check norm bound
	if ring.VecInfNorm(wit.S) > params.Beta {
		return errors.New("nizk: witness exceeds norm bound β")
	}
	// Check A·s = t
	as := ring.MatVecMul(stmt.A, wit.S)
	if !ring.VecEqual(as, stmt.T) {
		return errors.New("nizk: witness does not satisfy A·s = t")
	}
	return nil
}
