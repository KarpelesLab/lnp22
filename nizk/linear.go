package nizk

import (
	"errors"
	"io"

	"github.com/KarpelesLab/lnp22/sampler"
)

// ProveLinear generates a NIZK proof that the prover knows s such that
// A·s ≡ t (mod q) and ||s||_∞ ≤ β.
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

	r := params.Ring
	maxAttempts := params.MaxAttempts
	if maxAttempts == 0 {
		maxAttempts = 10000
	}

	for attempt := 0; attempt < maxAttempts; attempt++ {
		y := sampler.SampleGaussianVec(r, params.L, params.Sigma, rng)
		w := r.MatVecMul(stmt.A, y)
		c := HashToChallenge(params, stmt.A, stmt.T, w)

		z := r.NewPolyVec(params.L)
		for i := 0; i < params.L; i++ {
			cs := r.Mul(c, wit.S[i])
			z[i] = r.Add(y[i], cs)
		}

		if r.VecInfNorm(z) <= params.BoundZ {
			return &LinearProof{W: w, Z: z}, nil
		}
	}

	return nil, errors.New("nizk: rejection sampling exhausted MaxAttempts")
}

func validateStatement(params *Params, stmt *Statement) error {
	if len(stmt.A) != params.K {
		return errors.New("nizk: matrix A has wrong number of rows")
	}
	for _, row := range stmt.A {
		if len(row) != params.L {
			return errors.New("nizk: matrix A has wrong number of columns")
		}
	}
	if len(stmt.T) != params.K {
		return errors.New("nizk: target vector t has wrong length")
	}
	return nil
}

func validateWitness(params *Params, stmt *Statement, wit *Witness) error {
	r := params.Ring
	if len(wit.S) != params.L {
		return errors.New("nizk: witness has wrong length")
	}
	if r.VecInfNorm(wit.S) > params.Beta {
		return errors.New("nizk: witness exceeds norm bound β")
	}
	as := r.MatVecMul(stmt.A, wit.S)
	if !r.VecEqual(as, stmt.T) {
		return errors.New("nizk: witness does not satisfy A·s = t")
	}
	return nil
}
