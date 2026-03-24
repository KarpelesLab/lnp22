package nizk

import (
	"errors"
	"io"
	"math/bits"

	"github.com/KarpelesLab/lnp22/internal/modmath"
	"github.com/KarpelesLab/lnp22/ring"
)

// constantPoly returns a polynomial where every coefficient is v (mod Q).
func constantPoly(r *ring.Ring, v int64) ring.Poly {
	p := r.NewPoly()
	for i := range p {
		p[i] = modmath.Reduce(v, r.Q)
	}
	return p
}

// ProveRange generates a NIZK proof that the witness coefficients lie in [-β, β].
func ProveRange(params *Params, stmt *RangeStatement, wit *Witness, rng io.Reader) (*RangeProof, error) {
	if err := params.Validate(); err != nil {
		return nil, err
	}

	r := params.Ring

	linearProof, err := ProveLinear(params, &stmt.Statement, wit, rng)
	if err != nil {
		return nil, err
	}

	numBits := stmt.NumBits
	if numBits == 0 {
		numBits = bits.Len(uint(2 * stmt.Beta))
	}

	// Binary decomposition of shifted witness
	bitPolys := make([][]ring.Poly, params.L)
	for j := 0; j < params.L; j++ {
		bitPolys[j] = make([]ring.Poly, numBits)
		for b := 0; b < numBits; b++ {
			bitPolys[j][b] = r.NewPoly()
		}
		for k := 0; k < r.N; k++ {
			coeff := modmath.CenterReduce(wit.S[j][k], r.Q)
			shifted := coeff + stmt.Beta
			if shifted < 0 || shifted > 2*stmt.Beta {
				return nil, errors.New("nizk: witness coefficient out of range for range proof")
			}
			for b := 0; b < numBits; b++ {
				bitPolys[j][b][k] = (shifted >> b) & 1
			}
		}
	}

	// Prove each bit polynomial has ||b||_∞ ≤ 1
	bitParams := bitProofParams(params)
	bitProofs := make([]LinearProof, numBits*params.L)
	for j := 0; j < params.L; j++ {
		for b := 0; b < numBits; b++ {
			bp := bitPolys[j][b]
			bitStmt := &Statement{
				A: ring.PolyMat{ring.PolyVec{r.One()}},
				T: ring.PolyVec{bp},
			}
			bitWit := &Witness{S: ring.PolyVec{bp}}
			proof, err := ProveLinear(bitParams, bitStmt, bitWit, rng)
			if err != nil {
				return nil, err
			}
			bitProofs[j*numBits+b] = *proof
		}
	}

	// Prove reconstruction: Σ_i 2^i · b_i = s + β for each witness polynomial
	reconParams := reconProofParams(params)
	reconProofs := make([]LinearProof, params.L)
	for j := 0; j < params.L; j++ {
		reconWitVec := r.NewPolyVec(numBits)
		for b := 0; b < numBits; b++ {
			reconWitVec[b] = bitPolys[j][b]
		}

		betaFill := constantPoly(r, stmt.Beta)
		reconTarget := r.Add(wit.S[j], betaFill)

		reconA := r.NewPolyMat(1, numBits)
		for b := 0; b < numBits; b++ {
			reconA[0][b] = r.ScalarMul(r.One(), 1<<b)
		}

		reconStmt := &Statement{A: reconA, T: ring.PolyVec{reconTarget}}
		reconWit := &Witness{S: reconWitVec}

		proof, err := ProveLinear(reconParams, reconStmt, reconWit, rng)
		if err != nil {
			return nil, err
		}
		reconProofs[j] = *proof
	}

	return &RangeProof{
		LinearPart: *linearProof,
		BitProofs:  bitProofs,
		ReconProof: reconProofs[0],
	}, nil
}

func bitProofParams(params *Params) *Params {
	return &Params{
		Ring:        params.Ring,
		K:           1,
		L:           1,
		Kappa:       params.Kappa,
		Beta:        1,
		Sigma:       params.Sigma,
		BoundZ:      params.BoundZ,
		MaxAttempts: params.MaxAttempts,
	}
}

func reconProofParams(params *Params) *Params {
	numBits := bits.Len(uint(2 * params.Beta))
	return &Params{
		Ring:        params.Ring,
		K:           1,
		L:           numBits,
		Kappa:       params.Kappa,
		Beta:        1,
		Sigma:       params.Sigma,
		BoundZ:      params.BoundZ,
		MaxAttempts: params.MaxAttempts,
	}
}
