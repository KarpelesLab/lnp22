package nizk

import (
	"errors"
	"io"
	"math/bits"

	"github.com/KarpelesLab/lnp22/internal/modmath"
	"github.com/KarpelesLab/lnp22/ring"
)

// constantPoly returns a polynomial where every coefficient is v (mod Q).
func constantPoly(v int64) *ring.Poly {
	var p ring.Poly
	for i := range p {
		p[i] = modmath.Reduce(v, ring.Q)
	}
	return &p
}

// ProveRange generates a NIZK proof that the witness coefficients lie in [-β, β].
//
// The approach:
//  1. Prove the main linear relation A·s = t
//  2. Binary decomposition: shift s by β so each coefficient is in [0, 2β],
//     then decompose into bits: (s_j[k] + β) = Σ_i b_{j,i}[k] · 2^i
//  3. For each bit polynomial b_i, prove knowledge of b with ||b||_∞ ≤ 1
//  4. Prove reconstruction: Σ_i 2^i · b_i = s + β (coefficient-wise shift)
func ProveRange(params *Params, stmt *RangeStatement, wit *Witness, rng io.Reader) (*RangeProof, error) {
	if err := params.Validate(); err != nil {
		return nil, err
	}

	// Step 1: Prove the main linear relation
	linearProof, err := ProveLinear(params, &stmt.Statement, wit, rng)
	if err != nil {
		return nil, err
	}

	numBits := stmt.NumBits
	if numBits == 0 {
		numBits = bits.Len(uint(2 * stmt.Beta))
	}

	// Step 2: Compute binary decomposition of shifted witness
	// For each witness polynomial s_j, compute b_{i,j} polynomials where
	// (s_j[k] + β) = Σ_i b_{i,j}[k] · 2^i for each coefficient k
	bitPolys := make([][]ring.Poly, params.L) // [witness_idx][bit_idx]
	for j := 0; j < params.L; j++ {
		bitPolys[j] = make([]ring.Poly, numBits)
		for k := 0; k < ring.N; k++ {
			coeff := modmath.CenterReduce(wit.S[j][k], ring.Q)
			shifted := coeff + stmt.Beta // now in [0, 2β]
			if shifted < 0 || shifted > 2*stmt.Beta {
				return nil, errors.New("nizk: witness coefficient out of range for range proof")
			}
			for i := 0; i < numBits; i++ {
				bitPolys[j][i][k] = (shifted >> i) & 1
			}
		}
	}

	// Step 3: For each bit polynomial, prove knowledge of b with ||b||_∞ ≤ 1
	bitParams := bitProofParams(params)
	bitProofs := make([]LinearProof, numBits*params.L)
	for j := 0; j < params.L; j++ {
		for i := 0; i < numBits; i++ {
			bp := &bitPolys[j][i]
			// Statement: I · b = b (the key constraint is ||b||_∞ ≤ 1)
			bitStmt := &Statement{
				A: ring.PolyMat{ring.PolyVec{*ring.One()}},
				T: ring.PolyVec{*bp},
			}
			bitWit := &Witness{S: ring.PolyVec{*bp}}
			proof, err := ProveLinear(bitParams, bitStmt, bitWit, rng)
			if err != nil {
				return nil, err
			}
			bitProofs[j*numBits+i] = *proof
		}
	}

	// Step 4: Prove reconstruction for each witness polynomial
	// Σ_i 2^i · b_{i,j} = s_j + β (coefficient-wise)
	// The target uses constantPoly(β) so that every coefficient is shifted by β,
	// not just the constant term.
	reconParams := reconProofParams(params)
	reconProofs := make([]LinearProof, params.L)
	for j := 0; j < params.L; j++ {
		reconWitVec := ring.NewPolyVec(numBits)
		for i := 0; i < numBits; i++ {
			reconWitVec[i] = bitPolys[j][i]
		}

		// Target = s_j + β (all coefficients shifted)
		betaFill := constantPoly(stmt.Beta)
		sj := wit.S[j]
		reconTarget := ring.Add(&sj, betaFill)

		// A = [[1, 2, 4, ...]] — powers of 2 as constant polynomials
		reconA := ring.NewPolyMat(1, numBits)
		for i := 0; i < numBits; i++ {
			reconA[0][i] = *ring.ScalarMul(ring.One(), 1<<i)
		}

		reconStmt := &Statement{A: reconA, T: ring.PolyVec{*reconTarget}}
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
		ReconProof: reconProofs[0], // primary reconstruction proof
	}, nil
}

// bitProofParams returns parameters for proving that a polynomial is binary.
func bitProofParams(params *Params) *Params {
	return &Params{
		K:           1,
		L:           1,
		Kappa:       params.Kappa,
		Beta:        1,
		Sigma:       params.Sigma,
		BoundZ:      params.BoundZ,
		MaxAttempts: params.MaxAttempts,
	}
}

// reconProofParams returns parameters for the reconstruction proof.
func reconProofParams(params *Params) *Params {
	numBits := bits.Len(uint(2 * params.Beta))
	return &Params{
		K:           1,
		L:           numBits,
		Kappa:       params.Kappa,
		Beta:        1,
		Sigma:       params.Sigma,
		BoundZ:      params.BoundZ,
		MaxAttempts: params.MaxAttempts,
	}
}
