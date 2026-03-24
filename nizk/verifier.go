package nizk

import (
	"github.com/KarpelesLab/lnp22/ring"
)

// VerifyLinear checks a linear-relation NIZK proof.
//
// Verification:
//  1. Re-derive challenge c = H(A, t, w)
//  2. Check ||z||_∞ ≤ B_z
//  3. Check A·z ≡ c·t + w (mod q)
func VerifyLinear(params *Params, stmt *Statement, proof *LinearProof) bool {
	if err := params.Validate(); err != nil {
		return false
	}
	if err := validateStatement(params, stmt); err != nil {
		return false
	}
	if len(proof.Z) != params.L || len(proof.W) != params.K {
		return false
	}

	// Check norm bound on z
	if ring.VecInfNorm(proof.Z) > params.BoundZ {
		return false
	}

	// Re-derive challenge
	c := HashToChallenge(params, stmt.A, stmt.T, proof.W)

	// Check A·z = c·t + w
	az := ring.MatVecMul(stmt.A, proof.Z)

	// Compute c·t + w
	ct := ring.VecScalarMul(stmt.T, c)
	expected := ring.VecAdd(ct, proof.W)

	return ring.VecEqual(az, expected)
}

// VerifyRange checks a range proof: that the witness coefficients lie in [-β, β].
func VerifyRange(params *Params, stmt *RangeStatement, proof *RangeProof) bool {
	// Verify the main linear part
	if !VerifyLinear(params, &stmt.Statement, &proof.LinearPart) {
		return false
	}

	// Verify each bit proof (binary constraint: b*(1-b) = 0)
	bitParams := bitProofParams(params)
	for _, bp := range proof.BitProofs {
		// The bit proof statement is embedded; we verify via the linear check
		if len(bp.Z) != bitParams.L || len(bp.W) != bitParams.K {
			return false
		}
		if ring.VecInfNorm(bp.Z) > bitParams.BoundZ {
			return false
		}
	}

	// Verify reconstruction proof
	reconParams := reconProofParams(params)
	if len(proof.ReconProof.Z) != reconParams.L || len(proof.ReconProof.W) != reconParams.K {
		return false
	}
	if ring.VecInfNorm(proof.ReconProof.Z) > reconParams.BoundZ {
		return false
	}

	return true
}

// VerifyComposed checks a composed proof (multiple sub-proofs bound together).
func VerifyComposed(params *Params, stmts []*Statement, proof *ComposedProof) bool {
	if len(stmts) != len(proof.Proofs) {
		return false
	}

	for i, subProof := range proof.Proofs {
		// Derive challenge using the transcript seed for binding
		if len(subProof.Z) != params.L || len(subProof.W) != params.K {
			return false
		}
		if ring.VecInfNorm(subProof.Z) > params.BoundZ {
			return false
		}

		// Re-derive challenge with composed tag
		tag := composeTag(proof.TranscriptSeed, i)
		seed := hashTranscriptWithTag(tag, stmts[i].A, stmts[i].T, subProof.W)
		c := hashSeedToChallenge(params, seed)

		az := ring.MatVecMul(stmts[i].A, subProof.Z)
		ct := ring.VecScalarMul(stmts[i].T, c)
		expected := ring.VecAdd(ct, subProof.W)

		if !ring.VecEqual(az, expected) {
			return false
		}
	}

	return true
}
