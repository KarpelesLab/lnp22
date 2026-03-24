package nizk

// VerifyLinear checks a linear-relation NIZK proof.
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

	r := params.Ring

	if r.VecInfNorm(proof.Z) > params.BoundZ {
		return false
	}

	c := HashToChallenge(params, stmt.A, stmt.T, proof.W)
	az := r.MatVecMul(stmt.A, proof.Z)
	ct := r.VecScalarMul(stmt.T, c)
	expected := r.VecAdd(ct, proof.W)

	return r.VecEqual(az, expected)
}

// VerifyRange checks a range proof.
func VerifyRange(params *Params, stmt *RangeStatement, proof *RangeProof) bool {
	if !VerifyLinear(params, &stmt.Statement, &proof.LinearPart) {
		return false
	}

	bitParams := bitProofParams(params)
	for _, bp := range proof.BitProofs {
		if len(bp.Z) != bitParams.L || len(bp.W) != bitParams.K {
			return false
		}
		if params.Ring.VecInfNorm(bp.Z) > bitParams.BoundZ {
			return false
		}
	}

	reconParams := reconProofParams(params)
	if len(proof.ReconProof.Z) != reconParams.L || len(proof.ReconProof.W) != reconParams.K {
		return false
	}
	if params.Ring.VecInfNorm(proof.ReconProof.Z) > reconParams.BoundZ {
		return false
	}

	return true
}

// VerifyComposed checks a composed proof.
func VerifyComposed(params *Params, stmts []*Statement, proof *ComposedProof) bool {
	if len(stmts) != len(proof.Proofs) {
		return false
	}

	r := params.Ring

	for i, subProof := range proof.Proofs {
		if len(subProof.Z) != params.L || len(subProof.W) != params.K {
			return false
		}
		if r.VecInfNorm(subProof.Z) > params.BoundZ {
			return false
		}

		tag := composeTag(proof.TranscriptSeed, i)
		seed := hashTranscriptWithTag(r, tag, stmts[i].A, stmts[i].T, subProof.W)
		c := hashSeedToChallenge(params, seed)

		az := r.MatVecMul(stmts[i].A, subProof.Z)
		ct := r.VecScalarMul(stmts[i].T, c)
		expected := r.VecAdd(ct, subProof.W)

		if !r.VecEqual(az, expected) {
			return false
		}
	}

	return true
}
