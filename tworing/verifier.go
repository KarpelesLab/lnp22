package tworing

// Verify checks a two-ring NIZK proof.
//
// Verification:
//  1. Re-derive challenge c from the transcript
//  2. Check ||z||_∞ ≤ BoundZ
//  3. For each linear statement: check A·z = c·t + w (mod q)
//  4. For each quadratic statement: check ⟨z, B·z⟩ = c²·v + c·wq + ... (linearized)
func Verify(params *Params, stmt *Statement, proof *Proof) bool {
	if err := params.Validate(); err != nil {
		return false
	}

	r := params.Ring

	// Check response dimensions
	if len(proof.Z) != params.L {
		return false
	}
	if len(stmt.Linear) > 0 && len(proof.W) != params.K {
		return false
	}

	// Check infinity norm bound
	infNorm := r.VecInfNorm(proof.Z)
	if infNorm.Cmp(params.BoundZ) > 0 {
		return false
	}

	// Re-derive challenge
	c := hashToChallenge(params, stmt, proof.W, proof.Wq)

	// Verify linear relations: A·z = c·t + w
	for _, ls := range stmt.Linear {
		az := r.MatVecMul(ls.A, proof.Z)
		ct := r.VecScalarMul(ls.T, c)
		expected := r.VecAdd(ct, proof.W)
		if !r.VecEqual(az, expected) {
			return false
		}
	}

	// Verify quadratic relations
	// ⟨z, B·z⟩ = c²·v + c·wq + wy
	// where z = y + c·s, wq = ⟨y,B·s⟩+⟨s,B·y⟩, wy = ⟨y,B·y⟩
	cSq := r.Mul(c, c) // c²
	for j, qs := range stmt.Quadratic {
		if j >= len(proof.Wq) || j >= len(proof.Wy) {
			return false
		}
		// LHS: ⟨z, B·z⟩
		bz := r.MatVecMul(qs.B, proof.Z)
		zbz := r.InnerProduct(proof.Z, bz)

		// RHS: c²·v + c·wq + wy
		c2v := r.Mul(cSq, qs.V)
		cwq := r.Mul(c, proof.Wq[j])
		expected := r.Add(c2v, cwq)
		expected = r.Add(expected, proof.Wy[j])

		if !r.Equal(zbz, expected) {
			return false
		}
	}

	return true
}
