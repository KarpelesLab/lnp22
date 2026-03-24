package ring

// PolyVec is a vector of polynomials in R_q.
type PolyVec []Poly

// PolyMat is a matrix of polynomials in R_q (row-major).
type PolyMat []PolyVec

// NewPolyVec creates a zero vector of n polynomials.
func (r *Ring) NewPolyVec(n int) PolyVec {
	v := make(PolyVec, n)
	for i := range v {
		v[i] = r.NewPoly()
	}
	return v
}

// NewPolyMat creates a zero k×l matrix of polynomials.
func (r *Ring) NewPolyMat(k, l int) PolyMat {
	m := make(PolyMat, k)
	for i := range m {
		m[i] = r.NewPolyVec(l)
	}
	return m
}

// VecAdd returns a + b (component-wise polynomial addition).
func (r *Ring) VecAdd(a, b PolyVec) PolyVec {
	if len(a) != len(b) {
		panic("ring: vector length mismatch")
	}
	c := r.NewPolyVec(len(a))
	for i := range a {
		c[i] = r.Add(a[i], b[i])
	}
	return c
}

// VecSub returns a - b (component-wise polynomial subtraction).
func (r *Ring) VecSub(a, b PolyVec) PolyVec {
	if len(a) != len(b) {
		panic("ring: vector length mismatch")
	}
	c := r.NewPolyVec(len(a))
	for i := range a {
		c[i] = r.Sub(a[i], b[i])
	}
	return c
}

// VecScalarMul returns c * a for a scalar polynomial c and vector a.
func (r *Ring) VecScalarMul(a PolyVec, c Poly) PolyVec {
	out := r.NewPolyVec(len(a))
	for i := range a {
		out[i] = r.Mul(c, a[i])
	}
	return out
}

// InnerProduct computes ⟨a, b⟩ = Σ a_i * b_i in R_q.
func (r *Ring) InnerProduct(a, b PolyVec) Poly {
	if len(a) != len(b) {
		panic("ring: vector length mismatch")
	}
	result := r.NewPoly()
	for i := range a {
		prod := r.Mul(a[i], b[i])
		result = r.Add(result, prod)
	}
	return result
}

// MatVecMul computes A * s where A is a k×l matrix and s is an l-vector.
func (r *Ring) MatVecMul(A PolyMat, s PolyVec) PolyVec {
	k := len(A)
	if k == 0 {
		return nil
	}
	l := len(A[0])
	if len(s) != l {
		panic("ring: dimension mismatch in MatVecMul")
	}
	result := r.NewPolyVec(k)
	for i := 0; i < k; i++ {
		result[i] = r.InnerProduct(A[i], s)
	}
	return result
}

// VecInfNorm returns the maximum infinity norm across all polynomials in the vector.
func (r *Ring) VecInfNorm(v PolyVec) int64 {
	maxN := int64(0)
	for i := range v {
		n := r.InfNorm(v[i])
		if n > maxN {
			maxN = n
		}
	}
	return maxN
}

// VecEqual returns true if a and b represent the same polynomial vector.
func (r *Ring) VecEqual(a, b PolyVec) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !r.Equal(a[i], b[i]) {
			return false
		}
	}
	return true
}
