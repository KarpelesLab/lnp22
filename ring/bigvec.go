package ring

import "math/big"

// BigPolyVec is a vector of big.Int polynomials.
type BigPolyVec []BigPoly

// BigPolyMat is a matrix of big.Int polynomials (row-major).
type BigPolyMat []BigPolyVec

// NewPolyVec creates a zero vector of n polynomials.
func (r *BigRing) NewPolyVec(n int) BigPolyVec {
	v := make(BigPolyVec, n)
	for i := range v {
		v[i] = r.NewPoly()
	}
	return v
}

// NewPolyMat creates a zero k×l matrix.
func (r *BigRing) NewPolyMat(k, l int) BigPolyMat {
	m := make(BigPolyMat, k)
	for i := range m {
		m[i] = r.NewPolyVec(l)
	}
	return m
}

// VecAdd returns a + b.
func (r *BigRing) VecAdd(a, b BigPolyVec) BigPolyVec {
	if len(a) != len(b) {
		panic("ring: vector length mismatch")
	}
	c := r.NewPolyVec(len(a))
	for i := range a {
		c[i] = r.Add(a[i], b[i])
	}
	return c
}

// VecSub returns a - b.
func (r *BigRing) VecSub(a, b BigPolyVec) BigPolyVec {
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
func (r *BigRing) VecScalarMul(a BigPolyVec, c BigPoly) BigPolyVec {
	out := r.NewPolyVec(len(a))
	for i := range a {
		out[i] = r.Mul(c, a[i])
	}
	return out
}

// InnerProduct computes ⟨a, b⟩ = Σ a_i * b_i.
func (r *BigRing) InnerProduct(a, b BigPolyVec) BigPoly {
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

// MatVecMul computes A * s.
func (r *BigRing) MatVecMul(A BigPolyMat, s BigPolyVec) BigPolyVec {
	k := len(A)
	if k == 0 {
		return nil
	}
	if len(s) != len(A[0]) {
		panic("ring: dimension mismatch in MatVecMul")
	}
	result := r.NewPolyVec(k)
	for i := 0; i < k; i++ {
		result[i] = r.InnerProduct(A[i], s)
	}
	return result
}

// VecInfNorm returns the maximum infinity norm across all polynomials.
func (r *BigRing) VecInfNorm(v BigPolyVec) *big.Int {
	maxN := new(big.Int)
	for i := range v {
		n := r.InfNorm(v[i])
		if n.Cmp(maxN) > 0 {
			maxN.Set(n)
		}
	}
	return maxN
}

// VecEqual returns true if a and b are the same.
func (r *BigRing) VecEqual(a, b BigPolyVec) bool {
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
