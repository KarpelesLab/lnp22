package ring

import "github.com/KarpelesLab/lnp22/internal/modmath"

// Poly represents a polynomial in R_q = Z_q[X]/(X^N+1).
// Coefficients are stored as int64 values.
type Poly []int64

// NewPoly returns the zero polynomial for this ring.
func (r *Ring) NewPoly() Poly {
	return make(Poly, r.N)
}

// One returns the polynomial 1.
func (r *Ring) One() Poly {
	p := r.NewPoly()
	p[0] = 1
	return p
}

// Add returns a + b mod Q (coefficient-wise).
func (r *Ring) Add(a, b Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.Reduce(a[i]+b[i], r.Q)
	}
	return c
}

// Sub returns a - b mod Q (coefficient-wise).
func (r *Ring) Sub(a, b Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.Reduce(a[i]-b[i], r.Q)
	}
	return c
}

// Neg returns -a mod Q (coefficient-wise).
func (r *Ring) Neg(a Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		v := modmath.Reduce(a[i], r.Q)
		if v != 0 {
			c[i] = r.Q - v
		}
	}
	return c
}

// ScalarMul returns s*a mod Q (multiply each coefficient by scalar).
func (r *Ring) ScalarMul(a Poly, s int64) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.ModMul(a[i], s, r.Q)
	}
	return c
}

// Reduce ensures all coefficients are in [0, Q).
func (r *Ring) Reduce(a Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.Reduce(a[i], r.Q)
	}
	return c
}

// CenterReduce returns a copy with coefficients in [-(Q-1)/2, (Q-1)/2].
func (r *Ring) CenterReduce(a Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.CenterReduce(a[i], r.Q)
	}
	return c
}

// InfNorm returns ||a||_∞ (the maximum absolute value of centered coefficients).
func (r *Ring) InfNorm(a Poly) int64 {
	maxVal := int64(0)
	for i := 0; i < r.N; i++ {
		v := modmath.CenterReduce(a[i], r.Q)
		if v < 0 {
			v = -v
		}
		if v > maxVal {
			maxVal = v
		}
	}
	return maxVal
}

// Equal returns true if a and b represent the same polynomial in R_q.
func (r *Ring) Equal(a, b Poly) bool {
	for i := 0; i < r.N; i++ {
		if modmath.Reduce(a[i], r.Q) != modmath.Reduce(b[i], r.Q) {
			return false
		}
	}
	return true
}

// Mul computes a*b in R_q = Z_q[X]/(X^N+1).
// Uses NTT if available, otherwise falls back to Karatsuba.
func (r *Ring) Mul(a, b Poly) Poly {
	if r.hasNTT {
		aNTT := r.NTT(a)
		bNTT := r.NTT(b)
		cNTT := r.PointwiseMul(aNTT, bNTT)
		return r.InvNTT(cNTT)
	}
	return r.karatsubaMul(a, b)
}

// NaiveMul computes a*b in R_q using the schoolbook algorithm.
// O(N²), used for testing.
func (r *Ring) NaiveMul(a, b Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		for j := 0; j < r.N; j++ {
			prod := modmath.ModMul(a[i], b[j], r.Q)
			idx := i + j
			if idx < r.N {
				c[idx] = modmath.Reduce(c[idx]+prod, r.Q)
			} else {
				// X^N ≡ -1, so X^(N+k) ≡ -X^k
				c[idx-r.N] = modmath.Reduce(c[idx-r.N]-prod, r.Q)
			}
		}
	}
	return c
}

// karatsubaMul computes a*b in Z_q[X] using Karatsuba, then reduces mod X^N+1.
func (r *Ring) karatsubaMul(a, b Poly) Poly {
	// Multiply as polynomials in Z_q[X] (degree up to 2N-2)
	product := r.karatsubaRaw(a, b)

	// Reduce mod X^N + 1: for terms X^(N+i), subtract from coefficient i
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.Reduce(product[i]-product[i+r.N], r.Q)
	}
	return c
}

// karatsubaRaw multiplies two coefficient slices using Karatsuba and returns
// a slice of length len(a)+len(b)-1. Falls back to schoolbook for small inputs.
func (r *Ring) karatsubaRaw(a, b Poly) []int64 {
	n := len(a)
	result := make([]int64, 2*n)

	if n <= 32 {
		// Schoolbook for small sizes
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				result[i+j] = modmath.Reduce(result[i+j]+modmath.ModMul(a[i], b[j], r.Q), r.Q)
			}
		}
		return result
	}

	half := n / 2
	a0, a1 := Poly(a[:half]), Poly(a[half:])
	b0, b1 := Poly(b[:half]), Poly(b[half:])

	// z0 = a0 * b0
	z0 := r.karatsubaRaw(a0, b0)
	// z2 = a1 * b1
	z2 := r.karatsubaRaw(a1, b1)
	// z1 = (a0+a1)*(b0+b1) - z0 - z2
	aSum := make(Poly, half)
	bSum := make(Poly, half)
	for i := 0; i < half; i++ {
		aSum[i] = modmath.Reduce(a0[i]+a1[i], r.Q)
		bSum[i] = modmath.Reduce(b0[i]+b1[i], r.Q)
	}
	z1 := r.karatsubaRaw(aSum, bSum)

	// result = z0 + z1*X^half + z2*X^n
	for i := range z0 {
		result[i] = modmath.Reduce(result[i]+z0[i], r.Q)
	}
	for i := range z2 {
		result[i+n] = modmath.Reduce(result[i+n]+z2[i], r.Q)
	}
	for i := range z1 {
		v := modmath.Reduce(z1[i]-z0[i]-z2[i], r.Q)
		result[i+half] = modmath.Reduce(result[i+half]+v, r.Q)
	}

	return result
}
