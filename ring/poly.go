package ring

import "github.com/KarpelesLab/lnp22/internal/modmath"

// Poly represents a polynomial in R_q = Z_q[X]/(X^N+1).
// Coefficients are stored in the range [0, Q) in standard (non-NTT) representation.
type Poly [N]int64

// Zero returns the zero polynomial.
func Zero() *Poly {
	return &Poly{}
}

// One returns the polynomial 1.
func One() *Poly {
	var p Poly
	p[0] = 1
	return &p
}

// Add returns a + b mod Q (coefficient-wise).
func Add(a, b *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		c[i] = modmath.Reduce(a[i]+b[i], Q)
	}
	return &c
}

// Sub returns a - b mod Q (coefficient-wise).
func Sub(a, b *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		c[i] = modmath.Reduce(a[i]-b[i], Q)
	}
	return &c
}

// Neg returns -a mod Q (coefficient-wise).
func Neg(a *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		if a[i] != 0 {
			c[i] = Q - modmath.Reduce(a[i], Q)
		}
	}
	return &c
}

// ScalarMul returns s*a mod Q (multiply each coefficient by scalar).
func ScalarMul(a *Poly, s int64) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		c[i] = modmath.ModMul(a[i], s, Q)
	}
	return &c
}

// Reduce ensures all coefficients are in [0, Q).
func Reduce(a *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		c[i] = modmath.Reduce(a[i], Q)
	}
	return &c
}

// CenterReduce returns a copy with coefficients in [-(Q-1)/2, (Q-1)/2].
func CenterReduce(a *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		c[i] = modmath.CenterReduce(a[i], Q)
	}
	return &c
}

// InfNorm returns ||a||_∞ (the maximum absolute value of centered coefficients).
func InfNorm(a *Poly) int64 {
	maxVal := int64(0)
	for i := 0; i < N; i++ {
		v := modmath.CenterReduce(a[i], Q)
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
func Equal(a, b *Poly) bool {
	for i := 0; i < N; i++ {
		if modmath.Reduce(a[i], Q) != modmath.Reduce(b[i], Q) {
			return false
		}
	}
	return true
}

// NaiveMul computes a*b in R_q = Z_q[X]/(X^N+1) using the schoolbook algorithm.
// This is O(N²) and used only for testing the NTT-based multiplication.
func NaiveMul(a, b *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			idx := i + j
			prod := modmath.ModMul(a[i], b[j], Q)
			if idx < N {
				c[idx] = modmath.Reduce(c[idx]+prod, Q)
			} else {
				// X^N ≡ -1, so X^(N+k) ≡ -X^k
				c[idx-N] = modmath.Reduce(c[idx-N]-prod, Q)
			}
		}
	}
	return &c
}
