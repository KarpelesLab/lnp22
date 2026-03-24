package ring

import "github.com/KarpelesLab/lnp22/internal/modmath"

// NTT computes the Number Theoretic Transform of polynomial a (Cooley-Tukey butterfly).
// The result is in the NTT domain. Coefficients are reduced to [0, Q).
func NTT(a *Poly) *Poly {
	var out Poly = *a
	k := 0
	for length := N / 2; length >= 1; length >>= 1 {
		for start := 0; start < N; start += 2 * length {
			k++
			z := zetas[k]
			for j := start; j < start+length; j++ {
				t := modmath.ModMul(z, out[j+length], Q)
				out[j+length] = modmath.Reduce(out[j]-t, Q)
				out[j] = modmath.Reduce(out[j]+t, Q)
			}
		}
	}
	return &out
}

// InvNTT computes the inverse NTT (Gentleman-Sande butterfly).
// Transforms from NTT domain back to coefficient domain.
func InvNTT(a *Poly) *Poly {
	var out Poly = *a
	k := N
	for length := 1; length <= N/2; length <<= 1 {
		for start := 0; start < N; start += 2 * length {
			k--
			z := modmath.Reduce(-zetas[k], Q)
			for j := start; j < start+length; j++ {
				t := out[j]
				out[j] = modmath.Reduce(t+out[j+length], Q)
				out[j+length] = modmath.ModMul(z, modmath.Reduce(t-out[j+length], Q), Q)
			}
		}
	}
	// Multiply by N^{-1} mod Q
	for i := 0; i < N; i++ {
		out[i] = modmath.ModMul(out[i], NInv, Q)
	}
	return &out
}

// PointwiseMul computes the coefficient-wise product of two polynomials in NTT domain.
// Both inputs must already be in NTT representation.
func PointwiseMul(a, b *Poly) *Poly {
	var c Poly
	for i := 0; i < N; i++ {
		c[i] = modmath.ModMul(a[i], b[i], Q)
	}
	return &c
}

// Mul computes a*b in R_q using NTT: Mul(a,b) = InvNTT(NTT(a) ⊙ NTT(b)).
func Mul(a, b *Poly) *Poly {
	aNTT := NTT(a)
	bNTT := NTT(b)
	cNTT := PointwiseMul(aNTT, bNTT)
	return InvNTT(cNTT)
}
