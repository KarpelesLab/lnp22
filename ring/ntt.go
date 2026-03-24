package ring

import "github.com/KarpelesLab/lnp22/internal/modmath"

// NTT computes the Number Theoretic Transform (Cooley-Tukey butterfly).
// Panics if NTT is not available for this ring.
func (r *Ring) NTT(a Poly) Poly {
	if !r.hasNTT {
		panic("ring: NTT not available (Q ≢ 1 mod 2N)")
	}
	out := make(Poly, r.N)
	copy(out, a)
	k := 0
	for length := r.N / 2; length >= 1; length >>= 1 {
		for start := 0; start < r.N; start += 2 * length {
			k++
			z := r.zetas[k]
			for j := start; j < start+length; j++ {
				t := modmath.ModMul(z, out[j+length], r.Q)
				out[j+length] = modmath.Reduce(out[j]-t, r.Q)
				out[j] = modmath.Reduce(out[j]+t, r.Q)
			}
		}
	}
	return out
}

// InvNTT computes the inverse NTT (Gentleman-Sande butterfly).
// Panics if NTT is not available for this ring.
func (r *Ring) InvNTT(a Poly) Poly {
	if !r.hasNTT {
		panic("ring: NTT not available (Q ≢ 1 mod 2N)")
	}
	out := make(Poly, r.N)
	copy(out, a)
	k := r.N
	for length := 1; length <= r.N/2; length <<= 1 {
		for start := 0; start < r.N; start += 2 * length {
			k--
			z := modmath.Reduce(-r.zetas[k], r.Q)
			for j := start; j < start+length; j++ {
				t := out[j]
				out[j] = modmath.Reduce(t+out[j+length], r.Q)
				out[j+length] = modmath.ModMul(z, modmath.Reduce(t-out[j+length], r.Q), r.Q)
			}
		}
	}
	for i := 0; i < r.N; i++ {
		out[i] = modmath.ModMul(out[i], r.NInv, r.Q)
	}
	return out
}

// PointwiseMul computes the coefficient-wise product of two polynomials in NTT domain.
func (r *Ring) PointwiseMul(a, b Poly) Poly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = modmath.ModMul(a[i], b[i], r.Q)
	}
	return c
}
