package ring

import (
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
)

// NTT computes the forward NTT (Cooley-Tukey) for BigRing.
// Panics if NTT is not available.
func (r *BigRing) NTT(a BigPoly) BigPoly {
	if !r.hasNTT {
		panic("ring: NTT not available for this BigRing")
	}
	out := r.NewPoly()
	for i := 0; i < r.N; i++ {
		out[i].Set(a[i])
	}

	tmp := new(big.Int)
	k := 0
	for length := r.N / 2; length >= 1; length >>= 1 {
		for start := 0; start < r.N; start += 2 * length {
			k++
			z := r.zetas[k]
			for j := start; j < start+length; j++ {
				// t = z * out[j+length] mod Q
				tmp.Mul(z, out[j+length])
				tmp.Mod(tmp, r.Q)
				// out[j+length] = out[j] - t
				out[j+length].Sub(out[j], tmp)
				out[j+length].Mod(out[j+length], r.Q)
				// out[j] = out[j] + t
				out[j].Add(out[j], tmp)
				out[j].Mod(out[j], r.Q)
			}
		}
	}
	return out
}

// InvNTT computes the inverse NTT (Gentleman-Sande) for BigRing.
func (r *BigRing) InvNTT(a BigPoly) BigPoly {
	if !r.hasNTT {
		panic("ring: NTT not available for this BigRing")
	}
	out := r.NewPoly()
	for i := 0; i < r.N; i++ {
		out[i].Set(a[i])
	}

	negZ := new(big.Int)
	tmp := new(big.Int)
	k := r.N
	for length := 1; length <= r.N/2; length <<= 1 {
		for start := 0; start < r.N; start += 2 * length {
			k--
			negZ.Neg(r.zetas[k])
			negZ.Mod(negZ, r.Q)
			for j := start; j < start+length; j++ {
				tmp.Set(out[j])
				// out[j] = t + out[j+length]
				out[j].Add(tmp, out[j+length])
				out[j].Mod(out[j], r.Q)
				// out[j+length] = negZ * (t - out[j+length])
				out[j+length].Sub(tmp, out[j+length])
				bigmod.ModMulInPlace(out[j+length], out[j+length], negZ, r.Q)
			}
		}
	}
	// Multiply by N^{-1}
	for i := 0; i < r.N; i++ {
		bigmod.ModMulInPlace(out[i], out[i], r.NInv, r.Q)
	}
	return out
}

// PointwiseMul computes coefficient-wise product in NTT domain.
func (r *BigRing) PointwiseMul(a, b BigPoly) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		bigmod.ModMulInPlace(c[i], a[i], b[i], r.Q)
	}
	return c
}
