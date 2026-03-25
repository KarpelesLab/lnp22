package ring

import (
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
)

// BigPoly represents a polynomial in R_q with big.Int coefficients.
type BigPoly []*big.Int

// NewPoly returns the zero polynomial.
func (r *BigRing) NewPoly() BigPoly {
	p := make(BigPoly, r.N)
	for i := range p {
		p[i] = new(big.Int)
	}
	return p
}

// One returns the polynomial 1.
func (r *BigRing) One() BigPoly {
	p := r.NewPoly()
	p[0].SetInt64(1)
	return p
}

// Add returns a + b mod Q.
func (r *BigRing) Add(a, b BigPoly) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		bigmod.ModAddInPlace(c[i], a[i], b[i], r.Q)
	}
	return c
}

// Sub returns a - b mod Q.
func (r *BigRing) Sub(a, b BigPoly) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		bigmod.ModSubInPlace(c[i], a[i], b[i], r.Q)
	}
	return c
}

// Neg returns -a mod Q.
func (r *BigRing) Neg(a BigPoly) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		if a[i].Sign() != 0 {
			c[i].Sub(r.Q, bigmod.Reduce(a[i], r.Q))
		}
	}
	return c
}

// ScalarMul returns s*a mod Q.
func (r *BigRing) ScalarMul(a BigPoly, s *big.Int) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		bigmod.ModMulInPlace(c[i], a[i], s, r.Q)
	}
	return c
}

// Reduce ensures all coefficients are in [0, Q).
func (r *BigRing) Reduce(a BigPoly) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		bigmod.ReduceInPlace(c[i], a[i], r.Q)
	}
	return c
}

// CenterReduce returns a copy with coefficients in [-(Q-1)/2, (Q-1)/2].
func (r *BigRing) CenterReduce(a BigPoly) BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = bigmod.CenterReduce(a[i], r.Q)
	}
	return c
}

// InfNorm returns ||a||_∞ (max absolute centered coefficient).
func (r *BigRing) InfNorm(a BigPoly) *big.Int {
	maxVal := new(big.Int)
	for i := 0; i < r.N; i++ {
		v := bigmod.CenterReduce(a[i], r.Q)
		av := bigmod.Abs(v)
		if av.Cmp(maxVal) > 0 {
			maxVal.Set(av)
		}
	}
	return maxVal
}

// Equal returns true if a and b are the same polynomial in R_q.
func (r *BigRing) Equal(a, b BigPoly) bool {
	for i := 0; i < r.N; i++ {
		if bigmod.Reduce(a[i], r.Q).Cmp(bigmod.Reduce(b[i], r.Q)) != 0 {
			return false
		}
	}
	return true
}

// Mul computes a*b in R_q. Uses NTT if available, otherwise Karatsuba.
func (r *BigRing) Mul(a, b BigPoly) BigPoly {
	if r.hasNTT {
		aNTT := r.NTT(a)
		bNTT := r.NTT(b)
		cNTT := r.PointwiseMul(aNTT, bNTT)
		return r.InvNTT(cNTT)
	}
	return r.karatsubaMul(a, b)
}

// NaiveMul computes a*b using the schoolbook algorithm. For testing.
func (r *BigRing) NaiveMul(a, b BigPoly) BigPoly {
	c := r.NewPoly()
	tmp := new(big.Int)
	for i := 0; i < r.N; i++ {
		for j := 0; j < r.N; j++ {
			tmp.Mul(a[i], b[j])
			tmp.Mod(tmp, r.Q)
			idx := i + j
			if idx < r.N {
				c[idx].Add(c[idx], tmp)
				c[idx].Mod(c[idx], r.Q)
			} else {
				c[idx-r.N].Sub(c[idx-r.N], tmp)
				c[idx-r.N].Mod(c[idx-r.N], r.Q)
			}
		}
	}
	return c
}

// karatsubaMul computes a*b via Karatsuba then reduces mod X^N+1.
func (r *BigRing) karatsubaMul(a, b BigPoly) BigPoly {
	product := r.karatsubaRaw(a, b)
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i].Sub(product[i], product[i+r.N])
		c[i].Mod(c[i], r.Q)
	}
	return c
}

func (r *BigRing) karatsubaRaw(a, b BigPoly) []*big.Int {
	n := len(a)
	result := make([]*big.Int, 2*n)
	for i := range result {
		result[i] = new(big.Int)
	}

	if n <= 32 {
		tmp := new(big.Int)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				tmp.Mul(a[i], b[j])
				tmp.Mod(tmp, r.Q)
				result[i+j].Add(result[i+j], tmp)
				result[i+j].Mod(result[i+j], r.Q)
			}
		}
		return result
	}

	half := n / 2
	a0, a1 := BigPoly(a[:half]), BigPoly(a[half:])
	b0, b1 := BigPoly(b[:half]), BigPoly(b[half:])

	z0 := r.karatsubaRaw(a0, b0)
	z2 := r.karatsubaRaw(a1, b1)

	aSum := make(BigPoly, half)
	bSum := make(BigPoly, half)
	for i := 0; i < half; i++ {
		aSum[i] = bigmod.ModAdd(a0[i], a1[i], r.Q)
		bSum[i] = bigmod.ModAdd(b0[i], b1[i], r.Q)
	}
	z1 := r.karatsubaRaw(aSum, bSum)

	for i := range z0 {
		result[i].Add(result[i], z0[i])
		result[i].Mod(result[i], r.Q)
	}
	for i := range z2 {
		result[i+n].Add(result[i+n], z2[i])
		result[i+n].Mod(result[i+n], r.Q)
	}
	for i := range z1 {
		v := new(big.Int).Sub(z1[i], z0[i])
		v.Sub(v, z2[i])
		v.Mod(v, r.Q)
		result[i+half].Add(result[i+half], v)
		result[i+half].Mod(result[i+half], r.Q)
	}

	return result
}
