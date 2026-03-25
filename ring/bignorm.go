package ring

import (
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
)

// L2NormSq returns ||a||_2^2 = Σ (center-reduced coefficient)^2.
func (r *BigRing) L2NormSq(a BigPoly) *big.Int {
	sum := new(big.Int)
	tmp := new(big.Int)
	for i := 0; i < r.N; i++ {
		v := bigmod.CenterReduce(a[i], r.Q)
		tmp.Mul(v, v)
		sum.Add(sum, tmp)
	}
	return sum
}

// VecL2NormSq returns ||v||_2^2 = Σ_i ||v[i]||_2^2.
func (r *BigRing) VecL2NormSq(v BigPolyVec) *big.Int {
	sum := new(big.Int)
	for i := range v {
		sum.Add(sum, r.L2NormSq(v[i]))
	}
	return sum
}
