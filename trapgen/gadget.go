// Package trapgen implements NTRU trapdoor generation (TrapGen from DLP14)
// and the gadget matrix toolkit for lattice-based cryptography.
package trapgen

import (
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
)

// GadgetVec returns g = (1, b, b², ..., b^{k-1}) mod q
// where k = ceil(log_b(q)).
func GadgetVec(q *big.Int, base int) []*big.Int {
	b := big.NewInt(int64(base))
	k := GadgetLen(q, base)
	g := make([]*big.Int, k)
	power := big.NewInt(1)
	for i := 0; i < k; i++ {
		g[i] = new(big.Int).Set(power)
		power = bigmod.ModMul(power, b, q)
	}
	return g
}

// GadgetLen returns k = ceil(log_b(q)).
func GadgetLen(q *big.Int, base int) int {
	k := 0
	v := new(big.Int).Set(q)
	for v.Sign() > 0 {
		k++
		v.Div(v, big.NewInt(int64(base)))
	}
	return k
}

// Decompose returns the base-b signed digit decomposition of x.
// Returns v of length k such that ⟨g, v⟩ = x (mod q) and ||v||_∞ ≤ ⌈b/2⌉.
func Decompose(x, q *big.Int, base int) []*big.Int {
	k := GadgetLen(q, base)
	b := big.NewInt(int64(base))
	halfB := new(big.Int).Rsh(b, 1)

	v := make([]*big.Int, k)
	rem := bigmod.Reduce(new(big.Int).Set(x), q)

	for i := 0; i < k; i++ {
		digit := new(big.Int).Mod(rem, b)
		// Signed balanced: if digit > b/2, subtract b (carry)
		if digit.Cmp(halfB) > 0 {
			digit.Sub(digit, b)
			rem.Add(rem, b) // carry
		}
		v[i] = digit
		rem.Sub(rem, digit)
		rem.Div(rem, b)
	}
	return v
}

// DecomposePoly decomposes each coefficient of a polynomial using base-b.
// Returns k polynomials v[0], ..., v[k-1] such that
// Σ_i b^i · v[i] = p (mod q) coefficient-wise.
func DecomposePoly(r *ring.BigRing, p ring.BigPoly, base int) ring.BigPolyVec {
	k := GadgetLen(r.Q, base)
	result := r.NewPolyVec(k)
	for j := 0; j < r.N; j++ {
		digits := Decompose(p[j], r.Q, base)
		for i := 0; i < k; i++ {
			result[i][j] = bigmod.Reduce(digits[i], r.Q)
		}
	}
	return result
}

// RecomposePoly reconstructs a polynomial from its base-b decomposition.
// result = Σ_i b^i · v[i] (mod q) coefficient-wise.
func RecomposePoly(r *ring.BigRing, v ring.BigPolyVec, base int) ring.BigPoly {
	g := GadgetVec(r.Q, base)
	result := r.NewPoly()
	for i := range v {
		scaled := r.ScalarMul(v[i], g[i])
		result = r.Add(result, scaled)
	}
	return result
}

// Trapdoor represents the lattice trapdoor structure.
// The public matrix is [A | A·R + G] where R is the short trapdoor matrix
// and G is the gadget row.
type Trapdoor struct {
	R      ring.BigPolyVec // short trapdoor polynomials (l entries)
	Base   int             // gadget base b
	Params *Params
}
