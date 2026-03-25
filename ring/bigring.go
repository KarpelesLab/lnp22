package ring

import (
	"fmt"
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
)

// BigRing holds parameters for polynomial arithmetic in R_q = Z_q[X]/(X^N+1)
// with arbitrary-precision modulus Q (supporting q ≈ 2^152 and beyond).
type BigRing struct {
	N    int
	Q    *big.Int
	NInv *big.Int // N^{-1} mod Q

	hasNTT bool
	zetas  []*big.Int // length N, bit-reversed powers of primitive 2N-th root
}

// NewBig creates a BigRing with degree n and modulus q.
// n must be a power of 2. NTT is enabled if q ≡ 1 (mod 2n).
func NewBig(n int, q *big.Int) (*BigRing, error) {
	if n < 2 || n&(n-1) != 0 {
		return nil, fmt.Errorf("ring: N=%d must be a power of 2", n)
	}
	if q.Cmp(big.NewInt(3)) < 0 {
		return nil, fmt.Errorf("ring: Q must be >= 3")
	}

	r := &BigRing{
		N:    n,
		Q:    new(big.Int).Set(q),
		NInv: bigmod.ModInverse(big.NewInt(int64(n)), q),
	}

	// Try to enable NTT: need q ≡ 1 (mod 2n)
	twoN := big.NewInt(int64(2 * n))
	rem := new(big.Int).Mod(new(big.Int).Sub(q, big.NewInt(1)), twoN)
	if rem.Sign() == 0 {
		root, ok := findBigPrimitive2NthRoot(n, q)
		if ok {
			r.hasNTT = true
			r.zetas = make([]*big.Int, n)
			bits := log2(n)
			for i := 0; i < n; i++ {
				r.zetas[i] = bigmod.ModExp(root, big.NewInt(int64(bitrev(uint32(i), bits))), q)
			}
		}
	}

	return r, nil
}

// HasNTT reports whether NTT-based multiplication is available.
func (r *BigRing) HasNTT() bool {
	return r.hasNTT
}

// LNP22Big returns the BigRing for the full LNP22 proof system:
// N=1024, Q = 152-bit prime with Q ≡ 1 (mod 2048).
func LNP22Big() *BigRing {
	q, ok := new(big.Int).SetString("5708990770823839524233143877797980545530982401", 10)
	if !ok {
		panic("ring: failed to parse LNP22 big Q")
	}
	r, err := NewBig(1024, q)
	if err != nil {
		panic("ring: failed to create LNP22 big ring: " + err.Error())
	}
	return r
}

// findBigPrimitive2NthRoot searches for a primitive 2n-th root of unity mod q.
func findBigPrimitive2NthRoot(n int, q *big.Int) (*big.Int, bool) {
	// g^((q-1)/(2n)) is a 2n-th root for any generator g.
	// We need it to be primitive: root^n ≡ -1 (mod q).
	qm1 := new(big.Int).Sub(q, big.NewInt(1))
	exp := new(big.Int).Div(qm1, big.NewInt(int64(2*n)))
	nBig := big.NewInt(int64(n))
	qMinus1 := new(big.Int).Sub(q, big.NewInt(1))

	for candidate := int64(2); candidate < 10000; candidate++ {
		root := bigmod.ModExp(big.NewInt(candidate), exp, q)
		// Check root^n ≡ -1 (mod q)
		if bigmod.ModExp(root, nBig, q).Cmp(qMinus1) == 0 {
			return root, true
		}
	}
	return nil, false
}
