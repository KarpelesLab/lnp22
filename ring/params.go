// Package ring implements polynomial arithmetic over the ring R_q = Z_q[X]/(X^N+1)
// with configurable degree N and modulus Q.
package ring

import (
	"fmt"

	"github.com/KarpelesLab/lnp22/internal/modmath"
)

// Ring holds the parameters for polynomial arithmetic in R_q = Z_q[X]/(X^N+1).
type Ring struct {
	N    int   // polynomial degree (must be a power of 2)
	Q    int64 // prime modulus
	NInv int64 // N^{-1} mod Q

	// NTT support — nil if Q does not admit a 2N-th root of unity.
	hasNTT bool
	zetas  []int64 // length N, bit-reversed powers of the primitive 2N-th root
}

// New creates a Ring with degree n and modulus q.
// n must be a power of 2 and q must be prime.
// NTT is enabled automatically if q ≡ 1 (mod 2n).
func New(n int, q int64) (*Ring, error) {
	if n < 2 || n&(n-1) != 0 {
		return nil, fmt.Errorf("ring: N=%d must be a power of 2", n)
	}
	if q < 3 {
		return nil, fmt.Errorf("ring: Q=%d must be >= 3", q)
	}

	r := &Ring{
		N:    n,
		Q:    q,
		NInv: modmath.ModInverse(int64(n), q),
	}

	// Try to enable NTT: need q ≡ 1 (mod 2n), i.e., (q-1) divisible by 2n.
	if (q-1)%int64(2*n) == 0 {
		root, ok := findPrimitive2NthRoot(n, q)
		if ok {
			r.hasNTT = true
			r.zetas = make([]int64, n)
			for i := 0; i < n; i++ {
				r.zetas[i] = modmath.ModExp(root, int64(bitrev(uint32(i), log2(n))), q)
			}
		}
	}

	return r, nil
}

// HasNTT reports whether NTT-based multiplication is available for this ring.
func (r *Ring) HasNTT() bool {
	return r.hasNTT
}

// Dilithium returns the Ring used by CRYSTALS-Dilithium: N=256, Q=8380417.
func Dilithium() *Ring {
	r, err := New(256, 8380417)
	if err != nil {
		panic("ring: failed to create Dilithium ring: " + err.Error())
	}
	return r
}

// findPrimitive2NthRoot searches for a primitive 2n-th root of unity mod q.
// Returns (root, true) on success.
func findPrimitive2NthRoot(n int, q int64) (int64, bool) {
	// We need g such that g^(2n) ≡ 1 and g^n ≡ -1 (mod q).
	// Strategy: find a generator of Z_q* and compute g = generator^((q-1)/(2n)).
	exp := (q - 1) / int64(2*n)

	// Try small candidates as generators
	for candidate := int64(2); candidate < q && candidate < 10000; candidate++ {
		root := modmath.ModExp(candidate, exp, q)
		// Check it's a primitive 2n-th root: root^n ≡ -1 (mod q)
		if modmath.ModExp(root, int64(n), q) == q-1 {
			return root, true
		}
	}
	return 0, false
}

// log2 returns the base-2 log of n (n must be a power of 2).
func log2(n int) int {
	r := 0
	for n > 1 {
		n >>= 1
		r++
	}
	return r
}

// bitrev reverses the bottom `bits` bits of x.
func bitrev(x uint32, bits int) uint32 {
	var r uint32
	for i := 0; i < bits; i++ {
		r = (r << 1) | (x & 1)
		x >>= 1
	}
	return r
}
