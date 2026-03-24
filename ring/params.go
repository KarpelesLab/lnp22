// Package ring implements polynomial arithmetic over the ring R_q = Z_q[X]/(X^N+1).
package ring

import "github.com/KarpelesLab/lnp22/internal/modmath"

const (
	// N is the degree of the cyclotomic polynomial X^N+1.
	N = 256

	// Q is the prime modulus. Q ≡ 1 (mod 2N) which enables NTT.
	// This is the same prime used in CRYSTALS-Dilithium.
	Q int64 = 8380417

	// RootOfUnity is a primitive 512th root of unity mod Q.
	// ζ^256 ≡ -1 (mod Q) and ζ^512 ≡ 1 (mod Q).
	// This is the same root used in CRYSTALS-Dilithium.
	RootOfUnity int64 = 1753

	// NInv is the modular inverse of N modulo Q: N * NInv ≡ 1 (mod Q).
	NInv int64 = 8347681
)

// zetas contains precomputed powers of the primitive 2N-th root of unity
// in bit-reversed order, used for the NTT butterfly operations.
// zetas[i] = ζ^{bitrev8(i)} mod Q where ζ is a primitive 512th root of unity.
var zetas [N]int64

func init() {
	for i := 0; i < N; i++ {
		zetas[i] = modmath.ModExp(RootOfUnity, int64(bitrev8(uint8(i))), Q)
	}
}

// bitrev8 returns the 8-bit reversal of x.
func bitrev8(x uint8) uint8 {
	x = (x&0xAA)>>1 | (x&0x55)<<1
	x = (x&0xCC)>>2 | (x&0x33)<<2
	x = (x&0xF0)>>4 | (x&0x0F)<<4
	return x
}

// Zetas returns a copy of the precomputed zeta table (for testing).
func Zetas() [N]int64 {
	return zetas
}
