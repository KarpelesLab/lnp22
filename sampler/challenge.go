package sampler

import (
	"encoding/binary"

	"github.com/KarpelesLab/lnp22/ring"
	"golang.org/x/crypto/sha3"
)

// SampleChallenge produces a challenge polynomial with exactly kappa non-zero
// coefficients, each in {-1, +1}. The polynomial is deterministically derived
// from the given seed using SHAKE256.
func SampleChallenge(seed []byte, kappa int) *ring.Poly {
	if kappa < 0 || kappa > ring.N {
		panic("sampler: invalid kappa")
	}

	h := sha3.NewShake256()
	h.Write(seed)

	var p ring.Poly

	// Fisher-Yates shuffle to select kappa positions from [0, N).
	// We place non-zero values in the last kappa positions of the permutation.
	signs := readBytes(h, (kappa+7)/8) // sign bits
	signIdx := 0
	signBit := 0

	for i := ring.N - kappa; i < ring.N; i++ {
		// Sample j uniformly from [0, i]
		j := sampleIndex(h, i+1)
		p[i] = p[j]

		// Assign ±1 based on sign bit
		bit := (signs[signIdx] >> signBit) & 1
		signBit++
		if signBit == 8 {
			signBit = 0
			signIdx++
		}
		if bit == 0 {
			p[j] = 1
		} else {
			p[j] = ring.Q - 1 // -1 mod Q
		}
	}

	return &p
}

// sampleIndex returns a uniform random value in [0, bound) using rejection sampling.
func sampleIndex(h sha3.ShakeHash, bound int) int {
	// Find smallest power of 2 >= bound
	mask := 1
	for mask < bound {
		mask <<= 1
	}
	mask--

	var buf [4]byte
	for {
		h.Read(buf[:])
		v := int(binary.LittleEndian.Uint32(buf[:])) & mask
		if v < bound {
			return v
		}
	}
}

// readBytes reads n bytes from the SHAKE hash.
func readBytes(h sha3.ShakeHash, n int) []byte {
	buf := make([]byte, n)
	h.Read(buf)
	return buf
}
