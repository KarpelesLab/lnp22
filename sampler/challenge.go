package sampler

import (
	"encoding/binary"

	"github.com/KarpelesLab/lnp22/ring"
	"golang.org/x/crypto/sha3"
)

// SampleChallenge produces a challenge polynomial with exactly kappa non-zero
// coefficients, each in {-1, +1}. The polynomial is deterministically derived
// from the given seed using SHAKE256.
func SampleChallenge(r *ring.Ring, seed []byte, kappa int) ring.Poly {
	if kappa < 0 || kappa > r.N {
		panic("sampler: invalid kappa")
	}

	h := sha3.NewShake256()
	h.Write(seed)

	p := r.NewPoly()

	signs := readBytes(h, (kappa+7)/8)
	signIdx := 0
	signBit := 0

	for i := r.N - kappa; i < r.N; i++ {
		j := sampleIndex(h, i+1)
		p[i] = p[j]

		bit := (signs[signIdx] >> signBit) & 1
		signBit++
		if signBit == 8 {
			signBit = 0
			signIdx++
		}
		if bit == 0 {
			p[j] = 1
		} else {
			p[j] = r.Q - 1
		}
	}

	return p
}

func sampleIndex(h sha3.ShakeHash, bound int) int {
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

func readBytes(h sha3.ShakeHash, n int) []byte {
	buf := make([]byte, n)
	h.Read(buf)
	return buf
}
