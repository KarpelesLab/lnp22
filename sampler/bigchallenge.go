package sampler

import (
	"encoding/binary"
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
	"golang.org/x/crypto/sha3"
)

// SampleBigChallenge produces a challenge polynomial for a BigRing with exactly
// kappa non-zero coefficients in {-1, +1}, deterministically from seed.
func SampleBigChallenge(r *ring.BigRing, seed []byte, kappa int) ring.BigPoly {
	if kappa < 0 || kappa > r.N {
		panic("sampler: invalid kappa")
	}

	h := sha3.NewShake256()
	h.Write(seed)

	p := r.NewPoly()
	qm1 := new(big.Int).Sub(r.Q, big.NewInt(1))

	signs := readBytes(h, (kappa+7)/8)
	signIdx := 0
	signBit := 0

	for i := r.N - kappa; i < r.N; i++ {
		j := sampleBigIndex(h, i+1)
		p[i].Set(p[j])

		bit := (signs[signIdx] >> signBit) & 1
		signBit++
		if signBit == 8 {
			signBit = 0
			signIdx++
		}
		if bit == 0 {
			p[j].SetInt64(1)
		} else {
			p[j].Set(qm1)
		}
	}

	return p
}

func sampleBigIndex(h sha3.ShakeHash, bound int) int {
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
