package sampler

import (
	"encoding/binary"
	"io"

	"github.com/KarpelesLab/lnp22/ring"
)

// SampleUniformPoly samples a polynomial with coefficients uniform in [0, Q).
func SampleUniformPoly(rng io.Reader) *ring.Poly {
	var p ring.Poly
	var buf [4]byte
	for i := 0; i < ring.N; i++ {
		// Rejection sampling to get uniform in [0, Q)
		for {
			if _, err := io.ReadFull(rng, buf[:]); err != nil {
				panic("sampler: failed to read randomness: " + err.Error())
			}
			// Use 23 bits (since Q < 2^23)
			v := int64(binary.LittleEndian.Uint32(buf[:]) & 0x7FFFFF)
			if v < ring.Q {
				p[i] = v
				break
			}
		}
	}
	return &p
}

// SampleUniformVec samples a vector of l uniform polynomials.
func SampleUniformVec(l int, rng io.Reader) ring.PolyVec {
	v := ring.NewPolyVec(l)
	for i := 0; i < l; i++ {
		v[i] = *SampleUniformPoly(rng)
	}
	return v
}

// SampleUniformMat samples a k×l matrix of uniform polynomials.
func SampleUniformMat(k, l int, rng io.Reader) ring.PolyMat {
	m := ring.NewPolyMat(k, l)
	for i := 0; i < k; i++ {
		for j := 0; j < l; j++ {
			m[i][j] = *SampleUniformPoly(rng)
		}
	}
	return m
}

// SampleTernaryPoly samples a polynomial with coefficients uniform in {-1, 0, 1}.
func SampleTernaryPoly(rng io.Reader) *ring.Poly {
	var p ring.Poly
	var buf [1]byte
	for i := 0; i < ring.N; i++ {
		// Rejection sample to get uniform in {0, 1, 2}
		for {
			if _, err := io.ReadFull(rng, buf[:]); err != nil {
				panic("sampler: failed to read randomness: " + err.Error())
			}
			v := buf[0]
			if v < 252 { // 252 = 84*3, largest multiple of 3 < 256
				switch v % 3 {
				case 0:
					p[i] = ring.Q - 1 // -1 mod Q
				case 1:
					p[i] = 0
				case 2:
					p[i] = 1
				}
				break
			}
		}
	}
	return &p
}

// SampleTernaryVec samples a vector of l ternary polynomials.
func SampleTernaryVec(l int, rng io.Reader) ring.PolyVec {
	v := ring.NewPolyVec(l)
	for i := 0; i < l; i++ {
		v[i] = *SampleTernaryPoly(rng)
	}
	return v
}
