package sampler

import (
	"encoding/binary"
	"io"

	"github.com/KarpelesLab/lnp22/ring"
)

// SampleUniformPoly samples a polynomial with coefficients uniform in [0, Q).
func SampleUniformPoly(r *ring.Ring, rng io.Reader) ring.Poly {
	p := r.NewPoly()
	var buf [8]byte

	// Find the smallest bitmask >= Q
	mask := int64(1)
	for mask < r.Q {
		mask <<= 1
	}
	mask--

	for i := 0; i < r.N; i++ {
		for {
			if _, err := io.ReadFull(rng, buf[:]); err != nil {
				panic("sampler: failed to read randomness: " + err.Error())
			}
			v := int64(binary.LittleEndian.Uint64(buf[:])) & mask
			if v >= 0 && v < r.Q {
				p[i] = v
				break
			}
		}
	}
	return p
}

// SampleUniformVec samples a vector of l uniform polynomials.
func SampleUniformVec(r *ring.Ring, l int, rng io.Reader) ring.PolyVec {
	v := r.NewPolyVec(l)
	for i := 0; i < l; i++ {
		v[i] = SampleUniformPoly(r, rng)
	}
	return v
}

// SampleUniformMat samples a k×l matrix of uniform polynomials.
func SampleUniformMat(r *ring.Ring, k, l int, rng io.Reader) ring.PolyMat {
	m := r.NewPolyMat(k, l)
	for i := 0; i < k; i++ {
		for j := 0; j < l; j++ {
			m[i][j] = SampleUniformPoly(r, rng)
		}
	}
	return m
}

// SampleTernaryPoly samples a polynomial with coefficients uniform in {-1, 0, 1}.
func SampleTernaryPoly(r *ring.Ring, rng io.Reader) ring.Poly {
	p := r.NewPoly()
	var buf [1]byte
	for i := 0; i < r.N; i++ {
		for {
			if _, err := io.ReadFull(rng, buf[:]); err != nil {
				panic("sampler: failed to read randomness: " + err.Error())
			}
			if buf[0] < 252 {
				switch buf[0] % 3 {
				case 0:
					p[i] = r.Q - 1
				case 1:
					p[i] = 0
				case 2:
					p[i] = 1
				}
				break
			}
		}
	}
	return p
}

// SampleTernaryVec samples a vector of l ternary polynomials.
func SampleTernaryVec(r *ring.Ring, l int, rng io.Reader) ring.PolyVec {
	v := r.NewPolyVec(l)
	for i := 0; i < l; i++ {
		v[i] = SampleTernaryPoly(r, rng)
	}
	return v
}
