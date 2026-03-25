package sampler

import (
	"io"
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
)

// SampleBigUniformPoly samples a polynomial with coefficients uniform in [0, Q).
func SampleBigUniformPoly(r *ring.BigRing, rng io.Reader) ring.BigPoly {
	p := r.NewPoly()
	byteLen := (r.Q.BitLen() + 7) / 8
	buf := make([]byte, byteLen)
	mask := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), uint(r.Q.BitLen())), big.NewInt(1))

	for i := 0; i < r.N; i++ {
		for {
			if _, err := io.ReadFull(rng, buf); err != nil {
				panic("sampler: failed to read randomness: " + err.Error())
			}
			v := new(big.Int).SetBytes(buf)
			v.And(v, mask)
			if v.Cmp(r.Q) < 0 {
				p[i] = v
				break
			}
		}
	}
	return p
}

// SampleBigUniformVec samples a vector of l uniform BigPolynomials.
func SampleBigUniformVec(r *ring.BigRing, l int, rng io.Reader) ring.BigPolyVec {
	v := r.NewPolyVec(l)
	for i := 0; i < l; i++ {
		v[i] = SampleBigUniformPoly(r, rng)
	}
	return v
}

// SampleBigUniformMat samples a k×l matrix of uniform BigPolynomials.
func SampleBigUniformMat(r *ring.BigRing, k, l int, rng io.Reader) ring.BigPolyMat {
	m := r.NewPolyMat(k, l)
	for i := 0; i < k; i++ {
		for j := 0; j < l; j++ {
			m[i][j] = SampleBigUniformPoly(r, rng)
		}
	}
	return m
}

// SampleBigTernaryPoly samples a polynomial with coefficients in {-1, 0, 1}.
func SampleBigTernaryPoly(r *ring.BigRing, rng io.Reader) ring.BigPoly {
	p := r.NewPoly()
	var buf [1]byte
	qm1 := new(big.Int).Sub(r.Q, big.NewInt(1))
	for i := 0; i < r.N; i++ {
		for {
			if _, err := io.ReadFull(rng, buf[:]); err != nil {
				panic("sampler: failed to read randomness: " + err.Error())
			}
			if buf[0] < 252 {
				switch buf[0] % 3 {
				case 0:
					p[i].Set(qm1) // -1 mod Q
				case 1:
					// 0 (already set)
				case 2:
					p[i].SetInt64(1)
				}
				break
			}
		}
	}
	return p
}

// SampleBigTernaryVec samples a vector of l ternary BigPolynomials.
func SampleBigTernaryVec(r *ring.BigRing, l int, rng io.Reader) ring.BigPolyVec {
	v := r.NewPolyVec(l)
	for i := 0; i < l; i++ {
		v[i] = SampleBigTernaryPoly(r, rng)
	}
	return v
}
