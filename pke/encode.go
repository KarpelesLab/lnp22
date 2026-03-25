// Package pke implements NTRU-style lattice-based public key encryption
// with rounding-based encoding for ciphertext compression.
package pke

import (
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
)

// Encode lifts a plaintext polynomial from Z_p to Z_q.
// Each coefficient m_i ∈ [0, p) is mapped to round(q * m_i / p) ∈ Z_q.
func Encode(r *ring.BigRing, m ring.BigPoly, p *big.Int) ring.BigPoly {
	c := r.NewPoly()
	for i := 0; i < r.N; i++ {
		c[i] = bigmod.ScaleUp(m[i], p, r.Q)
	}
	return c
}

// Decode maps a Z_q polynomial to Z_p via rounding.
// Each coefficient x_i ∈ [0, q) is mapped to round(p * x_i / q) ∈ [0, p).
func Decode(r *ring.BigRing, c ring.BigPoly, p *big.Int) ring.BigPoly {
	m := r.NewPoly()
	for i := 0; i < r.N; i++ {
		m[i] = bigmod.RoundP(bigmod.Reduce(c[i], r.Q), p, r.Q)
	}
	return m
}
