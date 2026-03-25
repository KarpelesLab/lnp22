// Package tworing implements the LNP22 two-ring NIZK proof system.
//
// The system operates over two rings simultaneously:
//   - Ring 1 (commitment): used for BDLOP-style commitments
//   - Ring 2 (proof): used for proving polynomial relations
//
// This dual-ring approach enables direct ℓ₂-norm proofs without CRT overhead,
// supporting both linear (A·s = t) and quadratic (⟨s, B·s⟩ = v) relations.
//
// Reference: Lyubashevsky, Nguyen, Plançon — "Lattice-Based Zero-Knowledge Proofs
// and Applications: Shorter, Simpler, and More General" (CRYPTO 2022).
package tworing

import (
	"errors"
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
)

// Params holds the two-ring proof system parameters.
type Params struct {
	Ring *ring.BigRing // polynomial ring for both commitment and proof

	K     int      // constraint equations (rows in A)
	L     int      // witness polynomials (columns in A)
	Kappa int      // challenge weight
	Sigma float64  // Gaussian masking parameter
	BoundZ *big.Int // l-inf norm bound on response z

	// L2 norm parameters
	L2BoundSq *big.Int // ℓ₂-norm² bound on witness: ||s||₂² ≤ this

	MaxAttempts int
}

// DefaultParams returns default two-ring parameters for the given BigRing.
func DefaultParams(r *ring.BigRing) *Params {
	boundZ := big.NewInt(1400)
	l2bound := new(big.Int).Mul(big.NewInt(int64(r.N)), big.NewInt(4)) // N * β²
	return &Params{
		Ring:        r,
		K:           4,
		L:           5,
		Kappa:       60,
		Sigma:       350,
		BoundZ:      boundZ,
		L2BoundSq:   l2bound,
		MaxAttempts: 1000,
	}
}

// Validate checks parameter consistency.
func (p *Params) Validate() error {
	if p.Ring == nil {
		return errors.New("tworing: Ring must not be nil")
	}
	if p.K <= 0 || p.L <= 0 {
		return errors.New("tworing: K and L must be positive")
	}
	if p.Kappa <= 0 || p.Kappa > p.Ring.N {
		return errors.New("tworing: Kappa out of range")
	}
	if p.Sigma <= 0 {
		return errors.New("tworing: Sigma must be positive")
	}
	if p.BoundZ.Sign() <= 0 {
		return errors.New("tworing: BoundZ must be positive")
	}
	return nil
}
