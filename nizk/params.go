// Package nizk implements the LNP22 lattice-based non-interactive zero-knowledge proof system.
//
// The system proves knowledge of a short vector s satisfying A·s ≡ t (mod q)
// in the ring R_q = Z_q[X]/(X^N+1), using the Fiat-Shamir transform applied
// to a sigma protocol with rejection sampling.
//
// Reference: Lyubashevsky, Nguyen, Plançon — "Lattice-Based Zero-Knowledge Proofs
// and Applications: Shorter, Simpler, and More General" (CRYPTO 2022).
package nizk

import (
	"errors"
	"fmt"

	"github.com/KarpelesLab/lnp22/ring"
)

// Params holds the proof system parameters.
type Params struct {
	K      int     // number of constraint equations (rows in A)
	L      int     // number of witness polynomials (columns in A)
	Kappa  int     // challenge weight: number of non-zero coefficients in {-1,+1}
	Beta   int64   // witness infinity norm bound: ||s||_∞ ≤ β
	Sigma  float64 // Gaussian std dev for masking vectors
	BoundZ int64   // response infinity norm bound: ||z||_∞ ≤ B_z

	// MaxAttempts limits the rejection sampling loop. 0 means no limit.
	MaxAttempts int
}

// DefaultParams returns parameters suitable for proving knowledge of a ternary
// secret (β=1) in a 4×5 module with 128-bit security.
func DefaultParams() *Params {
	return &Params{
		K:           4,
		L:           5,
		Kappa:       60,
		Beta:        1,
		Sigma:       350,
		BoundZ:      1400,
		MaxAttempts: 1000,
	}
}

// Validate checks that the parameters are self-consistent.
func (p *Params) Validate() error {
	if p.K <= 0 {
		return errors.New("nizk: K must be positive")
	}
	if p.L <= 0 {
		return errors.New("nizk: L must be positive")
	}
	if p.Kappa <= 0 || p.Kappa > ring.N {
		return fmt.Errorf("nizk: Kappa must be in [1, %d]", ring.N)
	}
	if p.Beta <= 0 {
		return errors.New("nizk: Beta must be positive")
	}
	if p.Sigma <= 0 {
		return errors.New("nizk: Sigma must be positive")
	}
	if p.BoundZ <= 0 {
		return errors.New("nizk: BoundZ must be positive")
	}
	// Sigma should be large enough relative to Kappa*Beta for rejection sampling
	// to succeed with reasonable probability.
	if p.Sigma < 2*float64(p.Kappa)*float64(p.Beta) {
		return fmt.Errorf("nizk: Sigma=%f too small relative to Kappa*Beta=%d (need Sigma ≥ 2*κ*β)",
			p.Sigma, int64(p.Kappa)*p.Beta)
	}
	return nil
}
