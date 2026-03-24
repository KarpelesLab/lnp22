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
	Ring   *ring.Ring // underlying polynomial ring
	K      int        // number of constraint equations (rows in A)
	L      int        // number of witness polynomials (columns in A)
	Kappa  int        // challenge weight: number of non-zero coefficients in {-1,+1}
	Beta   int64      // witness infinity norm bound: ||s||_∞ ≤ β
	Sigma  float64    // Gaussian std dev for masking vectors
	BoundZ int64      // response infinity norm bound: ||z||_∞ ≤ B_z

	// MaxAttempts limits the rejection sampling loop. 0 means no limit.
	MaxAttempts int
}

// DefaultParams returns parameters suitable for the Dilithium ring (N=256, Q=8380417).
func DefaultParams() *Params {
	return &Params{
		Ring:        ring.Dilithium(),
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
	if p.Ring == nil {
		return errors.New("nizk: Ring must not be nil")
	}
	if p.K <= 0 {
		return errors.New("nizk: K must be positive")
	}
	if p.L <= 0 {
		return errors.New("nizk: L must be positive")
	}
	if p.Kappa <= 0 || p.Kappa > p.Ring.N {
		return fmt.Errorf("nizk: Kappa must be in [1, %d]", p.Ring.N)
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
	if p.Sigma < 2*float64(p.Kappa)*float64(p.Beta) {
		return fmt.Errorf("nizk: Sigma=%f too small relative to Kappa*Beta=%d",
			p.Sigma, int64(p.Kappa)*p.Beta)
	}
	return nil
}
