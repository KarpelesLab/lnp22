# LNP22 — Lattice-Based NIZK Framework in Go

[![Tests](https://github.com/KarpelesLab/lnp22/actions/workflows/test.yml/badge.svg)](https://github.com/KarpelesLab/lnp22/actions/workflows/test.yml)
[![Coverage Status](https://coveralls.io/repos/github/KarpelesLab/lnp22/badge.svg?branch=master)](https://coveralls.io/github/KarpelesLab/lnp22?branch=master)

A self-contained Go implementation of the LNP22 non-interactive zero-knowledge proof system from [Lyubashevsky, Nguyen, and Plançon (CRYPTO 2022)](https://eprint.iacr.org/2022/284).

The framework proves knowledge of short vectors satisfying linear relations over polynomial rings, with applications to range proofs and proof composition. No external cryptography dependencies beyond `golang.org/x/crypto/sha3` for SHAKE256.

## Features

- **Ring arithmetic** over R_q = Z_q[X]/(X^256+1) with NTT-accelerated polynomial multiplication
- **BDLOP commitment scheme** — additively homomorphic, with scalar multiplication support
- **Linear-relation proofs** — prove knowledge of short s such that A·s ≡ t (mod q)
- **Range proofs** — prove witness coefficients lie in [-β, β] via binary decomposition
- **Proof composition** — bind multiple sub-proofs under a single Fiat-Shamir transcript
- **Discrete Gaussian sampling** via cumulative distribution tables
- **Configurable parameters** — security level, module dimensions, norm bounds

## Packages

```
ring/           Polynomial ring R_q = Z_q[X]/(X^N+1), NTT, vectors, matrices
sampler/        Discrete Gaussian, uniform, ternary, and challenge sampling
commitment/     BDLOP lattice-based commitment scheme
nizk/           LNP22 NIZK proof system (linear, range, composed)
internal/       Modular arithmetic utilities
```

## Quick Start

```go
package main

import (
    "crypto/rand"
    "fmt"

    "github.com/KarpelesLab/lnp22/nizk"
    "github.com/KarpelesLab/lnp22/ring"
    "github.com/KarpelesLab/lnp22/sampler"
)

func main() {
    params := nizk.DefaultParams() // K=4, L=5, κ=60, β=1, σ=350

    // Generate a random statement: public matrix A, ternary secret s, target t = A·s
    A := sampler.SampleUniformMat(params.K, params.L, rand.Reader)
    s := sampler.SampleTernaryVec(params.L, rand.Reader)
    t := ring.MatVecMul(A, s)

    stmt := &nizk.Statement{A: A, T: t}
    wit := &nizk.Witness{S: s}

    // Prove
    proof, err := nizk.ProveLinear(params, stmt, wit, rand.Reader)
    if err != nil {
        panic(err)
    }

    // Verify
    ok := nizk.VerifyLinear(params, stmt, proof)
    fmt.Println("Proof valid:", ok) // true
}
```

## Cryptographic Parameters

Default parameters target 128-bit security:

| Parameter | Value | Description |
|-----------|-------|-------------|
| N | 256 | Ring degree (X^N+1 cyclotomic) |
| Q | 8380417 | Prime modulus (same as CRYSTALS-Dilithium) |
| K | 4 | Constraint equations (rows in A) |
| L | 5 | Witness polynomials (columns in A) |
| κ | 60 | Challenge weight (non-zero ternary coefficients) |
| β | 1 | Witness infinity norm bound |
| σ | 350 | Gaussian standard deviation for masking |
| B_z | 1400 | Response infinity norm bound |

Parameters are configurable via the `nizk.Params` struct. The `Validate()` method checks consistency (e.g., σ ≥ 2κβ for rejection sampling convergence).

## Protocol Overview

### Linear-Relation Proof

Proves knowledge of s with ||s||∞ ≤ β satisfying A·s ≡ t (mod q):

1. **Mask:** Sample y ← D_σ^l (Gaussian masking vector)
2. **Commit:** w = A·y (mod q)
3. **Challenge:** c = SHAKE256(A ‖ t ‖ w) → polynomial with κ non-zero ±1 coefficients
4. **Respond:** z = y + c·s
5. **Reject/accept:** If ||z||∞ > B_z, restart; otherwise output proof (w, z)

Verification checks: A·z ≡ c·t + w (mod q) and ||z||∞ ≤ B_z.

### Range Proof

Proves each coefficient of s lies in [-β, β]:

1. Prove the linear relation A·s = t
2. Shift and decompose: (s[k] + β) = Σ_i b_i[k] · 2^i
3. Prove each bit polynomial b_i has ||b_i||∞ ≤ 1
4. Prove reconstruction: Σ_i 2^i · b_i = s + β

### Proof Composition

Binds multiple independent proofs under a chained Fiat-Shamir transcript with a random seed, preventing proof reuse across contexts.

### BDLOP Commitment

Additively homomorphic commitment: Com(m; r) = (B₀·r, ⟨b_i, r⟩ + m_i). Supports:
- `AddCommitments(a, b)` — homomorphic addition
- `ScalarMulCommitment(c, s)` — multiplication by a ring element
- `Verify(key, commitment, opening)` — opening verification

## Performance

Benchmarked on Intel Core i9-14900K:

| Operation | Time |
|-----------|------|
| NTT (256-point) | 8 μs |
| Polynomial multiply | 28 μs |
| ProveLinear | 2.0 ms |
| VerifyLinear | 1.0 ms |
| ProveRange | 10.3 ms |
| ComposeProofs (×3) | 6.4 ms |
| BDLOP Commit | 2.4 ms |

## Testing

```sh
go test ./...
go test ./... -bench=. -benchtime=1s
```

47 tests covering:
- Ring arithmetic correctness (NTT roundtrip, naive multiplication equivalence, ring axioms, X^N ≡ -1)
- Sampling distributions (Gaussian mean/variance/tail bounds, uniform range, ternary balance, challenge determinism)
- Commitment scheme (open/verify roundtrip, homomorphic properties, tampering detection)
- NIZK proofs (prove/verify roundtrips, rejection of tampered proofs/statements/seeds, norm bound enforcement, cross-parameter-set validation)

## References

- V. Lyubashevsky, N. K. Nguyen, M. Plançon. *Lattice-Based Zero-Knowledge Proofs and Applications: Shorter, Simpler, and More General.* CRYPTO 2022. [ePrint 2022/284](https://eprint.iacr.org/2022/284)
- CRYSTALS-Dilithium (FIPS 204) — ring parameters (N=256, Q=8380417, ζ=1753)
- C. Baum, I. Damgård, V. Lyubashevsky, S. Oechsner, C. Peikert. *More Efficient Commitments from Structured Lattice Assumptions.* SCN 2018. — BDLOP commitment scheme

## License

See [LICENSE](LICENSE) for details.
