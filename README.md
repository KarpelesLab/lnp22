# LNP22 — Lattice-Based NIZK Framework in Go

[![Tests](https://github.com/KarpelesLab/lnp22/actions/workflows/test.yml/badge.svg)](https://github.com/KarpelesLab/lnp22/actions/workflows/test.yml)
[![Coverage Status](https://coveralls.io/repos/github/KarpelesLab/lnp22/badge.svg?branch=master)](https://coveralls.io/github/KarpelesLab/lnp22?branch=master)

A self-contained Go implementation of the LNP22 non-interactive zero-knowledge proof system from [Lyubashevsky, Nguyen, and Plançon (CRYPTO 2022)](https://eprint.iacr.org/2022/284).

The framework proves knowledge of short vectors satisfying linear relations over polynomial rings, with applications to range proofs and proof composition. No external cryptography dependencies beyond `golang.org/x/crypto/sha3` for SHAKE256.

## Features

- **Configurable ring** R_q = Z_q[X]/(X^N+1) — any power-of-2 degree N and prime Q
- **NTT or Karatsuba** — automatic NTT when Q ≡ 1 (mod 2N), Karatsuba fallback otherwise
- **BDLOP commitment scheme** — additively homomorphic, with scalar multiplication support
- **Linear-relation proofs** — prove knowledge of short s such that A·s ≡ t (mod q)
- **Range proofs** — prove witness coefficients lie in [-β, β] via binary decomposition
- **Proof composition** — bind multiple sub-proofs under a single Fiat-Shamir transcript
- **Discrete Gaussian sampling** via cumulative distribution tables
- **152-bit modular arithmetic** — `big.Int`-based `BigRing`/`BigPoly` types for large moduli
- **NTRU trapdoor generation** — TrapGen from DLP14 with gadget matrix toolkit
- **Gaussian pre-image sampling** — SamplePre via the GPV framework
- **Public key encryption** — NTRU-style PKE with rounding-based encoding/decoding
- **Two-ring NIZK** — linear + quadratic relation proofs with ℓ₂-norm bounds
- **Rounding function** ⌈·⌉_p — encode/decode between Z_p and Z_q
- **Configurable parameters** — security level, module dimensions, norm bounds

## Packages

```
ring/           Polynomial ring R_q, NTT, Karatsuba — both int64 and big.Int types
sampler/        Discrete Gaussian, uniform, ternary, challenge — int64 and big.Int
commitment/     BDLOP lattice-based commitment scheme
nizk/           NIZK proof system (linear, range, composed) — int64 ring
tworing/        Two-ring NIZK (linear + quadratic + ℓ₂-norm) — big.Int ring
trapgen/        NTRU trapdoor generation (TrapGen from DLP14)
preimage/       Gaussian pre-image sampling (SamplePre, GPV framework)
pke/            NTRU-style public key encryption with rounding
internal/       Modular arithmetic (int64 + big.Int) and rounding
```

## Quick Start

```go
package main

import (
    "crypto/rand"
    "fmt"

    "github.com/KarpelesLab/lnp22/nizk"
    "github.com/KarpelesLab/lnp22/sampler"
)

func main() {
    params := nizk.DefaultParams() // Dilithium ring (N=256, Q=8380417)

    // Generate a random statement: public matrix A, ternary secret s, target t = A·s
    r := params.Ring
    A := sampler.SampleUniformMat(r, params.K, params.L, rand.Reader)
    s := sampler.SampleTernaryVec(r, params.L, rand.Reader)
    t := r.MatVecMul(A, s)

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

Parameters are fully configurable. To use a different ring (e.g., Q=7933, N=512):

```go
r, _ := ring.New(512, 7933) // Karatsuba multiplication (no NTT for this Q)
params := &nizk.Params{
    Ring: r, K: 3, L: 4, Kappa: 40, Beta: 1,
    Sigma: 200, BoundZ: 800, MaxAttempts: 1000,
}
```

The `Validate()` method checks consistency (e.g., σ ≥ 2κβ for rejection sampling convergence). NTT is enabled automatically when Q ≡ 1 (mod 2N); otherwise Karatsuba multiplication is used.

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

148 tests across 10 packages covering:
- Ring arithmetic — NTT roundtrip, Karatsuba correctness, naive equivalence, ring axioms, X^N ≡ -1
- Multiple ring configurations — Dilithium (N=256, Q=8380417) and Q=7933/N=512
- Sampling distributions — Gaussian mean/variance/tail bounds, uniform range, ternary balance, challenge determinism
- Commitment scheme — open/verify roundtrip, homomorphic properties, tampering detection
- NIZK proofs — prove/verify roundtrips, rejection of tampered proofs/statements/seeds, norm bound enforcement, cross-parameter-set and cross-ring validation
- 152-bit big.Int ring — NTT roundtrip, Karatsuba, ring axioms at full 152-bit modulus
- TrapGen — NTRU key generation, gadget decomposition roundtrip, polynomial inverse
- Pre-image sampling — dimension checks, shortness bounds, multiple targets
- PKE — encode/decode roundtrip, encrypt/decrypt, out-of-range rejection
- Two-ring NIZK — linear, quadratic, ℓ₂-norm proofs, tamper rejection

## References

- V. Lyubashevsky, N. K. Nguyen, M. Plançon. *Lattice-Based Zero-Knowledge Proofs and Applications: Shorter, Simpler, and More General.* CRYPTO 2022. [ePrint 2022/284](https://eprint.iacr.org/2022/284)
- CRYSTALS-Dilithium (FIPS 204) — ring parameters (N=256, Q=8380417, ζ=1753)
- C. Baum, I. Damgård, V. Lyubashevsky, S. Oechsner, C. Peikert. *More Efficient Commitments from Structured Lattice Assumptions.* SCN 2018. — BDLOP commitment scheme

## License

See [LICENSE](LICENSE) for details.
