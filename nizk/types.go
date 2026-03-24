package nizk

import "github.com/KarpelesLab/lnp22/ring"

// Statement is the public part of a linear-relation proof: A·s ≡ t (mod q).
type Statement struct {
	A ring.PolyMat // k×l public matrix
	T ring.PolyVec // k-vector target
}

// Witness is the prover's secret: a short vector s with ||s||_∞ ≤ β.
type Witness struct {
	S ring.PolyVec // l-vector secret
}

// LinearProof is a NIZK proof of knowledge of s such that A·s ≡ t (mod q) and ||s||_∞ ≤ β.
type LinearProof struct {
	W ring.PolyVec // commitment w = A·y (k-vector)
	Z ring.PolyVec // response z = y + c·s (l-vector)
}

// RangeStatement extends a linear statement with a commitment to the witness
// and a bound on the coefficient range.
type RangeStatement struct {
	Statement            // linear relation A·s = t
	Beta      int64      // range bound: coefficients of s lie in [-β, β]
	NumBits   int        // number of bits for binary decomposition: ⌈log2(2β+1)⌉
}

// RangeProof proves that the witness coefficients lie in [-β, β].
// It uses binary decomposition: s_j + β = Σ_i b_{j,i} · 2^i
// and proves each b_{j,i} is binary via the quadratic relation b·(1-b) = 0.
type RangeProof struct {
	// LinearPart proves the main relation A·s = t
	LinearPart LinearProof

	// BitProofs proves that each decomposition bit is binary
	BitProofs []LinearProof

	// ReconProof proves the reconstruction: Σ 2^i · b_i = s + β
	ReconProof LinearProof
}

// ComposedProof binds multiple sub-proofs under a single Fiat-Shamir chain.
type ComposedProof struct {
	Proofs []LinearProof
	// Seed used to derive all challenges (ensures binding)
	TranscriptSeed []byte
}
