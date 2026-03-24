package nizk

import "github.com/KarpelesLab/lnp22/ring"

// Statement is the public part of a linear-relation proof: A·s ≡ t (mod q).
type Statement struct {
	A ring.PolyMat
	T ring.PolyVec
}

// Witness is the prover's secret: a short vector s with ||s||_∞ ≤ β.
type Witness struct {
	S ring.PolyVec
}

// LinearProof is a NIZK proof of knowledge of s such that A·s ≡ t (mod q) and ||s||_∞ ≤ β.
type LinearProof struct {
	W ring.PolyVec // commitment w = A·y
	Z ring.PolyVec // response z = y + c·s
}

// RangeStatement extends a linear statement with a bound on the coefficient range.
type RangeStatement struct {
	Statement
	Beta    int64
	NumBits int
}

// RangeProof proves that the witness coefficients lie in [-β, β].
type RangeProof struct {
	LinearPart LinearProof
	BitProofs  []LinearProof
	ReconProof LinearProof
}

// ComposedProof binds multiple sub-proofs under a single Fiat-Shamir chain.
type ComposedProof struct {
	Proofs         []LinearProof
	TranscriptSeed []byte
}
