package tworing

import (
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
)

// LinearStatement: prove knowledge of s with A·s = t.
type LinearStatement struct {
	A ring.BigPolyMat
	T ring.BigPolyVec
}

// QuadraticStatement: prove ⟨s, B·s⟩ = v (quadratic form).
type QuadraticStatement struct {
	B ring.BigPolyMat
	V ring.BigPoly
}

// NormStatement: prove ||s||₂² ≤ bound.
type NormStatement struct {
	L2BoundSq *big.Int
}

// Statement combines linear, quadratic, and norm constraints.
type Statement struct {
	Linear    []LinearStatement
	Quadratic []QuadraticStatement
	Norm      *NormStatement // optional ℓ₂-norm bound
}

// Witness is the prover's secret.
type Witness struct {
	S ring.BigPolyVec
}

// Proof is the complete two-ring NIZK proof.
type Proof struct {
	W  ring.BigPolyVec // commitment w = A·y for each linear statement
	Z  ring.BigPolyVec // response z = y + c·s
	Wq []ring.BigPoly  // quadratic cross-terms: ⟨y, B·s⟩ + ⟨s, B·y⟩
	Wy []ring.BigPoly  // quadratic self-terms: ⟨y, B·y⟩
}
