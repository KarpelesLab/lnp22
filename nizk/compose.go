package nizk

import (
	"encoding/binary"
	"errors"
	"io"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
	"golang.org/x/crypto/sha3"
)

// ComposeProofs generates a composed proof binding multiple linear-relation
// sub-proofs under a single Fiat-Shamir transcript chain.
//
// Each statement/witness pair is proved independently, but the challenges
// are chained so that all proofs are bound together and cannot be reused
// in a different context.
func ComposeProofs(params *Params, stmts []*Statement, wits []*Witness, rng io.Reader) (*ComposedProof, error) {
	if len(stmts) != len(wits) {
		return nil, errors.New("nizk: statements and witnesses count mismatch")
	}
	if err := params.Validate(); err != nil {
		return nil, err
	}

	// Generate a random transcript seed for binding
	transcriptSeed := make([]byte, 32)
	if _, err := io.ReadFull(rng, transcriptSeed); err != nil {
		return nil, err
	}

	proofs := make([]LinearProof, len(stmts))
	for i := range stmts {
		proof, err := proveLinearWithTag(params, stmts[i], wits[i], rng, composeTag(transcriptSeed, i))
		if err != nil {
			return nil, err
		}
		proofs[i] = *proof
	}

	return &ComposedProof{
		Proofs:         proofs,
		TranscriptSeed: transcriptSeed,
	}, nil
}

// proveLinearWithTag is like ProveLinear but uses a tagged hash for the challenge.
func proveLinearWithTag(params *Params, stmt *Statement, wit *Witness, rng io.Reader, tag []byte) (*LinearProof, error) {
	if err := validateStatement(params, stmt); err != nil {
		return nil, err
	}
	if err := validateWitness(params, stmt, wit); err != nil {
		return nil, err
	}

	maxAttempts := params.MaxAttempts
	if maxAttempts == 0 {
		maxAttempts = 10000
	}

	for attempt := 0; attempt < maxAttempts; attempt++ {
		y := sampler.SampleGaussianVec(params.L, params.Sigma, rng)
		w := ring.MatVecMul(stmt.A, y)

		// Challenge with tag
		seed := hashTranscriptWithTag(tag, stmt.A, stmt.T, w)
		c := hashSeedToChallenge(params, seed)

		z := ring.NewPolyVec(params.L)
		for i := 0; i < params.L; i++ {
			cs := ring.Mul(c, &wit.S[i])
			z[i] = *ring.Add(&y[i], cs)
		}

		if ring.VecInfNorm(z) <= params.BoundZ {
			return &LinearProof{W: w, Z: z}, nil
		}
	}

	return nil, errors.New("nizk: rejection sampling exhausted MaxAttempts")
}

// hashSeedToChallenge converts a seed to a challenge polynomial.
func hashSeedToChallenge(params *Params, seed []byte) *ring.Poly {
	return sampler.SampleChallenge(seed, params.Kappa)
}

// composeTag creates a unique tag for the i-th sub-proof in a composed proof.
func composeTag(transcriptSeed []byte, index int) []byte {
	h := sha3.NewShake256()
	h.Write([]byte("lnp22-compose-tag"))
	h.Write(transcriptSeed)
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], uint64(index))
	h.Write(buf[:])

	tag := make([]byte, 32)
	h.Read(tag)
	return tag
}
