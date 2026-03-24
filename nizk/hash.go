package nizk

import (
	"encoding/binary"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
	"golang.org/x/crypto/sha3"
)

const domainSep = "lnp22-challenge-v1"

// HashToChallenge derives a deterministic challenge polynomial from the
// proof transcript using SHAKE256 with domain separation.
// The challenge has exactly kappa non-zero coefficients in {-1, +1}.
func HashToChallenge(params *Params, A ring.PolyMat, t, w ring.PolyVec) *ring.Poly {
	seed := hashTranscript(A, t, w)
	return sampler.SampleChallenge(seed, params.Kappa)
}

// hashTranscript computes SHAKE256(domainSep || serialize(A, t, w)).
func hashTranscript(A ring.PolyMat, t, w ring.PolyVec) []byte {
	h := sha3.NewShake256()
	h.Write([]byte(domainSep))

	// Serialize matrix A
	writeInt(h, len(A))
	for i := range A {
		writeVec(h, A[i])
	}

	// Serialize vectors t and w
	writeVec(h, t)
	writeVec(h, w)

	seed := make([]byte, 64)
	h.Read(seed)
	return seed
}

// hashTranscriptWithTag extends the hash with an additional tag for composed proofs.
func hashTranscriptWithTag(tag []byte, A ring.PolyMat, t, w ring.PolyVec) []byte {
	h := sha3.NewShake256()
	h.Write([]byte(domainSep))
	h.Write(tag)

	writeInt(h, len(A))
	for i := range A {
		writeVec(h, A[i])
	}
	writeVec(h, t)
	writeVec(h, w)

	seed := make([]byte, 64)
	h.Read(seed)
	return seed
}

func writeVec(h sha3.ShakeHash, v ring.PolyVec) {
	writeInt(h, len(v))
	var buf [8]byte
	for i := range v {
		for j := 0; j < ring.N; j++ {
			binary.LittleEndian.PutUint64(buf[:], uint64(v[i][j]))
			h.Write(buf[:])
		}
	}
}

func writeInt(h sha3.ShakeHash, v int) {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], uint64(v))
	h.Write(buf[:])
}
