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
func HashToChallenge(params *Params, A ring.PolyMat, t, w ring.PolyVec) ring.Poly {
	seed := hashTranscript(params.Ring, A, t, w)
	return sampler.SampleChallenge(params.Ring, seed, params.Kappa)
}

func hashTranscript(r *ring.Ring, A ring.PolyMat, t, w ring.PolyVec) []byte {
	h := sha3.NewShake256()
	h.Write([]byte(domainSep))
	writeInt(h, len(A))
	for i := range A {
		writeVec(h, r, A[i])
	}
	writeVec(h, r, t)
	writeVec(h, r, w)
	seed := make([]byte, 64)
	h.Read(seed)
	return seed
}

func hashTranscriptWithTag(r *ring.Ring, tag []byte, A ring.PolyMat, t, w ring.PolyVec) []byte {
	h := sha3.NewShake256()
	h.Write([]byte(domainSep))
	h.Write(tag)
	writeInt(h, len(A))
	for i := range A {
		writeVec(h, r, A[i])
	}
	writeVec(h, r, t)
	writeVec(h, r, w)
	seed := make([]byte, 64)
	h.Read(seed)
	return seed
}

func writeVec(h sha3.ShakeHash, r *ring.Ring, v ring.PolyVec) {
	writeInt(h, len(v))
	var buf [8]byte
	for i := range v {
		for j := 0; j < r.N; j++ {
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
