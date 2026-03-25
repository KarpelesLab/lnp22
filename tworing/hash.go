package tworing

import (
	"encoding/binary"
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
	"golang.org/x/crypto/sha3"
)

const domainSep = "lnp22-tworing-v1"

// hashToChallenge derives a challenge polynomial from the two-ring transcript.
func hashToChallenge(params *Params, stmt *Statement, w ring.BigPolyVec, wq []ring.BigPoly) ring.BigPoly {
	h := sha3.NewShake256()
	h.Write([]byte(domainSep))

	// Hash linear statements
	writeInt(h, len(stmt.Linear))
	for _, ls := range stmt.Linear {
		writeInt(h, len(ls.A))
		for _, row := range ls.A {
			writeBigVec(h, params.Ring, row)
		}
		writeBigVec(h, params.Ring, ls.T)
	}

	// Hash quadratic statements
	writeInt(h, len(stmt.Quadratic))
	for _, qs := range stmt.Quadratic {
		writeInt(h, len(qs.B))
		for _, row := range qs.B {
			writeBigVec(h, params.Ring, row)
		}
		writeBigPoly(h, params.Ring, qs.V)
	}

	// Hash norm bound
	if stmt.Norm != nil {
		writeBigInt(h, stmt.Norm.L2BoundSq)
	}

	// Hash commitments
	writeBigVec(h, params.Ring, w)
	for _, wqi := range wq {
		writeBigPoly(h, params.Ring, wqi)
	}

	seed := make([]byte, 64)
	h.Read(seed)
	return sampler.SampleBigChallenge(params.Ring, seed, params.Kappa)
}

func writeBigVec(h sha3.ShakeHash, r *ring.BigRing, v ring.BigPolyVec) {
	writeInt(h, len(v))
	for i := range v {
		writeBigPoly(h, r, v[i])
	}
}

func writeBigPoly(h sha3.ShakeHash, r *ring.BigRing, p ring.BigPoly) {
	for j := 0; j < r.N; j++ {
		writeBigInt(h, p[j])
	}
}

func writeBigInt(h sha3.ShakeHash, v *big.Int) {
	b := v.Bytes()
	writeInt(h, len(b))
	h.Write(b)
}

func writeInt(h sha3.ShakeHash, v int) {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], uint64(v))
	h.Write(buf[:])
}
