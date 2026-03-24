package commitment

import (
	"crypto/rand"
	"testing"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

var dilithium = ring.Dilithium()

func TestCommitVerifyRoundtrip(t *testing.T) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	if !Verify(ck, com, opening) {
		t.Fatal("Valid commitment failed verification")
	}
}

func TestVerifyRejectsBadRandomness(t *testing.T) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	opening.R[0][0] = (opening.R[0][0] + 1) % dilithium.Q
	if Verify(ck, com, opening) {
		t.Fatal("Verification should fail with tampered randomness")
	}
}

func TestVerifyRejectsBadMessage(t *testing.T) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	opening.Msg[0][0] = (opening.Msg[0][0] + 1) % dilithium.Q
	if Verify(ck, com, opening) {
		t.Fatal("Verification should fail with tampered message")
	}
}

func TestVerifyRejectsBadCommitment(t *testing.T) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	com.T0[0][0] = (com.T0[0][0] + 1) % dilithium.Q
	if Verify(ck, com, opening) {
		t.Fatal("Verification should fail with tampered commitment")
	}
}

func TestHomomorphicAdd(t *testing.T) {
	r := dilithium
	params := DefaultParams(r)
	ck := KeyGen(params, rand.Reader)

	msg1 := sampler.SampleTernaryVec(r, params.Ell, rand.Reader)
	msg2 := sampler.SampleTernaryVec(r, params.Ell, rand.Reader)
	com1, op1 := Commit(ck, msg1, rand.Reader)
	com2, op2 := Commit(ck, msg2, rand.Reader)

	comSum := AddCommitments(r, com1, com2)
	opSum := AddOpenings(r, op1, op2)

	if !Verify(ck, comSum, opSum) {
		t.Fatal("Homomorphic addition: combined commitment failed verification")
	}
	if !r.VecEqual(opSum.Msg, r.VecAdd(msg1, msg2)) {
		t.Fatal("Homomorphic addition: combined message mismatch")
	}
}

func TestScalarMulCommitment(t *testing.T) {
	r := dilithium
	params := DefaultParams(r)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(r, params.Ell, rand.Reader)
	com, _ := Commit(ck, msg, rand.Reader)

	comScaled := ScalarMulCommitment(r, com, r.One())
	if !r.VecEqual(comScaled.T0, com.T0) {
		t.Fatal("ScalarMul by 1: T0 changed")
	}

	two := r.ScalarMul(r.One(), 2)
	comDouble := ScalarMulCommitment(r, com, two)
	comAddSelf := AddCommitments(r, com, com)
	if !r.VecEqual(comDouble.T0, comAddSelf.T0) {
		t.Fatal("2*C != C + C (T0)")
	}
	if !r.VecEqual(comDouble.Tm, comAddSelf.Tm) {
		t.Fatal("2*C != C + C (Tm)")
	}
}

func TestMultipleCommitmentsIndependent(t *testing.T) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)

	msg1 := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	msg2 := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	com1, op1 := Commit(ck, msg1, rand.Reader)
	com2, op2 := Commit(ck, msg2, rand.Reader)

	if !Verify(ck, com1, op1) {
		t.Fatal("First commitment failed")
	}
	if !Verify(ck, com2, op2) {
		t.Fatal("Second commitment failed")
	}
	if Verify(ck, com1, op2) {
		t.Fatal("Cross-verification should fail")
	}
	if Verify(ck, com2, op1) {
		t.Fatal("Cross-verification should fail")
	}
}

func TestCommitSmallRing(t *testing.T) {
	r, _ := ring.New(64, 7933)
	params := DefaultParams(r)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(r, params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	if !Verify(ck, com, opening) {
		t.Fatal("Commitment on small ring failed verification")
	}
}

func BenchmarkCommit(b *testing.B) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Commit(ck, msg, rand.Reader)
	}
}

func BenchmarkVerify(b *testing.B) {
	params := DefaultParams(dilithium)
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(dilithium, params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Verify(ck, com, opening)
	}
}
