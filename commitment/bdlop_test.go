package commitment

import (
	"crypto/rand"
	"testing"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

func TestCommitVerifyRoundtrip(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)

	if !Verify(ck, com, opening) {
		t.Fatal("Valid commitment failed verification")
	}
}

func TestVerifyRejectsBadRandomness(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)

	// Tamper with randomness
	opening.R[0][0] = (opening.R[0][0] + 1) % ring.Q
	if Verify(ck, com, opening) {
		t.Fatal("Verification should fail with tampered randomness")
	}
}

func TestVerifyRejectsBadMessage(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)

	// Tamper with message
	opening.Msg[0][0] = (opening.Msg[0][0] + 1) % ring.Q
	if Verify(ck, com, opening) {
		t.Fatal("Verification should fail with tampered message")
	}
}

func TestVerifyRejectsBadCommitment(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)

	// Tamper with commitment
	com.T0[0][0] = (com.T0[0][0] + 1) % ring.Q
	if Verify(ck, com, opening) {
		t.Fatal("Verification should fail with tampered commitment")
	}
}

func TestHomomorphicAdd(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg1 := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	msg2 := sampler.SampleTernaryVec(params.Ell, rand.Reader)

	com1, op1 := Commit(ck, msg1, rand.Reader)
	com2, op2 := Commit(ck, msg2, rand.Reader)

	// Homomorphic addition
	comSum := AddCommitments(com1, com2)
	opSum := AddOpenings(op1, op2)

	if !Verify(ck, comSum, opSum) {
		t.Fatal("Homomorphic addition: combined commitment failed verification")
	}

	// Check that the combined message is msg1 + msg2
	expectedMsg := ring.VecAdd(msg1, msg2)
	if !ring.VecEqual(opSum.Msg, expectedMsg) {
		t.Fatal("Homomorphic addition: combined message mismatch")
	}
}

func TestScalarMulCommitment(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)

	// Scalar multiply by a small polynomial
	s := ring.One()
	comScaled := ScalarMulCommitment(com, s)

	// Verify that scaling by 1 preserves the commitment
	if !ring.VecEqual(comScaled.T0, com.T0) {
		t.Fatal("ScalarMul by 1: T0 changed")
	}
	if !ring.VecEqual(comScaled.Tm, com.Tm) {
		t.Fatal("ScalarMul by 1: Tm changed")
	}

	// Scale by 2 and check consistency
	two := ring.ScalarMul(ring.One(), 2)
	comDouble := ScalarMulCommitment(com, two)
	comAddSelf := AddCommitments(com, com)

	if !ring.VecEqual(comDouble.T0, comAddSelf.T0) {
		t.Fatal("2*C != C + C (T0)")
	}
	if !ring.VecEqual(comDouble.Tm, comAddSelf.Tm) {
		t.Fatal("2*C != C + C (Tm)")
	}

	_ = opening // keep the linter happy
}

func TestMultipleCommitmentsIndependent(t *testing.T) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)

	msg1 := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	msg2 := sampler.SampleTernaryVec(params.Ell, rand.Reader)

	com1, op1 := Commit(ck, msg1, rand.Reader)
	com2, op2 := Commit(ck, msg2, rand.Reader)

	// Each commitment verifies with its own opening
	if !Verify(ck, com1, op1) {
		t.Fatal("First commitment failed verification")
	}
	if !Verify(ck, com2, op2) {
		t.Fatal("Second commitment failed verification")
	}

	// Cross-verify should fail
	if Verify(ck, com1, op2) {
		t.Fatal("Cross-verification should fail")
	}
	if Verify(ck, com2, op1) {
		t.Fatal("Cross-verification should fail")
	}
}

func BenchmarkCommit(b *testing.B) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Commit(ck, msg, rand.Reader)
	}
}

func BenchmarkVerify(b *testing.B) {
	params := DefaultParams()
	ck := KeyGen(params, rand.Reader)
	msg := sampler.SampleTernaryVec(params.Ell, rand.Reader)
	com, opening := Commit(ck, msg, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Verify(ck, com, opening)
	}
}
