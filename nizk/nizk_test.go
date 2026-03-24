package nizk

import (
	"crypto/rand"
	"testing"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

func makeValidInstance(params *Params) (*Statement, *Witness) {
	r := params.Ring
	A := sampler.SampleUniformMat(r, params.K, params.L, rand.Reader)
	s := sampler.SampleTernaryVec(r, params.L, rand.Reader)
	t := r.MatVecMul(A, s)
	return &Statement{A: A, T: t}, &Witness{S: s}
}

// --- Linear Proof Tests ---

func TestLinearProveVerify(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveLinear failed: %v", err)
	}
	if !VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear returned false for valid proof")
	}
}

func TestLinearProveVerifyMultiple(t *testing.T) {
	params := DefaultParams()
	for trial := 0; trial < 10; trial++ {
		stmt, wit := makeValidInstance(params)
		proof, err := ProveLinear(params, stmt, wit, rand.Reader)
		if err != nil {
			t.Fatalf("Trial %d: ProveLinear failed: %v", trial, err)
		}
		if !VerifyLinear(params, stmt, proof) {
			t.Fatalf("Trial %d: VerifyLinear returned false", trial)
		}
	}
}

func TestLinearRejectTamperedZ(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.Z[0][0] = (proof.Z[0][0] + 1) % params.Ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear should reject tampered z")
	}
}

func TestLinearRejectTamperedW(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.W[0][0] = (proof.W[0][0] + 1) % params.Ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear should reject tampered w")
	}
}

func TestLinearRejectTamperedStatement(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	stmt.T[0][0] = (stmt.T[0][0] + 1) % params.Ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear should reject tampered statement")
	}
}

func TestLinearRejectWrongWitness(t *testing.T) {
	params := DefaultParams()
	stmt, _ := makeValidInstance(params)
	wrongS := sampler.SampleTernaryVec(params.Ring, params.L, rand.Reader)
	_, err := ProveLinear(params, stmt, &Witness{S: wrongS}, rand.Reader)
	if err == nil {
		t.Fatal("ProveLinear should fail with invalid witness")
	}
}

func TestLinearRejectOversizedWitness(t *testing.T) {
	params := DefaultParams()
	r := params.Ring
	A := sampler.SampleUniformMat(r, params.K, params.L, rand.Reader)
	s := r.NewPolyVec(params.L)
	s[0][0] = params.Beta + 1
	target := r.MatVecMul(A, s)
	_, err := ProveLinear(params, &Statement{A: A, T: target}, &Witness{S: s}, rand.Reader)
	if err == nil {
		t.Fatal("ProveLinear should reject witness exceeding norm bound")
	}
}

func TestLinearNormBoundEnforced(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	if params.Ring.VecInfNorm(proof.Z) > params.BoundZ {
		t.Errorf("Proof z norm %d exceeds BoundZ %d", params.Ring.VecInfNorm(proof.Z), params.BoundZ)
	}
}

func TestLinearDifferentParams(t *testing.T) {
	paramSets := []*Params{
		{Ring: ring.Dilithium(), K: 2, L: 3, Kappa: 40, Beta: 1, Sigma: 200, BoundZ: 800, MaxAttempts: 1000},
		{Ring: ring.Dilithium(), K: 6, L: 8, Kappa: 60, Beta: 1, Sigma: 400, BoundZ: 1600, MaxAttempts: 1000},
	}
	for i, params := range paramSets {
		stmt, wit := makeValidInstance(params)
		proof, err := ProveLinear(params, stmt, wit, rand.Reader)
		if err != nil {
			t.Fatalf("ParamSet %d: ProveLinear failed: %v", i, err)
		}
		if !VerifyLinear(params, stmt, proof) {
			t.Fatalf("ParamSet %d: VerifyLinear returned false", i)
		}
	}
}

// --- Tests with q=7933, d=512 (no NTT, Karatsuba multiplication) ---

func smallRingParams() *Params {
	r, _ := ring.New(512, 7933)
	return &Params{
		Ring:        r,
		K:           3,
		L:           4,
		Kappa:       40,
		Beta:        1,
		Sigma:       200,
		BoundZ:      800,
		MaxAttempts: 1000,
	}
}

func TestLinearProveVerifySmallRing(t *testing.T) {
	params := smallRingParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveLinear (q=7933 d=512) failed: %v", err)
	}
	if !VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear (q=7933 d=512) returned false")
	}
}

func TestLinearRejectTamperedSmallRing(t *testing.T) {
	params := smallRingParams()
	stmt, wit := makeValidInstance(params)
	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.Z[0][0] = (proof.Z[0][0] + 1) % params.Ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear (q=7933) should reject tampered z")
	}
}

func TestLinearMultipleSmallRing(t *testing.T) {
	params := smallRingParams()
	for trial := 0; trial < 5; trial++ {
		stmt, wit := makeValidInstance(params)
		proof, err := ProveLinear(params, stmt, wit, rand.Reader)
		if err != nil {
			t.Fatalf("Trial %d: %v", trial, err)
		}
		if !VerifyLinear(params, stmt, proof) {
			t.Fatalf("Trial %d: verification failed", trial)
		}
	}
}

// --- Range Proof Tests ---

func TestRangeProveVerify(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	rangeStmt := &RangeStatement{Statement: *stmt, Beta: params.Beta}
	proof, err := ProveRange(params, rangeStmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveRange failed: %v", err)
	}
	if !VerifyRange(params, rangeStmt, proof) {
		t.Fatal("VerifyRange returned false for valid proof")
	}
}

func TestRangeProveVerifyLargerBeta(t *testing.T) {
	r := ring.Dilithium()
	params := &Params{Ring: r, K: 3, L: 4, Kappa: 40, Beta: 3, Sigma: 350, BoundZ: 1400, MaxAttempts: 1000}

	A := sampler.SampleUniformMat(r, params.K, params.L, rand.Reader)
	s := r.NewPolyVec(params.L)
	for i := range s {
		for j := 0; j < r.N; j++ {
			var buf [1]byte
			rand.Read(buf[:])
			v := int64(buf[0]%7) - 3
			s[i][j] = (v + r.Q) % r.Q
		}
	}
	target := r.MatVecMul(A, s)

	stmt := &RangeStatement{Statement: Statement{A: A, T: target}, Beta: 3}
	proof, err := ProveRange(params, stmt, &Witness{S: s}, rand.Reader)
	if err != nil {
		t.Fatalf("ProveRange failed: %v", err)
	}
	if !VerifyRange(params, stmt, proof) {
		t.Fatal("VerifyRange returned false for valid proof")
	}
}

// --- Composed Proof Tests ---

func TestComposeProveVerify(t *testing.T) {
	params := DefaultParams()
	stmts := make([]*Statement, 3)
	wits := make([]*Witness, 3)
	for i := 0; i < 3; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}
	proof, err := ComposeProofs(params, stmts, wits, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	if !VerifyComposed(params, stmts, proof) {
		t.Fatal("VerifyComposed returned false for valid composed proof")
	}
}

func TestComposeRejectTampered(t *testing.T) {
	params := DefaultParams()
	stmts := make([]*Statement, 2)
	wits := make([]*Witness, 2)
	for i := 0; i < 2; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}
	proof, err := ComposeProofs(params, stmts, wits, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.Proofs[0].Z[0][0] = (proof.Proofs[0].Z[0][0] + 1) % params.Ring.Q
	if VerifyComposed(params, stmts, proof) {
		t.Fatal("VerifyComposed should reject tampered sub-proof")
	}
}

func TestComposeRejectSwappedStatements(t *testing.T) {
	params := DefaultParams()
	stmts := make([]*Statement, 2)
	wits := make([]*Witness, 2)
	for i := 0; i < 2; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}
	proof, err := ComposeProofs(params, stmts, wits, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	stmts[0], stmts[1] = stmts[1], stmts[0]
	if VerifyComposed(params, stmts, proof) {
		t.Fatal("VerifyComposed should reject swapped statements")
	}
}

func TestComposeRejectTamperedSeed(t *testing.T) {
	params := DefaultParams()
	stmts := make([]*Statement, 2)
	wits := make([]*Witness, 2)
	for i := 0; i < 2; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}
	proof, err := ComposeProofs(params, stmts, wits, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.TranscriptSeed[0] ^= 0xFF
	if VerifyComposed(params, stmts, proof) {
		t.Fatal("VerifyComposed should reject tampered transcript seed")
	}
}

func TestComposeSmallRing(t *testing.T) {
	params := smallRingParams()
	stmts := make([]*Statement, 2)
	wits := make([]*Witness, 2)
	for i := 0; i < 2; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}
	proof, err := ComposeProofs(params, stmts, wits, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	if !VerifyComposed(params, stmts, proof) {
		t.Fatal("VerifyComposed (q=7933 d=512) returned false")
	}
}

// --- Benchmarks ---

func BenchmarkProveLinear(b *testing.B) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ProveLinear(params, stmt, wit, rand.Reader)
	}
}

func BenchmarkVerifyLinear(b *testing.B) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	proof, _ := ProveLinear(params, stmt, wit, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		VerifyLinear(params, stmt, proof)
	}
}

func BenchmarkProveLinearSmallRing(b *testing.B) {
	params := smallRingParams()
	stmt, wit := makeValidInstance(params)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ProveLinear(params, stmt, wit, rand.Reader)
	}
}

func BenchmarkProveRange(b *testing.B) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)
	rangeStmt := &RangeStatement{Statement: *stmt, Beta: params.Beta}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ProveRange(params, rangeStmt, wit, rand.Reader)
	}
}

func BenchmarkComposeProofs(b *testing.B) {
	params := DefaultParams()
	stmts := make([]*Statement, 3)
	wits := make([]*Witness, 3)
	for i := 0; i < 3; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ComposeProofs(params, stmts, wits, rand.Reader)
	}
}
