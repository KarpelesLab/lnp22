package nizk

import (
	"crypto/rand"
	"testing"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

// makeValidInstance generates a random statement and matching witness.
// A is sampled uniformly, s is ternary (||s||_∞ = 1), t = A·s.
func makeValidInstance(params *Params) (*Statement, *Witness) {
	A := sampler.SampleUniformMat(params.K, params.L, rand.Reader)
	s := sampler.SampleTernaryVec(params.L, rand.Reader)
	t := ring.MatVecMul(A, s)
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
		t.Fatalf("ProveLinear failed: %v", err)
	}

	// Tamper with z
	proof.Z[0][0] = (proof.Z[0][0] + 1) % ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear should reject tampered z")
	}
}

func TestLinearRejectTamperedW(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)

	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveLinear failed: %v", err)
	}

	// Tamper with w — this changes the derived challenge
	proof.W[0][0] = (proof.W[0][0] + 1) % ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear should reject tampered w")
	}
}

func TestLinearRejectTamperedStatement(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)

	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveLinear failed: %v", err)
	}

	// Tamper with the statement (change t)
	stmt.T[0][0] = (stmt.T[0][0] + 1) % ring.Q
	if VerifyLinear(params, stmt, proof) {
		t.Fatal("VerifyLinear should reject tampered statement")
	}
}

func TestLinearRejectWrongWitness(t *testing.T) {
	params := DefaultParams()
	stmt, _ := makeValidInstance(params)

	// Create a different witness that doesn't satisfy A·s = t
	wrongS := sampler.SampleTernaryVec(params.L, rand.Reader)
	wrongWit := &Witness{S: wrongS}

	_, err := ProveLinear(params, stmt, wrongWit, rand.Reader)
	if err == nil {
		t.Fatal("ProveLinear should fail with invalid witness")
	}
}

func TestLinearRejectOversizedWitness(t *testing.T) {
	params := DefaultParams()
	A := sampler.SampleUniformMat(params.K, params.L, rand.Reader)

	// Create a witness with ||s||_∞ > β
	s := ring.NewPolyVec(params.L)
	s[0][0] = params.Beta + 1
	target := ring.MatVecMul(A, s)

	stmt := &Statement{A: A, T: target}
	wit := &Witness{S: s}

	_, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err == nil {
		t.Fatal("ProveLinear should reject witness exceeding norm bound")
	}
}

func TestLinearNormBoundEnforced(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)

	proof, err := ProveLinear(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveLinear failed: %v", err)
	}

	// The proof's z should respect the norm bound
	if ring.VecInfNorm(proof.Z) > params.BoundZ {
		t.Errorf("Proof z norm %d exceeds BoundZ %d", ring.VecInfNorm(proof.Z), params.BoundZ)
	}
}

func TestLinearDifferentParams(t *testing.T) {
	// Test with different parameter sets
	paramSets := []*Params{
		{K: 2, L: 3, Kappa: 40, Beta: 1, Sigma: 200, BoundZ: 800, MaxAttempts: 1000},
		{K: 6, L: 8, Kappa: 60, Beta: 1, Sigma: 400, BoundZ: 1600, MaxAttempts: 1000},
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

// --- Range Proof Tests ---

func TestRangeProveVerify(t *testing.T) {
	params := DefaultParams()
	stmt, wit := makeValidInstance(params)

	rangeStmt := &RangeStatement{
		Statement: *stmt,
		Beta:      params.Beta,
	}

	proof, err := ProveRange(params, rangeStmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("ProveRange failed: %v", err)
	}

	if !VerifyRange(params, rangeStmt, proof) {
		t.Fatal("VerifyRange returned false for valid proof")
	}
}

func TestRangeProveVerifyLargerBeta(t *testing.T) {
	// Test with β=3 (coefficients in [-3, 3])
	params := &Params{
		K:           3,
		L:           4,
		Kappa:       40,
		Beta:        3,
		Sigma:       350,
		BoundZ:      1400,
		MaxAttempts: 1000,
	}

	A := sampler.SampleUniformMat(params.K, params.L, rand.Reader)
	// Create witness with small coefficients in [-3, 3]
	s := ring.NewPolyVec(params.L)
	for i := range s {
		for j := 0; j < ring.N; j++ {
			// Use values in {-3, -2, -1, 0, 1, 2, 3}
			var buf [1]byte
			for {
				rand.Read(buf[:])
				v := int64(buf[0] % 7)
				s[i][j] = (v - 3 + ring.Q) % ring.Q
				break
			}
		}
	}
	target := ring.MatVecMul(A, s)

	stmt := &RangeStatement{
		Statement: Statement{A: A, T: target},
		Beta:      3,
	}
	wit := &Witness{S: s}

	proof, err := ProveRange(params, stmt, wit, rand.Reader)
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

	numProofs := 3
	stmts := make([]*Statement, numProofs)
	wits := make([]*Witness, numProofs)
	for i := 0; i < numProofs; i++ {
		stmts[i], wits[i] = makeValidInstance(params)
	}

	proof, err := ComposeProofs(params, stmts, wits, rand.Reader)
	if err != nil {
		t.Fatalf("ComposeProofs failed: %v", err)
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
		t.Fatalf("ComposeProofs failed: %v", err)
	}

	// Tamper with one sub-proof's z
	proof.Proofs[0].Z[0][0] = (proof.Proofs[0].Z[0][0] + 1) % ring.Q
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
		t.Fatalf("ComposeProofs failed: %v", err)
	}

	// Swap statements
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
		t.Fatalf("ComposeProofs failed: %v", err)
	}

	// Tamper with transcript seed
	proof.TranscriptSeed[0] ^= 0xFF
	if VerifyComposed(params, stmts, proof) {
		t.Fatal("VerifyComposed should reject tampered transcript seed")
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
