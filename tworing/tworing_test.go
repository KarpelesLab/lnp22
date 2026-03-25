package tworing

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
)

var testRing, _ = ring.NewBig(64, big.NewInt(8380417))

func testParams() *Params {
	return &Params{
		Ring:        testRing,
		K:           3,
		L:           4,
		Kappa:       30,
		Sigma:       200,
		BoundZ:      big.NewInt(800),
		L2BoundSq:   big.NewInt(256), // N * 1² = 64 per poly, 4 polys = 256
		MaxAttempts: 1000,
	}
}

func makeLinearInstance(params *Params) (*Statement, *Witness) {
	r := params.Ring
	A := sampler.SampleBigUniformMat(r, params.K, params.L, rand.Reader)
	s := sampler.SampleBigTernaryVec(r, params.L, rand.Reader)
	t := r.MatVecMul(A, s)
	stmt := &Statement{
		Linear: []LinearStatement{{A: A, T: t}},
	}
	return stmt, &Witness{S: s}
}

func makeQuadraticInstance(params *Params) (*Statement, *Witness) {
	r := params.Ring
	A := sampler.SampleBigUniformMat(r, params.K, params.L, rand.Reader)
	s := sampler.SampleBigTernaryVec(r, params.L, rand.Reader)
	t := r.MatVecMul(A, s)

	// B = identity-like: B[i][j] = delta_{ij}
	B := r.NewPolyMat(params.L, params.L)
	for i := 0; i < params.L; i++ {
		B[i][i] = r.One()
	}
	// v = ⟨s, B·s⟩ = ⟨s, s⟩ = ||s||² in the ring
	bs := r.MatVecMul(B, s)
	v := r.InnerProduct(s, bs)

	stmt := &Statement{
		Linear:    []LinearStatement{{A: A, T: t}},
		Quadratic: []QuadraticStatement{{B: B, V: v}},
	}
	return stmt, &Witness{S: s}
}

// --- Linear-only tests ---

func TestLinearProveVerify(t *testing.T) {
	params := testParams()
	stmt, wit := makeLinearInstance(params)
	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("Prove failed: %v", err)
	}
	if !Verify(params, stmt, proof) {
		t.Fatal("Verify returned false for valid proof")
	}
}

func TestLinearMultiple(t *testing.T) {
	params := testParams()
	for trial := 0; trial < 5; trial++ {
		stmt, wit := makeLinearInstance(params)
		proof, err := Prove(params, stmt, wit, rand.Reader)
		if err != nil {
			t.Fatalf("Trial %d: %v", trial, err)
		}
		if !Verify(params, stmt, proof) {
			t.Fatalf("Trial %d: verification failed", trial)
		}
	}
}

func TestLinearRejectTamperedZ(t *testing.T) {
	params := testParams()
	stmt, wit := makeLinearInstance(params)
	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.Z[0][0].Add(proof.Z[0][0], big.NewInt(1))
	proof.Z[0][0].Mod(proof.Z[0][0], params.Ring.Q)
	if Verify(params, stmt, proof) {
		t.Fatal("Should reject tampered z")
	}
}

func TestLinearRejectTamperedW(t *testing.T) {
	params := testParams()
	stmt, wit := makeLinearInstance(params)
	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	proof.W[0][0].Add(proof.W[0][0], big.NewInt(1))
	proof.W[0][0].Mod(proof.W[0][0], params.Ring.Q)
	if Verify(params, stmt, proof) {
		t.Fatal("Should reject tampered w")
	}
}

func TestLinearRejectTamperedStatement(t *testing.T) {
	params := testParams()
	stmt, wit := makeLinearInstance(params)
	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	stmt.Linear[0].T[0][0].Add(stmt.Linear[0].T[0][0], big.NewInt(1))
	stmt.Linear[0].T[0][0].Mod(stmt.Linear[0].T[0][0], params.Ring.Q)
	if Verify(params, stmt, proof) {
		t.Fatal("Should reject tampered statement")
	}
}

func TestLinearRejectWrongWitness(t *testing.T) {
	params := testParams()
	stmt, _ := makeLinearInstance(params)
	wrongS := sampler.SampleBigTernaryVec(params.Ring, params.L, rand.Reader)
	_, err := Prove(params, stmt, &Witness{S: wrongS}, rand.Reader)
	if err == nil {
		t.Fatal("Should reject wrong witness")
	}
}

// --- Quadratic tests ---

func TestQuadraticProveVerify(t *testing.T) {
	params := testParams()
	stmt, wit := makeQuadraticInstance(params)
	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("Prove failed: %v", err)
	}
	if !Verify(params, stmt, proof) {
		t.Fatal("Verify returned false for valid quadratic proof")
	}
}

func TestQuadraticRejectTampered(t *testing.T) {
	params := testParams()
	stmt, wit := makeQuadraticInstance(params)
	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	// Tamper with quadratic commitment
	proof.Wq[0][0].Add(proof.Wq[0][0], big.NewInt(1))
	proof.Wq[0][0].Mod(proof.Wq[0][0], params.Ring.Q)
	if Verify(params, stmt, proof) {
		t.Fatal("Should reject tampered quadratic commitment")
	}
}

// --- Norm-bound tests ---

func TestNormBoundProveVerify(t *testing.T) {
	params := testParams()
	r := params.Ring
	A := sampler.SampleBigUniformMat(r, params.K, params.L, rand.Reader)
	s := sampler.SampleBigTernaryVec(r, params.L, rand.Reader)
	t_vec := r.MatVecMul(A, s)

	// Compute actual ℓ₂-norm² of witness
	l2sq := r.VecL2NormSq(s)
	// Set bound slightly above actual norm
	bound := new(big.Int).Add(l2sq, big.NewInt(10))

	stmt := &Statement{
		Linear: []LinearStatement{{A: A, T: t_vec}},
		Norm:   &NormStatement{L2BoundSq: bound},
	}
	wit := &Witness{S: s}

	proof, err := Prove(params, stmt, wit, rand.Reader)
	if err != nil {
		t.Fatalf("Prove with norm bound failed: %v", err)
	}
	if !Verify(params, stmt, proof) {
		t.Fatal("Verify returned false for valid norm-bounded proof")
	}
}

func TestNormBoundRejectExceeded(t *testing.T) {
	params := testParams()
	r := params.Ring
	A := sampler.SampleBigUniformMat(r, params.K, params.L, rand.Reader)
	s := sampler.SampleBigTernaryVec(r, params.L, rand.Reader)
	t_vec := r.MatVecMul(A, s)

	// Set bound too small — witness violates it
	stmt := &Statement{
		Linear: []LinearStatement{{A: A, T: t_vec}},
		Norm:   &NormStatement{L2BoundSq: big.NewInt(1)}, // too small for ternary
	}
	wit := &Witness{S: s}

	_, err := Prove(params, stmt, wit, rand.Reader)
	if err == nil {
		t.Fatal("Should reject witness exceeding norm bound")
	}
}

// --- Benchmarks ---

func BenchmarkProveLinear(b *testing.B) {
	params := testParams()
	stmt, wit := makeLinearInstance(params)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Prove(params, stmt, wit, rand.Reader)
	}
}

func BenchmarkVerifyLinear(b *testing.B) {
	params := testParams()
	stmt, wit := makeLinearInstance(params)
	proof, _ := Prove(params, stmt, wit, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Verify(params, stmt, proof)
	}
}

func BenchmarkProveQuadratic(b *testing.B) {
	params := testParams()
	stmt, wit := makeQuadraticInstance(params)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Prove(params, stmt, wit, rand.Reader)
	}
}
