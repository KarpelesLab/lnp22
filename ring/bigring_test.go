package ring

import (
	"math/big"
	"math/rand/v2"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
)

var bigDilithium, _ = NewBig(256, big.NewInt(8380417))

func randBigPoly(rng *rand.Rand, r *BigRing) BigPoly {
	p := r.NewPoly()
	for i := range p {
		p[i] = bigmod.Reduce(big.NewInt(rng.Int64()), r.Q)
	}
	return p
}

// --- BigRing construction ---

func TestNewBigDilithium(t *testing.T) {
	r := bigDilithium
	if r.N != 256 {
		t.Fatalf("N = %d, want 256", r.N)
	}
	if !r.HasNTT() {
		t.Fatal("Dilithium BigRing should have NTT")
	}
}

func TestNewBigLNP22(t *testing.T) {
	r := LNP22Big()
	if r.N != 1024 {
		t.Fatalf("N = %d, want 1024", r.N)
	}
	if !r.HasNTT() {
		t.Fatal("LNP22 BigRing should have NTT (q ≡ 1 mod 2048)")
	}
	if r.Q.BitLen() != 152 {
		t.Fatalf("Q bit length = %d, want 152", r.Q.BitLen())
	}
}

func TestNewBigNoNTT(t *testing.T) {
	r, err := NewBig(512, big.NewInt(7933))
	if err != nil {
		t.Fatal(err)
	}
	if r.HasNTT() {
		t.Fatal("q=7933 d=512 should NOT have NTT")
	}
}

// --- BigPoly basic ops ---

func TestBigAddSub(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	b := randBigPoly(rng, r)
	sum := r.Add(a, b)
	diff := r.Sub(sum, b)
	if !r.Equal(diff, a) {
		t.Error("(a + b) - b != a")
	}
}

func TestBigNeg(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	negA := r.Neg(a)
	sum := r.Add(a, negA)
	if !r.Equal(sum, r.NewPoly()) {
		t.Error("a + (-a) != 0")
	}
}

func TestBigScalarMul(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	double := r.ScalarMul(a, big.NewInt(2))
	sum := r.Add(a, a)
	if !r.Equal(double, sum) {
		t.Error("2*a != a + a")
	}
}

func TestBigInfNorm(t *testing.T) {
	r := bigDilithium
	p := r.NewPoly()
	p[0].SetInt64(10)
	p[1].Sub(r.Q, big.NewInt(5)) // -5
	got := r.InfNorm(p)
	if got.Cmp(big.NewInt(10)) != 0 {
		t.Errorf("InfNorm = %s, want 10", got)
	}
}

func TestBigEqual(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	b := r.NewPoly()
	for i := range a {
		b[i].Set(a[i])
	}
	if !r.Equal(a, b) {
		t.Error("Equal failed on identical polys")
	}
	b[0].Add(b[0], big.NewInt(1))
	b[0].Mod(b[0], r.Q)
	if r.Equal(a, b) {
		t.Error("Equal returned true for different polys")
	}
}

// --- BigRing NTT tests ---

func TestBigNTTRoundtrip(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 20; trial++ {
		a := randBigPoly(rng, r)
		aRed := r.Reduce(a)
		recovered := r.InvNTT(r.NTT(a))
		if !r.Equal(recovered, aRed) {
			t.Fatalf("InvNTT(NTT(a)) != a on trial %d", trial)
		}
	}
}

func TestBigNTTMulMatchesNaive(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 10; trial++ {
		a := randBigPoly(rng, r)
		b := randBigPoly(rng, r)
		nttResult := r.Mul(a, b)
		naiveResult := r.NaiveMul(a, b)
		if !r.Equal(nttResult, naiveResult) {
			t.Fatalf("NTT mul != naive mul on trial %d", trial)
		}
	}
}

func TestBigMulOne(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	if !r.Equal(r.Mul(a, r.One()), r.Reduce(a)) {
		t.Error("a * 1 != a")
	}
}

func TestBigMulCommutative(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 10; trial++ {
		a := randBigPoly(rng, r)
		b := randBigPoly(rng, r)
		if !r.Equal(r.Mul(a, b), r.Mul(b, a)) {
			t.Fatalf("a*b != b*a on trial %d", trial)
		}
	}
}

func TestBigMulDistributive(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 10; trial++ {
		a := randBigPoly(rng, r)
		b := randBigPoly(rng, r)
		c := randBigPoly(rng, r)
		lhs := r.Mul(a, r.Add(b, c))
		rhs := r.Add(r.Mul(a, b), r.Mul(a, c))
		if !r.Equal(lhs, rhs) {
			t.Fatalf("distributivity failed on trial %d", trial)
		}
	}
}

// --- Karatsuba with big.Int (small q, no NTT) ---

func TestBigKaratsubaMulMatchesNaive(t *testing.T) {
	r, _ := NewBig(64, big.NewInt(7933))
	rng := rand.New(rand.NewPCG(99, 0))
	for trial := 0; trial < 10; trial++ {
		a := randBigPoly(rng, r)
		b := randBigPoly(rng, r)
		if !r.Equal(r.Mul(a, b), r.NaiveMul(a, b)) {
			t.Fatalf("Karatsuba != naive on trial %d", trial)
		}
	}
}

// --- 152-bit NTT tests ---

func TestLNP22BigNTTRoundtrip(t *testing.T) {
	r := LNP22Big()
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	aRed := r.Reduce(a)
	recovered := r.InvNTT(r.NTT(a))
	if !r.Equal(recovered, aRed) {
		t.Fatal("152-bit NTT roundtrip failed")
	}
}

func TestLNP22BigMulOne(t *testing.T) {
	r := LNP22Big()
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	if !r.Equal(r.Mul(a, r.One()), r.Reduce(a)) {
		t.Error("152-bit: a * 1 != a")
	}
}

func TestLNP22BigMulMatchesNaiveSmall(t *testing.T) {
	// Use small degree to make naive tractable for 152-bit q
	q, _ := new(big.Int).SetString("5708990770823839524233143877797980545530982401", 10)
	r, _ := NewBig(16, q)
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 5; trial++ {
		a := randBigPoly(rng, r)
		b := randBigPoly(rng, r)
		if !r.Equal(r.Mul(a, b), r.NaiveMul(a, b)) {
			t.Fatalf("152-bit NTT mul != naive on trial %d", trial)
		}
	}
}

// --- BigVec tests ---

func TestBigVecAddSub(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := BigPolyVec{randBigPoly(rng, r), randBigPoly(rng, r)}
	b := BigPolyVec{randBigPoly(rng, r), randBigPoly(rng, r)}
	sum := r.VecAdd(a, b)
	diff := r.VecSub(sum, b)
	if !r.VecEqual(diff, a) {
		t.Error("VecAdd/VecSub roundtrip failed")
	}
}

func TestBigMatVecMul(t *testing.T) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	A := r.NewPolyMat(2, 3)
	s := r.NewPolyVec(3)
	for i := range A {
		for j := range A[i] {
			A[i][j] = randBigPoly(rng, r)
		}
	}
	for i := range s {
		s[i] = randBigPoly(rng, r)
	}
	result := r.MatVecMul(A, s)
	if len(result) != 2 {
		t.Fatalf("MatVecMul result length = %d, want 2", len(result))
	}
	expected := r.Add(r.Add(r.Mul(A[0][0], s[0]), r.Mul(A[0][1], s[1])), r.Mul(A[0][2], s[2]))
	if !r.Equal(result[0], expected) {
		t.Error("MatVecMul result[0] mismatch")
	}
}

// --- L2 norm tests ---

func TestBigL2NormSq(t *testing.T) {
	r := bigDilithium
	p := r.NewPoly()
	p[0].SetInt64(3)
	p[1].SetInt64(4)
	// ||p||_2^2 = 9 + 16 = 25
	got := r.L2NormSq(p)
	if got.Cmp(big.NewInt(25)) != 0 {
		t.Errorf("L2NormSq = %s, want 25", got)
	}
}

func TestBigVecL2NormSq(t *testing.T) {
	r := bigDilithium
	p1 := r.NewPoly()
	p2 := r.NewPoly()
	p1[0].SetInt64(3) // contributes 9
	p2[0].SetInt64(4) // contributes 16
	v := BigPolyVec{p1, p2}
	got := r.VecL2NormSq(v)
	if got.Cmp(big.NewInt(25)) != 0 {
		t.Errorf("VecL2NormSq = %s, want 25", got)
	}
}

// --- Benchmarks ---

func BenchmarkBigNTT(b *testing.B) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	p := randBigPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.NTT(p)
	}
}

func BenchmarkBigMulNTT(b *testing.B) {
	r := bigDilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	p := randBigPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.Mul(a, p)
	}
}

func BenchmarkBigMul152bit(b *testing.B) {
	r := LNP22Big()
	rng := rand.New(rand.NewPCG(42, 0))
	a := randBigPoly(rng, r)
	p := randBigPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.Mul(a, p)
	}
}
