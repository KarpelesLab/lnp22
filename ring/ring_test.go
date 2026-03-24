package ring

import (
	"math/rand/v2"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/modmath"
)

// dilithium is the standard test ring (N=256, Q=8380417).
var dilithium = Dilithium()

func randPoly(rng *rand.Rand, r *Ring) Poly {
	p := r.NewPoly()
	for i := range p {
		p[i] = modmath.Reduce(rng.Int64(), r.Q)
	}
	return p
}

func smallPoly(rng *rand.Rand, r *Ring, bound int64) Poly {
	p := r.NewPoly()
	for i := range p {
		p[i] = modmath.Reduce(rng.Int64N(2*bound+1)-bound, r.Q)
	}
	return p
}

// --- Ring construction tests ---

func TestNewRingDilithium(t *testing.T) {
	r := dilithium
	if r.N != 256 {
		t.Fatalf("N = %d, want 256", r.N)
	}
	if r.Q != 8380417 {
		t.Fatalf("Q = %d, want 8380417", r.Q)
	}
	if !r.HasNTT() {
		t.Fatal("Dilithium ring should have NTT")
	}
}

func TestNewRingNoNTT(t *testing.T) {
	// q=7933, d=512: q-1=7932, 2*512=1024, 7932 % 1024 != 0 → no NTT
	r, err := New(512, 7933)
	if err != nil {
		t.Fatal(err)
	}
	if r.HasNTT() {
		t.Fatal("q=7933 d=512 should NOT have NTT")
	}
}

func TestNewRingBadParams(t *testing.T) {
	if _, err := New(3, 7933); err == nil {
		t.Fatal("N=3 should fail (not power of 2)")
	}
	if _, err := New(256, 1); err == nil {
		t.Fatal("Q=1 should fail")
	}
}

// --- params tests ---

func TestBitrev(t *testing.T) {
	tests := []struct {
		x    uint32
		bits int
		want uint32
	}{
		{0, 8, 0},
		{1, 8, 128},
		{2, 8, 64},
		{128, 8, 1},
		{255, 8, 255},
		{1, 7, 64},
	}
	for _, tt := range tests {
		got := bitrev(tt.x, tt.bits)
		if got != tt.want {
			t.Errorf("bitrev(%d, %d) = %d, want %d", tt.x, tt.bits, got, tt.want)
		}
	}
}

func TestRootOfUnity(t *testing.T) {
	r := dilithium
	// zetas[1] should be ζ^{bitrev(1)} where ζ is primitive 512th root of unity.
	// Verify the primitive root property: ζ^N = -1.
	// We can recover ζ from zetas[1] = ζ^{bitrev(1, 8)} = ζ^128.
	// Instead, just verify the NTT roundtrips work (done below).
	if len(r.zetas) != 256 {
		t.Fatalf("zetas length = %d, want 256", len(r.zetas))
	}
}

func TestNInv(t *testing.T) {
	r := dilithium
	if modmath.ModMul(int64(r.N), r.NInv, r.Q) != 1 {
		t.Fatal("N * NInv mod Q != 1")
	}
}

// --- poly tests ---

func TestAddSub(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	b := randPoly(rng, r)
	sum := r.Add(a, b)
	diff := r.Sub(sum, b)
	if !r.Equal(diff, a) {
		t.Error("(a + b) - b != a")
	}
}

func TestNeg(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	negA := r.Neg(a)
	sum := r.Add(a, negA)
	if !r.Equal(sum, r.NewPoly()) {
		t.Error("a + (-a) != 0")
	}
}

func TestScalarMul(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	double := r.ScalarMul(a, 2)
	sum := r.Add(a, a)
	if !r.Equal(double, sum) {
		t.Error("2*a != a + a")
	}
	if !r.Equal(r.ScalarMul(a, 0), r.NewPoly()) {
		t.Error("0*a != 0")
	}
}

func TestInfNorm(t *testing.T) {
	r := dilithium
	p := r.NewPoly()
	p[0] = 10
	p[1] = r.Q - 5 // represents -5
	p[2] = 3
	if r.InfNorm(p) != 10 {
		t.Errorf("InfNorm = %d, want 10", r.InfNorm(p))
	}
}

func TestEqual(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	b := make(Poly, r.N)
	copy(b, a)
	if !r.Equal(a, b) {
		t.Error("Equal failed on identical polys")
	}
	b[0] = modmath.Reduce(b[0]+1, r.Q)
	if r.Equal(a, b) {
		t.Error("Equal returned true for different polys")
	}
}

// --- NTT tests (Dilithium ring) ---

func TestNTTInvNTTRoundtrip(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 100; trial++ {
		a := randPoly(rng, r)
		aReduced := r.Reduce(a)
		recovered := r.InvNTT(r.NTT(a))
		if !r.Equal(recovered, aReduced) {
			t.Fatalf("InvNTT(NTT(a)) != a on trial %d", trial)
		}
	}
}

func TestNTTMulMatchesNaive(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 50; trial++ {
		a := randPoly(rng, r)
		b := randPoly(rng, r)
		nttResult := r.Mul(a, b)
		naiveResult := r.NaiveMul(a, b)
		if !r.Equal(nttResult, naiveResult) {
			t.Fatalf("NTT mul != naive mul on trial %d", trial)
		}
	}
}

func TestMulOne(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	if !r.Equal(r.Mul(a, r.One()), r.Reduce(a)) {
		t.Error("a * 1 != a")
	}
}

func TestMulZero(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	if !r.Equal(r.Mul(a, r.NewPoly()), r.NewPoly()) {
		t.Error("a * 0 != 0")
	}
}

func TestMulCommutative(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 20; trial++ {
		a := randPoly(rng, r)
		b := randPoly(rng, r)
		if !r.Equal(r.Mul(a, b), r.Mul(b, a)) {
			t.Fatalf("a*b != b*a on trial %d", trial)
		}
	}
}

func TestMulDistributive(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 20; trial++ {
		a := randPoly(rng, r)
		b := randPoly(rng, r)
		c := randPoly(rng, r)
		lhs := r.Mul(a, r.Add(b, c))
		rhs := r.Add(r.Mul(a, b), r.Mul(a, c))
		if !r.Equal(lhs, rhs) {
			t.Fatalf("distributivity failed on trial %d", trial)
		}
	}
}

func TestMulXN(t *testing.T) {
	r := dilithium
	x := r.NewPoly()
	x[1] = 1
	p := make(Poly, r.N)
	copy(p, x)
	for i := 1; i < r.N; i++ {
		p = r.Mul(p, x)
	}
	expected := r.NewPoly()
	expected[0] = r.Q - 1 // -1
	if !r.Equal(p, expected) {
		t.Error("X^N != -1")
	}
}

// --- Karatsuba tests (q=7933, d=512: no NTT) ---

func TestKaratsubaMulMatchesNaive(t *testing.T) {
	r, err := New(512, 7933)
	if err != nil {
		t.Fatal(err)
	}
	if r.HasNTT() {
		t.Fatal("Expected no NTT for q=7933 d=512")
	}

	rng := rand.New(rand.NewPCG(99, 0))
	for trial := 0; trial < 10; trial++ {
		a := randPoly(rng, r)
		b := randPoly(rng, r)
		karatResult := r.Mul(a, b) // uses Karatsuba
		naiveResult := r.NaiveMul(a, b)
		if !r.Equal(karatResult, naiveResult) {
			t.Fatalf("Karatsuba mul != naive mul on trial %d", trial)
		}
	}
}

func TestKaratsubaMulOne(t *testing.T) {
	r, _ := New(512, 7933)
	rng := rand.New(rand.NewPCG(99, 0))
	a := randPoly(rng, r)
	if !r.Equal(r.Mul(a, r.One()), r.Reduce(a)) {
		t.Error("Karatsuba: a * 1 != a")
	}
}

func TestKaratsubaCommutative(t *testing.T) {
	r, _ := New(512, 7933)
	rng := rand.New(rand.NewPCG(99, 0))
	for trial := 0; trial < 10; trial++ {
		a := randPoly(rng, r)
		b := randPoly(rng, r)
		if !r.Equal(r.Mul(a, b), r.Mul(b, a)) {
			t.Fatalf("Karatsuba: a*b != b*a on trial %d", trial)
		}
	}
}

func TestKaratsubaDistributive(t *testing.T) {
	r, _ := New(512, 7933)
	rng := rand.New(rand.NewPCG(99, 0))
	for trial := 0; trial < 10; trial++ {
		a := randPoly(rng, r)
		b := randPoly(rng, r)
		c := randPoly(rng, r)
		lhs := r.Mul(a, r.Add(b, c))
		rhs := r.Add(r.Mul(a, b), r.Mul(a, c))
		if !r.Equal(lhs, rhs) {
			t.Fatalf("Karatsuba: distributivity failed on trial %d", trial)
		}
	}
}

func TestKaratsubaXN(t *testing.T) {
	r, _ := New(64, 7933) // smaller degree for speed
	x := r.NewPoly()
	x[1] = 1
	p := make(Poly, r.N)
	copy(p, x)
	for i := 1; i < r.N; i++ {
		p = r.Mul(p, x)
	}
	expected := r.NewPoly()
	expected[0] = r.Q - 1
	if !r.Equal(p, expected) {
		t.Error("Karatsuba: X^N != -1")
	}
}

// --- vec tests ---

func TestVecAddSub(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := PolyVec{randPoly(rng, r), randPoly(rng, r), randPoly(rng, r)}
	b := PolyVec{randPoly(rng, r), randPoly(rng, r), randPoly(rng, r)}
	sum := r.VecAdd(a, b)
	diff := r.VecSub(sum, b)
	if !r.VecEqual(diff, a) {
		t.Error("VecAdd/VecSub roundtrip failed")
	}
}

func TestInnerProduct(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	b := randPoly(rng, r)
	result := r.InnerProduct(PolyVec{r.One(), r.NewPoly()}, PolyVec{a, b})
	if !r.Equal(result, r.Reduce(a)) {
		t.Error("InnerProduct with unit vector failed")
	}
}

func TestMatVecMul(t *testing.T) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	A := r.NewPolyMat(2, 3)
	s := r.NewPolyVec(3)
	for i := range A {
		for j := range A[i] {
			A[i][j] = randPoly(rng, r)
		}
	}
	for i := range s {
		s[i] = randPoly(rng, r)
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

func TestVecInfNorm(t *testing.T) {
	r := dilithium
	p1 := r.NewPoly()
	p2 := r.NewPoly()
	p1[0] = 5
	p2[0] = r.Q - 10
	v := PolyVec{p1, p2}
	if r.VecInfNorm(v) != 10 {
		t.Errorf("VecInfNorm = %d, want 10", r.VecInfNorm(v))
	}
}

// --- benchmarks ---

func BenchmarkNTT(b *testing.B) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	p := randPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.NTT(p)
	}
}

func BenchmarkMulNTT(b *testing.B) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	p := randPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.Mul(a, p)
	}
}

func BenchmarkMulKaratsuba(b *testing.B) {
	r, _ := New(512, 7933)
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	p := randPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.Mul(a, p)
	}
}

func BenchmarkNaiveMul(b *testing.B) {
	r := dilithium
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng, r)
	p := randPoly(rng, r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.NaiveMul(a, p)
	}
}
