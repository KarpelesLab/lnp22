package ring

import (
	"math/rand/v2"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/modmath"
)

func randPoly(rng *rand.Rand) *Poly {
	var p Poly
	for i := range p {
		p[i] = modmath.Reduce(rng.Int64(), Q)
	}
	return &p
}

func smallPoly(rng *rand.Rand, bound int64) *Poly {
	var p Poly
	for i := range p {
		p[i] = modmath.Reduce(rng.Int64N(2*bound+1)-bound, Q)
	}
	return &p
}

// --- params tests ---

func TestBitrev8(t *testing.T) {
	tests := []struct {
		in, out uint8
	}{
		{0, 0},
		{1, 128},
		{2, 64},
		{128, 1},
		{255, 255},
	}
	for _, tt := range tests {
		got := bitrev8(tt.in)
		if got != tt.out {
			t.Errorf("bitrev8(%d) = %d, want %d", tt.in, got, tt.out)
		}
	}
	// bitrev8 should be an involution
	for i := 0; i < 256; i++ {
		if bitrev8(bitrev8(uint8(i))) != uint8(i) {
			t.Errorf("bitrev8 is not an involution at %d", i)
		}
	}
}

func TestRootOfUnity(t *testing.T) {
	zeta := int64(RootOfUnity)
	// ζ^(2N) should be 1
	if modmath.ModExp(zeta, 2*N, Q) != 1 {
		t.Fatal("ζ^(2N) != 1")
	}
	// ζ^N should be Q-1 (i.e., -1 mod Q)
	if modmath.ModExp(zeta, N, Q) != Q-1 {
		t.Fatal("ζ^N != -1")
	}
	// ζ^k != 1 for 0 < k < 2N (primitive)
	for k := int64(1); k < 2*N; k++ {
		if modmath.ModExp(zeta, k, Q) == 1 {
			t.Fatalf("ζ^%d == 1, so ζ is not primitive", k)
		}
	}
}

func TestNInv(t *testing.T) {
	if modmath.ModMul(N, NInv, Q) != 1 {
		t.Fatalf("N * NInv mod Q = %d, want 1", modmath.ModMul(N, NInv, Q))
	}
}

// --- poly tests ---

func TestAddSub(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	b := randPoly(rng)

	sum := Add(a, b)
	diff := Sub(sum, b)
	if !Equal(diff, a) {
		t.Error("(a + b) - b != a")
	}
}

func TestNeg(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	negA := Neg(a)
	sum := Add(a, negA)
	if !Equal(sum, Zero()) {
		t.Error("a + (-a) != 0")
	}
}

func TestScalarMul(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)

	double := ScalarMul(a, 2)
	sum := Add(a, a)
	if !Equal(double, sum) {
		t.Error("2*a != a + a")
	}

	zero := ScalarMul(a, 0)
	if !Equal(zero, Zero()) {
		t.Error("0*a != 0")
	}
}

func TestInfNorm(t *testing.T) {
	var p Poly
	p[0] = 10
	p[1] = Q - 5 // represents -5
	p[2] = 3
	if InfNorm(&p) != 10 {
		t.Errorf("InfNorm = %d, want 10", InfNorm(&p))
	}
}

func TestEqual(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	b := *a
	if !Equal(a, &b) {
		t.Error("Equal failed on identical polys")
	}
	b[0] = modmath.Reduce(b[0]+1, Q)
	if Equal(a, &b) {
		t.Error("Equal returned true for different polys")
	}
}

// --- NTT tests ---

func TestNTTInvNTTRoundtrip(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 100; trial++ {
		a := randPoly(rng)
		aReduced := Reduce(a)
		recovered := InvNTT(NTT(a))
		if !Equal(recovered, aReduced) {
			t.Fatalf("InvNTT(NTT(a)) != a on trial %d", trial)
		}
	}
}

func TestNTTMulMatchesNaive(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 50; trial++ {
		a := randPoly(rng)
		b := randPoly(rng)
		nttResult := Mul(a, b)
		naiveResult := NaiveMul(a, b)
		if !Equal(nttResult, naiveResult) {
			t.Fatalf("NTT mul != naive mul on trial %d", trial)
		}
	}
}

func TestMulOne(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	one := One()
	result := Mul(a, one)
	if !Equal(result, Reduce(a)) {
		t.Error("a * 1 != a")
	}
}

func TestMulZero(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	result := Mul(a, Zero())
	if !Equal(result, Zero()) {
		t.Error("a * 0 != 0")
	}
}

func TestMulCommutative(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 20; trial++ {
		a := randPoly(rng)
		b := randPoly(rng)
		ab := Mul(a, b)
		ba := Mul(b, a)
		if !Equal(ab, ba) {
			t.Fatalf("a*b != b*a on trial %d", trial)
		}
	}
}

func TestMulDistributive(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	for trial := 0; trial < 20; trial++ {
		a := randPoly(rng)
		b := randPoly(rng)
		c := randPoly(rng)
		// a*(b+c) = a*b + a*c
		lhs := Mul(a, Add(b, c))
		rhs := Add(Mul(a, b), Mul(a, c))
		if !Equal(lhs, rhs) {
			t.Fatalf("distributivity failed on trial %d", trial)
		}
	}
}

func TestMulXN(t *testing.T) {
	// X^N ≡ -1 in R_q, so multiplying by X^N should negate all coefficients.
	// We can test: X * X^(N-1) should give X^N ≡ -1.
	var x Poly
	x[1] = 1 // X
	// X^2 = X * X
	p := &x
	for i := 1; i < N; i++ {
		p = Mul(p, &x)
	}
	// p should be X^N ≡ -1 mod (X^N+1), i.e., coefficients: [-1, 0, 0, ...]
	expected := Zero()
	expected[0] = Q - 1 // -1 mod Q
	if !Equal(p, expected) {
		t.Errorf("X^N != -1; got coeff[0]=%d", p[0])
	}
}

// --- vec tests ---

func TestVecAddSub(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := PolyVec{*randPoly(rng), *randPoly(rng), *randPoly(rng)}
	b := PolyVec{*randPoly(rng), *randPoly(rng), *randPoly(rng)}

	sum := VecAdd(a, b)
	diff := VecSub(sum, b)
	if !VecEqual(diff, a) {
		t.Error("VecAdd/VecSub roundtrip failed")
	}
}

func TestInnerProduct(t *testing.T) {
	// ⟨[1, 0], [a, b]⟩ = a
	rng := rand.New(rand.NewPCG(42, 0))
	a := *randPoly(rng)
	b := *randPoly(rng)
	one := *One()
	zero := Poly{}
	result := InnerProduct(PolyVec{one, zero}, PolyVec{a, b})
	if !Equal(result, Reduce(&a)) {
		t.Error("InnerProduct with unit vector failed")
	}
}

func TestMatVecMul(t *testing.T) {
	rng := rand.New(rand.NewPCG(42, 0))
	// 2×3 matrix times 3-vector = 2-vector
	A := NewPolyMat(2, 3)
	s := NewPolyVec(3)
	for i := range A {
		for j := range A[i] {
			A[i][j] = *randPoly(rng)
		}
	}
	for i := range s {
		s[i] = *randPoly(rng)
	}

	result := MatVecMul(A, s)
	if len(result) != 2 {
		t.Fatalf("MatVecMul result length = %d, want 2", len(result))
	}

	// Verify manually: result[0] = A[0][0]*s[0] + A[0][1]*s[1] + A[0][2]*s[2]
	expected := Add(Add(Mul(&A[0][0], &s[0]), Mul(&A[0][1], &s[1])), Mul(&A[0][2], &s[2]))
	if !Equal(&result[0], expected) {
		t.Error("MatVecMul result[0] mismatch")
	}
}

func TestVecInfNorm(t *testing.T) {
	var p1, p2 Poly
	p1[0] = 5
	p2[0] = Q - 10 // represents -10
	v := PolyVec{p1, p2}
	if VecInfNorm(v) != 10 {
		t.Errorf("VecInfNorm = %d, want 10", VecInfNorm(v))
	}
}

// --- benchmarks ---

func BenchmarkNTT(b *testing.B) {
	rng := rand.New(rand.NewPCG(42, 0))
	p := randPoly(rng)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		NTT(p)
	}
}

func BenchmarkMul(b *testing.B) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	p := randPoly(rng)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Mul(a, p)
	}
}

func BenchmarkNaiveMul(b *testing.B) {
	rng := rand.New(rand.NewPCG(42, 0))
	a := randPoly(rng)
	p := randPoly(rng)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		NaiveMul(a, p)
	}
}
