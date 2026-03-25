package bigmod

import (
	"math/big"
	"testing"
)

var (
	smallQ = big.NewInt(8380417)
	// 152-bit prime with q ≡ 1 (mod 2048), enabling NTT for N=1024
	bigQ, _ = new(big.Int).SetString("5708990770823839524233143877797980545530982401", 10)
)

func TestReduce(t *testing.T) {
	q := smallQ
	tests := []struct {
		a    int64
		want int64
	}{
		{0, 0}, {1, 1}, {-1, 8380416}, {8380417, 0}, {8380418, 1},
	}
	for _, tt := range tests {
		got := Reduce(big.NewInt(tt.a), q)
		if got.Cmp(big.NewInt(tt.want)) != 0 {
			t.Errorf("Reduce(%d) = %s, want %d", tt.a, got, tt.want)
		}
	}
}

func TestCenterReduce(t *testing.T) {
	q := smallQ
	half := new(big.Int).Rsh(q, 1) // 4190208

	tests := []struct {
		a    int64
		want int64
	}{
		{0, 0}, {1, 1}, {-1, -1},
	}
	for _, tt := range tests {
		got := CenterReduce(big.NewInt(tt.a), q)
		if got.Cmp(big.NewInt(tt.want)) != 0 {
			t.Errorf("CenterReduce(%d) = %s, want %d", tt.a, got, tt.want)
		}
	}

	// half+1 should map to negative
	got := CenterReduce(new(big.Int).Add(half, bigOne), q)
	if got.Sign() >= 0 {
		t.Errorf("CenterReduce(half+1) = %s, expected negative", got)
	}
}

func TestModMul(t *testing.T) {
	q := smallQ
	// (-1) * (-1) = 1
	a := new(big.Int).Sub(q, bigOne)
	got := ModMul(a, a, q)
	if got.Cmp(bigOne) != 0 {
		t.Errorf("(-1)*(-1) mod q = %s, want 1", got)
	}
}

func TestModExp(t *testing.T) {
	q := smallQ
	// Fermat's little theorem: a^(q-1) ≡ 1 (mod q) for prime q
	qm1 := new(big.Int).Sub(q, bigOne)
	for _, base := range []int64{2, 3, 7, 1000} {
		got := ModExp(big.NewInt(base), qm1, q)
		if got.Cmp(bigOne) != 0 {
			t.Errorf("%d^(q-1) mod q = %s, want 1", base, got)
		}
	}
}

func TestModInverse(t *testing.T) {
	q := smallQ
	for _, a := range []int64{1, 2, 3, 7, 1000, 8380416} {
		inv := ModInverse(big.NewInt(a), q)
		prod := ModMul(big.NewInt(a), inv, q)
		if prod.Cmp(bigOne) != 0 {
			t.Errorf("ModInverse(%d): %d * %s mod q = %s, want 1", a, a, inv, prod)
		}
	}
}

func TestModInversePanics(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Error("ModInverse(0, q) should panic")
		}
	}()
	ModInverse(bigZero, smallQ)
}

// --- 152-bit tests ---

func TestBigQOperations(t *testing.T) {
	q := bigQ
	if q == nil {
		t.Fatal("Failed to parse bigQ")
	}

	// Basic reduce
	neg := new(big.Int).Neg(bigOne)
	r := Reduce(neg, q)
	expected := new(big.Int).Sub(q, bigOne)
	if r.Cmp(expected) != 0 {
		t.Errorf("Reduce(-1, bigQ) = %s, want %s", r, expected)
	}

	// Fermat's little theorem with big q
	qm1 := new(big.Int).Sub(q, bigOne)
	got := ModExp(big.NewInt(2), qm1, q)
	if got.Cmp(bigOne) != 0 {
		t.Error("2^(bigQ-1) mod bigQ != 1 — bigQ may not be prime")
	}

	// Inverse
	inv := ModInverse(big.NewInt(42), q)
	prod := ModMul(big.NewInt(42), inv, q)
	if prod.Cmp(bigOne) != 0 {
		t.Errorf("42 * 42^{-1} mod bigQ = %s, want 1", prod)
	}
}

func TestBigQCenterReduce(t *testing.T) {
	q := bigQ
	half := new(big.Int).Rsh(q, 1)

	// q-1 should center-reduce to -1
	qm1 := new(big.Int).Sub(q, bigOne)
	got := CenterReduce(qm1, q)
	if got.Cmp(big.NewInt(-1)) != 0 {
		t.Errorf("CenterReduce(q-1, bigQ) = %s, want -1", got)
	}

	// half should stay positive
	got = CenterReduce(half, q)
	if got.Sign() < 0 {
		t.Errorf("CenterReduce(half) = %s, expected non-negative", got)
	}
}

// --- Rounding tests ---

func TestRoundPBasic(t *testing.T) {
	q := big.NewInt(100)
	p := big.NewInt(10)

	tests := []struct {
		x    int64
		want int64
	}{
		{0, 0},
		{10, 1},
		{50, 5},
		{94, 9}, // round(10*94/100) = floor((940+50)/100) = 9
	}
	for _, tt := range tests {
		got := RoundP(big.NewInt(tt.x), p, q)
		if got.Cmp(big.NewInt(tt.want)) != 0 {
			t.Errorf("RoundP(%d, 10, 100) = %s, want %d", tt.x, got, tt.want)
		}
	}
}

func TestRoundPBigQ(t *testing.T) {
	q := bigQ
	p := big.NewInt(257)

	// RoundP(0) = 0
	got := RoundP(bigZero, p, q)
	if got.Sign() != 0 {
		t.Errorf("RoundP(0) = %s, want 0", got)
	}

	// RoundP(q-1) wraps to 0 with nearest-integer rounding (since q-1 ≈ q → p*q/q ≈ p → mod p = 0)
	qm1 := new(big.Int).Sub(q, bigOne)
	got = RoundP(qm1, p, q)
	if got.Sign() != 0 {
		t.Errorf("RoundP(q-1, 257, bigQ) = %s, want 0", got)
	}

	// Test a value that should decode to p-1=256: x ≈ q*(p-1)/p
	target := ScaleUp(new(big.Int).Sub(p, bigOne), p, q)
	got = RoundP(target, p, q)
	pm1 := new(big.Int).Sub(p, bigOne)
	if got.Cmp(pm1) != 0 {
		t.Errorf("RoundP(ScaleUp(p-1)) = %s, want %s", got, pm1)
	}
}

func TestScaleUpRoundTrip(t *testing.T) {
	q := big.NewInt(1000003)
	p := big.NewInt(257)

	// For m in [0, p), ScaleUp then RoundP should recover m (approximately).
	for m := int64(0); m < 257; m++ {
		scaled := ScaleUp(big.NewInt(m), p, q)
		recovered := RoundP(scaled, p, q)
		if recovered.Cmp(big.NewInt(m)) != 0 {
			t.Errorf("ScaleUp/RoundP roundtrip: m=%d, recovered=%s", m, recovered)
		}
	}
}

func TestScaleUpRoundTripBigQ(t *testing.T) {
	q := bigQ
	p := big.NewInt(257)

	for m := int64(0); m < 257; m++ {
		scaled := ScaleUp(big.NewInt(m), p, q)
		recovered := RoundP(scaled, p, q)
		if recovered.Cmp(big.NewInt(m)) != 0 {
			t.Errorf("Big roundtrip: m=%d, recovered=%s", m, recovered)
		}
	}
}

// --- In-place variant tests ---

func TestInPlaceVariants(t *testing.T) {
	q := smallQ
	a := big.NewInt(1234567)
	b := big.NewInt(7654321)

	dst := new(big.Int)

	ModAddInPlace(dst, a, b, q)
	expected := ModAdd(a, b, q)
	if dst.Cmp(expected) != 0 {
		t.Errorf("ModAddInPlace: %s != %s", dst, expected)
	}

	ModSubInPlace(dst, a, b, q)
	expected = ModSub(a, b, q)
	if dst.Cmp(expected) != 0 {
		t.Errorf("ModSubInPlace: %s != %s", dst, expected)
	}

	ModMulInPlace(dst, a, b, q)
	expected = ModMul(a, b, q)
	if dst.Cmp(expected) != 0 {
		t.Errorf("ModMulInPlace: %s != %s", dst, expected)
	}
}

func BenchmarkModMulSmallQ(b *testing.B) {
	q := smallQ
	x := big.NewInt(4190208)
	y := big.NewInt(7654321)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ModMul(x, y, q)
	}
}

func BenchmarkModMulBigQ(b *testing.B) {
	q := bigQ
	x := new(big.Int).Rsh(q, 1)
	y := new(big.Int).Sub(q, big.NewInt(42))
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ModMul(x, y, q)
	}
}

func BenchmarkModMulInPlaceBigQ(b *testing.B) {
	q := bigQ
	x := new(big.Int).Rsh(q, 1)
	y := new(big.Int).Sub(q, big.NewInt(42))
	dst := new(big.Int)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ModMulInPlace(dst, x, y, q)
	}
}
