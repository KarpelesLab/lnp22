package trapgen

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
)

var testRing, _ = ring.NewBig(64, big.NewInt(8380417))

func TestGadgetVec(t *testing.T) {
	q := big.NewInt(8380417)
	g := GadgetVec(q, 2)
	// g[0] = 1, g[1] = 2, g[2] = 4, ...
	for i, v := range g {
		expected := new(big.Int).Lsh(big.NewInt(1), uint(i))
		expected.Mod(expected, q)
		if v.Cmp(expected) != 0 {
			t.Errorf("GadgetVec[%d] = %s, want %s", i, v, expected)
		}
	}
}

func TestGadgetLen(t *testing.T) {
	// log2(8380417) ≈ 23, so GadgetLen should be 23
	q := big.NewInt(8380417)
	k := GadgetLen(q, 2)
	if k != 23 {
		t.Errorf("GadgetLen(8380417, 2) = %d, want 23", k)
	}
}

func TestDecomposeRoundtrip(t *testing.T) {
	q := big.NewInt(8380417)
	g := GadgetVec(q, 2)

	for _, x := range []*big.Int{big.NewInt(0), big.NewInt(1), big.NewInt(42), big.NewInt(8380416)} {
		v := Decompose(x, q, 2)
		// Recompose: Σ g[i] * v[i] mod q
		sum := new(big.Int)
		for i := range v {
			term := new(big.Int).Mul(g[i], v[i])
			sum.Add(sum, term)
		}
		sum.Mod(sum, q)
		xr := bigmod.Reduce(x, q)
		if sum.Cmp(xr) != 0 {
			t.Errorf("Decompose roundtrip failed for x=%s: got %s", x, sum)
		}
	}
}

func TestDecomposeSmallDigits(t *testing.T) {
	q := big.NewInt(8380417)
	x := big.NewInt(123456)
	v := Decompose(x, q, 2)
	// All digits should be in {-1, 0, 1} for base 2
	for i, d := range v {
		if d.CmpAbs(big.NewInt(1)) > 0 {
			t.Errorf("Decompose digit[%d] = %s, expected |d| ≤ 1 for base 2", i, d)
		}
	}
}

func TestDecomposePolyRecompose(t *testing.T) {
	r := testRing
	// Create a random-ish polynomial
	p := r.NewPoly()
	for i := 0; i < r.N; i++ {
		p[i] = big.NewInt(int64(i * 17 % 1000))
	}
	p = r.Reduce(p)

	decomp := DecomposePoly(r, p, 2)
	recomp := RecomposePoly(r, decomp, 2)
	if !r.Equal(recomp, p) {
		t.Error("DecomposePoly/RecomposePoly roundtrip failed")
	}
}

func TestPolyInverse(t *testing.T) {
	r := testRing
	// Create a simple polynomial and check its inverse
	one := r.One()

	// Try to invert 1 (trivially invertible)
	inv, err := polyInverse(r, one)
	if err != nil {
		t.Fatalf("polyInverse(1) failed: %v", err)
	}
	check := r.Mul(one, inv)
	if !r.Equal(check, one) {
		t.Error("1 * 1^{-1} != 1")
	}
}

func TestPolyInverseRandom(t *testing.T) {
	r := testRing
	// A Gaussian polynomial with small coefficients is likely invertible
	for trial := 0; trial < 10; trial++ {
		f := ring.BigPoly(r.NewPoly())
		// Small random coefficients
		for i := 0; i < r.N; i++ {
			var buf [1]byte
			rand.Read(buf[:])
			f[i] = big.NewInt(int64(buf[0]%5) - 2)
			f[i].Mod(f[i], r.Q)
		}
		// Ensure f[0] != 0
		f[0] = big.NewInt(1)

		inv, err := polyInverse(r, f)
		if err != nil {
			continue // ok, some polynomials aren't invertible
		}

		check := r.Mul(f, inv)
		one := r.One()
		if !r.Equal(check, one) {
			t.Fatalf("f * f^{-1} != 1 on trial %d", trial)
		}
	}
}

func TestGenerateNTRUKey(t *testing.T) {
	r := testRing
	params := DefaultParams(r)

	key, err := GenerateNTRUKey(params, rand.Reader)
	if err != nil {
		t.Fatalf("GenerateNTRUKey failed: %v", err)
	}

	// Verify h = f^{-1} * g: f * h = g
	fh := r.Mul(key.F, key.H)
	if !r.Equal(fh, r.Reduce(key.G)) {
		t.Error("f * h != g (NTRU relation violated)")
	}

	// Check that f and g are short
	fNorm := r.InfNorm(key.F)
	gNorm := r.InfNorm(key.G)
	bound := big.NewInt(100) // Gaussian with sigma=1.17 should be very small
	if fNorm.Cmp(bound) > 0 {
		t.Errorf("f norm %s exceeds expected bound", fNorm)
	}
	if gNorm.Cmp(bound) > 0 {
		t.Errorf("g norm %s exceeds expected bound", gNorm)
	}
}

func TestTrapGen(t *testing.T) {
	r := testRing
	params := DefaultParams(r)

	pk, td, err := TrapGen(params, rand.Reader)
	if err != nil {
		t.Fatalf("TrapGen failed: %v", err)
	}

	if pk == nil || td == nil {
		t.Fatal("TrapGen returned nil")
	}

	// Check public key dimensions
	k := GadgetLen(r.Q, params.Base)
	if len(pk.A) != k+1 {
		t.Errorf("Public key A length = %d, want %d", len(pk.A), k+1)
	}

	// Check trapdoor dimensions
	if len(td.R) != k+1 {
		t.Errorf("Trapdoor R length = %d, want %d", len(td.R), k+1)
	}

	// Verify trapdoor polynomials are short
	for i := range td.R {
		norm := r.InfNorm(td.R[i])
		if norm.Cmp(big.NewInt(100)) > 0 {
			t.Errorf("Trapdoor R[%d] norm %s exceeds expected bound", i, norm)
		}
	}
}

func TestSolveNTRUBase(t *testing.T) {
	q := big.NewInt(101)
	f := big.NewInt(3)
	g := big.NewInt(7)

	F, G, err := solveNTRUBase(f, g, q)
	if err != nil {
		t.Fatalf("solveNTRUBase failed: %v", err)
	}

	// Check: f*G[0] - g*F[0] = q
	fG := new(big.Int).Mul(f, G[0])
	gF := new(big.Int).Mul(g, F[0])
	diff := new(big.Int).Sub(fG, gF)
	if diff.Cmp(q) != 0 {
		t.Errorf("fG - gF = %s, want %s", diff, q)
	}
}

func BenchmarkGenerateNTRUKey(b *testing.B) {
	r := testRing
	params := DefaultParams(r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		GenerateNTRUKey(params, rand.Reader)
	}
}

func BenchmarkTrapGen(b *testing.B) {
	r := testRing
	params := DefaultParams(r)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		TrapGen(params, rand.Reader)
	}
}
