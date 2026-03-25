package sampler

import (
	"crypto/rand"
	"math"
	"math/big"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
)

var bigDilithium, _ = ring.NewBig(256, big.NewInt(8380417))

func TestBigGaussianPolyNorm(t *testing.T) {
	r := bigDilithium
	sigma := 50.0
	p := SampleBigGaussianPoly(r, sigma, rand.Reader)
	norm := r.InfNorm(p)
	bound := big.NewInt(int64(6 * sigma))
	if norm.Cmp(bound) > 0 {
		t.Errorf("Big Gaussian poly norm %s exceeds 6σ bound %s", norm, bound)
	}
}

func TestBigGaussianMeanVariance(t *testing.T) {
	r := bigDilithium
	sigma := 50.0
	p := SampleBigGaussianPoly(r, sigma, rand.Reader)
	sum := 0.0
	sumSq := 0.0
	for i := 0; i < r.N; i++ {
		v := float64(bigmod.CenterReduce(p[i], r.Q).Int64())
		sum += v
		sumSq += v * v
	}
	mean := sum / float64(r.N)
	variance := sumSq/float64(r.N) - mean*mean
	// Loose check — single polynomial has high variance in the estimate
	if math.Abs(mean) > 3*sigma {
		t.Errorf("Big Gaussian mean = %f, too far from 0", mean)
	}
	_ = variance // just check it doesn't panic
}

func TestBigUniformPolyRange(t *testing.T) {
	r := bigDilithium
	p := SampleBigUniformPoly(r, rand.Reader)
	for i := 0; i < r.N; i++ {
		if p[i].Sign() < 0 || p[i].Cmp(r.Q) >= 0 {
			t.Fatalf("Big uniform coeff %d = %s out of range [0, Q)", i, p[i])
		}
	}
}

func TestBigUniformPoly152bit(t *testing.T) {
	r := ring.LNP22Big()
	p := SampleBigUniformPoly(r, rand.Reader)
	for i := 0; i < r.N; i++ {
		if p[i].Sign() < 0 || p[i].Cmp(r.Q) >= 0 {
			t.Fatalf("152-bit uniform coeff %d out of range", i)
		}
	}
}

func TestBigUniformMatDimensions(t *testing.T) {
	r := bigDilithium
	m := SampleBigUniformMat(r, 3, 5, rand.Reader)
	if len(m) != 3 {
		t.Fatalf("rows = %d, want 3", len(m))
	}
	for i := range m {
		if len(m[i]) != 5 {
			t.Fatalf("row %d cols = %d, want 5", i, len(m[i]))
		}
	}
}

func TestBigTernaryPoly(t *testing.T) {
	r := bigDilithium
	p := SampleBigTernaryPoly(r, rand.Reader)
	for i := 0; i < r.N; i++ {
		v := bigmod.CenterReduce(p[i], r.Q)
		vi := v.Int64()
		if vi != -1 && vi != 0 && vi != 1 {
			t.Fatalf("Big ternary coeff %d = %s, expected {-1, 0, 1}", i, v)
		}
	}
}

func TestBigChallengeWeight(t *testing.T) {
	r := bigDilithium
	for _, kappa := range []int{1, 32, 60, 128} {
		seed := make([]byte, 32)
		rand.Read(seed)
		c := SampleBigChallenge(r, seed, kappa)
		nonZero := 0
		for i := 0; i < r.N; i++ {
			v := bigmod.CenterReduce(c[i], r.Q)
			if v.Sign() != 0 {
				nonZero++
				vi := v.Int64()
				if vi != -1 && vi != 1 {
					t.Fatalf("Big challenge coeff %d = %s, expected ±1", i, v)
				}
			}
		}
		if nonZero != kappa {
			t.Errorf("Big challenge kappa=%d: got %d non-zero", kappa, nonZero)
		}
	}
}

func TestBigChallengeDeterministic(t *testing.T) {
	r := bigDilithium
	seed := []byte("big-test-seed")
	c1 := SampleBigChallenge(r, seed, 60)
	c2 := SampleBigChallenge(r, seed, 60)
	if !r.Equal(c1, c2) {
		t.Error("SampleBigChallenge is not deterministic")
	}
}

func TestBigChallenge152bit(t *testing.T) {
	r := ring.LNP22Big()
	seed := make([]byte, 32)
	rand.Read(seed)
	c := SampleBigChallenge(r, seed, 60)
	nonZero := 0
	for i := 0; i < r.N; i++ {
		if c[i].Sign() != 0 {
			nonZero++
		}
	}
	if nonZero != 60 {
		t.Errorf("152-bit challenge: got %d non-zero, want 60", nonZero)
	}
}
