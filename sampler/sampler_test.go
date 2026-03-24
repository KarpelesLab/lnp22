package sampler

import (
	"crypto/rand"
	"math"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/modmath"
	"github.com/KarpelesLab/lnp22/ring"
)

var dilithium = ring.Dilithium()

func TestGaussianMeanVariance(t *testing.T) {
	sigma := 50.0
	numSamples := 100000
	cdt := newGaussianCDT(sigma)
	sum := 0.0
	sumSq := 0.0
	for i := 0; i < numSamples; i++ {
		v := float64(cdt.sample(rand.Reader))
		sum += v
		sumSq += v * v
	}
	mean := sum / float64(numSamples)
	variance := sumSq/float64(numSamples) - mean*mean
	if math.Abs(mean) > 3*sigma/math.Sqrt(float64(numSamples)) {
		t.Errorf("Gaussian mean = %f, expected ≈ 0", mean)
	}
	expectedVar := sigma * sigma
	if math.Abs(variance-expectedVar)/expectedVar > 0.05 {
		t.Errorf("Gaussian variance = %f, expected ≈ %f", variance, expectedVar)
	}
}

func TestGaussianTailBound(t *testing.T) {
	sigma := 50.0
	cdt := newGaussianCDT(sigma)
	tailBound := int64(6 * sigma)
	for i := 0; i < 100000; i++ {
		v := cdt.sample(rand.Reader)
		if v > tailBound || v < -tailBound {
			t.Fatalf("Gaussian sample %d exceeds 6σ bound %d", v, tailBound)
		}
	}
}

func TestGaussianPolyNorm(t *testing.T) {
	r := dilithium
	sigma := 50.0
	p := SampleGaussianPoly(r, sigma, rand.Reader)
	norm := r.InfNorm(p)
	bound := int64(6 * sigma)
	if norm > bound {
		t.Errorf("Gaussian poly infinity norm %d exceeds expected bound %d", norm, bound)
	}
}

func TestUniformPolyRange(t *testing.T) {
	r := dilithium
	p := SampleUniformPoly(r, rand.Reader)
	for i := 0; i < r.N; i++ {
		if p[i] < 0 || p[i] >= r.Q {
			t.Fatalf("Uniform coefficient %d = %d out of range [0, Q)", i, p[i])
		}
	}
}

func TestUniformPolySmallQ(t *testing.T) {
	r, _ := ring.New(64, 7933)
	p := SampleUniformPoly(r, rand.Reader)
	for i := 0; i < r.N; i++ {
		if p[i] < 0 || p[i] >= r.Q {
			t.Fatalf("Uniform coefficient %d = %d out of range [0, Q=%d)", i, p[i], r.Q)
		}
	}
}

func TestUniformMatDimensions(t *testing.T) {
	r := dilithium
	m := SampleUniformMat(r, 3, 5, rand.Reader)
	if len(m) != 3 {
		t.Fatalf("Matrix rows = %d, want 3", len(m))
	}
	for i := range m {
		if len(m[i]) != 5 {
			t.Fatalf("Matrix row %d cols = %d, want 5", i, len(m[i]))
		}
	}
}

func TestTernaryPoly(t *testing.T) {
	r := dilithium
	p := SampleTernaryPoly(r, rand.Reader)
	for i := 0; i < r.N; i++ {
		v := modmath.CenterReduce(p[i], r.Q)
		if v != -1 && v != 0 && v != 1 {
			t.Fatalf("Ternary coefficient %d = %d, expected {-1, 0, 1}", i, v)
		}
	}
}

func TestTernaryDistribution(t *testing.T) {
	r := dilithium
	counts := [3]int{}
	numPolys := 100
	for k := 0; k < numPolys; k++ {
		p := SampleTernaryPoly(r, rand.Reader)
		for i := 0; i < r.N; i++ {
			v := modmath.CenterReduce(p[i], r.Q)
			counts[v+1]++
		}
	}
	total := numPolys * r.N
	for i, c := range counts {
		frac := float64(c) / float64(total)
		if math.Abs(frac-1.0/3.0) > 0.02 {
			t.Errorf("Ternary value %d: fraction = %f, expected ≈ 1/3", i-1, frac)
		}
	}
}

func TestChallengeWeight(t *testing.T) {
	r := dilithium
	for _, kappa := range []int{1, 32, 60, 128, 256} {
		seed := make([]byte, 32)
		rand.Read(seed)
		c := SampleChallenge(r, seed, kappa)
		nonZero := 0
		for i := 0; i < r.N; i++ {
			v := modmath.CenterReduce(c[i], r.Q)
			if v != 0 {
				nonZero++
				if v != -1 && v != 1 {
					t.Fatalf("Challenge coeff %d = %d, expected {-1, 0, 1}", i, v)
				}
			}
		}
		if nonZero != kappa {
			t.Errorf("Challenge kappa=%d: got %d non-zero coefficients", kappa, nonZero)
		}
	}
}

func TestChallengeWeightSmallRing(t *testing.T) {
	r, _ := ring.New(64, 7933)
	seed := make([]byte, 32)
	rand.Read(seed)
	c := SampleChallenge(r, seed, 20)
	nonZero := 0
	for i := 0; i < r.N; i++ {
		v := modmath.CenterReduce(c[i], r.Q)
		if v != 0 {
			nonZero++
		}
	}
	if nonZero != 20 {
		t.Errorf("Challenge kappa=20 on small ring: got %d non-zero", nonZero)
	}
}

func TestChallengeDeterministic(t *testing.T) {
	r := dilithium
	seed := []byte("test-seed-for-deterministic-challenge")
	c1 := SampleChallenge(r, seed, 60)
	c2 := SampleChallenge(r, seed, 60)
	if !r.Equal(c1, c2) {
		t.Error("SampleChallenge is not deterministic with same seed")
	}
}

func TestChallengeDifferentSeeds(t *testing.T) {
	r := dilithium
	c1 := SampleChallenge(r, []byte("seed-a"), 60)
	c2 := SampleChallenge(r, []byte("seed-b"), 60)
	if r.Equal(c1, c2) {
		t.Error("Different seeds produced the same challenge")
	}
}

func BenchmarkSampleGaussianPoly(b *testing.B) {
	r := dilithium
	for i := 0; i < b.N; i++ {
		SampleGaussianPoly(r, 350, rand.Reader)
	}
}

func BenchmarkSampleUniformPoly(b *testing.B) {
	r := dilithium
	for i := 0; i < b.N; i++ {
		SampleUniformPoly(r, rand.Reader)
	}
}

func BenchmarkSampleChallenge(b *testing.B) {
	r := dilithium
	seed := make([]byte, 32)
	rand.Read(seed)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		SampleChallenge(r, seed, 60)
	}
}
