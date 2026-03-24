package sampler

import (
	"crypto/rand"
	"math"
	"testing"

	"github.com/KarpelesLab/lnp22/internal/modmath"
	"github.com/KarpelesLab/lnp22/ring"
)

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

	// Mean should be close to 0
	if math.Abs(mean) > 3*sigma/math.Sqrt(float64(numSamples)) {
		t.Errorf("Gaussian mean = %f, expected ≈ 0 (3σ/√n = %f)", mean, 3*sigma/math.Sqrt(float64(numSamples)))
	}

	// Variance should be close to σ²
	expectedVar := sigma * sigma
	relError := math.Abs(variance-expectedVar) / expectedVar
	if relError > 0.05 {
		t.Errorf("Gaussian variance = %f, expected ≈ %f (rel error = %.2f%%)", variance, expectedVar, relError*100)
	}
}

func TestGaussianTailBound(t *testing.T) {
	sigma := 50.0
	cdt := newGaussianCDT(sigma)
	tailBound := int64(6 * sigma) // 6σ tail — should almost never exceed

	for i := 0; i < 100000; i++ {
		v := cdt.sample(rand.Reader)
		if v > tailBound || v < -tailBound {
			t.Fatalf("Gaussian sample %d exceeds 6σ bound %d", v, tailBound)
		}
	}
}

func TestGaussianPolyNorm(t *testing.T) {
	sigma := 50.0
	p := SampleGaussianPoly(sigma, rand.Reader)

	// The infinity norm should be bounded with high probability
	norm := ring.InfNorm(p)
	bound := int64(6 * sigma)
	if norm > bound {
		t.Errorf("Gaussian poly infinity norm %d exceeds expected bound %d", norm, bound)
	}
}

func TestUniformPolyRange(t *testing.T) {
	p := SampleUniformPoly(rand.Reader)
	for i := 0; i < ring.N; i++ {
		if p[i] < 0 || p[i] >= ring.Q {
			t.Fatalf("Uniform coefficient %d = %d out of range [0, Q)", i, p[i])
		}
	}
}

func TestUniformMatDimensions(t *testing.T) {
	m := SampleUniformMat(3, 5, rand.Reader)
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
	p := SampleTernaryPoly(rand.Reader)
	for i := 0; i < ring.N; i++ {
		v := modmath.CenterReduce(p[i], ring.Q)
		if v != -1 && v != 0 && v != 1 {
			t.Fatalf("Ternary coefficient %d = %d, expected {-1, 0, 1}", i, v)
		}
	}
}

func TestTernaryDistribution(t *testing.T) {
	counts := [3]int{} // -1, 0, 1
	numPolys := 100
	for k := 0; k < numPolys; k++ {
		p := SampleTernaryPoly(rand.Reader)
		for i := 0; i < ring.N; i++ {
			v := modmath.CenterReduce(p[i], ring.Q)
			counts[v+1]++
		}
	}
	total := numPolys * ring.N
	for i, c := range counts {
		frac := float64(c) / float64(total)
		if math.Abs(frac-1.0/3.0) > 0.02 {
			t.Errorf("Ternary value %d: fraction = %f, expected ≈ 1/3", i-1, frac)
		}
	}
}

func TestChallengeWeight(t *testing.T) {
	for _, kappa := range []int{1, 32, 60, 128, 256} {
		seed := make([]byte, 32)
		rand.Read(seed)
		c := SampleChallenge(seed, kappa)

		nonZero := 0
		for i := 0; i < ring.N; i++ {
			v := modmath.CenterReduce(c[i], ring.Q)
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

func TestChallengeDeterministic(t *testing.T) {
	seed := []byte("test-seed-for-deterministic-challenge")
	c1 := SampleChallenge(seed, 60)
	c2 := SampleChallenge(seed, 60)
	if !ring.Equal(c1, c2) {
		t.Error("SampleChallenge is not deterministic with same seed")
	}
}

func TestChallengeDifferentSeeds(t *testing.T) {
	c1 := SampleChallenge([]byte("seed-a"), 60)
	c2 := SampleChallenge([]byte("seed-b"), 60)
	if ring.Equal(c1, c2) {
		t.Error("Different seeds produced the same challenge")
	}
}

func BenchmarkSampleGaussianPoly(b *testing.B) {
	for i := 0; i < b.N; i++ {
		SampleGaussianPoly(350, rand.Reader)
	}
}

func BenchmarkSampleUniformPoly(b *testing.B) {
	for i := 0; i < b.N; i++ {
		SampleUniformPoly(rand.Reader)
	}
}

func BenchmarkSampleChallenge(b *testing.B) {
	seed := make([]byte, 32)
	rand.Read(seed)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		SampleChallenge(seed, 60)
	}
}
