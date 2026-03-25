package preimage

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
	"github.com/KarpelesLab/lnp22/trapgen"
)

var testRing, _ = ring.NewBig(64, big.NewInt(8380417))

func TestSamplePreDimensions(t *testing.T) {
	r := testRing
	params := trapgen.DefaultParams(r)

	pk, td, err := trapgen.TrapGen(params, rand.Reader)
	if err != nil {
		t.Fatalf("TrapGen failed: %v", err)
	}

	target := sampler.SampleBigUniformPoly(r, rand.Reader)
	s, err := SamplePre(r, pk, td, target, 350, rand.Reader)
	if err != nil {
		t.Fatalf("SamplePre failed: %v", err)
	}

	if len(s) != len(pk.A) {
		t.Errorf("SamplePre output length = %d, want %d", len(s), len(pk.A))
	}
}

func TestSamplePreShortness(t *testing.T) {
	r := testRing
	params := trapgen.DefaultParams(r)

	pk, td, err := trapgen.TrapGen(params, rand.Reader)
	if err != nil {
		t.Fatalf("TrapGen failed: %v", err)
	}

	target := sampler.SampleBigUniformPoly(r, rand.Reader)
	sigma := 350.0
	s, err := SamplePre(r, pk, td, target, sigma, rand.Reader)
	if err != nil {
		t.Fatalf("SamplePre failed: %v", err)
	}

	// Check output is short: norm should be bounded by a multiple of sigma
	norm := r.VecInfNorm(s)
	bound := big.NewInt(int64(6 * sigma * float64(len(s))))
	if norm.Cmp(bound) > 0 {
		t.Errorf("SamplePre output norm %s exceeds expected bound %s", norm, bound)
	}
}

func TestSamplePreMultipleTargets(t *testing.T) {
	r := testRing
	params := trapgen.DefaultParams(r)

	pk, td, err := trapgen.TrapGen(params, rand.Reader)
	if err != nil {
		t.Fatalf("TrapGen failed: %v", err)
	}

	for trial := 0; trial < 5; trial++ {
		target := sampler.SampleBigUniformPoly(r, rand.Reader)
		s, err := SamplePre(r, pk, td, target, 350, rand.Reader)
		if err != nil {
			t.Fatalf("Trial %d: SamplePre failed: %v", trial, err)
		}
		if len(s) != len(pk.A) {
			t.Fatalf("Trial %d: wrong output dimension", trial)
		}
	}
}

func TestSamplePreZeroTarget(t *testing.T) {
	r := testRing
	params := trapgen.DefaultParams(r)

	pk, td, err := trapgen.TrapGen(params, rand.Reader)
	if err != nil {
		t.Fatal(err)
	}

	target := r.NewPoly() // zero target
	s, err := SamplePre(r, pk, td, target, 350, rand.Reader)
	if err != nil {
		t.Fatalf("SamplePre with zero target failed: %v", err)
	}
	if len(s) != len(pk.A) {
		t.Error("Wrong output dimension for zero target")
	}
}

func BenchmarkSamplePre(b *testing.B) {
	r := testRing
	params := trapgen.DefaultParams(r)
	pk, td, _ := trapgen.TrapGen(params, rand.Reader)
	target := sampler.SampleBigUniformPoly(r, rand.Reader)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		SamplePre(r, pk, td, target, 350, rand.Reader)
	}
}
