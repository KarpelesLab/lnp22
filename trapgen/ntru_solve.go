package trapgen

import (
	"errors"
	"math/big"

	"github.com/KarpelesLab/lnp22/internal/bigmod"
	"github.com/KarpelesLab/lnp22/ring"
)

// SolveNTRU finds F, G such that f·G - g·F = q in Z[X]/(X^N+1).
// Uses the tower-based solver (Pornin-Prest, ASIACRYPT 2019).
//
// Base case (N=1): extended GCD on integers f, g to get F, G with fG - gF = q.
// Recursive case: project f, g to half degree via f'(x) = f(√x), solve,
// then lift back.
func SolveNTRU(r *ring.BigRing, f, g ring.BigPoly, q *big.Int) (F, G ring.BigPoly, err error) {
	if r.N == 1 {
		return solveNTRUBase(f[0], g[0], q)
	}

	// Project to half degree: f'(x) = Norm_{K/K'}(f) = f(√x) * f(-√x)
	// For f ∈ Z[x]/(x^N+1), the field norm to the sub-field Z[x]/(x^{N/2}+1)
	// is N_f(x) = f(√x) * f(-√x) = f_even(x)^2 - x * f_odd(x)^2
	halfN := r.N / 2
	rHalf, err := ring.NewBig(halfN, q)
	if err != nil {
		return nil, nil, err
	}

	fNorm := fieldNorm(r, rHalf, f)
	gNorm := fieldNorm(r, rHalf, g)

	// Recursively solve for the projected polynomials
	Fp, Gp, err := SolveNTRU(rHalf, fNorm, gNorm, q)
	if err != nil {
		return nil, nil, err
	}

	// Lift: F = F' * adj(f), G = G' * adj(g)
	// where adj(f)(x) = f(-√x) lifted back to degree N
	fAdj := adjoint(r, f)
	gAdj := adjoint(r, g)

	// Lift F', G' from half-degree to full degree (embed in even coefficients)
	FpLift := liftPoly(r, rHalf, Fp)
	GpLift := liftPoly(r, rHalf, Gp)

	F = r.Mul(FpLift, gAdj)
	G = r.Mul(GpLift, fAdj)

	// Reduce: make F, G small using Babai-style reduction
	// Normalize: F and G may need sign correction
	// Verify: f*G - g*F = q
	fG := r.Mul(f, G)
	gF := r.Mul(g, F)
	diff := r.Sub(fG, gF)

	qPoly := r.NewPoly()
	qPoly[0].Set(q)

	if !r.Equal(diff, qPoly) {
		return nil, nil, errors.New("trapgen: NTRU equation verification failed")
	}

	return F, G, nil
}

// solveNTRUBase handles the base case N=1: find F, G ∈ Z with fG - gF = q.
func solveNTRUBase(f, g, q *big.Int) (ring.BigPoly, ring.BigPoly, error) {
	// Extended GCD: find u, v such that f*u + g*v = gcd(f, g)
	fv := bigmod.Reduce(f, q)
	gv := bigmod.Reduce(g, q)
	if fv.Sign() == 0 && gv.Sign() == 0 {
		return nil, nil, errors.New("trapgen: f and g are both zero mod q")
	}

	// Work in Z (not Z_q) for the NTRU equation
	u := new(big.Int)
	v := new(big.Int)
	d := new(big.Int).GCD(u, v, f, g)

	// We need fG - gF = q, i.e., f*G + g*(-F) = q.
	// From f*u + g*v = d, multiply by q/d: f*(u*q/d) + g*(v*q/d) = q
	if new(big.Int).Mod(q, d).Sign() != 0 {
		return nil, nil, errors.New("trapgen: gcd(f,g) does not divide q")
	}

	scale := new(big.Int).Div(q, d)
	G := make(ring.BigPoly, 1)
	F := make(ring.BigPoly, 1)
	G[0] = new(big.Int).Mul(u, scale)                   // G = u * q/d
	F[0] = new(big.Int).Neg(new(big.Int).Mul(v, scale)) // F = -(v * q/d), so fG - gF = f*u*q/d + g*v*q/d = q

	return F, G, nil
}

// fieldNorm computes the field norm from R = Z[x]/(x^N+1) to R' = Z[x]/(x^{N/2}+1).
// N_f(x) = f_even(x)^2 - x * f_odd(x)^2
// where f(y) = f_even(y^2) + y * f_odd(y^2).
func fieldNorm(r *ring.BigRing, rHalf *ring.BigRing, f ring.BigPoly) ring.BigPoly {
	halfN := rHalf.N
	fEven := rHalf.NewPoly()
	fOdd := rHalf.NewPoly()

	for i := 0; i < halfN; i++ {
		fEven[i].Set(f[2*i])
		if 2*i+1 < r.N {
			fOdd[i].Set(f[2*i+1])
		}
	}

	evenSq := rHalf.Mul(fEven, fEven)
	oddSq := rHalf.Mul(fOdd, fOdd)

	// x * oddSq in Z[x]/(x^{N/2}+1): multiply by x means shift by 1 with wraparound
	xOddSq := rHalf.NewPoly()
	for i := 0; i < halfN; i++ {
		if i == 0 {
			// x * x^{N/2-1} = x^{N/2} ≡ -1, so coefficient 0 gets -last
			xOddSq[0].Neg(oddSq[halfN-1])
			xOddSq[0].Mod(xOddSq[0], rHalf.Q)
		} else {
			xOddSq[i].Set(oddSq[i-1])
		}
	}

	return rHalf.Sub(evenSq, xOddSq)
}

// adjoint computes the "conjugate" adj(f)(x) where if f(y) = f_even(y^2) + y*f_odd(y^2),
// then adj(f)(y) = f_even(y^2) - y*f_odd(y^2).
// This corresponds to f(-√x) in the tower decomposition.
func adjoint(r *ring.BigRing, f ring.BigPoly) ring.BigPoly {
	result := r.NewPoly()
	for i := 0; i < r.N; i++ {
		if i%2 == 0 {
			result[i].Set(f[i])
		} else {
			result[i].Neg(f[i])
			result[i].Mod(result[i], r.Q)
		}
	}
	return result
}

// liftPoly embeds a polynomial from R' = Z[x]/(x^{N/2}+1) into R = Z[x]/(x^N+1)
// by placing coefficients at even indices.
func liftPoly(r *ring.BigRing, rHalf *ring.BigRing, p ring.BigPoly) ring.BigPoly {
	result := r.NewPoly()
	for i := 0; i < rHalf.N; i++ {
		result[2*i].Set(p[i])
	}
	return result
}
