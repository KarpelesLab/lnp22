// Package bigmod provides modular arithmetic over arbitrary-precision integers
// for lattice cryptography with large moduli (e.g., 152-bit).
package bigmod

import "math/big"

var (
	bigZero = new(big.Int)
	bigOne  = big.NewInt(1)
	bigTwo  = big.NewInt(2)
)

// Reduce returns a mod q in [0, q).
func Reduce(a, q *big.Int) *big.Int {
	r := new(big.Int).Mod(a, q) // Go's Mod always returns [0, q)
	return r
}

// CenterReduce returns a mod q in [-(q-1)/2, (q-1)/2].
func CenterReduce(a, q *big.Int) *big.Int {
	r := Reduce(a, q)
	half := new(big.Int).Rsh(q, 1) // q/2
	if r.Cmp(half) > 0 {
		r.Sub(r, q)
	}
	return r
}

// ModAdd returns (a + b) mod q.
func ModAdd(a, b, q *big.Int) *big.Int {
	r := new(big.Int).Add(a, b)
	return r.Mod(r, q)
}

// ModSub returns (a - b) mod q.
func ModSub(a, b, q *big.Int) *big.Int {
	r := new(big.Int).Sub(a, b)
	return r.Mod(r, q)
}

// ModMul returns (a * b) mod q.
func ModMul(a, b, q *big.Int) *big.Int {
	r := new(big.Int).Mul(a, b)
	return r.Mod(r, q)
}

// ModExp returns (base^exp) mod q.
func ModExp(base, exp, q *big.Int) *big.Int {
	return new(big.Int).Exp(base, exp, q)
}

// ModInverse returns a^{-1} mod q. Panics if gcd(a,q) != 1.
func ModInverse(a, q *big.Int) *big.Int {
	r := new(big.Int).ModInverse(a, q)
	if r == nil {
		panic("bigmod: element not invertible")
	}
	return r
}

// Abs returns |a| as a new big.Int.
func Abs(a *big.Int) *big.Int {
	return new(big.Int).Abs(a)
}

// --- In-place variants for performance-critical loops ---

// ReduceInPlace sets dst = a mod q in [0, q). Returns dst.
func ReduceInPlace(dst, a, q *big.Int) *big.Int {
	return dst.Mod(a, q)
}

// ModAddInPlace sets dst = (a + b) mod q. Returns dst.
func ModAddInPlace(dst, a, b, q *big.Int) *big.Int {
	dst.Add(a, b)
	return dst.Mod(dst, q)
}

// ModSubInPlace sets dst = (a - b) mod q. Returns dst.
func ModSubInPlace(dst, a, b, q *big.Int) *big.Int {
	dst.Sub(a, b)
	return dst.Mod(dst, q)
}

// ModMulInPlace sets dst = (a * b) mod q. Returns dst.
func ModMulInPlace(dst, a, b, q *big.Int) *big.Int {
	dst.Mul(a, b)
	return dst.Mod(dst, q)
}

// --- Rounding ---

// RoundP computes the rounding function ⌈x⌉_p = round(p·x / q).
// Maps x ∈ Z_q to Z_p. Result is in [0, p).
func RoundP(x, p, q *big.Int) *big.Int {
	// result = floor((p * x + q/2) / q)  [nearest integer rounding]
	num := new(big.Int).Mul(p, x)
	halfQ := new(big.Int).Rsh(q, 1)
	num.Add(num, halfQ)
	r := num.Div(num, q)
	// Clamp to [0, p) in case of boundary rounding
	if r.Cmp(p) >= 0 {
		r.Mod(r, p)
	}
	return r
}

// ScaleUp computes round(q * m / p), lifting m ∈ Z_p to Z_q.
// This is the encoding counterpart of RoundP: RoundP(ScaleUp(m)) = m.
func ScaleUp(m, p, q *big.Int) *big.Int {
	// result = floor((q * m + p/2) / p)  [nearest integer rounding]
	num := new(big.Int).Mul(q, m)
	half := new(big.Int).Rsh(p, 1)
	num.Add(num, half)
	return num.Div(num, p)
}

// NewInt is a convenience wrapper for big.NewInt.
func NewInt(v int64) *big.Int {
	return big.NewInt(v)
}
