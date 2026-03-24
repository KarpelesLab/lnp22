// Package modmath provides modular arithmetic utilities for lattice cryptography.
package modmath

// Reduce returns a mod q in the range [0, q).
func Reduce(a, q int64) int64 {
	r := a % q
	if r < 0 {
		r += q
	}
	return r
}

// CenterReduce returns a mod q in the range [-(q-1)/2, (q-1)/2].
func CenterReduce(a, q int64) int64 {
	r := Reduce(a, q)
	if r > q/2 {
		r -= q
	}
	return r
}

// ModMul returns (a * b) mod q without overflow for values where |a*b| < 2^63.
func ModMul(a, b, q int64) int64 {
	return Reduce(a%q*(b%q), q)
}

// ModAdd returns (a + b) mod q.
func ModAdd(a, b, q int64) int64 {
	return Reduce(a+b, q)
}

// ModSub returns (a - b) mod q.
func ModSub(a, b, q int64) int64 {
	return Reduce(a-b, q)
}

// ModExp returns (base^exp) mod q using square-and-multiply.
func ModExp(base, exp, q int64) int64 {
	base = Reduce(base, q)
	result := int64(1)
	for exp > 0 {
		if exp&1 == 1 {
			result = ModMul(result, base, q)
		}
		base = ModMul(base, base, q)
		exp >>= 1
	}
	return result
}

// ModInverse returns the modular inverse of a mod q using the extended Euclidean algorithm.
// Panics if gcd(a, q) != 1.
func ModInverse(a, q int64) int64 {
	a = Reduce(a, q)
	if a == 0 {
		panic("modmath: inverse of zero")
	}
	// Extended Euclidean algorithm
	old_r, r := a, q
	old_s, s := int64(1), int64(0)
	for r != 0 {
		quotient := old_r / r
		old_r, r = r, old_r-quotient*r
		old_s, s = s, old_s-quotient*s
	}
	if old_r != 1 {
		panic("modmath: element not invertible")
	}
	return Reduce(old_s, q)
}

// Abs returns the absolute value of a.
func Abs(a int64) int64 {
	if a < 0 {
		return -a
	}
	return a
}
