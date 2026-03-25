package pke

import (
	"errors"
	"io"
	"math/big"

	"github.com/KarpelesLab/lnp22/ring"
	"github.com/KarpelesLab/lnp22/sampler"
	"github.com/KarpelesLab/lnp22/trapgen"
)

// Params holds PKE configuration.
type Params struct {
	Ring       *ring.BigRing
	P          *big.Int        // plaintext modulus
	TrapParams *trapgen.Params // parameters for NTRU key generation
	ErrorSigma float64         // Gaussian std dev for encryption error
}

// DefaultParams returns default PKE parameters for the given ring.
func DefaultParams(r *ring.BigRing) *Params {
	return &Params{
		Ring:       r,
		P:          big.NewInt(257),
		TrapParams: trapgen.DefaultParams(r),
		ErrorSigma: 3.2,
	}
}

// PublicKey is the NTRU public key h = f^{-1} · g mod q.
type PublicKey struct {
	H ring.BigPoly
}

// SecretKey contains the short NTRU polynomial f for decryption.
type SecretKey struct {
	F ring.BigPoly // short polynomial f
}

// Ciphertext is an encrypted message.
type Ciphertext struct {
	C ring.BigPoly
}

// KeyGen generates a PKE keypair using NTRU key generation.
func KeyGen(params *Params, rng io.Reader) (*PublicKey, *SecretKey, error) {
	key, err := trapgen.GenerateNTRUKey(params.TrapParams, rng)
	if err != nil {
		return nil, nil, err
	}
	return &PublicKey{H: key.H}, &SecretKey{F: key.F}, nil
}

// Encrypt encrypts a plaintext polynomial m with coefficients in [0, p).
// c = h · e + Encode(m) mod q
func Encrypt(params *Params, pk *PublicKey, m ring.BigPoly, rng io.Reader) (*Ciphertext, error) {
	r := params.Ring

	// Validate plaintext range
	for i := 0; i < r.N; i++ {
		if m[i].Sign() < 0 || m[i].Cmp(params.P) >= 0 {
			return nil, errors.New("pke: plaintext coefficient out of range [0, p)")
		}
	}

	// Sample small error e
	e := sampler.SampleBigGaussianPoly(r, params.ErrorSigma, rng)

	// Encode plaintext: lift m from Z_p to Z_q
	encoded := Encode(r, m, params.P)

	// c = h · e + encoded mod q
	he := r.Mul(pk.H, e)
	c := r.Add(he, encoded)

	return &Ciphertext{C: c}, nil
}

// Decrypt decrypts a ciphertext using the secret key.
//
// Compute u = f · c mod q. Since c = h·e + Encode(m) and h = f^{-1}·g:
//
//	f·c = g·e + f·Encode(m)
//
// The g·e term is small noise. After rounding from Z_q to Z_p,
// the noise is absorbed and the original m is recovered.
func Decrypt(params *Params, sk *SecretKey, ct *Ciphertext) (ring.BigPoly, error) {
	r := params.Ring
	u := r.Mul(sk.F, ct.C)
	m := Decode(r, u, params.P)
	return m, nil
}
