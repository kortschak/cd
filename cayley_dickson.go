// Copyright ©2023 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package cd implements Cayley–Dickson construction.
package cd

import (
	"math"
	"unsafe"
)

type Field interface {
	~float32 | ~float64
}

// Value is a Cayley-Dickson algebra value.
type Value[A any, F Field] interface {
	// Real returns the real part of the value.
	Real() F
	// Imag returns the imaginary vector part of
	// the value.
	Imag() A

	// Scale returns the value element-wise
	// scaled by f.
	Scale(f F) A

	// Neg returns the negation of the value.
	Neg() A
	// Conj returns the Cayley-Dickson conjugate
	// of the value.
	Conj() A
	// Add returns the value with a added.
	Add(a A) A
	// Mul returns the Cayley-Dickson product
	// of the value and a.
	Mul(a A) A

	// Elems returns the field elements of the
	// value.
	Elems() []F

	comparable
}

// Abs returns the absolute value (also called the modulus) of x.
func Abs[A Value[A, F], F Field](x A) F {
	return F(math.Sqrt(float64(x.Mul(x.Conj()).Real())))
}

// Exp returns e**x, the base-e exponential of x.
func Exp[A Value[A, F], F Field](x A) A {
	w := x.Real()
	uv := x.Imag()
	var zero A
	if uv == zero {
		return Lift[A](F(math.Exp(float64(w))))
	}
	v := float64(Abs(uv))
	e := math.Exp(float64(w))
	s, c := math.Sincos(v)
	return Lift[A](F(e * c)).Add(uv.Scale(F(e * s / v)))
}

// Inf returns an infinity for the algebra, with all elements positive infinity.
func Inf[A Value[A, F], F Field]() A {
	var zero A
	e := zero.Elems()
	for i := range e {
		e[i] = F(math.Inf(1))
	}
	return *(*A)(unsafe.Pointer(&e))
}

// Inv returns the inverse of x.
func Inv[A Value[A, F], F Field](x A) A {
	xc := x.Conj()
	return xc.Scale(1 / x.Mul(xc).Real())
}

// Lift returns an element of the algebra with the real part equal to f.
func Lift[A Value[A, F], F Field](f F) A {
	var zero A
	e := zero.Elems()
	e[0] = f
	return *(*A)(unsafe.Pointer(&e[0]))
}

// Log returns the natural logarithm of x.
func Log[A Value[A, F], F Field](x A) A {
	w := float64(x.Real())
	uv := x.Imag()
	var zero A
	if uv == zero {
		return Lift[A](F(math.Log(w)))
	}
	v := float64(Abs(uv))
	return Lift[A](F(math.Log(float64(Abs(x))))).Add(uv.Scale(F(math.Atan2(v, w) / v)))
}

// Pow returns x**r, the base-x exponential of r.
// For generalized compatibility with math.Pow:
//
//	Pow(0, ±0) returns Lift(1)
//	Pow(0, c) for real(c)<0 returns Lift(Inf) if imag(c) is zero,
//	    otherwise Inf.
func Pow[A Value[A, F], F Field](x, r A) A {
	var zero A
	if x == zero {
		w := r.Real()
		uv := r.Imag()
		switch {
		case w == 0:
			return Lift[A](1)
		case w < 0:
			if uv == zero {
				return Lift[A](F(math.Inf(1)))
			}
			return Inf[A]()
		case w < 0:
			return zero
		}
	}
	return Exp(Log(x).Mul(r))
}

// PowFloat returns x**r, the base-x exponential of r.
// For generalized compatibility with math.Pow:
//
//	PowFloat(0, ±0) returns Lift(1)
//	PowFloat(0, c) for c<0 returns Lift(inf).
func PowFloat[A Value[A, F], F Field](x A, r F) A {
	var zero A
	if x == zero {
		switch {
		case r == 0:
			return Lift[A](1)
		case r < 0:
			return Inf[A]()
		case r > 0:
			return zero
		}
	}
	return Exp(Log(x).Scale(r))
}

// Sqrt returns the square root of x.
func Sqrt[A Value[A, F], F Field](x A) A {
	var zero A
	if x == zero {
		return zero
	}
	return PowFloat(x, 0.5)
}

// R is a real algebra.
type R float64

func (x R) Neg() R            { return -x }
func (x R) Conj() R           { return x }
func (x R) Add(y R) R         { return x + y }
func (x R) Mul(y R) R         { return x * y }
func (x R) Scale(f float64) R { return x * R(f) }
func (x R) Real() float64     { return float64(x) }
func (x R) Imag() R           { return 0 }
func (x R) Elems() []float64  { return []float64{float64(x)} }

// C is a complex algebra value.
type C = Construction[R, float64]

// NewC returns a new float64-based complex number.
func NewC(a, b float64) C {
	return cons(
		R(a),
		R(b),
	)
}

// H is a quaternion algebra value.
type H = Construction[C, float64]

// NewH returns an float64-based quaternion.
func NewH(a, b, c, d float64) H {
	return cons(
		cons(
			R(a),
			R(b),
		),
		cons(
			R(c),
			R(d),
		),
	)
}

// O is an octonion algebra value.
type O = Construction[H, float64]

// NewO returns an float64-based octonion.
func NewO(a, b, c, d, e, f, g, h float64) O {
	return cons(
		cons(
			cons(
				R(a),
				R(b),
			),
			cons(
				R(c),
				R(d),
			),
		),
		cons(
			cons(
				R(e),
				R(f),
			),
			cons(
				R(g),
				R(h),
			),
		),
	)
}

// O is a sedonion algebra value.
type S = Construction[O, float64]

// NewS returns an float64-based sedonion.
func NewS(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p float64) S {
	return cons(
		cons(
			cons(
				cons(
					R(a),
					R(b),
				),
				cons(
					R(c),
					R(d),
				),
			),
			cons(
				cons(
					R(e),
					R(f),
				),
				cons(
					R(g),
					R(h),
				),
			),
		),
		cons(
			cons(
				cons(
					R(i),
					R(j),
				),
				cons(
					R(k),
					R(l),
				),
			),
			cons(
				cons(
					R(m),
					R(n),
				),
				cons(
					R(o),
					R(p),
				),
			),
		),
	)
}

// Constructions is a Cayley-Dickson algebra construction.
type Construction[A Value[A, F], F Field] [2]A

// Real returns the real parts of x.
func (x Construction[A, E]) Real() E {
	return x[0].Real()
}

// Imag returns the imaginary parts of x.
func (x Construction[A, F]) Imag() Construction[A, F] {
	return cons(x[0].Imag(), x[1])
}

// Neg returns a negated copy of x. It is equivalent to x.Scale(-1).
func (x Construction[A, F]) Neg() Construction[A, F] {
	return cons(x[0].Neg(), x[1].Neg())
}

// Conj returns the conjugate of x.
func (x Construction[A, F]) Conj() Construction[A, F] {
	return cons(x[0].Conj(), x[1].Neg())
}

// Add returns the result of adding x and y element wise.
func (x Construction[A, F]) Add(y Construction[A, F]) Construction[A, F] {
	return cons(x[0].Add(y[0]), x[1].Add(y[1]))
}

// Mul returns the product of x and y.
func (x Construction[A, F]) Mul(y Construction[A, F]) Construction[A, F] {
	a, b := x[0], x[1]
	c, d := y[0], y[1]
	return cons(
		a.Mul(c).Add((d.Conj().Mul(b)).Neg()),
		d.Mul(a).Add(b.Mul(c.Conj())),
	)
}

// Scale returns the result of scaling each element of x by f.
func (x Construction[A, F]) Scale(f F) Construction[A, F] {
	return cons(x[0].Scale(f), x[1].Scale(f))
}

// Elems returns the elements of x.
func (x Construction[A, F]) Elems() []F {
	var zero F
	return unsafe.Slice((*F)(unsafe.Pointer(&x)), unsafe.Sizeof(x)/unsafe.Sizeof(zero))
}

func cons[A Value[A, F], F Field](a, b A) Construction[A, F] {
	return Construction[A, F]{a, b}
}
