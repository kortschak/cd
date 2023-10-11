// Copyright ©2023 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Copyright ©2018 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Original example at https://pkg.go.dev/gonum.org/v1/gonum/num/quat#example-package-Rotate

package cd_test

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/floats/scalar"

	"github.com/kortschak/cd"
)

// point is a 3-dimensional point/vector.
type point struct {
	x, y, z float64
}

// raise raises the dimensionality of a point to a quaternion.
func raise(p point) cd.H {
	return cd.NewH(0, p.x, p.y, p.z)
}

// rotate performs the quaternion rotation of p by the given quaternion
// and scaling by the scale factor.
func rotate(p point, by cd.H, scale float64) point {
	// Ensure the modulus of by is correctly scaled.
	if len := cd.Abs(by); len != scale {
		by = by.Scale(math.Sqrt(scale) / len)
	}

	// Perform the rotation/scaling.
	pp := by.Mul(raise(p)).Mul(by.Conj())

	// Extract the point.
	e := pp.Elems()
	return point{x: e[1], y: e[2], z: e[3]}
}

// Rotate a cube 120° around the diagonal vector [1, 1, 1].
func Example_rotate() {
	alpha := 2 * math.Pi / 3
	q := raise(point{1, 1, 1})
	scale := 1.0

	q = q.Scale(math.Sin(alpha/2) / cd.Abs(q))
	q = q.Add(cd.Lift[cd.H](math.Cos(alpha / 2)))

	for i, p := range []point{
		{x: 0, y: 0, z: 0},
		{x: 0, y: 0, z: 1},
		{x: 0, y: 1, z: 0},
		{x: 0, y: 1, z: 1},
		{x: 1, y: 0, z: 0},
		{x: 1, y: 0, z: 1},
		{x: 1, y: 1, z: 0},
		{x: 1, y: 1, z: 1},
	} {
		pp := rotate(p, q, scale)

		// Clean up floating point error for clarity.
		pp.x = scalar.Round(pp.x, 2)
		pp.y = scalar.Round(pp.y, 2)
		pp.z = scalar.Round(pp.z, 2)

		fmt.Printf("%d %+v -> %+v\n", i, p, pp)
	}

	// Output:
	//
	// 0 {x:0 y:0 z:0} -> {x:0 y:0 z:0}
	// 1 {x:0 y:0 z:1} -> {x:1 y:0 z:0}
	// 2 {x:0 y:1 z:0} -> {x:0 y:0 z:1}
	// 3 {x:0 y:1 z:1} -> {x:1 y:0 z:1}
	// 4 {x:1 y:0 z:0} -> {x:0 y:1 z:0}
	// 5 {x:1 y:0 z:1} -> {x:1 y:1 z:0}
	// 6 {x:1 y:1 z:0} -> {x:0 y:1 z:1}
	// 7 {x:1 y:1 z:1} -> {x:1 y:1 z:1}
}
