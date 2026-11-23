package geom

import (
	"fmt"
	"math"
)

// Moments 面积矩汇总：一阶矩定质心，二阶中心矩定主轴方向。
// 质心与主轴角在平移、旋转下随面域做同一变换。
type Moments struct {
	Area      float64 // 有向面积
	Mx, My    float64 // 一阶矩 ∫x dA、∫y dA
	Centroid  Point   // 质心
	Sxx, Syy  float64 // 二阶中心矩 ∫(x-xc)²dA、∫(y-yc)²dA
	Sxy       float64 // 二阶中心矩 ∫(x-xc)(y-yc)dA
	AxisAngle float64 // 主轴方向角，度，[0,180)
}

func (m Moments) String() string {
	return fmt.Sprintf("面积=%.6g Mx=%.6g My=%.6g 质心=(%.6g,%.6g) Sxx=%.6g Syy=%.6g Sxy=%.6g 主轴角=%.4f°",
		m.Area, m.Mx, m.My, m.Centroid.X, m.Centroid.Y, m.Sxx, m.Syy, m.Sxy, m.AxisAngle)
}

// ringRaw 单环有向积分：面积、一阶矩、关于原点的二阶矩。
func ringRaw(p []Point) (a, mx, my, ixx, iyy, ixy float64) {
	n := len(p)
	for i := 0; i < n; i++ {
		u, v := p[i], p[(i+1)%n]
		c := u.X*v.Y - v.X*u.Y
		a += c
		mx += c * (u.X + v.X)
		my += c * (u.Y + v.Y)
		ixx += c * (u.X*u.X + u.X*v.X + v.X*v.X)
		iyy += c * (u.Y*u.Y + u.Y*v.Y + v.Y*v.Y)
		ixy += c * (u.X*v.Y + 2*u.X*u.Y + 2*v.X*v.Y + v.X*u.Y)
	}
	return a / 2, mx / 6, my / 6, ixx / 12, iyy / 12, ixy / 24
}

// buildMoments 由有向积分推出质心、中心矩与主轴角。
func buildMoments(a, mx, my, ixx, iyy, ixy float64) Moments {
	cx, cy := mx/a, my/a
	sxx := ixx - a*cx*cx
	syy := iyy - a*cy*cy
	sxy := ixy - a*cx*cy
	ang := 0.5 * math.Atan2(2*sxy, sxx-syy) * 180 / math.Pi
	if ang = math.Mod(ang, 180); ang < 0 {
		ang += 180
	}
	return Moments{Area: a, Mx: mx, My: my, Centroid: Point{cx, cy},
		Sxx: sxx, Syy: syy, Sxy: sxy, AxisAngle: ang}
}

func scaleOf(p []Point) float64 {
	s := 1.0
	for _, q := range p {
		s = math.Max(s, math.Max(math.Abs(q.X), math.Abs(q.Y)))
	}
	return s
}

// extremumVertex 找鞋带项最负（sign<0）或最正（sign>0）的顶点。
func extremumVertex(p []Point, sign float64) int {
	best, bi := 0.0, 0
	for i := range p {
		u, v := p[i], p[(i+1)%len(p)]
		if c := sign * (u.X*v.Y - v.X*u.Y); c < best {
			best, bi = c, i
		}
	}
	return bi
}

// maxTurnVertex 找拐弯最厉害的顶点，用于共线/零面积定位。
func maxTurnVertex(p []Point) (idx int, turn float64) {
	n := len(p)
	for i := 0; i < n; i++ {
		a, b, c := p[(i+n-1)%n], p[i], p[(i+1)%n]
		t := math.Abs((b.X-a.X)*(c.Y-b.Y) - (c.X-b.X)*(b.Y-a.Y))
		if t > turn {
			turn, idx = t, i
		}
	}
	return
}

func degenerateErr(ring int, p []Point) error {
	idx, turn := maxTurnVertex(p)
	kind := "零面积（自交或回折）"
	if turn <= 1e-9*scaleOf(p)*scaleOf(p) {
		kind = "退化共线"
	}
	return &VertexError{Ring: ring, Index: idx, Err: fmt.Errorf("%w: %s", ErrDegenerate, kind)}
}

// PolygonMoments 单环面积矩；环需逆时针（有向面积为正）。
// 面积为负、顶点不足、退化共线/零面积均返回带顶点号的 *VertexError。
func PolygonMoments(poly []Point) (Moments, error) {
	if len(poly) < 3 {
		return Moments{}, &VertexError{Index: len(poly), Err: ErrTooFewPts}
	}
	tol := 1e-12 * scaleOf(poly) * scaleOf(poly)
	a, mx, my, ixx, iyy, ixy := ringRaw(poly)
	switch {
	case a < -tol:
		return Moments{}, &VertexError{Index: extremumVertex(poly, 1), Err: ErrNegativeArea}
	case a <= tol:
		return Moments{}, degenerateErr(0, poly)
	}
	return buildMoments(a, mx, my, ixx, iyy, ixy), nil
}

// RegionMoments 多环面域面积矩：rings[0] 为外环（逆时针），其余为洞（顺时针）。
// 洞的有向面积为负，自动从外环矩中扣除。
func RegionMoments(rings ...[]Point) (Moments, error) {
	if len(rings) == 0 {
		return Moments{}, &VertexError{Err: ErrTooFewPts}
	}
	scale := 1.0
	for _, r := range rings {
		scale = math.Max(scale, scaleOf(r))
	}
	tol := 1e-12 * scale * scale
	var A, Mx, My, Ixx, Iyy, Ixy float64
	for ri, ring := range rings {
		if len(ring) < 3 {
			return Moments{}, &VertexError{Ring: ri, Index: len(ring), Err: ErrTooFewPts}
		}
		a, mx, my, ixx, iyy, ixy := ringRaw(ring)
		if ri == 0 {
			switch {
			case a < -tol:
				return Moments{}, &VertexError{Index: extremumVertex(ring, 1), Err: ErrNegativeArea}
			case a <= tol:
				return Moments{}, degenerateErr(0, ring)
			}
		} else {
			switch {
			case a > tol:
				return Moments{}, &VertexError{Ring: ri, Index: extremumVertex(ring, -1), Err: ErrBadWinding}
			case a >= -tol:
				return Moments{}, degenerateErr(ri, ring)
			}
		}
		A += a
		Mx += mx
		My += my
		Ixx += ixx
		Iyy += iyy
		Ixy += ixy
	}
	if A <= tol {
		return Moments{}, &VertexError{Err: fmt.Errorf("%w: 洞口抵消外环，总面积非正", ErrDegenerate)}
	}
	return buildMoments(A, Mx, My, Ixx, Iyy, Ixy), nil
}
