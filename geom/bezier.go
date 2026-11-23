package geom

import "math"

type Point struct{ X, Y float64 }

func Bezier(p0, p1, p2, p3 Point, t float64) Point {
	u := 1 - t
	a, b := u*u*u, 3*u*u*t
	c, d := 3*u*t*t, t*t*t
	return Point{a*p0.X + b*p1.X + c*p2.X + d*p3.X, a*p0.Y + b*p1.Y + c*p2.Y + d*p3.Y}
}

func Flatten(p0, p1, p2, p3 Point, n int) ([]Point, error) {
	if n < 1 {
		return nil, ErrBadSegments
	}
	out := make([]Point, n+1)
	for i := 0; i <= n; i++ {
		out[i] = Bezier(p0, p1, p2, p3, float64(i)/float64(n))
	}
	return out, nil
}

func Length(pts []Point) float64 {
	var s float64
	for i := 0; i+1 < len(pts); i++ {
		s += math.Hypot(pts[i+1].X-pts[i].X, pts[i+1].Y-pts[i].Y)
	}
	return s
}
