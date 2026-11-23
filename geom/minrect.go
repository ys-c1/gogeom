package geom

import (
	"fmt"
	"math"
	"sort"
)

// ConvexHull 单调链凸包，逆时针，剔除共线点，不含重复端点。
func ConvexHull(pts []Point) []Point {
	p := make([]Point, len(pts))
	copy(p, pts)
	sort.Slice(p, func(i, j int) bool {
		if p[i].X != p[j].X {
			return p[i].X < p[j].X
		}
		return p[i].Y < p[j].Y
	})
	cross := func(o, a, b Point) float64 {
		return (a.X-o.X)*(b.Y-o.Y) - (a.Y-o.Y)*(b.X-o.X)
	}
	H := make([]Point, 2*len(p)+1)
	k := 0
	for _, q := range p { // 下壳
		for k >= 2 && cross(H[k-2], H[k-1], q) <= 0 {
			k--
		}
		H[k] = q
		k++
	}
	for i, t := len(p)-2, k+1; i >= 0; i-- { // 上壳
		for k >= t && cross(H[k-2], H[k-1], p[i]) <= 0 {
			k--
		}
		H[k] = p[i]
		k++
	}
	if k <= 1 {
		return H[:k]
	}
	return H[:k-1]
}

// MinAreaRect 最小面积外接矩形（旋转卡壳）。
// 返回中心、长边、短边与长边方向角（度，[0,180)）。
func MinAreaRect(poly []Point) (center Point, long, short, angleDeg float64, err error) {
	if len(poly) < 3 {
		return Point{}, 0, 0, 0, &VertexError{Index: len(poly), Err: ErrTooFewPts}
	}
	hull := ConvexHull(poly)
	if len(hull) < 3 {
		idx, _ := maxTurnVertex(poly)
		return Point{}, 0, 0, 0, &VertexError{Index: idx, Err: fmt.Errorf("%w: 退化共线", ErrDegenerate)}
	}
	best := math.Inf(1)
	n := len(hull)
	for i := 0; i < n; i++ {
		a, b := hull[i], hull[(i+1)%n]
		l := math.Hypot(b.X-a.X, b.Y-a.Y)
		ux, uy := (b.X-a.X)/l, (b.Y-a.Y)/l // 边单位方向
		minU, maxU := math.Inf(1), math.Inf(-1)
		minV, maxV := math.Inf(1), math.Inf(-1)
		for _, q := range hull {
			du := q.X*ux + q.Y*uy
			dv := -q.X*uy + q.Y*ux
			minU, maxU = math.Min(minU, du), math.Max(maxU, du)
			minV, maxV = math.Min(minV, dv), math.Max(maxV, dv)
		}
		w, h := maxU-minU, maxV-minV
		if w*h < best {
			best = w * h
			long, short = math.Max(w, h), math.Min(w, h)
			ang := math.Atan2(uy, ux) * 180 / math.Pi
			if h > w {
				ang += 90 // 长边在垂直方向
			}
			if ang = math.Mod(ang, 180); ang < 0 {
				ang += 180
			}
			angleDeg = ang
			mu, mv := (minU+maxU)/2, (minV+maxV)/2
			center = Point{mu*ux - mv*uy, mu*uy + mv*ux}
		}
	}
	return center, long, short, angleDeg, nil
}
