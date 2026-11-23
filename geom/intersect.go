package geom

import (
	"math"
	"sort"
)

// seg 是边界被交点切开后的小段。
type seg struct{ from, to Point }

type xMark struct {
	t  float64 // 在边上的参数
	id int     // 交点编号
}

// Validate 检查多边形是否可作为布尔运算输入：
// 顶点不足、顶点重合、自交未拆分三种情况拒绝并报出序号。
func Validate(poly []Point) error {
	if len(poly) < 3 {
		return &VertexError{Err: ErrTooFewPts, Index: len(poly)}
	}
	for i := range poly {
		for j := i + 1; j < len(poly); j++ {
			if poly[i] == poly[j] {
				return &VertexError{Err: ErrDupVertex, Index: j}
			}
		}
	}
	n := len(poly)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if j == i+1 || (i == 0 && j == n-1) {
				continue // 邻边共享顶点，不算自交
			}
			if edgesTouch(poly[i], poly[(i+1)%n], poly[j], poly[(j+1)%n]) {
				return &VertexError{Err: ErrSelfIntersect, Index: j}
			}
		}
	}
	return nil
}

// Intersect 计算两个简单多边形的布尔交，返回组成交集的环（统一为逆时针）。
// 输入非法时返回带序号的 *VertexError。
func Intersect(a, b []Point) ([][]Point, error) {
	if err := Validate(a); err != nil {
		return nil, err
	}
	if err := Validate(b); err != nil {
		return nil, err
	}
	a, b = ccw(a), ccw(b)
	scale := bboxScale(a, b)
	tol := 1e-9 * scale

	// 1. 求全部边界交点（含共线重叠端点），按位置去重。
	var xpts []Point
	addX := func(p Point) int {
		for i, q := range xpts {
			if math.Hypot(p.X-q.X, p.Y-q.Y) <= tol {
				return i
			}
		}
		xpts = append(xpts, p)
		return len(xpts) - 1
	}
	aMarks := make([][]xMark, len(a))
	bMarks := make([][]xMark, len(b))
	addMark := func(ea int, ta float64, eb int, tb float64, p Point) {
		id := addX(p)
		aMarks[ea] = append(aMarks[ea], xMark{ta, id})
		bMarks[eb] = append(bMarks[eb], xMark{tb, id})
	}
	for i := range a {
		for j := range b {
			edgeX(a[i], a[(i+1)%len(a)], b[j], b[(j+1)%len(b)], tol, i, j, addMark)
		}
	}

	// 2. 无交点：只能是包含或相离。
	if len(xpts) == 0 {
		switch {
		case Inside(b, a[0].X, a[0].Y):
			return [][]Point{a}, nil
		case Inside(a, b[0].X, b[0].Y):
			return [][]Point{b}, nil
		default:
			return nil, nil
		}
	}

	// 3. 交点按参数排序，把两条边界切成小段。
	aSegs := splitEdges(a, aMarks, xpts)
	bSegs := splitEdges(b, bMarks, xpts)

	// 4. 保留落在对方内部的小段（A、B 重合段只留一份）。
	seen := map[[2]Point]bool{}
	var kept []seg
	keep := func(s seg) {
		k := segKey(s)
		if !seen[k] {
			seen[k] = true
			kept = append(kept, s)
		}
	}
	for _, s := range aSegs {
		if insideOther(s, b, tol, scale) {
			keep(s)
		}
	}
	for _, s := range bSegs {
		if insideOther(s, a, tol, scale) {
			keep(s)
		}
	}

	// 5. 把保留的小段首尾相接串成环。
	return chainRings(kept), nil
}

// edgeX 求两条边的交点：一般相交给出一个点，共线重叠给出重叠区间端点。
func edgeX(a1, a2, b1, b2 Point, tol float64, ea, eb int, add func(ea int, ta float64, eb int, tb float64, p Point)) {
	d1 := Point{a2.X - a1.X, a2.Y - a1.Y}
	d2 := Point{b2.X - b1.X, b2.Y - b1.Y}
	l1 := math.Hypot(d1.X, d1.Y)
	l2 := math.Hypot(d2.X, d2.Y)
	den := d1.X*d2.Y - d1.Y*d2.X
	w := Point{b1.X - a1.X, b1.Y - a1.Y}

	// 参数恰好在端点处时，把交点吸附为顶点原值，保证两侧链上坐标一致。
	snap := func(t float64, p0, p1 Point) (float64, Point, bool) {
		if t <= 1e-9 {
			return 0, p0, true
		}
		if t >= 1-1e-9 {
			return 1, p1, true
		}
		return t, Point{}, false
	}

	if math.Abs(den) > 1e-12*l1*l2 {
		t := (w.X*d2.Y - w.Y*d2.X) / den
		u := (w.X*d1.Y - w.Y*d1.X) / den
		if t < -1e-9 || t > 1+1e-9 || u < -1e-9 || u > 1+1e-9 {
			return
		}
		var p Point
		if tt, pp, ok := snap(t, a1, a2); ok {
			t, p = tt, pp
		} else if uu, pp, ok := snap(u, b1, b2); ok {
			u, p = uu, pp
		} else {
			p = Point{a1.X + d1.X*t, a1.Y + d1.Y*t}
		}
		add(ea, t, eb, u, p)
		return
	}

	// 平行：不共线则无交。
	if math.Abs(w.X*d1.Y-w.Y*d1.X) > tol*l1 {
		return
	}
	// 共线：投影到边 ea 上求重叠区间 [lo, hi]。
	tt1 := (w.X*d1.X + w.Y*d1.Y) / (l1 * l1)
	tt2 := tt1 + (d2.X*d1.X+d2.Y*d1.Y)/(l1*l1)
	lo, hi := tt1, tt2
	if lo > hi {
		lo, hi = hi, lo
	}
	if lo > 1+1e-9 || hi < -1e-9 {
		return
	}
	lo = math.Max(lo, 0)
	hi = math.Min(hi, 1)
	addEnd := func(t float64) {
		var p Point
		if tt, pp, ok := snap(t, a1, a2); ok {
			t, p = tt, pp
		} else {
			p = Point{a1.X + d1.X*t, a1.Y + d1.Y*t}
		}
		var u float64
		if math.Abs(d2.X) >= math.Abs(d2.Y) && d2.X != 0 {
			u = (p.X - b1.X) / d2.X
		} else if d2.Y != 0 {
			u = (p.Y - b1.Y) / d2.Y
		}
		add(ea, t, eb, u, p)
	}
	addEnd(lo)
	if hi-lo > 1e-9 {
		addEnd(hi)
	}
}

// splitEdges 按边上交点参数排序后把多边形切成小段。
func splitEdges(poly []Point, marks [][]xMark, xpts []Point) []seg {
	var out []seg
	for i := range poly {
		ms := marks[i]
		sort.Slice(ms, func(x, y int) bool { return ms[x].t < ms[y].t })
		cur := poly[i]
		for _, m := range ms {
			p := xpts[m.id]
			if p != cur {
				out = append(out, seg{cur, p})
				cur = p
			}
		}
		next := poly[(i+1)%len(poly)]
		if cur != next {
			out = append(out, seg{cur, next})
		}
	}
	return out
}

// insideOther 判定小段是否属于交集边界：中点在对方内部，或贴在对方
// 边界上且本多边形（已规范为逆时针，内部在行进方向左侧）的内侧在对方内。
func insideOther(s seg, other []Point, tol, scale float64) bool {
	m := Point{(s.from.X + s.to.X) / 2, (s.from.Y + s.to.Y) / 2}
	if onBoundary(other, m, tol) {
		d := Point{s.to.X - s.from.X, s.to.Y - s.from.Y}
		l := math.Hypot(d.X, d.Y)
		delta := 1e-6 * scale
		return Inside(other, m.X-d.Y/l*delta, m.Y+d.X/l*delta)
	}
	return Inside(other, m.X, m.Y)
}

// chainRings 把小段按端点串成闭合环。
func chainRings(segs []seg) [][]Point {
	adj := map[Point][]int{}
	for i, s := range segs {
		adj[s.from] = append(adj[s.from], i)
	}
	used := make([]bool, len(segs))
	var rings [][]Point
	for i := range segs {
		if used[i] {
			continue
		}
		start := segs[i].from
		ring := []Point{}
		cur, k := segs[i], i
		for {
			used[k] = true
			ring = append(ring, cur.from)
			if cur.to == start {
				break
			}
			nxt := -1
			for _, c := range adj[cur.to] {
				if !used[c] {
					nxt = c
					break
				}
			}
			if nxt < 0 {
				break // 开口链，防御性退出
			}
			k, cur = nxt, segs[nxt]
		}
		if len(ring) >= 3 {
			rings = append(rings, ring)
		}
	}
	return rings
}

// ccw 用鞋带公式符号自推绕向，统一返回逆时针顶点序。
func ccw(p []Point) []Point {
	var s float64
	for i := range p {
		a, b := p[i], p[(i+1)%len(p)]
		s += a.X*b.Y - b.X*a.Y
	}
	if s < 0 {
		out := make([]Point, len(p))
		for i, q := range p {
			out[len(p)-1-i] = q
		}
		return out
	}
	return append([]Point(nil), p...)
}

func bboxScale(a, b []Point) float64 {
	minX, minY := math.Inf(1), math.Inf(1)
	maxX, maxY := math.Inf(-1), math.Inf(-1)
	for _, ps := range [][]Point{a, b} {
		for _, p := range ps {
			minX = math.Min(minX, p.X)
			minY = math.Min(minY, p.Y)
			maxX = math.Max(maxX, p.X)
			maxY = math.Max(maxY, p.Y)
		}
	}
	s := math.Max(maxX-minX, maxY-minY)
	if s < 1 {
		s = 1
	}
	return s
}

func segKey(s seg) [2]Point {
	if s.to.X < s.from.X || (s.to.X == s.from.X && s.to.Y < s.from.Y) {
		return [2]Point{s.to, s.from}
	}
	return [2]Point{s.from, s.to}
}

func onBoundary(poly []Point, p Point, tol float64) bool {
	for i := range poly {
		if distSeg(poly[i], poly[(i+1)%len(poly)], p) <= tol {
			return true
		}
	}
	return false
}

func distSeg(a, b, p Point) float64 {
	dx, dy := b.X-a.X, b.Y-a.Y
	l2 := dx*dx + dy*dy
	if l2 == 0 {
		return math.Hypot(p.X-a.X, p.Y-a.Y)
	}
	t := ((p.X-a.X)*dx + (p.Y-a.Y)*dy) / l2
	t = math.Max(0, math.Min(1, t))
	return math.Hypot(p.X-a.X-t*dx, p.Y-a.Y-t*dy)
}

// edgesTouch 判断两条边是否相交或相触（用于自交检测，精确谓词）。
func edgesTouch(p1, p2, p3, p4 Point) bool {
	x := func(o, a, b Point) float64 {
		return (a.X-o.X)*(b.Y-o.Y) - (a.Y-o.Y)*(b.X-o.X)
	}
	on := func(a, b, p Point) bool {
		return p.X >= math.Min(a.X, b.X) && p.X <= math.Max(a.X, b.X) &&
			p.Y >= math.Min(a.Y, b.Y) && p.Y <= math.Max(a.Y, b.Y)
	}
	d1, d2 := x(p3, p4, p1), x(p3, p4, p2)
	d3, d4 := x(p1, p2, p3), x(p1, p2, p4)
	if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
		((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)) {
		return true
	}
	return (d1 == 0 && on(p3, p4, p1)) || (d2 == 0 && on(p3, p4, p2)) ||
		(d3 == 0 && on(p1, p2, p3)) || (d4 == 0 && on(p1, p2, p4))
}
