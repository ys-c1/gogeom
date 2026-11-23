package geom

import (
	"fmt"
	"math"
	"sort"
)

const strokeEps = 1e-9

// miter 极限：1+cos(夹角) 低于该值时 miter 过长，退化为 bevel 两点。
const miterLimitDenom = 0.125

// StrokeError 描边错误，携带相关段号。
type StrokeError struct {
	Op   string // "线宽为零" / "法向退化" / "自交未拆分"
	Seg  int    // 相关段号
	SegB int    // 第二个段号（仅自交），无为 -1
}

func (e *StrokeError) Error() string {
	if e.SegB >= 0 {
		return fmt.Sprintf("%s: 段 %d 与段 %d", e.Op, e.Seg, e.SegB)
	}
	return fmt.Sprintf("%s: 段 %d", e.Op, e.Seg)
}

// Stroke 原始直折描边：两侧法向偏移链直接首尾相接成多边形，不做自交裁剪。
// 保持朴素行为不变，作为 Outline 裁剪效果的对照基线。
func Stroke(pts []Point, width float64) ([]Point, error) {
	left, right, _, err := offsetSides(pts, width)
	if err != nil {
		return nil, err
	}
	return stitch(left, right), nil
}

// OutlineResult 自交裁剪后的描边结果。
type OutlineResult struct {
	Shells        [][]Point // 保留的外环
	Holes         [][]Point // 内孔（闭合折线的内侧边界）
	Discarded     int       // 丢弃的绕入环数
	Intersections int       // 由线宽决定的自交切分点个数
	Area          float64   // 净面积（鞋带公式，壳减孔）
}

// Outline 折线加粗成轮廓：按每段法向偏移出内外两条链，
// 在自交处切分并丢掉绕进去的环，返回保留环与统计。
func Outline(pts []Point, width float64) (*OutlineResult, error) {
	left, right, closed, err := offsetSides(pts, width)
	if err != nil {
		return nil, err
	}
	res := &OutlineResult{}
	if closed {
		// 闭合折线：两条链各自闭合成环，分别解缠；面积大者为壳、小者为孔。
		loopsL, nL := splitLoops(left)
		loopsR, nR := splitLoops(right)
		res.Intersections = nL + nR
		keepL, dropL := keepOuter(loopsL)
		keepR, dropR := keepOuter(loopsR)
		res.Discarded = dropL + dropR
		areaL, areaR := absSum(keepL), absSum(keepR)
		shells, holes := keepL, keepR
		if areaR > areaL {
			shells, holes, areaL, areaR = keepR, keepL, areaR, areaL
		}
		res.Shells, res.Holes = shells, holes
		res.Area = areaL - areaR
	} else {
		loops, nX := splitLoops(stitch(left, right))
		res.Intersections = nX
		keep, drop := keepOuter(loops)
		res.Discarded = drop
		res.Shells = keep
		res.Area = math.Abs(signedSum(keep))
	}
	if len(res.Shells) == 0 {
		return nil, ErrDegenerate
	}
	// 校验：保留环不得残留未拆分的自交，否则报出段号。
	for _, lp := range append(append([][]Point{}, res.Shells...), res.Holes...) {
		if e := checkSimple(lp); e != nil {
			return nil, e
		}
	}
	return res, nil
}

// offsetSides 计算折线两侧的法向偏移链。left 为左法向一侧，right 为右侧。
// 首末点重合视为闭合环，两条链各自闭合。线宽为零、段长为零（法向退化）报段号。
func offsetSides(pts []Point, width float64) (left, right []Point, closed bool, err error) {
	if width == 0 {
		return nil, nil, false, &StrokeError{Op: "线宽为零", Seg: 0, SegB: -1}
	}
	if width < 0 {
		width = -width
	}
	if len(pts) < 2 {
		return nil, nil, false, ErrTooFewPts
	}
	v := pts
	closed = len(pts) >= 4 && pts[0] == pts[len(pts)-1]
	if closed {
		v = pts[:len(pts)-1]
	}
	m := len(v)
	nseg := m - 1
	if closed {
		nseg = m
	}
	n := make([]Point, nseg)
	for i := 0; i < nseg; i++ {
		a, b := v[i], v[(i+1)%m]
		d := Point{b.X - a.X, b.Y - a.Y}
		l := math.Hypot(d.X, d.Y)
		if l < strokeEps {
			return nil, nil, false, &StrokeError{Op: "法向退化", Seg: i, SegB: -1}
		}
		n[i] = Point{-d.Y / l, d.X / l} // 单位左法向
	}
	h := width / 2
	return offsetChain(v, n, h, +1, closed), offsetChain(v, n, h, -1, closed), closed, nil
}

// offsetChain 沿各段单位左法向 n 把每个顶点偏移 s*h。
// 顶点处相邻法向用 miter 合并；夹角过锐（miter 超限）或反向时退化为 bevel 两点。
func offsetChain(v []Point, n []Point, h, s float64, closed bool) []Point {
	m := len(v)
	off := func(i int, nv Point) Point {
		return Point{v[i].X + s*h*nv.X, v[i].Y + s*h*nv.Y}
	}
	var out []Point
	join := func(i int) {
		var pin, pout Point
		if closed {
			pin = n[(i-1+m)%m]
			pout = n[i%m]
		} else {
			pin = n[i-1]
			pout = n[i]
		}
		denom := 1 + pin.X*pout.X + pin.Y*pout.Y
		if denom < miterLimitDenom {
			out = append(out, off(i, pin), off(i, pout)) // bevel
			return
		}
		mx := (pin.X + pout.X) / denom
		my := (pin.Y + pout.Y) / denom
		out = append(out, Point{v[i].X + s*h*mx, v[i].Y + s*h*my})
	}
	if closed {
		for i := 0; i < m; i++ {
			join(i)
		}
		return out
	}
	out = append(out, off(0, n[0]))
	for i := 1; i < m-1; i++ {
		join(i)
	}
	out = append(out, off(m-1, n[m-2]))
	return out
}

// stitch 左链正向接右链反向，闭合成原始描边多边形（平头端帽）。
func stitch(left, right []Point) []Point {
	out := make([]Point, 0, len(left)+len(right))
	out = append(out, left...)
	for i := len(right) - 1; i >= 0; i-- {
		out = append(out, right[i])
	}
	return out
}

type isect struct {
	a, b   int // 两条边的段号
	ta, tb float64
	p      Point
}

// selfIntersections 求闭合多边形全部非相邻边的真交点。
func selfIntersections(poly []Point) []isect {
	n := len(poly)
	var out []isect
	for i := 0; i < n; i++ {
		a1, a2 := poly[i], poly[(i+1)%n]
		for j := i + 1; j < n; j++ {
			if j == i+1 || (i == 0 && j == n-1) {
				continue // 相邻边
			}
			b1, b2 := poly[j], poly[(j+1)%n]
			if t, u, p, ok := segCross(a1, a2, b1, b2); ok {
				out = append(out, isect{i, j, t, u, p})
			}
		}
	}
	return out
}

// segCross 两线段真交点（端点接触不算）。
func segCross(a1, a2, b1, b2 Point) (t, u float64, p Point, ok bool) {
	d1x, d1y := a2.X-a1.X, a2.Y-a1.Y
	d2x, d2y := b2.X-b1.X, b2.Y-b1.Y
	den := d1x*d2y - d1y*d2x
	if math.Abs(den) < strokeEps {
		return 0, 0, Point{}, false
	}
	rx, ry := b1.X-a1.X, b1.Y-a1.Y
	t = (rx*d2y - ry*d2x) / den
	u = (rx*d1y - ry*d1x) / den
	if t <= strokeEps || t >= 1-strokeEps || u <= strokeEps || u >= 1-strokeEps {
		return 0, 0, Point{}, false
	}
	return t, u, Point{a1.X + t*d1x, a1.Y + t*d1y}, true
}

// splitLoops 在自交处把闭合多边形切分为若干简单环：交点处换边对接（解缠），
// 返回全部环（含绕入环）与切分的交点个数。所有环的有向面积之和等于原多边形有向面积。
func splitLoops(poly []Point) ([][]Point, int) {
	xs := selfIntersections(poly)
	if len(xs) == 0 {
		return [][]Point{poly}, 0
	}
	n := len(poly)
	pos := make([]Point, n, n+len(xs))
	copy(pos, poly)
	for _, x := range xs {
		pos = append(pos, x.p)
	}
	type cut struct {
		t   float64
		nid int
	}
	cuts := make([][]cut, n)
	for k, x := range xs {
		nid := n + k
		cuts[x.a] = append(cuts[x.a], cut{x.ta, nid})
		cuts[x.b] = append(cuts[x.b], cut{x.tb, nid})
	}
	type sub struct{ from, to, edge int }
	var subs []sub
	inAt := map[int][]int{}
	outAt := map[int][]int{}
	addSub := func(from, to, edge int) {
		idx := len(subs)
		subs = append(subs, sub{from, to, edge})
		outAt[from] = append(outAt[from], idx)
		inAt[to] = append(inAt[to], idx)
	}
	for i := 0; i < n; i++ {
		cs := cuts[i]
		sort.Slice(cs, func(p, q int) bool { return cs[p].t < cs[q].t })
		prev := i
		for _, c := range cs {
			addSub(prev, c.nid, i)
			prev = c.nid
		}
		addSub(prev, (i+1)%n, i)
	}
	// 配对：原顶点处直通；交点处换边（a 边进接 b 边出），把绕进去的结解开。
	pair := make([]int, len(subs))
	for node, ins := range inAt {
		outs := outAt[node]
		if len(ins) == 1 && len(outs) == 1 {
			pair[ins[0]] = outs[0]
			continue
		}
		for _, in := range ins {
			for _, o := range outs {
				if subs[o].edge != subs[in].edge {
					pair[in] = o
				}
			}
		}
	}
	visited := make([]bool, len(subs))
	var loops [][]Point
	for s := range subs {
		if visited[s] {
			continue
		}
		var loop []Point
		for cur := s; !visited[cur]; cur = pair[cur] {
			visited[cur] = true
			loop = append(loop, pos[subs[cur].from])
		}
		loops = append(loops, loop)
	}
	return loops, len(xs)
}

// keepOuter 只保留最外环（|面积|最大者）。
// 开折线的加粗轮廓拓扑上必为单连通盘：解缠后其余环都绕进了外环内部，
// 正是要丢掉的结；闭合折线的每条链同理各自取最外环。
func keepOuter(loops [][]Point) (kept [][]Point, discarded int) {
	best, bestAbs := -1, 0.0
	for i, l := range loops {
		a, err := Area(l)
		if err != nil {
			a = 0
		}
		if math.Abs(a) > bestAbs {
			bestAbs, best = math.Abs(a), i
		}
	}
	if best < 0 {
		return nil, len(loops)
	}
	return loops[best : best+1 : best+1], len(loops) - 1
}

// checkSimple 校验环已无残留自交；若有，报出未拆分的两段段号。
func checkSimple(loop []Point) *StrokeError {
	if xs := selfIntersections(loop); len(xs) > 0 {
		return &StrokeError{Op: "自交未拆分", Seg: xs[0].a, SegB: xs[0].b}
	}
	return nil
}

func signedSum(loops [][]Point) float64 {
	var s float64
	for _, l := range loops {
		if a, err := Area(l); err == nil {
			s += a
		}
	}
	return s
}

func absSum(loops [][]Point) float64 {
	var s float64
	for _, l := range loops {
		if a, err := Area(l); err == nil {
			s += math.Abs(a)
		}
	}
	return s
}
