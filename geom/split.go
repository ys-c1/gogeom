package geom

import (
	"fmt"
	"math"
)

// Polygon 表示可带孔的多边形：Outer 为外环，Holes 为孔环。
// 各环方向任意，SplitEqualArea 内部会统一规范化
// （外环逆时针、孔环顺时针）。
type Polygon struct {
	Outer []Point
	Holes [][]Point
}

// Piece 为切分得到的一块：Area 为其面积，Rings 为构成它的环
// （首个为外环，其余为孔环；块落在原多边形零宽部分时可能为空）。
type Piece struct {
	Area  float64
	Rings [][]Point
}

// SplitResult 为等面积切分结果。
// 切线为一组互相平行的竖直线 x = Cuts[i]，因此两两不相交。
type SplitResult struct {
	Parts   int       // 目标块数
	Cuts    []float64 // 各切线位置（x 坐标），len == Parts-1，严格递增
	Pieces  []Piece   // 各块，len == Parts
	Total   float64   // 原多边形面积（外环减孔）
	Spread  float64   // 各块面积极差（最大 - 最小）
	SumDiff float64   // 各块面积之和 - 原面积
}

// SelfIntersectionError 报告多边形自交（或环与环相交）的位置：
// 环 RingA 的边 EdgeA 与环 RingB 的边 EdgeB 相交。
// 环编号 0 为外环，1 起为各孔环；边号即起点下标。
type SelfIntersectionError struct {
	RingA, EdgeA int
	RingB, EdgeB int
}

func (e *SelfIntersectionError) Error() string {
	if e.RingA == e.RingB {
		return fmt.Sprintf("多边形自交：环 %d 的边 %d 与边 %d 相交", e.RingA, e.EdgeA, e.EdgeB)
	}
	return fmt.Sprintf("多边形自交：环 %d 的边 %d 与环 %d 的边 %d 相交", e.RingA, e.EdgeA, e.RingB, e.EdgeB)
}

// SplitEqualArea 把多边形（可带孔）等面积切成 n 块。
// 对累计面积函数 F(c) = “x <= c 部分的面积” 逐条二分反解切线位置，
// 使 F(Cuts[i]) = (i+1)/n * 总面积；切线互相平行，保证不相交。
// 块数为零、面积为零、多边形自交分别返回
// ErrZeroParts、ErrZeroArea、*SelfIntersectionError（含边号）。
func SplitEqualArea(p Polygon, n int) (*SplitResult, error) {
	if n <= 0 {
		return nil, ErrZeroParts
	}
	if len(p.Outer) < 3 {
		return nil, ErrTooFewPts
	}
	if err := validateTopology(p); err != nil {
		return nil, err
	}
	q := Polygon{Outer: ensureOrientation(p.Outer, true)}
	for _, h := range p.Holes {
		q.Holes = append(q.Holes, ensureOrientation(h, false))
	}
	total := q.areaLeftOf(math.Inf(1))
	if total == 0 {
		return nil, ErrZeroArea
	}
	if total < 0 {
		return nil, ErrDegenerate // 孔面积超过外环
	}

	minX, maxX := q.Outer[0].X, q.Outer[0].X
	for _, pt := range q.Outer[1:] {
		minX = math.Min(minX, pt.X)
		maxX = math.Max(maxX, pt.X)
	}

	// 反解每条切线：F(c) = i/n * total
	cuts := make([]float64, 0, n-1)
	for i := 1; i < n; i++ {
		target := total * float64(i) / float64(n)
		lo, hi := minX, maxX
		for iter := 0; iter < 100; iter++ {
			mid := 0.5 * (lo + hi)
			if q.areaLeftOf(mid) < target {
				lo = mid
			} else {
				hi = mid
			}
		}
		cuts = append(cuts, 0.5*(lo+hi))
	}

	// 逐块裁剪并量算面积
	pieces := make([]Piece, n)
	for i := 0; i < n; i++ {
		lo, hi := minX, maxX
		if i > 0 {
			lo = cuts[i-1]
		}
		if i < n-1 {
			hi = cuts[i]
		}
		var rings [][]Point
		var area float64
		if r := clipStrip(q.Outer, lo, hi); len(r) >= 3 {
			rings = append(rings, r)
			area += signedArea(r)
		}
		for _, h := range q.Holes {
			if r := clipStrip(h, lo, hi); len(r) >= 3 {
				rings = append(rings, r)
				area += signedArea(r) // 孔环为顺时针，面积为负
			}
		}
		pieces[i] = Piece{Area: area, Rings: rings}
	}

	res := &SplitResult{Parts: n, Cuts: cuts, Pieces: pieces, Total: total}
	mn, mx := math.Inf(1), math.Inf(-1)
	var sum float64
	for _, pc := range pieces {
		sum += pc.Area
		mn = math.Min(mn, pc.Area)
		mx = math.Max(mx, pc.Area)
	}
	res.Spread = mx - mn
	res.SumDiff = sum - total
	return res, nil
}

// Verify 逐块核验切分结果：
//  1. 每块面积与目标值（Total/Parts）的偏差不超过 tol；
//  2. 各块面积极差不超过 2*tol（上界）；
//  3. 各块面积之和与原面积之差不超过 tol（守恒）；
//  4. 切线位置严格递增（平行竖线，互不相交）。
func (r *SplitResult) Verify(tol float64) error {
	target := r.Total / float64(r.Parts)
	for i, pc := range r.Pieces {
		if d := math.Abs(pc.Area - target); d > tol {
			return fmt.Errorf("第 %d 块面积 %.10g 偏离目标 %.10g 达 %.3g，超过上界 %.3g", i, pc.Area, target, d, tol)
		}
	}
	if r.Spread > 2*tol {
		return fmt.Errorf("各块面积极差 %.3g 超过上界 %.3g", r.Spread, 2*tol)
	}
	if d := math.Abs(r.SumDiff); d > tol {
		return fmt.Errorf("各块面积之和与原面积之差 %.3g 超过上界 %.3g", d, tol)
	}
	for i := 1; i < len(r.Cuts); i++ {
		if !(r.Cuts[i] > r.Cuts[i-1]) {
			return fmt.Errorf("切线 %d 与 %d 位置重合或次序颠倒，切线可能相交", i-1, i)
		}
	}
	return nil
}

// areaLeftOf 返回多边形区域在竖直半平面 x <= c 内的面积。
func (p Polygon) areaLeftOf(c float64) float64 {
	a := signedArea(clipX(p.Outer, c, true))
	for _, h := range p.Holes {
		a += signedArea(clipX(h, c, true))
	}
	return a
}

// signedArea 为鞋带公式有符号面积（逆时针为正）。
func signedArea(r []Point) float64 {
	var s float64
	for i := range r {
		a, b := r[i], r[(i+1)%len(r)]
		s += a.X*b.Y - b.X*a.Y
	}
	return s / 2
}

// clipX 用半平面裁剪环：keepLE 为真保留 x <= c，否则保留 x >= c。
// Sutherland–Hodgman 算法，保持顶点顺序与环方向。
func clipX(ring []Point, c float64, keepLE bool) []Point {
	n := len(ring)
	if n == 0 {
		return nil
	}
	inside := func(p Point) bool {
		if keepLE {
			return p.X <= c
		}
		return p.X >= c
	}
	var out []Point
	for i := 0; i < n; i++ {
		a, b := ring[i], ring[(i+1)%n]
		ain, bin := inside(a), inside(b)
		if ain {
			out = append(out, a)
		}
		if ain != bin {
			t := (c - a.X) / (b.X - a.X)
			out = append(out, Point{X: c, Y: a.Y + t*(b.Y-a.Y)})
		}
	}
	return out
}

// clipStrip 把环裁剪到竖条 lo <= x <= hi。
func clipStrip(ring []Point, lo, hi float64) []Point {
	return clipX(clipX(ring, hi, true), lo, false)
}

// ensureOrientation 调整环方向：ccw 为真保证逆时针，否则保证顺时针。
func ensureOrientation(r []Point, ccw bool) []Point {
	if len(r) < 3 {
		return r
	}
	a := signedArea(r)
	if a != 0 && (a > 0) != ccw {
		out := make([]Point, len(r))
		for i, p := range r {
			out[len(r)-1-i] = p
		}
		return out
	}
	return r
}

// cross 为向量 (a-o) x (b-o) 的 z 分量。
func cross(o, a, b Point) float64 {
	return (a.X-o.X)*(b.Y-o.Y) - (a.Y-o.Y)*(b.X-o.X)
}

// edgesCross 判断线段 p1p2 与 p3p4 是否严格相交（不含共端点）。
func edgesCross(p1, p2, p3, p4 Point) bool {
	return cross(p1, p2, p3)*cross(p1, p2, p4) < 0 &&
		cross(p3, p4, p1)*cross(p3, p4, p2) < 0
}

// validateTopology 检查每个环是否自交、环与环是否相交，
// 发现即返回含环号与边号的 *SelfIntersectionError。
func validateTopology(p Polygon) error {
	rings := make([][]Point, 0, 1+len(p.Holes))
	rings = append(rings, p.Outer)
	rings = append(rings, p.Holes...)
	for ri, r := range rings {
		n := len(r)
		for i := 0; i < n; i++ {
			for j := i + 2; j < n; j++ {
				if i == 0 && j == n-1 {
					continue // 首尾边相邻
				}
				if edgesCross(r[i], r[(i+1)%n], r[j], r[(j+1)%n]) {
					return &SelfIntersectionError{RingA: ri, EdgeA: i, RingB: ri, EdgeB: j}
				}
			}
		}
	}
	for a := 0; a+1 < len(rings); a++ {
		for b := a + 1; b < len(rings); b++ {
			ra, rb := rings[a], rings[b]
			for i := 0; i < len(ra); i++ {
				for j := 0; j < len(rb); j++ {
					if edgesCross(ra[i], ra[(i+1)%len(ra)], rb[j], rb[(j+1)%len(rb)]) {
						return &SelfIntersectionError{RingA: a, EdgeA: i, RingB: b, EdgeB: j}
					}
				}
			}
		}
	}
	return nil
}
