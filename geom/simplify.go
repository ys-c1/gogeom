package geom

import (
	"math"
	"time"
)

// IndexRange 表示闭区间的点号范围 [Start, End]。
type IndexRange struct {
	Start int
	End   int
}

// Len 返回区间内的点号个数。
func (r IndexRange) Len() int { return r.End - r.Start + 1 }

// Contains 判断点号 i 是否落在区间内。
func (r IndexRange) Contains(i int) bool { return i >= r.Start && i <= r.End }

// DeletedRanges 根据保留点号（须严格递增）计算被删点号的连续区间。
// 返回的区间互不重叠、按点号升序排列，且每个区间的长度至少为 1。
func DeletedRanges(n int, kept []int) []IndexRange {
	var out []IndexRange
	prev := -1
	for _, k := range kept {
		if k <= prev || k >= n {
			continue
		}
		if k > prev+1 {
			out = append(out, IndexRange{Start: prev + 1, End: k - 1})
		}
		prev = k
	}
	if prev < n-1 {
		out = append(out, IndexRange{Start: prev + 1, End: n - 1})
	}
	return out
}

// Result 是一次路径简化的完整结果。Points 为简化后的新切片，输入不会被修改。
type Result struct {
	// Points 是简化后的路径；首尾点与原始路径相同（闭合环首尾值也相同）。
	Points []Point
	// Kept 是保留下来的原始点号，严格递增。
	Kept []int
	// Deleted 是被删除点号的连续区间，互不重叠、可逐个点号查询。
	Deleted []IndexRange

	Tolerance      float64
	OriginalPoints int
	RetainedPoints int
	// RetentionRatio 是保留点数 / 原始点数。
	RetentionRatio float64
	// MaxDeviation 是原始路径每个点到简化路径的最小距离的最大值。
	MaxDeviation float64
	// WorstPoint 是取得 MaxDeviation 的原始点号。
	WorstPoint int
	// Closed 表示输入是否为闭合环（首点与末点完全相同）。
	Closed bool
	// ClosureGap 是简化后首尾点的直线距离；闭合环恒为 0。
	ClosureGap float64
	// Elapsed 是抽稀本身的耗时（不含逐点核对）。
	Elapsed time.Duration
	// VerifyElapsed 是逐点核对偏差的耗时。
	VerifyElapsed time.Duration
}

// DeletedContains 判断原始点号 i 是否被删除。
func (r *Result) DeletedContains(i int) bool {
	for _, rg := range r.Deleted {
		if rg.Contains(i) {
			return true
		}
	}
	return false
}

// IsClosed 判断路径是否为首尾点完全相同的闭合环。
func IsClosed(pts []Point) bool {
	return len(pts) >= 4 && pts[0] == pts[len(pts)-1]
}

// Simplify 用 Douglas–Peucker 按最大偏差抽稀路径，容差可配。
//
// 输入必须至少 3 个点、容差非负且不能全共线。首尾点始终保留；
// 若输入是闭合环（首点与末点完全相同），简化结果的首尾点也完全相同。
// 抽稀后会把每个原始点与简化路径逐点核对，最大偏差不得超过容差上界，
// 否则返回带点号的 ErrToleranceExceeded。
func Simplify(pts []Point, tolerance float64) (*Result, error) {
	if math.IsNaN(tolerance) || tolerance < 0 {
		return nil, &PointError{Err: ErrNegativeTolerance, Index: 0}
	}
	if len(pts) < 3 {
		return nil, &PointError{Err: ErrTooFewPts, Index: len(pts)}
	}
	if err := checkCollinear(pts); err != nil {
		return nil, err
	}

	start := time.Now()
	n := len(pts)
	keep := make([]bool, n)
	keep[0], keep[n-1] = true, true
	closed := IsClosed(pts)

	if closed {
		// 首尾重合时不能直接把整圈交给 DP（两端是同一个点，弦退化为点）。
		// 选离首点最远的点把环切成两条开链，分别抽稀后再拼回。
		k := farthestFrom(pts, 0, 1, n-2)
		keep[k] = true
		dpMark(pts, 0, k, tolerance, keep)
		dpMark(pts, k, n-1, tolerance, keep)
	} else {
		dpMark(pts, 0, n-1, tolerance, keep)
	}

	kept := make([]int, 0, n)
	for i, ok := range keep {
		if ok {
			kept = append(kept, i)
		}
	}
	out := make([]Point, len(kept))
	for i, idx := range kept {
		out[i] = pts[idx]
	}
	elapsed := time.Since(start)

	verifyStart := time.Now()
	maxDev, worst := MaxDeviation(pts, out)
	verifyElapsed := time.Since(verifyStart)

	// DP 在精确算术下保证偏差不超过容差；这里留出浮点舍入余量。
	if maxDev > tolerance+slack(pts, tolerance) {
		return nil, &PointError{Err: ErrToleranceExceeded, Index: worst}
	}

	res := &Result{
		Points:         out,
		Kept:           kept,
		Deleted:        DeletedRanges(n, kept),
		Tolerance:      tolerance,
		OriginalPoints: n,
		RetainedPoints: len(out),
		RetentionRatio: float64(len(out)) / float64(n),
		MaxDeviation:   maxDev,
		WorstPoint:     worst,
		Closed:         closed,
		Elapsed:        elapsed,
		VerifyElapsed:  verifyElapsed,
	}
	if closed {
		res.ClosureGap = math.Hypot(out[0].X-out[len(out)-1].X, out[0].Y-out[len(out)-1].Y)
	}
	return res, nil
}

// dpMark 用显式栈做 Douglas–Peucker，把需要保留的中间点号标记到 keep。
// 端点 first、last 必须已标记为保留。
func dpMark(pts []Point, first, last int, tolerance float64, keep []bool) {
	type span struct{ a, b int }
	stack := []span{{first, last}}
	for len(stack) > 0 {
		s := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if s.b-s.a < 2 {
			continue
		}
		maxD, idx := -1.0, -1
		for i := s.a + 1; i < s.b; i++ {
			d := pointSegmentDist(pts[i], pts[s.a], pts[s.b])
			if d > maxD {
				maxD, idx = d, i
			}
		}
		if idx >= 0 && maxD > tolerance {
			keep[idx] = true
			stack = append(stack, span{s.a, idx}, span{idx, s.b})
		}
	}
}

// MaxDeviation 返回原始路径逐点对照简化路径的最大最小距离，以及取得该值的原始点号。
// 简化路径为空或退化时返回 +Inf（没有可对照的路径）。
func MaxDeviation(orig, simplified []Point) (float64, int) {
	if len(simplified) == 0 {
		return math.Inf(1), 0
	}
	maxD, worst := 0.0, 0
	for i, p := range orig {
		d := distToPolyline(p, simplified)
		if d > maxD {
			maxD, worst = d, i
		}
	}
	return maxD, worst
}

// Verify 逐点核对简化路径与原始路径，最大偏差不得超上界。
// 超限时返回带点号的 ErrToleranceExceeded。
func Verify(orig, simplified []Point, bound float64) (float64, error) {
	if math.IsNaN(bound) || bound < 0 {
		return 0, &PointError{Err: ErrNegativeTolerance, Index: 0}
	}
	if len(orig) < 3 {
		return 0, &PointError{Err: ErrTooFewPts, Index: len(orig)}
	}
	if len(simplified) < 2 {
		return 0, &PointError{Err: ErrTooFewPts, Index: len(simplified)}
	}
	maxD, worst := MaxDeviation(orig, simplified)
	if maxD > bound+slack(orig, bound) {
		return maxD, &PointError{Err: ErrToleranceExceeded, Index: worst}
	}
	return maxD, nil
}

// slack 是逐点核对时允许的浮点舍入余量，随坐标量级放大。
func slack(pts []Point, tolerance float64) float64 {
	return 1e-12 * (1 + tolerance + coordScale(pts))
}

func distToPolyline(p Point, line []Point) float64 {
	if len(line) == 1 {
		return math.Hypot(p.X-line[0].X, p.Y-line[0].Y)
	}
	best := math.Inf(1)
	for i := 0; i+1 < len(line); i++ {
		if d := pointSegmentDist(p, line[i], line[i+1]); d < best {
			best = d
			if best == 0 {
				return 0
			}
		}
	}
	return best
}

func pointSegmentDist(p, a, b Point) float64 {
	dx, dy := b.X-a.X, b.Y-a.Y
	if dx == 0 && dy == 0 {
		return math.Hypot(p.X-a.X, p.Y-a.Y)
	}
	t := ((p.X-a.X)*dx + (p.Y-a.Y)*dy) / (dx*dx + dy*dy)
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	return math.Hypot(p.X-(a.X+t*dx), p.Y-(a.Y+t*dy))
}

// checkCollinear 检查所有点是否共线（含全部重合的情况）。
func checkCollinear(pts []Point) error {
	p0 := pts[0]
	second := Point{}
	found := false
	for _, p := range pts[1:] {
		if p != p0 {
			second, found = p, true
			break
		}
	}
	if !found {
		return &PointError{Err: ErrCollinear, Index: 0}
	}
	dx, dy := second.X-p0.X, second.Y-p0.Y
	span := math.Hypot(dx, dy)
	eps := 1e-12 * (1 + coordScale(pts))
	for _, p := range pts {
		// 点到直线 p0-second 的距离。
		cross := math.Abs((p.X-p0.X)*dy - (p.Y-p0.Y)*dx)
		if cross/span > eps {
			return nil
		}
	}
	return &PointError{Err: ErrCollinear, Index: 0}
}

func coordScale(pts []Point) float64 {
	var s float64
	for _, p := range pts {
		s = math.Max(s, math.Max(math.Abs(p.X), math.Abs(p.Y)))
	}
	return s
}

func farthestFrom(pts []Point, from, lo, hi int) int {
	best, bestD := from, -1.0
	for i := lo; i <= hi; i++ {
		d := math.Hypot(pts[i].X-pts[from].X, pts[i].Y-pts[from].Y)
		if d > bestD {
			best, bestD = i, d
		}
	}
	return best
}
