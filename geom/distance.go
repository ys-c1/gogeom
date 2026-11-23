package geom

import (
	"fmt"
	"math"
)

// Cubic 是一段三次贝塞尔曲线。
type Cubic struct{ P0, P1, P2, P3 Point }

// At 返回参数 t∈[0,1] 处的点。
func (c Cubic) At(t float64) Point { return Bezier(c.P0, c.P1, c.P2, c.P3, t) }

// Deriv 返回参数 t 处的一阶导数（切向量）。
func (c Cubic) Deriv(t float64) Point {
	u := 1 - t
	a, b, k := 3*u*u, 6*u*t, 3*t*t
	return Point{
		a*(c.P1.X-c.P0.X) + b*(c.P2.X-c.P1.X) + k*(c.P3.X-c.P2.X),
		a*(c.P1.Y-c.P0.Y) + b*(c.P2.Y-c.P1.Y) + k*(c.P3.Y-c.P2.Y),
	}
}

// degenerate 判定：四个控制点两两距离都小于 1e-12，曲线缩成一个点，
// 距离函数失去参数化意义，一律拦下。
func (c Cubic) degenerate() bool {
	pts := [4]Point{c.P0, c.P1, c.P2, c.P3}
	for i := 0; i < 4; i++ {
		for j := i + 1; j < 4; j++ {
			if math.Hypot(pts[i].X-pts[j].X, pts[i].Y-pts[j].Y) >= 1e-12 {
				return false
			}
		}
	}
	return true
}

// DistanceResult 是两条曲线间距离的结果。
type DistanceResult struct {
	Dist     float64 // 最近距离
	PA, PB   Point   // 最近点对（PA 在曲线 A 上，PB 在曲线 B 上）
	S, T     float64 // 最近点各自所在段的参数
	SegA     int     // 最近点在曲线 A 的第几段
	SegB     int     // 最近点在曲线 B 的第几段
	Steps    int     // 交替投影收敛所用步数
}

// ---- 点到曲线：距离写成参数函数，求导解驻点 ----
//
// 点 p 到曲线 C 的距离平方为 f(s) = |C(s)-p|²，s∈[0,1]。
// 驻点满足 f'(s) = 2(C(s)-p)·C'(s) = 0，即 g(s) = (C(s)-p)·C'(s) = 0。
// C(s)-p 是三次、C'(s) 是二次，g(s) 是五次多项式。
// 直接在 Bernstein 基下做二项加权卷积得到 g 的五次 Bernstein 系数，
// 再用凸包性质递归二分隔离全部实根，最后比较各驻点与端点 s=0,1 的 f 值取最小。

var (
	bin3 = [4]float64{1, 3, 3, 1}
	bin2 = [3]float64{1, 2, 1}
	bin5 = [6]float64{1, 5, 10, 10, 5, 1}
)

// distDerivBernstein 返回 g(s) = (C(s)-p)·C'(s) 的五次 Bernstein 系数。
func distDerivBernstein(p Point, c Cubic) [6]float64 {
	f := [4]Point{
		{c.P0.X - p.X, c.P0.Y - p.Y},
		{c.P1.X - p.X, c.P1.Y - p.Y},
		{c.P2.X - p.X, c.P2.Y - p.Y},
		{c.P3.X - p.X, c.P3.Y - p.Y},
	}
	g := [3]Point{
		{3 * (c.P1.X - c.P0.X), 3 * (c.P1.Y - c.P0.Y)},
		{3 * (c.P2.X - c.P1.X), 3 * (c.P2.Y - c.P1.Y)},
		{3 * (c.P3.X - c.P2.X), 3 * (c.P3.Y - c.P2.Y)},
	}
	var h [6]float64
	for k := 0; k <= 5; k++ {
		for i := 0; i <= 3; i++ {
			j := k - i
			if j < 0 || j > 2 {
				continue
			}
			w := bin3[i] * bin2[j] / bin5[k]
			h[k] += w * (f[i].X*g[j].X + f[i].Y*g[j].Y)
		}
	}
	return h
}

// splitBernstein 在区间中点做 de Casteljau 分裂。
func splitBernstein(c []float64) (left, right []float64) {
	n := len(c)
	left, right = make([]float64, n), make([]float64, n)
	tmp := append([]float64(nil), c...)
	for k := 0; k < n; k++ {
		left[k] = tmp[0]
		right[n-1-k] = tmp[n-1-k]
		for i := 0; i+1 < len(tmp); i++ {
			tmp[i] = 0.5 * (tmp[i] + tmp[i+1])
		}
		tmp = tmp[:len(tmp)-1]
	}
	return left, right
}

// isolateRoots 用 Bernstein 系数的凸包性质隔离 (a,b) 内的所有实根。
// 系数全正或全负则区间内无根；否则二分，直到区间宽度小于 1e-12。
func isolateRoots(h []float64, a, b float64, out *[]float64) {
	allPos, allNeg := true, true
	for _, v := range h {
		if v <= 0 {
			allPos = false
		}
		if v >= 0 {
			allNeg = false
		}
	}
	if allPos || allNeg {
		return
	}
	if b-a < 1e-12 {
		*out = append(*out, 0.5*(a+b))
		return
	}
	l, r := splitBernstein(h)
	m := 0.5 * (a + b)
	isolateRoots(l, a, m, out)
	isolateRoots(r, m, b, out)
}

// pointSegDist 返回点 p 到曲线段 c 的最短距离平方与对应参数。
// 候选集 = g(s) 在 (0,1) 内的全部驻点 ∪ 端点 {0,1}。
func pointSegDist(p Point, c Cubic) (d2, s float64) {
	h := distDerivBernstein(p, c)
	var roots []float64
	isolateRoots(h[:], 0, 1, &roots)
	cands := append(roots, 0, 1)
	d2, s = math.Inf(1), 0
	for _, v := range cands {
		if v < 0 {
			v = 0
		} else if v > 1 {
			v = 1
		}
		q := c.At(v)
		d := (q.X-p.X)*(q.X-p.X) + (q.Y-p.Y)*(q.Y-p.Y)
		if d < d2 {
			d2, s = d, v
		}
	}
	return d2, s
}

// ---- 段到段：交替投影 ----
//
// 在 (s,t) 上交替固定一个参数、对另一个参数做全局投影：
//   s_{k+1} = argmin_s |A(s) - B(t_k)|
//   t_{k+1} = argmin_t |A(s_{k+1}) - B(t)|
// 每个半步都取全局最小，所以距离序列单调不增且有下界 0，必然收敛。
//
// 收敛判据（自定）：同时满足
//   1) 距离相对下降量 prev-d <= 1e-12*max(1,prev)
//   2) 参数步长 max(|Δs|,|Δt|) <= 1e-9
// 为什么够用：单调不增序列的相邻下降量趋于 0 是收敛的必要条件；
// 再加上参数本身也停止移动，说明 (s,t) 已落在驻点（或边界）上，
// 继续迭代只能在该驻点的浮点误差内抖动，不可能再获得超过容差的改进。
// 另设 d <= 1e-13 直接判交、200 步上限兜底。
func segDist(a, b Cubic) (d, s, t float64, hist []float64) {
	// 粗网格播种，避免交替投影落入非全局的局部极小。
	s, t = 0.5, 0.5
	best := math.Inf(1)
	for i := 0; i <= 8; i++ {
		for j := 0; j <= 8; j++ {
			si, tj := float64(i)/8, float64(j)/8
			pa, pb := a.At(si), b.At(tj)
			dd := math.Hypot(pa.X-pb.X, pa.Y-pb.Y)
			if dd < best {
				best, s, t = dd, si, tj
			}
		}
	}
	prev := math.Inf(1)
	ps, pt := s, t
	for iter := 0; iter < 200; iter++ {
		_, s = pointSegDist(b.At(t), a)
		d2, t2 := pointSegDist(a.At(s), b)
		t = t2
		d = math.Sqrt(d2)
		hist = append(hist, d)
		if d <= 1e-13 {
			break
		}
		if iter > 0 && prev-d <= 1e-12*math.Max(1, prev) &&
			math.Max(math.Abs(s-ps), math.Abs(t-pt)) <= 1e-9 {
			break
		}
		prev, ps, pt = d, s, t
	}
	return d, s, t, hist
}

// pointCurveDist 返回点 p 到整条曲线（多段）的最短距离。
func pointCurveDist(p Point, curve []Cubic) float64 {
	best := math.Inf(1)
	for _, c := range curve {
		if d2, _ := pointSegDist(p, c); d2 < best {
			best = d2
		}
	}
	return math.Sqrt(best)
}

// CurveDistance 计算两条（可多段的）三次贝塞尔曲线之间的最近距离。
// 曲线重合、控制点退化、段数为零三种情形一律报错，并指出是第几段。
func CurveDistance(a, b []Cubic) (DistanceResult, error) {
	if len(a) == 0 {
		return DistanceResult{}, fmt.Errorf("曲线A段数为0: %w", ErrBadSegments)
	}
	if len(b) == 0 {
		return DistanceResult{}, fmt.Errorf("曲线B段数为0: %w", ErrBadSegments)
	}
	for i, c := range a {
		if c.degenerate() {
			return DistanceResult{}, fmt.Errorf("曲线A第%d段控制点退化: %w", i, ErrDegenerate)
		}
	}
	for j, c := range b {
		if c.degenerate() {
			return DistanceResult{}, fmt.Errorf("曲线B第%d段控制点退化: %w", j, ErrDegenerate)
		}
	}

	var res DistanceResult
	res.Dist = math.Inf(1)
	for i, ca := range a {
		for j, cb := range b {
			d, s, t, hist := segDist(ca, cb)
			if d < res.Dist {
				res = DistanceResult{
					Dist: d, S: s, T: t, SegA: i, SegB: j, Steps: len(hist),
					PA: ca.At(s), PB: cb.At(t),
				}
			}
		}
	}

	// 重合判定：最近距离为零时，沿 A 的最近段均匀取 7 个内点，
	// 若它们到 B 的距离也全部为零，则两段有公共弧而非孤立交点，拦下。
	if res.Dist < 1e-9 {
		ca := a[res.SegA]
		coincident := true
		for k := 1; k <= 7; k++ {
			p := ca.At(float64(k) / 8)
			if pointCurveDist(p, b) > 1e-9 {
				coincident = false
				break
			}
		}
		if coincident {
			return DistanceResult{}, fmt.Errorf(
				"曲线A第%d段与曲线B第%d段重合: %w", res.SegA, res.SegB, ErrCoincident)
		}
	}
	return res, nil
}
