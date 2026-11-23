// Package geom 提供三次 Bezier 曲线的几何工具。
//
// 本文件新增弧长参数化: 先对每段曲线建立累计长度表, 再由"等弧长"目标反解
// 单调的重参数函数 t = R(s)。原有 Bezier/Flatten 接口保持不变。
package geom

import (
	"fmt"
	"math"
	"sort"
)

// Cubic 是一段三次 Bezier 曲线。
type Cubic struct{ P0, P1, P2, P3 Point }

// At 按曲线自身参数 t∈[0,1] 取点, 与 Bezier(p0,p1,p2,p3,t) 完全一致。
func (c Cubic) At(t float64) Point { return Bezier(c.P0, c.P1, c.P2, c.P3, t) }

// Deriv 返回一阶导 dC/dt。
func (c Cubic) Deriv(t float64) Point {
	u := 1 - t
	b0 := 3 * u * u
	b1 := 6 * u * t
	b2 := 3 * t * t
	return Point{
		X: b0*(c.P1.X-c.P0.X) + b1*(c.P2.X-c.P1.X) + b2*(c.P3.X-c.P2.X),
		Y: b0*(c.P1.Y-c.P0.Y) + b1*(c.P2.Y-c.P1.Y) + b2*(c.P3.Y-c.P2.Y),
	}
}

// Speed 返回参数 t 处的速率 |dC/dt|。
func (c Cubic) Speed(t float64) float64 {
	d := c.Deriv(t)
	return math.Hypot(d.X, d.Y)
}

func finitePoint(p Point) bool {
	return !math.IsNaN(p.X) && !math.IsNaN(p.Y) && !math.IsInf(p.X, 0) && !math.IsInf(p.Y, 0)
}

func (c Cubic) validate(idx int) error {
	for _, p := range [...]Point{c.P0, c.P1, c.P2, c.P3} {
		if !finitePoint(p) {
			return fmt.Errorf("段 %d: %w", idx, ErrDegenerate)
		}
	}
	// 四个控制点全部重合: 曲线不移动, 属于"长度为零"。
	if c.P0 == c.P1 && c.P1 == c.P2 && c.P2 == c.P3 {
		return fmt.Errorf("段 %d: %w", idx, ErrZeroLength)
	}
	// 端点处控制点重合 => 端点速率为零, 参数化退化 (切向/反解病态)。
	if c.P0 == c.P1 || c.P2 == c.P3 {
		return fmt.Errorf("段 %d: %w", idx, ErrDegenerate)
	}
	return nil
}

// CubicsFromPoints 把每段 4 个控制点转换成 Cubic。每段点数不足会带段号报错。
func CubicsFromPoints(ctrl [][]Point) ([]Cubic, error) {
	if len(ctrl) == 0 {
		return nil, fmt.Errorf("路径: %w", ErrTooFewPts)
	}
	out := make([]Cubic, len(ctrl))
	for i, pts := range ctrl {
		if len(pts) < 4 {
			return nil, fmt.Errorf("段 %d: 需要 4 个控制点, 实际 %d 个: %w", i, len(pts), ErrTooFewPts)
		}
		out[i] = Cubic{pts[0], pts[1], pts[2], pts[3]}
	}
	return out, nil
}

// ---------- 自适应 Gauss-Legendre 积分 ----------

var glNodes = [5]float64{
	0,
	-0.5384693101056831, 0.5384693101056831,
	-0.9061798459386640, 0.9061798459386640,
}
var glWeights = [5]float64{
	0.5688888888888889,
	0.4786286704993665, 0.4786286704993665,
	0.2369268850561891, 0.2369268850561891,
}

func gl5(f func(float64) float64, a, b float64) float64 {
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i := 0; i < 5; i++ {
		s += glWeights[i] * f(mid+half*glNodes[i])
	}
	return s * half
}

// integrate 在 [a,b] 上自适应积分, 绝对/相对容差均为 tol。
func integrate(f func(float64) float64, a, b, tol float64) float64 {
	return integrateRec(f, a, b, gl5(f, a, b), tol, 0)
}

func integrateRec(f func(float64) float64, a, b, whole, tol float64, depth int) float64 {
	mid := 0.5 * (a + b)
	left := gl5(f, a, mid)
	right := gl5(f, mid, b)
	est := left + right
	if depth >= 22 || math.Abs(est-whole) <= tol*(1+math.Abs(est)) {
		return est
	}
	return integrateRec(f, a, mid, left, tol, depth+1) +
		integrateRec(f, mid, b, right, tol, depth+1)
}

// ---------- 累计长度表 ----------

// ArcSegment 是一段曲线的累计长度表: T[k] 上升, L[k] = 从本段起点到 T[k] 的弧长。
type ArcSegment struct {
	Cubic  Cubic
	Length float64
	T, L   []float64
}

// ArcTable 是整条路径的累计长度表。
type ArcTable struct {
	Segments []ArcSegment
	Total    float64
}

const (
	nodeTol     = 1e-12 // 单段积分容差
	speedRatio  = 1.20  // 相邻节点速率比超过该值时加密
	maxNodeDep  = 20
	maxSegNodes = 8192
	minSegLen   = 1e-13
)

// BuildArcTable 逐段建立累计长度表。所有错误都带段号。
func BuildArcTable(segs []Cubic) (*ArcTable, error) {
	if len(segs) == 0 {
		return nil, fmt.Errorf("路径: %w", ErrTooFewPts)
	}
	tbl := &ArcTable{Segments: make([]ArcSegment, len(segs))}
	for i, c := range segs {
		if err := c.validate(i); err != nil {
			return nil, err
		}
		ts, ls := buildNodes(c)
		if len(ts) < 2 {
			return nil, fmt.Errorf("段 %d: %w", i, ErrZeroLength)
		}
		length := ls[len(ls)-1]
		if length <= minSegLen {
			return nil, fmt.Errorf("段 %d: %w", i, ErrZeroLength)
		}
		tbl.Segments[i] = ArcSegment{Cubic: c, Length: length, T: ts, L: ls}
		tbl.Total += length
	}
	return tbl, nil
}

// buildNodes 自适应细分参数区间, 在速率变化大的地方加密。
// 返回的 T 严格上升, L 严格上升 (对非退化曲线), L 为累计弧长。
func buildNodes(c Cubic) ([]float64, []float64) {
	ts := []float64{0}
	ls := []float64{0}
	var rec func(a, la, b float64, depth int)
	rec = func(a, la, b float64, depth int) {
		if len(ts) >= maxSegNodes {
			// 达到节点上限, 剩余区间一次积分后收尾。
			ts = append(ts, b)
			ls = append(ls, la+integrate(c.Speed, a, b, nodeTol))
			return
		}
		m := 0.5 * (a + b)
		lm := integrate(c.Speed, a, m, nodeTol)
		rm := integrate(c.Speed, m, b, nodeTol)
		sa, sm, sb := c.Speed(a), c.Speed(m), c.Speed(b)
		lo := math.Min(sa, math.Min(sm, sb))
		hi := math.Max(sa, math.Max(sm, sb))
		refine := depth < maxNodeDep && lo > 0 && hi > speedRatio*lo
		if !refine {
			// 同时用中点积分残差兜底。
			lab := integrate(c.Speed, a, b, nodeTol)
			if depth < maxNodeDep && math.Abs(lm+rm-lab) > 1e-9*math.Max(lab, 1e-15) {
				refine = true
			} else {
				ts = append(ts, b)
				ls = append(ls, la+lab)
				return
			}
		}
		rec(a, la, m, depth+1)
		rec(m, la+lm, b, depth+1)
	}
	rec(0, 0, 1, 0)
	return ts, ls
}

// LocalArc 返回段 seg 从起点到参数 t 的真实弧长 (自适应积分)。
func (t *ArcTable) LocalArc(seg int, param float64) float64 {
	s := &t.Segments[seg]
	param = clamp01(param)
	if param == 0 {
		return 0
	}
	if param == 1 {
		return s.Length
	}
	return integrate(s.Cubic.Speed, 0, param, nodeTol)
}

// ---------- 单调三次 Hermite 插值 ----------

// hermiteNodes 是 L(t) 的单调三次 Hermite 模型 (Fritsch-Carlson 限幅)。
type hermiteNodes struct {
	ts, ls, ms []float64
}

func newHermite(t *ArcSegment) *hermiteNodes {
	n := len(t.T)
	h := &hermiteNodes{ts: t.T, ls: t.L, ms: make([]float64, n)}
	for k := 0; k < n; k++ {
		// 端点及内部都用真实速率作为导数, 再限幅保证单调。
		h.ms[k] = t.Cubic.Speed(t.T[k])
	}
	for k := 0; k+1 < n; k++ {
		dt := h.ts[k+1] - h.ts[k]
		if dt <= 0 {
			continue
		}
		delta := (h.ls[k+1] - h.ls[k]) / dt
		if delta <= 0 {
			h.ms[k], h.ms[k+1] = 0, 0
			continue
		}
		a := h.ms[k] / delta
		b := h.ms[k+1] / delta
		if a < 0 {
			a, h.ms[k] = 0, 0
		}
		if b < 0 {
			b, h.ms[k+1] = 0, 0
		}
		if a*a+b*b > 9 {
			tau := 3 / math.Sqrt(a*a+b*b)
			h.ms[k] = tau * a * delta
			h.ms[k+1] = tau * b * delta
		}
	}
	return h
}

func (h *hermiteNodes) idx(t float64) int {
	// 返回满足 ts[k] <= t 的最大 k。
	k := sort.SearchFloat64s(h.ts, t)
	if k >= len(h.ts) {
		return len(h.ts) - 2
	}
	if h.ts[k] > t {
		k--
	}
	if k < 0 {
		k = 0
	}
	return k
}

// at 计算模型 L(t)。
func (h *hermiteNodes) at(t float64) float64 {
	if t <= h.ts[0] {
		return h.ls[0]
	}
	if t >= h.ts[len(h.ts)-1] {
		return h.ls[len(h.ts)-1]
	}
	k := h.idx(t)
	dt := h.ts[k+1] - h.ts[k]
	u := (t - h.ts[k]) / dt
	u2, u3 := u*u, u*u*u
	h00 := 2*u3 - 3*u2 + 1
	h10 := u3 - 2*u2 + u
	h01 := -2*u3 + 3*u2
	h11 := u3 - u2
	return h00*h.ls[k] + h10*dt*h.ms[k] + h01*h.ls[k+1] + h11*dt*h.ms[k+1]
}

// deriv 计算模型 L'(t)。
func (h *hermiteNodes) deriv(t float64) float64 {
	if t <= h.ts[0] || t >= h.ts[len(h.ts)-1] {
		return 0
	}
	k := h.idx(t)
	dt := h.ts[k+1] - h.ts[k]
	u := (t - h.ts[k]) / dt
	u2 := u * u
	d00 := 6*u2 - 6*u
	d10 := 3*u2 - 4*u + 1
	d01 := -6*u2 + 6*u
	d11 := 3*u2 - 2*u
	return (d00*h.ls[k]+d01*h.ls[k+1])/dt + d10*h.ms[k] + d11*h.ms[k+1]
}

// invert 反解 L(t)=a, 结果关于 a 严格单调。
func (h *hermiteNodes) invert(a float64) float64 {
	n := len(h.ts)
	if a <= h.ls[0] {
		return h.ts[0]
	}
	if a >= h.ls[n-1] {
		return h.ts[n-1]
	}
	k := sort.SearchFloat64s(h.ls, a)
	if k >= n {
		k = n - 1
	}
	if h.ls[k] > a {
		k--
	}
	if k < 0 {
		k = 0
	}
	lo, hi := h.ts[k], h.ts[k+1]
	t := lo + (hi-lo)*(a-h.ls[k])/(h.ls[k+1]-h.ls[k])
	for iter := 0; iter < 80; iter++ {
		f := h.at(t) - a
		if math.Abs(f) <= 1e-14*math.Max(1, math.Abs(a)) {
			return t
		}
		if f > 0 {
			hi = t
		} else {
			lo = t
		}
		d := h.deriv(t)
		next := t
		if d > 0 {
			next = t - f/d
		}
		if !(next > lo && next < hi) {
			next = 0.5 * (lo + hi)
		}
		if math.Abs(next-t) <= 1e-15 {
			return next
		}
		t = next
	}
	return t
}

// ---------- 重参数函数 ----------

type segReparam struct {
	c      Cubic
	length float64
	h      *hermiteNodes
}

// trueArc 返回段内从起点到参数 t 的真实弧长。
func (sg *segReparam) trueArc(t float64) float64 {
	if t <= 0 {
		return 0
	}
	if t >= 1 {
		return sg.length
	}
	return integrate(sg.c.Speed, 0, t, nodeTol)
}

// paramAt 反解真实弧长: 先由累计长度表的单调 Hermite 模型给出初值,
// 再对真实 L(t)-a=0 做带二分保护的 Newton 迭代。逆函数因此严格单调。
func (sg *segReparam) paramAt(a float64) float64 {
	if a <= 0 {
		return 0
	}
	if a >= sg.length {
		return 1
	}
	t := sg.h.invert(a)
	lo, hi := 0.0, 1.0
	for iter := 0; iter < 64; iter++ {
		f := sg.trueArc(t) - a
		if math.Abs(f) <= 1e-14*math.Max(1, math.Abs(a)) {
			return t
		}
		if f > 0 {
			hi = t
		} else {
			lo = t
		}
		d := sg.c.Speed(t)
		next := t
		if d > 0 {
			next = t - f/d
		}
		if !(next > lo && next < hi) {
			next = 0.5 * (lo + hi)
		}
		if math.Abs(next-t) <= 1e-16 {
			return next
		}
		t = next
	}
	return t
}

// Reparam 是整条路径的重参数函数: 归一化弧长 u∈[0,1] -> (段号, 段内参数 t)。
// 每段 t=R_i(s) 由累计长度表反解得到, 在 s 上严格单调。
type Reparam struct {
	segs   []segReparam
	prefix []float64
	total  float64
}

// BuildReparam 依次完成"建累计长度表 -> 逐段反解重参数函数"。
func BuildReparam(segs []Cubic) (*Reparam, error) {
	tbl, err := BuildArcTable(segs)
	if err != nil {
		return nil, err
	}
	return NewReparam(tbl)
}

// NewReparam 由已建好的累计长度表构造重参数函数。
func NewReparam(tbl *ArcTable) (*Reparam, error) {
	if tbl == nil || len(tbl.Segments) == 0 {
		return nil, fmt.Errorf("路径: %w", ErrTooFewPts)
	}
	r := &Reparam{segs: make([]segReparam, len(tbl.Segments)), prefix: make([]float64, len(tbl.Segments))}
	for i := range tbl.Segments {
		s := &tbl.Segments[i]
		if s.Length <= minSegLen {
			return nil, fmt.Errorf("段 %d: %w", i, ErrZeroLength)
		}
		r.prefix[i] = r.total
		r.segs[i] = segReparam{c: s.Cubic, length: s.Length, h: newHermite(s)}
		r.total += s.Length
	}
	return r, nil
}

// TotalLength 返回整条路径弧长。
func (r *Reparam) TotalLength() float64 { return r.total }

// NumSegments 返回段数。
func (r *Reparam) NumSegments() int { return len(r.segs) }

// SegmentLength 返回第 seg 段的弧长。
func (r *Reparam) SegmentLength(seg int) float64 { return r.segs[seg].length }

// SegmentAt 返回第 seg 段在其局部归一化弧长 s∈[0,1] 处的曲线参数 t=R_i(s)。
func (r *Reparam) SegmentAt(seg int, s float64) float64 {
	sg := &r.segs[seg]
	return sg.paramAt(clamp01(s) * sg.length)
}

// Eval 返回全局归一化弧长 u 对应的 (段号, 段内参数)。
// 等间隔的 u 对应等弧长。
func (r *Reparam) Eval(u float64) (seg int, t float64) {
	u = clamp01(u)
	if u == 0 {
		return 0, 0
	}
	if u >= 1 {
		last := len(r.segs) - 1
		return last, 1
	}
	target := u * r.total
	// 找到第一个累计末端 >= target 的段。
	seg = sort.Search(len(r.segs), func(i int) bool {
		return r.prefix[i]+r.segs[i].length >= target
	})
	if seg >= len(r.segs) {
		seg = len(r.segs) - 1
	}
	local := (target - r.prefix[seg]) / r.segs[seg].length
	return seg, r.SegmentAt(seg, local)
}

// PointAt 返回全局归一化弧长 u 处的点。
func (r *Reparam) PointAt(u float64) Point {
	seg, t := r.Eval(u)
	return r.segs[seg].c.At(t)
}

// GlobalArc 返回点 p(u) 的真实累计弧长 (用于核验)。
func (r *Reparam) GlobalArc(u float64) float64 {
	seg, t := r.Eval(u)
	return r.prefix[seg] + r.segs[seg].trueArc(t)
}

// ParamArc 返回第 seg 段参数 t 处整条路径的累计弧长 (按原参数取点的对照量)。
func (r *Reparam) ParamArc(seg int, t float64) float64 {
	return r.prefix[seg] + r.segs[seg].trueArc(clamp01(t))
}

// FlattenEqualArc 按等弧长取 n+1 个点 (含首尾)。旧的按参数接口不受影响。
func (r *Reparam) FlattenEqualArc(n int) ([]Point, error) {
	if n < 1 {
		return nil, ErrBadSegments
	}
	out := make([]Point, n+1)
	for i := 0; i <= n; i++ {
		out[i] = r.PointAt(float64(i) / float64(n))
	}
	return out, nil
}

// ---------- 核验 ----------

// SpeedDeviation 在 steps 个等 u 步长上比较真实弧长与目标步长,
// 返回带符号相对偏差的最小值、最大值与最大绝对值。
// 偏差 = 实际弧长/目标弧长 - 1。
func (r *Reparam) SpeedDeviation(steps int) (minDev, maxDev, maxAbs float64) {
	if steps < 1 {
		steps = 1
	}
	target := r.total / float64(steps)
	prev := 0.0
	minDev, maxDev = math.Inf(1), math.Inf(-1)
	for k := 1; k <= steps; k++ {
		cur := r.GlobalArc(float64(k) / float64(steps))
		d := (cur-prev)/target - 1
		if d < minDev {
			minDev = d
		}
		if d > maxDev {
			maxDev = d
		}
		if math.Abs(d) > maxAbs {
			maxAbs = math.Abs(d)
		}
		prev = cur
	}
	return
}

// CheckMonotone 逐点核验每段重参数函数严格单调, 并检查端点 R(0)=0, R(1)=1。
// 返回第一个不满足的位置; 全部通过返回 nil。
func (r *Reparam) CheckMonotone(steps int) error {
	if steps < 2 {
		steps = 2
	}
	for i := range r.segs {
		if t := r.SegmentAt(i, 0); t != 0 {
			return fmt.Errorf("段 %d: R(0)=%g, 应为 0", i, t)
		}
		if t := r.SegmentAt(i, 1); math.Abs(t-1) > 1e-9 {
			return fmt.Errorf("段 %d: R(1)=%g, 应为 1", i, t)
		}
		prev := -1.0
		for k := 0; k <= steps; k++ {
			s := float64(k) / float64(steps)
			t := r.SegmentAt(i, s)
			if t <= prev {
				return fmt.Errorf("段 %d: R 在 s=%g 处非严格单调 (t=%g <= %g)", i, s, t, prev)
			}
			prev = t
		}
	}
	return nil
}

func clamp01(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x > 1 {
		return 1
	}
	return x
}
