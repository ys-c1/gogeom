package geom

import "math"

// Extremum 描述曲率的一个局部极值点。
type Extremum struct {
	T    float64 // 参数位置，位于 (0,1)
	K    float64 // 有符号曲率；尖点处为 +Inf（曲率不存在，|κ|→∞）
	Max  bool    // true 为极大值，false 为极小值
	Cusp bool    // 是否为尖点（速度为零、曲率导数不存在的位置）
}

// minExtremaSegments 是 CurvatureExtrema 接受的最小采样段数。
const minExtremaSegments = 8

// derivCoeffs 返回三次贝塞尔一阶导数的多项式系数：
// B'(t) = d0 + d1·t + d2·t²。
func derivCoeffs(p0, p1, p2, p3 Point) (d0, d1, d2 Point) {
	a := Point{p1.X - p0.X, p1.Y - p0.Y}
	b := Point{p2.X - p1.X, p2.Y - p1.Y}
	c := Point{p3.X - p2.X, p3.Y - p2.Y}
	d0 = Point{3 * a.X, 3 * a.Y}
	d1 = Point{6 * (b.X - a.X), 6 * (b.Y - a.Y)}
	d2 = Point{3 * (a.X - 2*b.X + c.X), 3 * (a.Y - 2*b.Y + c.Y)}
	return
}

// derivs 返回 t 处的一阶、二阶导数。三阶导数为常数 2·d2。
func derivs(d0, d1, d2 Point, t float64) (v1, v2 Point) {
	v1 = Point{d0.X + d1.X*t + d2.X*t*t, d0.Y + d1.Y*t + d2.Y*t*t}
	v2 = Point{d1.X + 2*d2.X*t, d1.Y + 2*d2.Y*t}
	return
}

// Curvature 返回参数 t 处的有符号曲率
//
//	κ(t) = (x'y'' − y'x'') / (x'² + y'²)^{3/2}
//
// 速度为零处（尖点）曲率不存在，返回 NaN。
func Curvature(p0, p1, p2, p3 Point, t float64) float64 {
	d0, d1, d2 := derivCoeffs(p0, p1, p2, p3)
	v1, v2 := derivs(d0, d1, d2, t)
	s := v1.X*v1.X + v1.Y*v1.Y
	if s == 0 {
		return math.NaN()
	}
	return (v1.X*v2.Y - v1.Y*v2.X) / (s * math.Sqrt(s))
}

// kappaStationary 计算 κ'(t)=0 的分子（去掉恒正因子 S^{5/2}）：
//
//	κ = N·S^{-3/2}，N = x'y'' − y'x''，S = x'² + y'²
//	κ' = 0  ⟺  f(t) = N'·S − 3·N·(x'x'' + y'y'') = 0
//
// 其中 N' = x'·y‴ − y'·x‴。对三次贝塞尔 f 为六次多项式。
func kappaStationary(d0, d1, d2 Point, t float64) float64 {
	v1, v2 := derivs(d0, d1, d2, t)
	x3, y3 := 2*d2.X, 2*d2.Y // 三阶导数（常数）
	n := v1.X*v2.Y - v1.Y*v2.X
	np := v1.X*y3 - v1.Y*x3
	s := v1.X*v1.X + v1.Y*v1.Y
	return np*s - 3*n*(v1.X*v2.X+v1.Y*v2.Y)
}

// speed2 返回 |B'(t)|²。
func speed2(d0, d1, d2 Point, t float64) float64 {
	v1, _ := derivs(d0, d1, d2, t)
	return v1.X*v1.X + v1.Y*v1.Y
}

// validateCurve 校验控制点与采样段数，错误均带段号。
func validateCurve(p0, p1, p2, p3 Point, n int) error {
	if n < minExtremaSegments {
		return &CurveError{Segment: n, Err: ErrBadSegments}
	}
	if p0 == p1 && p1 == p2 && p2 == p3 {
		return &CurveError{Segment: 0, Err: ErrZeroDeriv}
	}
	pts := [4]Point{p0, p1, p2, p3}
	for i := 0; i < 3; i++ {
		if pts[i] == pts[i+1] {
			return &CurveError{Segment: i, Err: ErrCoincident}
		}
	}
	return nil
}

// CurvatureExtrema 在 (0,1) 内定位曲率的全部局部极值。
//
// 方法：先求 κ'(t)=0 的根（曲率驻点），再求速度为零、曲率导数不存在
// 的尖点；每个候选点都要求左右两侧曲率单调（极大值左升右降，极小值
// 反之），不满足的驻点（如曲率的驻留拐点）被丢弃。
//
// n 为采样段数，n 越大根的定位越稳；不同 n 下结果应收敛到同一位置。
// 端点 t=0、1 不作为极值点返回。
func CurvatureExtrema(p0, p1, p2, p3 Point, n int) ([]Extremum, error) {
	if err := validateCurve(p0, p1, p2, p3, n); err != nil {
		return nil, err
	}
	d0, d1, d2 := derivCoeffs(p0, p1, p2, p3)

	var out []Extremum

	// 1. 尖点：S(t)=|B'|² 在 (0,1) 内触及零的位置。
	out = append(out, findCusps(d0, d1, d2, n)...)

	// 2. 曲率驻点：f(t)=0 的符号变化区间做二分；f 恰好在格点上
	// 为零（对称曲线中点常见）时直接记为候选根。
	dt := 1 / float64(n)
	tryRoot := func(t float64) {
		if e, ok := classifyExtremum(p0, p1, p2, p3, t, dt); ok {
			out = append(out, e)
		}
	}
	prevT, prevF := 0.0, kappaStationary(d0, d1, d2, 0)
	for i := 1; i <= n; i++ {
		t := float64(i) * dt
		f := kappaStationary(d0, d1, d2, t)
		if f == 0 {
			if t < 1 {
				tryRoot(t)
			}
			continue // 不更新 prevF，保留最近一个非零符号
		}
		if prevF != 0 && (prevF < 0) != (f < 0) {
			tryRoot(bisect(func(x float64) float64 { return kappaStationary(d0, d1, d2, x) }, prevT, t))
		}
		prevT, prevF = t, f
	}

	// 3. 去重（尖点附近的驻点根是曲率发散的数值假象）并排序。
	out = dedupeExtrema(out, 2*dt)
	return out, nil
}

// bisect 在 [a,b] 上二分求 f 的根，要求 f(a)、f(b) 异号。
func bisect(f func(float64) float64, a, b float64) float64 {
	fa := f(a)
	for i := 0; i < 200; i++ {
		m := 0.5 * (a + b)
		fm := f(m)
		if fm == 0 || b-a < 1e-15 {
			return m
		}
		if (fa < 0) != (fm < 0) {
			b = m
		} else {
			a, fa = m, fm
		}
	}
	return 0.5 * (a + b)
}

// classifyExtremum 用根两侧各两个采样点核验单调性，判定极大/极小。
// 两侧曲率不单调（驻留拐点）时返回 ok=false。
func classifyExtremum(p0, p1, p2, p3 Point, t, dt float64) (Extremum, bool) {
	d := dt // 采样步长，保证在根所在的符号变化区间量级上
	tl2, tl1 := t-2*d, t-d
	tr1, tr2 := t+d, t+2*d
	if tl2 < 0 || tr2 > 1 {
		// 根太靠近端点，缩小步长再试
		d = math.Min(t, 1-t) / 3
		if d <= 0 {
			return Extremum{}, false
		}
		tl2, tl1 = t-2*d, t-d
		tr1, tr2 = t+d, t+2*d
	}
	kc := Curvature(p0, p1, p2, p3, t)
	kl1, kl2 := Curvature(p0, p1, p2, p3, tl1), Curvature(p0, p1, p2, p3, tl2)
	kr1, kr2 := Curvature(p0, p1, p2, p3, tr1), Curvature(p0, p1, p2, p3, tr2)
	if math.IsNaN(kc) || math.IsNaN(kl1) || math.IsNaN(kr1) {
		return Extremum{}, false // 尖点附近，由尖点分支处理
	}
	switch {
	case kl2 < kl1 && kl1 < kc && kc > kr1 && kr1 > kr2:
		return Extremum{T: t, K: kc, Max: true}, true
	case kl2 > kl1 && kl1 > kc && kc < kr1 && kr1 < kr2:
		return Extremum{T: t, K: kc}, true
	}
	return Extremum{}, false
}

// findCusps 在 (0,1) 内寻找速度为零的尖点：对 S(t)=|B'|² 找低于
// 阈值的局部最小值并用黄金分割细化。
func findCusps(d0, d1, d2 Point, n int) []Extremum {
	// 阈值按控制多边形尺度归一
	scale := math.Max(math.Hypot(d0.X, d0.Y),
		math.Max(math.Hypot(d0.X+d1.X+d2.X, d0.Y+d1.Y+d2.Y), 1))
	eps := 1e-10 * scale * scale

	dt := 1 / float64(n)
	var out []Extremum
	prev := speed2(d0, d1, d2, 0)
	cur := speed2(d0, d1, d2, dt)
	for i := 1; i < n; i++ {
		next := speed2(d0, d1, d2, float64(i+1)*dt)
		if cur < eps && cur <= prev && cur <= next {
			lo, hi := float64(i-1)*dt, float64(i+1)*dt
			t := goldenMin(func(x float64) float64 { return speed2(d0, d1, d2, x) }, lo, hi)
			if speed2(d0, d1, d2, t) < eps {
				out = append(out, Extremum{T: t, K: math.Inf(1), Max: true, Cusp: true})
			}
		}
		prev, cur = cur, next
	}
	return out
}

// goldenMin 黄金分割法求 f 在 [a,b] 上的最小值点。
func goldenMin(f func(float64) float64, a, b float64) float64 {
	const gr = 0.6180339887498949
	c, d := b-gr*(b-a), a+gr*(b-a)
	fc, fd := f(c), f(d)
	for i := 0; i < 200 && b-a > 1e-15; i++ {
		if fc < fd {
			b, d, fd = d, c, fc
			c = b - gr*(b-a)
			fc = f(c)
		} else {
			a, c, fc = c, d, fd
			d = a + gr*(b-a)
			fd = f(d)
		}
	}
	return 0.5 * (a + b)
}

// dedupeExtrema 合并间距小于 tol 的极值点（尖点优先），按 T 排序。
func dedupeExtrema(in []Extremum, tol float64) []Extremum {
	// 简单插入排序，极值点数量很小
	for i := 1; i < len(in); i++ {
		for j := i; j > 0 && in[j-1].T > in[j].T; j-- {
			in[j-1], in[j] = in[j], in[j-1]
		}
	}
	out := in[:0]
	for _, e := range in {
		if len(out) > 0 && e.T-out[len(out)-1].T < tol {
			if e.Cusp && !out[len(out)-1].Cusp {
				out[len(out)-1] = e
			}
			continue
		}
		out = append(out, e)
	}
	return out
}
