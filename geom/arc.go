package geom

// 圆弧逼近：三点定圆 → 弦高误差反解分段数 → 折线逼近。
//
// 误差模型：把圆弧等分成 n 段等长弦，每段扫角 θ = ∠/n，弦中点处
// 径向收缩 r·(1−cos(θ/2))，即该段相对解析圆的最大径向偏差。
// 因此 n ≥ θ / (2·acos(1−d/r)) 时，全弧径向偏差不超过 d。
// 该值为解析上界，实际采样点严格落在圆上，径向偏差量级为舍入误差。

import (
	"errors"
	"fmt"
	"math"
)

const (
	// MaxArcSegments 是分段数的硬上限，避免极苛刻公差下的超大分配。
	MaxArcSegments = 1 << 22
	// 判定三点共线与两点重合时，相对尺度下的无量纲阈值。
	collinearEps  = 1e-10
	coincidentEps = 1e-12
)

// PointError 把错误绑定到输入点号（0 起）。
type PointError struct {
	Err error
	// Pts 为涉及的点号；重合两点时长度为 2，其余为 3。
	Pts []int
}

func (e *PointError) Error() string {
	if len(e.Pts) == 2 {
		return fmt.Sprintf("点 %d 与点 %d 重合，无法定圆", e.Pts[0], e.Pts[1])
	}
	if len(e.Pts) == 3 {
		return fmt.Sprintf("%v：点 %d、%d、%d", e.Err, e.Pts[0], e.Pts[1], e.Pts[2])
	}
	return e.Err.Error()
}

func (e *PointError) Unwrap() error { return e.Err }

// Circle 是三点定圆的结果。
type Circle struct {
	Center Point
	Radius float64
}

// Circumcircle 按三点解出外接圆圆心与半径。
// 三点共线返回 *PointError{Err: ErrCollinear}；
// 任意两点重合返回 *PointError{Err: ErrCoincident}（Pts 为重合的两点）。
func Circumcircle(a, b, c Point) (Circle, error) {
	p := [3]Point{a, b, c}
	// 重合优先于共线：重合点使叉积也为零，不先查会被误报为共线。
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 3; j++ {
			if coincident(p[i], p[j]) {
				return Circle{}, &PointError{Err: ErrCoincident, Pts: []int{i, j}}
			}
		}
	}

	ax, ay := a.X, a.Y
	bx, by := b.X, b.Y
	cx, cy := c.X, c.Y
	d := 2 * (ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))

	// 尺度归一的共线判据：叉积 / 最长边的平方。
	ab := math.Hypot(bx-ax, by-ay)
	bc := math.Hypot(cx-bx, cy-by)
	ca := math.Hypot(ax-cx, ay-cy)
	longest := math.Max(ab, math.Max(bc, ca))
	area2 := math.Abs((bx-ax)*(cy-ay) - (by-ay)*(cx-ax))
	if longest == 0 || area2/(longest*longest) < collinearEps {
		return Circle{}, &PointError{Err: ErrCollinear, Pts: []int{0, 1, 2}}
	}

	a2 := ax*ax + ay*ay
	b2 := bx*bx + by*by
	c2 := cx*cx + cy*cy
	ux := (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by)) / d
	uy := (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax)) / d
	center := Point{ux, uy}
	return Circle{Center: center, Radius: math.Hypot(ax-ux, ay-uy)}, nil
}

// Canvas 是绘图画幅的轴对齐包围盒。
type Canvas struct {
	Min, Max Point
}

// Diagonal 返回画幅对角线长度，作为可绘制圆弧半径的自然上限。
func (c Canvas) Diagonal() float64 {
	return math.Hypot(c.Max.X-c.Min.X, c.Max.Y-c.Min.Y)
}

// CircumcircleInCanvas 在三点定圆之上加画幅检查：
// 半径超过画幅对角线时返回 *PointError{Err: ErrRadiusTooLarge}。
func CircumcircleInCanvas(a, b, c Point, canvas Canvas) (Circle, error) {
	circ, err := Circumcircle(a, b, c)
	if err != nil {
		return Circle{}, err
	}
	if circ.Radius > canvas.Diagonal() {
		return Circle{}, &PointError{Err: ErrRadiusTooLarge, Pts: []int{0, 1, 2}}
	}
	return circ, nil
}

// coincident 在同一坐标尺度下判定两点是否重合。
func coincident(a, b Point) bool {
	scale := math.Max(1, math.Max(
		math.Max(math.Abs(a.X), math.Abs(a.Y)),
		math.Max(math.Abs(b.X), math.Abs(b.Y)),
	))
	return math.Hypot(a.X-b.X, a.Y-b.Y) <= coincidentEps*scale
}

// SegmentsForSagitta 按弦高误差反解整圆等分所需段数。
// 返回满足「每段弦中点径向偏差 ≤ sagitta」的最小段数，至少为 1。
// 整圆是最保守的情形，任何圆弧的段数都不超过它。
func SegmentsForSagitta(r, sagitta float64) (int, error) {
	return SegmentsForSweep(r, 2*math.Pi, sagitta)
}

// SegmentsForSweep 按弦高误差反解扫角为 sweep（弧度）的圆弧所需段数。
func SegmentsForSweep(r, sweep, sagitta float64) (int, error) {
	if !(r > 0) || math.IsInf(r, 0) {
		return 0, ErrBadRadius
	}
	if !(sagitta > 0) || math.IsInf(sagitta, 0) {
		return 0, ErrBadSagitta
	}
	s := math.Abs(sweep)
	if !(s > 0) || s > 2*math.Pi+1e-12 {
		return 0, ErrBadSweep
	}
	if s > 2*math.Pi {
		s = 2 * math.Pi
	}
	// 允许偏差达到 2r（直径）时一段足够；再大也无从更粗。
	if sagitta >= 2*r {
		return 1, nil
	}
	// 段扫角 θ 需满足 r(1−cos(θ/2)) ≤ d，即 θ ≤ 2·acos(1−d/r)。
	halfAngle := math.Acos(1 - sagitta/r)
	if !(halfAngle > 0) {
		// 公差小到浮点无法表示，所需段数必然远超上限。
		return 0, ErrTooManySegments
	}
	n := int(math.Ceil(s / (2 * halfAngle)))
	if n < 1 {
		n = 1
	}
	if n > MaxArcSegments {
		return 0, ErrTooManySegments
	}
	return n, nil
}

// Arc 是一段圆弧的解析描述。
type Arc struct {
	Center Point
	Radius float64
	// Start 为起始极角（弧度），Sweep 为带符号扫角，正为逆时针。
	Start float64
	Sweep float64
}

// ArcThroughPoints 按三点构造圆弧：圆心为外接圆圆心，弧从 a 经 b 到 c。
// 角度按相对圆心的位置取，b 必须落在 a→c 的指定一侧，否则该弧不经过 b。
func ArcThroughPoints(a, b, c Point) (Arc, error) {
	circ, err := Circumcircle(a, b, c)
	if err != nil {
		return Arc{}, err
	}
	start := math.Atan2(a.Y-circ.Center.Y, a.X-circ.Center.X)
	mid := math.Atan2(b.Y-circ.Center.Y, b.X-circ.Center.X)
	end := math.Atan2(c.Y-circ.Center.Y, c.X-circ.Center.X)

	// 逆时针从 start 到 mid、再到 end 的扫角，均归一到 (0, 2π]。
	ccwSweep := func(from, to float64) float64 {
		s := math.Mod(to-from, 2*math.Pi)
		if s <= 0 {
			s += 2 * math.Pi
		}
		return s
	}
	toMid := ccwSweep(start, mid)
	toEnd := ccwSweep(start, end)

	var sweep float64
	if toMid <= toEnd {
		sweep = toEnd // 逆时针，b 在中间
	} else {
		sweep = toEnd - 2*math.Pi // 顺时针
	}
	return Arc{Center: circ.Center, Radius: circ.Radius, Start: start, Sweep: sweep}, nil
}

// Validate 检查扫角是否在 (0, 2π] 内（负值表示顺时针）。
func (a Arc) Validate() error {
	if !(a.Radius > 0) || math.IsInf(a.Radius, 0) {
		return ErrBadRadius
	}
	s := math.Abs(a.Sweep)
	if !(s > 0) || s > 2*math.Pi+1e-12 {
		return ErrBadSweep
	}
	return nil
}

// PointAt 取扫角参数 t∈[0,1] 处的弧上点。
func (a Arc) PointAt(t float64) Point {
	ang := a.Start + a.Sweep*t
	return Point{a.Center.X + a.Radius*math.Cos(ang), a.Center.Y + a.Radius*math.Sin(ang)}
}

// RadialDeviation 返回弧上采样点 q 相对解析圆的径向偏差 |‖q−Center‖ − r|。
// 折线顶点应给出舍入量级的值，弦中点应不超过弦高公差。
func (a Arc) RadialDeviation(q Point) float64 {
	return math.Abs(math.Hypot(q.X-a.Center.X, q.Y-a.Center.Y) - a.Radius)
}

// FlattenArc 按弦高误差 sagitta 把圆弧折线化，返回 n+1 个点及所用段数。
// 输出为逆时针（Sweep>0）或顺时针（Sweep<0）顺序，起止点分别在 t=0、t=1。
func FlattenArc(a Arc, sagitta float64) ([]Point, int, error) {
	if err := a.Validate(); err != nil {
		return nil, 0, err
	}
	n, err := SegmentsForSweep(a.Radius, a.Sweep, sagitta)
	if err != nil {
		return nil, 0, err
	}
	out := make([]Point, n+1)
	for i := 0; i <= n; i++ {
		out[i] = a.PointAt(float64(i) / float64(n))
	}
	// 端点用精确输入角，避免端点处累积舍入。
	out[0] = Point{a.Center.X + a.Radius*math.Cos(a.Start), a.Center.Y + a.Radius*math.Sin(a.Start)}
	endAng := a.Start + a.Sweep
	out[n] = Point{a.Center.X + a.Radius*math.Cos(endAng), a.Center.Y + a.Radius*math.Sin(endAng)}
	return out, n, nil
}

// ArcDeviationReport 汇总一段折线逼近的偏差统计与圆心反算。
type ArcDeviationReport struct {
	Segments int
	Sagitta  float64 // 要求的弦高误差上界
	// MaxDeviation 为折线顶点相对解析圆的最大径向偏差（应近似舍入量级）。
	MaxDeviation float64
	// MaxSagitta 为每段弦中点相对解析圆的最大径向偏差，应不超过 Sagitta。
	MaxSagitta float64
	// Refit 为由折线首、中、末三点反算出的圆。
	Refit Circle
	// RefitResidual 为反算圆心相对解析圆心的距离。
	RefitResidual float64
	// RadiusResidual 为反算半径相对解析半径的绝对差。
	RadiusResidual float64
	// RefitMaxPointDeviation 为折线各顶点到反算圆的最大径向偏差。
	RefitMaxPointDeviation float64
}

// Report 对一段已折线化的弧做偏差统计，并把首、中、末三点反算回圆验证。
// pts 须由 FlattenArc 依序产生（至少 3 点）。
func (a Arc) Report(pts []Point, sagitta float64) (ArcDeviationReport, error) {
	if err := a.Validate(); err != nil {
		return ArcDeviationReport{}, err
	}
	if len(pts) < 3 {
		return ArcDeviationReport{}, ErrTooFewPts
	}
	n := len(pts) - 1
	rep := ArcDeviationReport{Segments: n, Sagitta: sagitta}

	for _, q := range pts {
		if dv := a.RadialDeviation(q); dv > rep.MaxDeviation {
			rep.MaxDeviation = dv
		}
	}
	// 每段弦中点：解析圆上的真正弦高中点在 Center 与折线中点的连线上。
	for i := 0; i+1 < len(pts); i++ {
		mid := Point{(pts[i].X + pts[i+1].X) / 2, (pts[i].Y + pts[i+1].Y) / 2}
		if dv := a.RadialDeviation(mid); dv > rep.MaxSagitta {
			rep.MaxSagitta = dv
		}
	}

	refit, err := Circumcircle(pts[0], pts[n/2], pts[n])
	if err != nil {
		return rep, fmt.Errorf("反算圆心失败：%w", err)
	}
	rep.Refit = refit
	rep.RefitResidual = math.Hypot(refit.Center.X-a.Center.X, refit.Center.Y-a.Center.Y)
	rep.RadiusResidual = math.Abs(refit.Radius - a.Radius)
	for _, q := range pts {
		if dv := math.Abs(math.Hypot(q.X-refit.Center.X, q.Y-refit.Center.Y) - refit.Radius); dv > rep.RefitMaxPointDeviation {
			rep.RefitMaxPointDeviation = dv
		}
	}
	return rep, nil
}

// CircleFromPoints 是对 Circumcircle 的显式别名，便于调用侧表达「三点定圆」。
func CircleFromPoints(a, b, c Point) (Circle, error) { return Circumcircle(a, b, c) }

// IsPointError 判断 err 是否为带点号的输入错误，并抽出点号。
func IsPointError(err error) ([]int, bool) {
	var pe *PointError
	if errors.As(err, &pe) {
		return pe.Pts, true
	}
	return nil, false
}
