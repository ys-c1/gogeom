package geom

import (
	"fmt"
	"math"
	"strings"
	"time"
)

// Cubic 表示一段三次贝塞尔曲线。
type Cubic struct{ P0, P1, P2, P3 Point }

// JoinDepth 为拼接深度档位。
type JoinDepth int

const (
	// JoinTangent 一档：仅按端点导数反解缩放系数做切向对齐（G1）。
	JoinTangent JoinDepth = iota + 1
	// JoinCurvature 二档：在一档基础上追加一次曲率对齐修正（G2）。
	JoinCurvature
	// JoinCurvatureIter 三档：迭代切向+曲率修正，直到接缝切向夹角为零
	// 且曲率差收敛到上界以内。
	JoinCurvatureIter
)

// CurvatureTol 是曲率差上界的相对基准：上界 = CurvatureTol * max(1, |κ左|, |κ右|)。
const CurvatureTol = 1e-9

// moveTol 是控制点修正计入修正次数的最小位移。
const moveTol = 1e-12

// maxJoinIter 是三档迭代的最大扫描轮数。
const maxJoinIter = 100

// SeamInfo 记录一处接缝的修正结果。
type SeamInfo struct {
	Seam           int     // 接缝号：第 Seam 段与第 Seam+1 段之间
	Scale          float64 // 由端点导数反解出的第二段切向缩放系数 λ
	TangentAngle   float64 // 修正后接缝两侧切向夹角（弧度），须为零
	CurvatureDiff  float64 // 修正后接缝两侧曲率差绝对值
	CurvatureBound float64 // 曲率差上界
	Corrections    int     // 该接缝累计修正次数
}

// LevelStat 记录一档拼接深度的修正次数与耗时。
type LevelStat struct {
	Depth       JoinDepth
	Corrections int
	Elapsed     time.Duration
}

// JoinReport 为拼接报告。
type JoinReport struct {
	Seams  []SeamInfo
	Levels []LevelStat
}

// Table 输出接缝两项指标（切向夹角、曲率差）与各档修正次数表。
func (r *JoinReport) Table() string {
	var b strings.Builder
	b.WriteString("接缝 | 缩放系数λ | 切向夹角(rad) | 曲率差 | 曲率差上界 | 修正次数\n")
	for _, s := range r.Seams {
		fmt.Fprintf(&b, "%4d | %9.4f | %13.3e | %9.3e | %11.3e | %6d\n",
			s.Seam, s.Scale, s.TangentAngle, s.CurvatureDiff, s.CurvatureBound, s.Corrections)
	}
	for _, l := range r.Levels {
		fmt.Fprintf(&b, "深度%d: 修正 %d 次, 耗时 %s\n", l.Depth, l.Corrections, l.Elapsed)
	}
	return b.String()
}

// Join 将多段三次贝塞尔首尾相接并消除接缝折点：
// 先按相邻两段端点导数反解第二段控制点的切向缩放系数使切向对齐，
// 再按曲率对齐修正第二段 P2。depth 取三档之一，超出范围时收敛到最近档。
// 输入退化（段数不足、控制点重合、端点导数为零）时报错并给出段号。
// 不改动入参，返回修正后的新切片。
func Join(segs []Cubic, depth JoinDepth) ([]Cubic, *JoinReport, error) {
	if err := validateJoin(segs); err != nil {
		return nil, nil, err
	}
	if depth < JoinTangent {
		depth = JoinTangent
	}
	if depth > JoinCurvatureIter {
		depth = JoinCurvatureIter
	}
	out := append([]Cubic(nil), segs...)
	nSeam := len(out) - 1
	rep := &JoinReport{Seams: make([]SeamInfo, nSeam)}
	corr := make([]int, nSeam)
	scale := make([]float64, nSeam)

	// 一档：切向对齐扫描。
	t0 := time.Now()
	n1 := 0
	for i := 0; i < nSeam; i++ {
		lambda, c := alignTangent(&out[i], &out[i+1])
		scale[i] = lambda
		corr[i] += c
		n1 += c
	}
	rep.Levels = append(rep.Levels, LevelStat{JoinTangent, n1, time.Since(t0)})

	// 二档：逐接缝（切向复核 + 一次曲率修正）扫描。
	// 曲率修正会改变本段末端切向，故按接缝顺序交织处理，
	// 使下游接缝读到的是修正后的导数。
	if depth >= JoinCurvature {
		t1 := time.Now()
		n2 := 0
		for i := 0; i < nSeam; i++ {
			_, c := alignTangent(&out[i], &out[i+1])
			c += alignCurvature(&out[i], &out[i+1])
			corr[i] += c
			n2 += c
		}
		rep.Levels = append(rep.Levels, LevelStat{JoinCurvature, n2, time.Since(t1)})
	}

	// 三档：重复交织扫描，直到一整轮没有任何有效修正（收敛复核）。
	if depth >= JoinCurvatureIter {
		t2 := time.Now()
		n3 := 0
		for iter := 0; iter < maxJoinIter; iter++ {
			moved := 0
			for i := 0; i < nSeam; i++ {
				_, c := alignTangent(&out[i], &out[i+1])
				c += alignCurvature(&out[i], &out[i+1])
				corr[i] += c
				moved += c
			}
			n3 += moved
			if moved == 0 {
				break
			}
		}
		rep.Levels = append(rep.Levels, LevelStat{JoinCurvatureIter, n3, time.Since(t2)})
	}

	// 汇总测量：切向夹角、曲率差与上界。
	for i := 0; i < nSeam; i++ {
		ka := curvature(endTangent(out[i]), second1(out[i]))
		kb := curvature(startTangent(out[i+1]), second0(out[i+1]))
		bound := CurvatureTol * math.Max(1, math.Max(math.Abs(ka), math.Abs(kb)))
		rep.Seams[i] = SeamInfo{
			Seam:           i,
			Scale:          scale[i],
			TangentAngle:   seamAngle(out[i], out[i+1]),
			CurvatureDiff:  math.Abs(ka - kb),
			CurvatureBound: bound,
			Corrections:    corr[i],
		}
	}
	return out, rep, nil
}

// validateJoin 检查退化输入，错误信息携带段号。
func validateJoin(segs []Cubic) error {
	if len(segs) < 2 {
		return fmt.Errorf("仅 %d 段: %w", len(segs), ErrTooFewSegments)
	}
	for i, c := range segs {
		if c.P0 == c.P1 && c.P1 == c.P2 && c.P2 == c.P3 {
			return fmt.Errorf("第 %d 段: %w", i, ErrCoincidentControls)
		}
		if c.P0 == c.P1 || c.P2 == c.P3 {
			return fmt.Errorf("第 %d 段: %w", i, ErrZeroDerivative)
		}
	}
	return nil
}

// alignTangent 反解第二段切向缩放系数 λ = |dB|/|dA|，
// 将 b.P1 放到 a 末端切线方向上（保持原切向模长），并把 b.P0 吸附到 a.P3。
// 返回 λ 与修正次数（0 或 1）。
func alignTangent(a, b *Cubic) (float64, int) {
	moved := pdist(b.P0, a.P3)
	b.P0 = a.P3 // C0 吸附
	dA := endTangent(*a)
	dB := startTangent(*b)
	lambda := pnorm(dB) / pnorm(dA)
	// P1 - P0 是控制边向量（导数的 1/3），方向取 dA、模长保持 |dB|/3。
	np := padd(b.P0, pscale(dA, lambda/3))
	if d := pdist(np, b.P1); d > moved {
		moved = d
	}
	b.P1 = np
	if moved > moveTol {
		return lambda, 1
	}
	return lambda, 0
}

// alignCurvature 按曲率对齐修正 b.P2：令 b 起点曲率等于 a 终点曲率。
// 二阶导数的法向分量取 λ²·ddA 的法向分量，切向分量保持原值（改动最小）。
// 返回修正次数（0 或 1）。
func alignCurvature(a, b *Cubic) int {
	dA := endTangent(*a)
	ddA := second1(*a)
	dB := startTangent(*b)
	na, nb := pnorm(dA), pnorm(dB)
	lambda := nb / na
	tx, ty := dA.X/na, dA.Y/na // 单位切向
	nx, ny := -ty, tx          // 单位法向
	ddB := second0(*b)
	targetN := lambda * lambda * (ddA.X*nx + ddA.Y*ny)
	curT := ddB.X*tx + ddB.Y*ty
	ddNew := Point{tx*curT + nx*targetN, ty*curT + ny*targetN}
	// 由 B''(0) = 6(P0 - 2P1 + P2) 反解 P2 = 2P1 - P0 + B''/6。
	np := Point{
		2*b.P1.X - b.P0.X + ddNew.X/6,
		2*b.P1.Y - b.P0.Y + ddNew.Y/6,
	}
	moved := pdist(np, b.P2)
	b.P2 = np
	if moved > moveTol {
		return 1
	}
	return 0
}

// --- 导数与曲率 ---

func startTangent(c Cubic) Point { return pscale(psub(c.P1, c.P0), 3) }
func endTangent(c Cubic) Point   { return pscale(psub(c.P3, c.P2), 3) }

// second0 为 B”(0) = 6(P0 - 2P1 + P2)。
func second0(c Cubic) Point {
	return Point{6 * (c.P0.X - 2*c.P1.X + c.P2.X), 6 * (c.P0.Y - 2*c.P1.Y + c.P2.Y)}
}

// second1 为 B”(1) = 6(P1 - 2P2 + P3)。
func second1(c Cubic) Point {
	return Point{6 * (c.P1.X - 2*c.P2.X + c.P3.X), 6 * (c.P1.Y - 2*c.P2.Y + c.P3.Y)}
}

// curvature 为带符号曲率 κ = (d × dd) / |d|³。
func curvature(d, dd Point) float64 {
	n := pnorm(d)
	return pcross(d, dd) / (n * n * n)
}

// seamAngle 为接缝两侧切向夹角（弧度，非负）。
func seamAngle(a, b Cubic) float64 {
	dA, dB := endTangent(a), startTangent(b)
	return math.Atan2(math.Abs(pcross(dA, dB)), pdot(dA, dB))
}

// --- 二维向量小工具 ---

func padd(a, b Point) Point           { return Point{a.X + b.X, a.Y + b.Y} }
func psub(a, b Point) Point           { return Point{a.X - b.X, a.Y - b.Y} }
func pscale(a Point, s float64) Point { return Point{a.X * s, a.Y * s} }
func pdot(a, b Point) float64         { return a.X*b.X + a.Y*b.Y }
func pcross(a, b Point) float64 {
	return a.X*b.Y - a.Y*b.X
}
func pnorm(a Point) float64    { return math.Hypot(a.X, a.Y) }
func pdist(a, b Point) float64 { return pnorm(psub(a, b)) }
