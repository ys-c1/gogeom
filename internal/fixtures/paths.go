// Package paths 提供测试与示例共用的曲线用例: 直线、圆弧、S 形、尖角、近退化。
// 每个用例接受"曲率档" 1 或 2, 档 2 在同一几何上用更急的转折 / 更大的幅度。
package fixtures

import "math"

// Named 是一条命名路径及其按曲率档生成控制点的函数。
type Named struct {
	Name string
	Fn   func(tier float64) []CubicPts
}

// CubicPts 是一段三次 Bezier 的 4 个控制点 [P0,P1,P2,P3]。
type CubicPts [4][2]float64

// linePath 几何始终是 x∈[0,10] 的直线; 档 2 把控制点挤到一端,
// 使 |dC/dt| 沿直线剧烈变化 (参数化极不均匀), 检验重参数化是否只看弧长。
func linePath(tier float64) []CubicPts {
	if tier > 1 {
		return []CubicPts{{{0, 0}, {0.2, 0}, {0.4, 0}, {10, 0}}}
	}
	return []CubicPts{{{0, 0}, {2, 0}, {7, 0}, {10, 0}}}
}

// arcCubic 用三次 Bezier 逼近半径 r、圆心角 theta 的圆弧。
func arcCubic(r, theta float64) CubicPts {
	k := 4.0 / 3.0 * math.Tan(theta/4)
	ct, st := math.Cos(theta), math.Sin(theta)
	return CubicPts{
		{r, 0},
		{r, k * r},
		{r*ct + k*r*st, r*st - k*r*ct},
		{r * ct, r * st},
	}
}

func arcPath(tier float64) []CubicPts {
	sw := math.Pi / 2
	if tier > 1 {
		sw = math.Pi
	}
	return []CubicPts{arcCubic(10, sw)}
}

func sPath(tier float64) []CubicPts {
	a := 6.0
	if tier > 1 {
		a = 14.0
	}
	return []CubicPts{
		{{0, 0}, {5, a}, {10, a}, {15, 0}},
		{{15, 0}, {20, -a}, {25, -a}, {30, 0}},
	}
}

func cornerPath(tier float64) []CubicPts {
	turn := math.Pi / 4
	if tier > 1 {
		turn = math.Pi * 5 / 6
	}
	cx, cy := 10.0, 0.0
	dx, dy := math.Cos(turn), -math.Sin(turn)
	// 档 1: 45° 转折; 档 2: 150° 转折, 且拐点侧手柄收短到 1.2,
	// 曲率在拐角附近高度集中 (端点速率仍远离零, 不产生尖点/自交)。
	h1, h2 := 3.0, 7.0
	if tier > 1 {
		h1, h2 = 1.2, 8.8
	}
	return []CubicPts{
		{{0, 0}, {3, 0.6}, {7, 0.6}, {cx, cy}},
		{{cx, cy}, {cx + dx*h1, cy + dy*h1}, {cx + dx*h2, cy + dy*h2}, {cx + dx*10, cy + dy*10}},
	}
}

// nearPath 的近退化包含两层含义: 相邻控制点 P0/P1 几乎重合 (端点速率趋零),
// 且两段长度相差悬殊; 档 2 把 eps 压到 1e-7。
func nearPath(tier float64) []CubicPts {
	eps := 1e-4
	if tier > 1 {
		eps = 1e-7
	}
	return []CubicPts{
		{{0, 0}, {eps, 0}, {10, 10}, {10, 0}},
		{{10, 0}, {eps, 0}, {-10, 10}, {-10, 0}},
	}
}

// All 是五个用例的固定顺序。
var All = []Named{
	{"直线", linePath},
	{"圆弧", arcPath},
	{"S形", sPath},
	{"尖角", cornerPath},
	{"近退化", nearPath},
}
