package geom

import (
	"fmt"
	"math"
	"sort"
)

// HullLayer 是一层洋葱剥离的结果。
// Index 从 0 开始连续编号；Hull 为该层凸包拐角（逆时针，不含边上共线点）；
// Area 为凸包面积；Size 为该层剥去的点数（含落在边上的共线点）。
// Band 是该层相对下一层新增的环带面积：Area[k]-Area[k+1]，最内层为 Area[last]；
// 全部环带之和恰等于最外层面积 Area[0]。
type HullLayer struct {
	Index int
	Hull  []Point
	Area  float64
	Band  float64
	Size  int
}

// LayerSet 是完整分层剥离的结果。
// StoppedAt 是剥离停止的层号：此时剩余点不足三个，无法再成层；
// Remaining 是停止时的剩余点数。
type LayerSet struct {
	Layers    []HullLayer
	StoppedAt int
	Remaining int
}

// CoveragePlan 是面积覆盖目标反解结果。
// Keep 为需要保留的层数；LastLayer 为保留的最深层号（Keep-1）；
// Peeled 为相对于全部层需要剥掉的层数（全部层数 - Keep）；
// Coverage 为保留这些层时的实际面积覆盖比例。
type CoveragePlan struct {
	Target    float64
	Keep      int
	LastLayer int
	Peeled    int
	Coverage  float64
}

// ConvexHull 返回点集的单层凸包，逆时针排列，只保留拐角顶点。
// 少于三个互异顶点返回 ErrTooFewPts / ErrDegenerate；
// 全部共线时返回 ErrAllCollinear。
func ConvexHull(pts []Point) ([]Point, error) {
	if len(pts) < 3 {
		return nil, &LayerError{Layer: 0, Err: ErrTooFewPts}
	}
	uniq := sortedUnique(pts)
	if len(uniq) < 3 {
		return nil, &LayerError{Layer: 0, Err: ErrDegenerate}
	}
	hull := hullCorners(uniq)
	if len(hull) < 3 {
		return nil, &LayerError{Layer: 0, Err: ErrAllCollinear}
	}
	return hull, nil
}

// Peel 自外向内逐层求凸包并剥离，直到剩余点不足三个。
// 每层都基于剩余点重新求包，因此能拿到内层结构。
func Peel(pts []Point) (*LayerSet, error) {
	if len(pts) < 3 {
		return nil, &LayerError{Layer: 0, Err: ErrTooFewPts}
	}
	rem := sortedUnique(pts)
	if len(rem) < 3 {
		return nil, &LayerError{Layer: 0, Err: ErrDegenerate}
	}
	layers := make([]HullLayer, 0, 8)
	for len(rem) >= 3 {
		hull := hullCorners(rem)
		if len(hull) < 3 {
			return nil, &LayerError{Layer: len(layers), Err: ErrAllCollinear}
		}
		area := math.Abs(polygonArea(hull))
		if area == 0 {
			return nil, &LayerError{Layer: len(layers), Err: ErrAllCollinear}
		}
		next := make([]Point, 0, len(rem))
		for _, p := range rem {
			if !onHullBoundary(hull, p) {
				next = append(next, p)
			}
		}
		if len(next) == len(rem) {
			// 理论上不会发生；防御性返回，避免死循环。
			return nil, &LayerError{Layer: len(layers), Err: ErrDegenerate}
		}
		layers = append(layers, HullLayer{
			Index: len(layers),
			Hull:  hull,
			Area:  area,
			Size:  len(rem) - len(next),
		})
		rem = next
	}
	// 回填环带面积：第 k 层环带 = Area[k]-Area[k+1]，最内层取其自身面积。
	for i := range layers {
		if i+1 < len(layers) {
			layers[i].Band = layers[i].Area - layers[i+1].Area
		} else {
			layers[i].Band = layers[i].Area
		}
	}
	return &LayerSet{Layers: layers, StoppedAt: len(layers), Remaining: len(rem)}, nil
}

// HullLayers 返回分层剥离的层列表，不保留停止信息。
func HullLayers(pts []Point) ([]HullLayer, error) {
	set, err := Peel(pts)
	if err != nil {
		return nil, err
	}
	return set.Layers, nil
}

// ValidateLayers 校验层号连续、每层有效且面积随层号单调不增。
// 层号跳号时返回携带缺失层号的 ErrLayerGap。
func ValidateLayers(layers []HullLayer) error {
	for i := range layers {
		if layers[i].Index != i {
			return &LayerError{
				Layer: i,
				Err:   fmt.Errorf("%w: 期望 %d 实得 %d", ErrLayerGap, i, layers[i].Index),
			}
		}
	}
	for i := range layers {
		l := layers[i]
		if len(l.Hull) < 3 || !(l.Area > 0) {
			return &LayerError{Layer: i, Err: ErrDegenerate}
		}
		if i > 0 && l.Area > layers[i-1].Area+1e-12*layers[i-1].Area {
			return &LayerError{Layer: i, Err: ErrAreaNotMonotone}
		}
		if l.Band < -1e-9*l.Area {
			return &LayerError{Layer: i, Err: ErrAreaNotMonotone}
		}
	}
	return nil
}

// LayerAt 返回指定层号的层。
// 列表本身层号跳号时返回指向缺号处的 ErrLayerGap；
// 请求的层号超出已生成的层数（剩余点不足）时返回指向首个不可用层的 ErrTooFewRemaining。
func LayerAt(layers []HullLayer, idx int) (HullLayer, error) {
	if err := ValidateLayers(layers); err != nil {
		return HullLayer{}, err
	}
	if idx < 0 {
		return HullLayer{}, &LayerError{Layer: 0, Err: ErrTooFewRemaining}
	}
	if idx >= len(layers) {
		return HullLayer{}, &LayerError{Layer: len(layers), Err: ErrTooFewRemaining}
	}
	return layers[idx], nil
}

// TrimLayers 把分层结果裁剪到前 n 层。
// n 超过实际层数时返回携带首个不可用层号的 ErrTooFewRemaining。
func TrimLayers(layers []HullLayer, n int) ([]HullLayer, error) {
	if err := ValidateLayers(layers); err != nil {
		return nil, err
	}
	if n < 0 {
		return nil, &LayerError{Layer: 0, Err: ErrTooFewRemaining}
	}
	if n > len(layers) {
		return nil, &LayerError{Layer: len(layers), Err: ErrTooFewRemaining}
	}
	return append([]HullLayer(nil), layers[:n]...), nil
}

// CoverageAt 返回保留前 keep 层后的面积覆盖比例。
// 各层凸包互相嵌套，第 k 层实际新增覆盖的是它与下一层之间的环带 Band[k]。
// 全部环带之和恰为最外层面积，故覆盖比例 = 已保留环带之和 / 全部环带之和，
// 等价于 (Area[0]-Area[keep])/Area[0]（keep 为层数时含最内层环带，比例为 1；keep 为 0 时为 0）。
// 该口径与 SolveCoverage 一致，且随 keep 单调不减。
func CoverageAt(layers []HullLayer, keep int) (float64, error) {
	if err := ValidateLayers(layers); err != nil {
		return 0, err
	}
	if len(layers) == 0 {
		return 0, &LayerError{Layer: 0, Err: ErrTooFewRemaining}
	}
	if keep < 0 {
		return 0, &LayerError{Layer: 0, Err: ErrTooFewRemaining}
	}
	if keep > len(layers) {
		// 首个不可用的层号即已有层数。
		return 0, &LayerError{Layer: len(layers), Err: ErrTooFewRemaining}
	}
	total := 0.0
	covered := 0.0
	for i, l := range layers {
		total += l.Band
		if i < keep {
			covered += l.Band
		}
	}
	if !(total > 0) {
		return 0, &LayerError{Layer: 0, Err: ErrDegenerate}
	}
	return covered / total, nil
}

// SolveCoverage 由面积覆盖目标反解需要保留到第几层。
// 覆盖口径见 CoverageAt：随层号单调不减，目标须落在 (0,1]；
// 返回首个达到目标的方案，即所需的最小层数。
func SolveCoverage(layers []HullLayer, target float64) (CoveragePlan, error) {
	if err := ValidateLayers(layers); err != nil {
		return CoveragePlan{}, err
	}
	if !(target > 0 && target <= 1) || math.IsNaN(target) {
		return CoveragePlan{}, ErrBadCoverage
	}
	if len(layers) == 0 {
		return CoveragePlan{}, &LayerError{Layer: 0, Err: ErrTooFewRemaining}
	}
	total := 0.0
	for _, l := range layers {
		total += l.Band
	}
	if !(total > 0) {
		return CoveragePlan{}, &LayerError{Layer: 0, Err: ErrDegenerate}
	}
	last := len(layers) - 1
	cum := 0.0
	for keep := 1; keep <= len(layers); keep++ {
		cum += layers[keep-1].Band
		cov := cum / total
		if cov >= target-1e-12 {
			return CoveragePlan{
				Target: target, Keep: keep, LastLayer: keep - 1,
				Peeled: len(layers) - keep, Coverage: cov,
			}, nil
		}
	}
	return CoveragePlan{
		Target: target, Keep: len(layers), LastLayer: last,
		Peeled: 0, Coverage: 1,
	}, nil
}

// PeelToCoverage 先分层剥离，再按面积覆盖目标反解层号并裁剪输出。
func PeelToCoverage(pts []Point, target float64) ([]HullLayer, CoveragePlan, error) {
	set, err := Peel(pts)
	if err != nil {
		return nil, CoveragePlan{}, err
	}
	plan, err := SolveCoverage(set.Layers, target)
	if err != nil {
		return nil, plan, err
	}
	trimmed, err := TrimLayers(set.Layers, plan.Keep)
	if err != nil {
		return nil, plan, err
	}
	return trimmed, plan, nil
}

func cross(o, a, b Point) float64 {
	return (a.X-o.X)*(b.Y-o.Y) - (a.Y-o.Y)*(b.X-o.X)
}

func polygonArea(poly []Point) float64 {
	var s float64
	for i := range poly {
		a, b := poly[i], poly[(i+1)%len(poly)]
		s += a.X*b.Y - b.X*a.Y
	}
	return s / 2
}

func sortedUnique(pts []Point) []Point {
	out := append([]Point(nil), pts...)
	sort.Slice(out, func(i, j int) bool {
		if out[i].X != out[j].X {
			return out[i].X < out[j].X
		}
		return out[i].Y < out[j].Y
	})
	w := 0
	for _, p := range out {
		if w == 0 || out[w-1] != p {
			out[w] = p
			w++
		}
	}
	return out[:w]
}

// hullCorners 对按 (X,Y) 去重排序后的点求逆时针凸包，只保留严格拐角。
// 共线输入会退化为两个端点（长度 < 3）。
func hullCorners(pts []Point) []Point {
	n := len(pts)
	if n < 3 {
		return append([]Point(nil), pts...)
	}
	lower := make([]Point, 0, n)
	for _, p := range pts {
		for len(lower) >= 2 && cross(lower[len(lower)-2], lower[len(lower)-1], p) <= 0 {
			lower = lower[:len(lower)-1]
		}
		lower = append(lower, p)
	}
	upper := make([]Point, 0, n)
	for i := n - 1; i >= 0; i-- {
		p := pts[i]
		for len(upper) >= 2 && cross(upper[len(upper)-2], upper[len(upper)-1], p) <= 0 {
			upper = upper[:len(upper)-1]
		}
		upper = append(upper, p)
	}
	hull := append([]Point(nil), lower[:len(lower)-1]...)
	hull = append(hull, upper[:len(upper)-1]...)
	return hull
}

// onHullBoundary 判断点是否落在凸包边界上（含拐角与边上共线点）。
// 用点到线段的距离判定，避免把锐角附近、落在某条边延长线附近的内部点误判为边界。
func onHullBoundary(hull []Point, p Point) bool {
	for i := range hull {
		a, b := hull[i], hull[(i+1)%len(hull)]
		scale := 1 + math.Abs(p.X) + math.Abs(p.Y) + math.Abs(a.X) + math.Abs(a.Y) + math.Abs(b.X) + math.Abs(b.Y)
		if pointSegDist(p, a, b) <= 1e-12*scale {
			return true
		}
	}
	return false
}

// pointSegDist 返回点 p 到线段 ab 的欧氏距离。
func pointSegDist(p, a, b Point) float64 {
	dx, dy := b.X-a.X, b.Y-a.Y
	l2 := dx*dx + dy*dy
	if l2 == 0 {
		return math.Hypot(p.X-a.X, p.Y-a.Y)
	}
	t := ((p.X-a.X)*dx + (p.Y-a.Y)*dy) / l2
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	return math.Hypot(p.X-(a.X+t*dx), p.Y-(a.Y+t*dy))
}
