package geom

import (
	"math"
	"testing"
)

func absArea(p []Point) float64 {
	a, _ := Area(p)
	return math.Abs(a)
}

// TestJoinAnalyticArea 接头面积对照手算解析公式：直折风筝形 h²tan(α/2)、
// 斜切 h·r·tan(α/2)、圆角扇形 h²α/2（含离散上界），并验证左右转镜像对称。
func TestJoinAnalyticArea(t *testing.T) {
	h := 5.0
	v := Point{X: 10, Y: 0}
	t0 := Point{X: 1, Y: 0}
	for _, deg := range []int{30, 45, 60, 90, 120, 150, 170} {
		rad := float64(deg) * math.Pi / 180
		t1l := Point{X: math.Cos(rad), Y: math.Sin(rad)}
		mp := miterJoin(v, t0, t1l, h)
		rp := roundJoin(v, t0, t1l, h, h, 256)
		bp := bevelJoin(v, t0, t1l, h, h)
		bpHalf := bevelJoin(v, t0, t1l, h, h/2)

		wantMiter := h * h * math.Tan(rad/2)
		wantSector := 0.5 * h * h * rad
		wantBevelFull := h * h * math.Tan(rad/2)
		// 一般半径 r：鞋带公式化简为 h·r·tan(α/2)（与 r 线性）
		r := h / 2
		wantBevelHalf := h * r * math.Tan(rad/2)

		if math.Abs(absArea(mp)-wantMiter) > 1e-8*wantMiter {
			t.Fatalf("deg=%d miter 面积 %v want %v", deg, absArea(mp), wantMiter)
		}
		if math.Abs(absArea(bp)-wantBevelFull) > 1e-8*wantBevelFull {
			t.Fatalf("deg=%d bevel(r=h) 面积 %v want %v", deg, absArea(bp), wantBevelFull)
		}
		if math.Abs(absArea(bpHalf)-wantBevelHalf) > 1e-8*math.Max(1, wantBevelHalf) {
			t.Fatalf("deg=%d bevel(r=h/2) 面积 %v want %v", deg, absArea(bpHalf), wantBevelHalf)
		}
		// 圆弧离散近似误差随段数收敛：256 段下相对误差 < 1e-5
		if math.Abs(absArea(rp)-wantSector) > 5e-5*wantSector {
			t.Fatalf("deg=%d round 面积 %v want %v", deg, absArea(rp), wantSector)
		}
		// 弧上点（除 v 与 h 偏移点）到 v 距离必须为 r
		for _, q := range rp {
			d := math.Hypot(q.X-v.X, q.Y-v.Y)
			if d > 1e-9 && math.Abs(d-h) > 1e-9 {
				t.Fatalf("deg=%d 圆弧点 %v 距 v=%v 非 r", deg, q, d)
			}
		}
		// 右转镜像对称
		t1r := Point{X: math.Cos(rad), Y: -math.Sin(rad)}
		mr := miterJoin(v, t0, t1r, h)
		rr := roundJoin(v, t0, t1r, h, h, 256)
		br := bevelJoin(v, t0, t1r, h, h)
		if math.Abs(absArea(mr)-wantMiter) > 1e-8*wantMiter ||
			math.Abs(absArea(rr)-wantSector) > 5e-5*wantSector ||
			math.Abs(absArea(br)-wantBevelFull) > 1e-8*wantBevelFull {
			t.Fatalf("deg=%d 右转面积不对称", deg)
		}
	}
}
