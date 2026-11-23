package geom

import (
	"errors"
	"math"
	"testing"
)

// 两档半径：小半径 10、大半径 1000。
const (
	smallR = 10.0
	largeR = 1000.0
)

// arcCase 是一组三点定圆的测试输入与期望。
type arcCase struct {
	name    string
	a, b, c Point
	wantR   float64
	kind    string  // 大弧 / 小弧 / 半圆 / 近直线 / 退化
	wantArc float64 // 期望扫角（绝对值），退化例为 0
}

func arcOnCircle(r, deg float64) Point {
	t := deg * math.Pi / 180
	return Point{r * math.Cos(t), r * math.Sin(t)}
}

// cases 返回覆盖大弧、小弧、半圆、近直线四类的用例；退化例另见
// TestDegenerateErrorsCarryPointNumbers。
func cases(r float64) []arcCase {
	rise := r / 1000
	return []arcCase{
		{"大弧 300°", arcOnCircle(r, 20), arcOnCircle(r, 260), arcOnCircle(r, 320), r, "大弧", 300 * math.Pi / 180},
		{"小弧 30°", arcOnCircle(r, 5), arcOnCircle(r, 20), arcOnCircle(r, 35), r, "小弧", 30 * math.Pi / 180},
		{"半圆 180°", arcOnCircle(r, 180), arcOnCircle(r, 270), arcOnCircle(r, 360), r, "半圆", math.Pi},
		// 近直线：三点几乎共线，但仍有可解的外接圆，半径远大于 r。
		{"近直线", Point{0, 0}, Point{r, rise}, Point{2 * r, 0}, 0, "近直线", 0}}
}

func TestCircumcircleExactOnCircle(t *testing.T) {
	for _, r := range []float64{smallR, largeR} {
		for _, tc := range cases(r) {
			t.Run(tc.name, func(t *testing.T) {
				got, err := Circumcircle(tc.a, tc.b, tc.c)
				if err != nil {
					t.Fatalf("应可解：%v", err)
				}
				if tc.kind == "近直线" {
					// 近直线半径由拱高决定，由 TestNearLineRadiusLarge 校验。
					for i, p := range []Point{tc.a, tc.b, tc.c} {
						if d := math.Hypot(p.X-got.Center.X, p.Y-got.Center.Y); math.Abs(d-got.Radius) > 1e-9*got.Radius {
							t.Fatalf("点 %d 到圆心 %v，半径 %v", i, d, got.Radius)
						}
					}
					return
				}
				if math.Abs(got.Radius-tc.wantR) > 1e-9*tc.wantR {
					t.Fatalf("半径 %v，期望 %v", got.Radius, tc.wantR)
				}
				if math.Hypot(got.Center.X, got.Center.Y) > 1e-9*tc.wantR {
					t.Fatalf("圆心 %v 偏离原点", got.Center)
				}
				// 三点到圆心距离必须一致：这是外接圆的定义性验证。
				for i, p := range []Point{tc.a, tc.b, tc.c} {
					if d := math.Hypot(p.X-got.Center.X, p.Y-got.Center.Y); math.Abs(d-got.Radius) > 1e-9*tc.wantR {
						t.Fatalf("点 %d 到圆心 %v，半径 %v", i, d, got.Radius)
					}
				}
			})
		}
	}
}

func TestArcSweepMatchesCase(t *testing.T) {
	for _, r := range []float64{smallR, largeR} {
		for _, tc := range cases(r) {
			if tc.wantArc == 0 {
				continue // 近直线的扫角接近 0，不按期望扫角断言。
			}
			arc, err := ArcThroughPoints(tc.a, tc.b, tc.c)
			if err != nil {
				t.Fatal(err)
			}
			if math.Abs(math.Abs(arc.Sweep)-tc.wantArc) > 1e-9 {
				t.Fatalf("%s r=%v：扫角 %v，期望 %v", tc.name, r, arc.Sweep, tc.wantArc)
			}
			// 中间点必须落在弧内，而不是补弧上。
			midAng := math.Atan2(tc.b.Y-arc.Center.Y, tc.b.X-arc.Center.X)
			ccwOff := math.Mod(midAng-arc.Start, 2*math.Pi)
			if ccwOff <= 0 {
				ccwOff += 2 * math.Pi
			}
			if arc.Sweep > 0 && ccwOff > arc.Sweep {
				t.Fatalf("%s r=%v：中间点不在逆时针弧内", tc.name, r)
			}
			if arc.Sweep < 0 && (-arc.Sweep) < 2*math.Pi-ccwOff {
				t.Fatalf("%s r=%v：中间点不在顺时针弧内", tc.name, r)
			}
		}
	}
}

func TestNearLineRadiusLarge(t *testing.T) {
	r := 1000.0
	circ, err := Circumcircle(Point{0, 0}, Point{r, 1e-6}, Point{2 * r, 0})
	if err != nil {
		t.Fatalf("近直线仍应可解：%v", err)
	}
	// 解析近似：R ≈ (弦长/2)^2 / (2·拱高) = r^2/(2e-6)
	want := (r * r) / (2e-6)
	if math.Abs(circ.Radius-want)/want > 1e-6 {
		t.Fatalf("近直线半径 %v，解析近似 %v", circ.Radius, want)
	}
}

func TestDegenerateErrorsCarryPointNumbers(t *testing.T) {
	type test struct {
		name    string
		pts     [3]Point
		wantErr error
		wantPts []int
	}
	tests := []test{
		{"三点共线", [3]Point{{0, 0}, {5, 5}, {10, 10}}, ErrCollinear, []int{0, 1, 2}},
		{"点 0 与 1 重合", [3]Point{{3, 4}, {3, 4}, {10, 0}}, ErrCoincident, []int{0, 1}},
		{"点 1 与 2 重合", [3]Point{{0, 0}, {3, 4}, {3, 4}}, ErrCoincident, []int{1, 2}},
		{"点 0 与 2 重合", [3]Point{{3, 4}, {0, 0}, {3, 4}}, ErrCoincident, []int{0, 2}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			_, err := Circumcircle(tc.pts[0], tc.pts[1], tc.pts[2])
			if !errors.Is(err, tc.wantErr) {
				t.Fatalf("错误 %v，期望 %v", err, tc.wantErr)
			}
			pts, ok := IsPointError(err)
			if !ok || len(pts) != len(tc.wantPts) {
				t.Fatalf("点号 %v，期望 %v", pts, tc.wantPts)
			}
			for i := range pts {
				if pts[i] != tc.wantPts[i] {
					t.Fatalf("点号 %v，期望 %v", pts, tc.wantPts)
				}
			}
		})
	}
}

func TestRadiusExceedsCanvas(t *testing.T) {
	canvas := Canvas{Min: Point{0, 0}, Max: Point{100, 100}}
	// 半径 1000 的圆弧，画幅对角线约 141。
	_, err := CircumcircleInCanvas(arcOnCircle(1000, 0), arcOnCircle(1000, 30), arcOnCircle(1000, 60), canvas)
	if !errors.Is(err, ErrRadiusTooLarge) {
		t.Fatalf("应报半径超画幅：%v", err)
	}
	if pts, ok := IsPointError(err); !ok || len(pts) != 3 {
		t.Fatalf("点号缺失：%v", pts)
	}
	// 画幅足够时正常返回。
	if _, err := CircumcircleInCanvas(arcOnCircle(10, 0), arcOnCircle(10, 30), arcOnCircle(10, 60), canvas); err != nil {
		t.Fatalf("画幅内不应报错：%v", err)
	}
}

func TestSegmentsForSagittaBound(t *testing.T) {
	for _, r := range []float64{smallR, largeR} {
		for _, sag := range []float64{1e-3, 1e-2, 0.1, 1} {
			if sag >= 2*r {
				continue
			}
			n, err := SegmentsForSagitta(r, sag)
			if err != nil {
				t.Fatal(err)
			}
			// 用全圆扫角验证：n 段时实际弦高 <= 公差，n-1 段时应超出。
			check := func(n int) float64 {
				theta := 2 * math.Pi / float64(n)
				return r * (1 - math.Cos(theta/2))
			}
			if got := check(n); got > sag+1e-15 {
				t.Fatalf("r=%v sag=%v n=%d 实际弦高 %v 超界", r, sag, n, got)
			}
			if n > 1 {
				if got := check(n - 1); got <= sag {
					t.Fatalf("r=%v sag=%v 段数 %d 非最小", r, sag, n)
				}
			}
		}
	}
}

func TestSagittaBadInputs(t *testing.T) {
	if _, err := SegmentsForSagitta(0, 0.1); !errors.Is(err, ErrBadRadius) {
		t.Fatal("半径 0 应报错")
	}
	if _, err := SegmentsForSagitta(math.Inf(1), 0.1); !errors.Is(err, ErrBadRadius) {
		t.Fatal("半径无穷应报错")
	}
	for _, sag := range []float64{0, -1, math.NaN(), math.Inf(1)} {
		if _, err := SegmentsForSagitta(10, sag); !errors.Is(err, ErrBadSagitta) {
			t.Fatalf("弦高 %v 应报错", sag)
		}
	}
}

func TestSegmentsForSweepMinimal(t *testing.T) {
	for _, r := range []float64{smallR, largeR} {
		for _, sweep := range []float64{math.Pi / 6, math.Pi, 5 * math.Pi / 3, 2 * math.Pi} {
			for _, sag := range []float64{1e-3, 1e-2, 0.1} {
				n, err := SegmentsForSweep(r, sweep, sag)
				if err != nil {
					t.Fatal(err)
				}
				sagittaAt := func(n int) float64 {
					theta := sweep / float64(n)
					return r * (1 - math.Cos(theta/2))
				}
				if got := sagittaAt(n); got > sag+1e-15 {
					t.Fatalf("r=%v sweep=%v sag=%v n=%d 弦高 %v 超界", r, sweep, sag, n, got)
				}
				if n > 1 {
					if got := sagittaAt(n - 1); got <= sag {
						t.Fatalf("r=%v sweep=%v sag=%v 段数 %d 非最小", r, sweep, sag, n)
					}
				}
			}
		}
	}
	if _, err := SegmentsForSweep(1, 0, 1); err == nil {
		t.Fatal("零扫角应报错")
	}
	if _, err := SegmentsForSweep(1, -1, 1); err != nil {
		t.Fatalf("负扫角应接受（表示顺时针）：%v", err)
	}
}

func TestFlattenArcDeviationBound(t *testing.T) {
	sag := 0.01
	for _, r := range []float64{smallR, largeR} {
		for _, sweepDeg := range []float64{30, 180, 300} {
			arc, err := ArcThroughPoints(
				arcOnCircle(r, 0), arcOnCircle(r, sweepDeg/2), arcOnCircle(r, sweepDeg))
			if err != nil {
				t.Fatalf("r=%v sweep=%v：%v", r, sweepDeg, err)
			}
			pts, n, err := FlattenArc(arc, sag)
			if err != nil {
				t.Fatal(err)
			}
			rep, err := arc.Report(pts, sag)
			if err != nil {
				t.Fatal(err)
			}
			if n != rep.Segments {
				t.Fatalf("段数不一致 %d vs %d", n, rep.Segments)
			}
			// 顶点严格在圆上：径向偏差应为舍入量级。
			if rep.MaxDeviation > 1e-9*r {
				t.Fatalf("顶点径向偏差 %v 过大", rep.MaxDeviation)
			}
			// 弦中点径向偏差（即弦高误差）不超过公差。
			if rep.MaxSagitta > sag*(1+1e-9) {
				t.Fatalf("弦高 %v 超过公差 %v", rep.MaxSagitta, sag)
			}
			// 端点保真。
			for _, q := range []Point{pts[0], pts[n]} {
				if dv := arc.RadialDeviation(q); dv > 1e-9*r {
					t.Fatalf("端点偏差 %v", dv)
				}
			}
		}
	}
}

func TestArcDirectionFollowsMiddlePoint(t *testing.T) {
	// 逆时针：0° → 90° 经 45°。
	arc, err := ArcThroughPoints(arcOnCircle(1, 0), arcOnCircle(1, 45), arcOnCircle(1, 90))
	if err != nil {
		t.Fatal(err)
	}
	if arc.Sweep <= 0 || math.Abs(arc.Sweep-math.Pi/2) > 1e-12 {
		t.Fatalf("逆时针扫角 %v", arc.Sweep)
	}
	// 顺时针：0° → -90° 经 -45°，等价 0° → 270° 经 315°。
	arc, err = ArcThroughPoints(arcOnCircle(1, 0), arcOnCircle(1, -45), arcOnCircle(1, -90))
	if err != nil {
		t.Fatal(err)
	}
	if arc.Sweep >= 0 || math.Abs(arc.Sweep+math.Pi/2) > 1e-12 {
		t.Fatalf("顺时针扫角 %v", arc.Sweep)
	}
}

func TestFlattenUnchanged(t *testing.T) {
	pts, err := Flatten(Point{}, Point{Y: 10}, Point{X: 10, Y: 10}, Point{X: 10}, 8)
	if err != nil || len(pts) != 9 {
		t.Fatalf("细分异常 %v %d", err, len(pts))
	}
}

func TestArcValidate(t *testing.T) {
	for _, bad := range []Arc{
		{Radius: 0, Start: 0, Sweep: 1},
		{Radius: 1, Start: 0, Sweep: 0},
		{Radius: 1, Start: 0, Sweep: 3 * math.Pi},
	} {
		if err := bad.Validate(); err == nil {
			t.Fatalf("非法弧未报错 %+v", bad)
		}
	}
}
