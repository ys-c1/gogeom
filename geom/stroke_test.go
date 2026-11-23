package geom

import (
	"errors"
	"fmt"
	"math"
	"strings"
	"testing"
	"time"
)

var (
	rightAngle = []Point{{0, 0}, {10, 0}, {10, 10}}
	acute      = []Point{{0, 0}, {10, 0}, {4, 2}}
	foldback   = []Point{{0, 0}, {10, 0}, {4, 0.5}}
	spiral     = []Point{{0, 0}, {14, 0}, {14, 14}, {0, 14}, {0, 4}, {10, 4}, {10, 10}, {4, 10}, {4, 7}, {7, 7}}
	square     = []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}, {0, 0}}
)

// 原始直折描边行为不变：直角折线的朴素描边多边形逐点固定。
func TestStrokeBaselineUnchanged(t *testing.T) {
	got, err := Stroke(rightAngle, 2)
	if err != nil {
		t.Fatal(err)
	}
	want := []Point{{0, 1}, {9, 1}, {9, 10}, {11, 10}, {11, -1}, {0, -1}}
	if len(got) != len(want) {
		t.Fatalf("顶点数 %d != %d", len(got), len(want))
	}
	for i, p := range want {
		if got[i] != p {
			t.Fatalf("顶点 %d: got %v want %v", i, got[i], p)
		}
	}
}

// 直角：无自交，面积 = 线宽 × (L1+L2) = 20w，与手算一致。
func TestOutlineRightAngle(t *testing.T) {
	for _, w := range []float64{2, 6} {
		res, err := Outline(rightAngle, w)
		if err != nil {
			t.Fatal(err)
		}
		if res.Intersections != 0 || res.Discarded != 0 || len(res.Shells) != 1 {
			t.Fatalf("w=%v 不应有自交: %+v", w, res)
		}
		if want := 20 * w; math.Abs(res.Area-want) > 1e-9 {
			t.Fatalf("w=%v 面积 %v != 手算 %v", w, res.Area, want)
		}
	}
}

// 锐角：miter 超限退化 bevel，不得报错；面积守恒。
func TestOutlineAcute(t *testing.T) {
	for _, w := range []float64{2, 6} {
		res, err := Outline(acute, w)
		if err != nil {
			t.Fatal(err)
		}
		checkConservation(t, acute, w)
		if res.Area <= 0 {
			t.Fatalf("w=%v 面积应正: %v", w, res.Area)
		}
	}
}

// 回折：两档线宽都绕进去，必须切分并丢环；面积守恒。
func TestOutlineFoldback(t *testing.T) {
	for _, w := range []float64{2, 6} {
		res, err := Outline(foldback, w)
		if err != nil {
			t.Fatal(err)
		}
		if res.Intersections < 1 || res.Discarded < 1 {
			t.Fatalf("w=%v 回折应有绕入环: %+v", w, res)
		}
		checkConservation(t, foldback, w)
	}
}

// 螺旋：大线宽自交多于小线宽；面积守恒。
func TestOutlineSpiral(t *testing.T) {
	small, err := Outline(spiral, 2)
	if err != nil {
		t.Fatal(err)
	}
	big, err := Outline(spiral, 6)
	if err != nil {
		t.Fatal(err)
	}
	if big.Intersections <= small.Intersections {
		t.Fatalf("大线宽应交点更多: w=2 → %d, w=6 → %d", small.Intersections, big.Intersections)
	}
	if big.Discarded < 1 {
		t.Fatalf("w=6 应丢绕入环: %+v", big)
	}
	checkConservation(t, spiral, 2)
	checkConservation(t, spiral, 6)
}

// 闭合方环：面积 = 4·边长·线宽（外方减内方），与手算一致。
func TestOutlineClosedSquare(t *testing.T) {
	for _, w := range []float64{2, 6} {
		res, err := Outline(square, w)
		if err != nil {
			t.Fatal(err)
		}
		if len(res.Shells) != 1 || len(res.Holes) != 1 {
			t.Fatalf("w=%v 应为一壳一孔: %+v", w, res)
		}
		if want := 4 * 10 * w; math.Abs(res.Area-want) > 1e-9 {
			t.Fatalf("w=%v 面积 %v != 手算 %v", w, res.Area, want)
		}
	}
}

// 守恒：切分后全部环（含被丢的）有向面积和 == 原始描边多边形鞋带面积。
func checkConservation(t *testing.T, pts []Point, w float64) {
	t.Helper()
	raw, err := Stroke(pts, w)
	if err != nil {
		t.Fatal(err)
	}
	a0, err := Area(raw)
	if err != nil {
		t.Fatal(err)
	}
	loops, _ := splitLoops(raw)
	if sum := signedSum(loops); math.Abs(sum-a0) > 1e-6*math.Max(1, math.Abs(a0)) {
		t.Fatalf("w=%v 面积不守恒: 环和 %v != 原始 %v", w, sum, a0)
	}
}

// 线宽为零、法向退化、自交未拆分：各报段号。
func TestStrokeErrors(t *testing.T) {
	var se *StrokeError

	if _, err := Outline(rightAngle, 0); !errors.As(err, &se) || se.Op != "线宽为零" {
		t.Fatalf("线宽为零应报 StrokeError: %v", err)
	}

	degen := []Point{{0, 0}, {5, 0}, {5, 0}, {8, 0}}
	if _, err := Outline(degen, 2); !errors.As(err, &se) || se.Op != "法向退化" || se.Seg != 1 {
		t.Fatalf("法向退化应报段 1: %v", err)
	}

	bowtie := []Point{{0, 0}, {2, 2}, {0, 2}, {2, 0}}
	if e := checkSimple(bowtie); e == nil || e.Op != "自交未拆分" || e.Seg != 0 || e.SegB != 2 {
		t.Fatalf("自交未拆分应报段 0 与段 2: %v", e)
	}
}

// 两档线宽 × 五种用例：环数、交点数、面积与耗时对照。
func TestOutlineReport(t *testing.T) {
	cases := []struct {
		name string
		pts  []Point
		hand func(w float64) float64 // 手算面积，无则 nil
	}{
		{"直角", rightAngle, func(w float64) float64 { return 20 * w }},
		{"锐角", acute, nil},
		{"回折", foldback, nil},
		{"螺旋", spiral, nil},
		{"闭合", square, func(w float64) float64 { return 40 * w }},
	}
	t.Log("用例  线宽  交点  保留环  丢弃环  孔  面积(裁剪)  面积(原始)  手算值  耗时")
	for _, c := range cases {
		for _, w := range []float64{2, 6} {
			start := time.Now()
			res, err := Outline(c.pts, w)
			el := time.Since(start)
			if err != nil {
				t.Fatalf("%s w=%v: %v", c.name, w, err)
			}
			raw, _ := Stroke(c.pts, w)
			a0, _ := Area(raw)
			hand := "-"
			if c.hand != nil {
				h := c.hand(w)
				hand = f2(h)
				if math.Abs(res.Area-h) > 1e-9 {
					t.Errorf("%s w=%v 裁剪面积 %v != 手算 %v", c.name, w, res.Area, h)
				}
			}
			t.Logf("%s  %.1f  %4d  %5d  %5d  %2d  %9s  %9s  %7s  %v",
				c.name, w, res.Intersections, len(res.Shells), res.Discarded,
				len(res.Holes), f2(res.Area), f2(math.Abs(a0)), hand, el)
		}
	}
}

func f2(v float64) string {
	return strings.TrimRight(strings.TrimRight(fmt.Sprintf("%.2f", v), "0"), ".")
}
