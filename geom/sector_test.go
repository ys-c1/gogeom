package geom

import (
	"errors"
	"math"
	"strings"
	"testing"
)

// arcPts 生成圆心 c、半径 r、扫角 theta 的 n 段内接折线顶点(含端点共 n+1 个)。
func arcPts(c Point, r, theta float64, n int) []Point {
	pts := make([]Point, n+1)
	for i := 0; i <= n; i++ {
		a := theta * float64(i) / float64(n)
		pts[i] = Point{c.X + r*math.Cos(a), c.Y + r*math.Sin(a)}
	}
	return pts
}

// 折线近似:扇形 = 多边形 (c, 弧上点…);弓形 = 多边形 (弧上点…, 弦闭合)。复用现有 Area。
func polySectorArea(c Point, r, theta float64, n int) float64 {
	pts := append([]Point{c}, arcPts(c, r, theta, n)...)
	a, _ := Area(pts)
	return a
}

func polySegmentArea(c Point, r, theta float64, n int) float64 {
	a, _ := Area(arcPts(c, r, theta, n))
	return a
}

func TestSectorSegmentFormulas(t *testing.T) {
	c := Point{}
	r := 2.0
	cases := []struct {
		name  string
		theta float64
	}{
		{"小角", 0.1},
		{"直角", math.Pi / 2},
		{"半圆", math.Pi},
		{"大角", 3 * math.Pi / 2},
		{"近整圆", 2*math.Pi - 1e-6},
	}
	for _, tc := range cases {
		sec, err := SectorArea(r, tc.theta)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if want := r * r * tc.theta / 2; sec != want {
			t.Fatalf("%s: 扇形 %v != %v", tc.name, sec, want)
		}
		seg, err := SegmentArea(r, tc.theta)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		// 弓形 = 扇形 − 圆心三角形
		tri := r * r * math.Sin(tc.theta) / 2
		if got := sec - tri; math.Abs(seg-got) > 1e-12*math.Abs(seg) {
			t.Fatalf("%s: 弓形 %v != 扇形−三角形 %v", tc.name, seg, got)
		}
		// 与细折线近似对照 (n = 200000)
		if got := polySectorArea(c, r, tc.theta, 200000); math.Abs(got-sec) > 1e-7*math.Abs(sec) {
			t.Fatalf("%s: 扇形折线近似 %v 偏离解析值 %v", tc.name, got, sec)
		}
		if got := polySegmentArea(c, r, tc.theta, 200000); math.Abs(got-seg) > 1e-7*math.Abs(seg) {
			t.Fatalf("%s: 弓形折线近似 %v 偏离解析值 %v", tc.name, got, seg)
		}
	}
	// 已知值: r=2, 直角 → 扇形 π, 弓形 π−2
	if sec, _ := SectorArea(2, math.Pi/2); math.Abs(sec-math.Pi) > 1e-15 {
		t.Fatalf("直角扇形 %v != π", sec)
	}
	if seg, _ := SegmentArea(2, math.Pi/2); math.Abs(seg-(math.Pi-2)) > 1e-15 {
		t.Fatalf("直角弓形 %v != π−2", seg)
	}
}

func TestSignedArcArea(t *testing.T) {
	c := Point{1, -2}
	r, theta := 3.0, 2.4
	pts := arcPts(c, r, theta, 256)
	got, err := SignedArcArea(c, pts)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := SectorArea(r, theta)
	if math.Abs(got-want) > 1e-9*want {
		t.Fatalf("带符号面积 %v, 期望 %v", got, want)
	}
	// 反向顶点 → 符号翻转
	rev := make([]Point, len(pts))
	for i, p := range pts {
		rev[len(pts)-1-i] = p
	}
	neg, _ := SignedArcArea(c, rev)
	if math.Abs(neg+got) > 1e-9*want {
		t.Fatalf("反向未翻号: %v vs %v", got, neg)
	}
}

func TestRejectBadInput(t *testing.T) {
	if _, err := SectorArea(0, 1); !errors.Is(err, ErrBadRadius) || !strings.Contains(err.Error(), "radius") {
		t.Fatalf("半径为零应报 radius 字段: %v", err)
	}
	if _, err := SegmentArea(-1, 1); !errors.Is(err, ErrBadRadius) {
		t.Fatalf("半径为负应报 ErrBadRadius: %v", err)
	}
	if _, err := SectorArea(1, 3*math.Pi); !errors.Is(err, ErrBadAngle) || !strings.Contains(err.Error(), "theta") {
		t.Fatalf("角度超一周应报 theta 字段: %v", err)
	}
	if _, err := SegmentArea(1, -2*math.Pi-1e-9); !errors.Is(err, ErrBadAngle) {
		t.Fatalf("负向超一周应报 ErrBadAngle: %v", err)
	}
	// 顶点与圆心共线(圆心位于两端点之间)
	if _, err := SignedArcArea(Point{}, []Point{{1, 0}, {-1, 0}}); !errors.Is(err, ErrCollinear) || !strings.Contains(err.Error(), "pts") {
		t.Fatalf("顶点共线应报 pts 字段: %v", err)
	}
	// 顶点即圆心 → 半径为零
	if _, err := SignedArcArea(Point{}, []Point{{0, 0}, {1, 0}}); !errors.Is(err, ErrBadRadius) {
		t.Fatalf("顶点落在圆心应报 ErrBadRadius: %v", err)
	}
}

// 折线近似误差: E(n) = ½r²[θ − n·sin(θ/n)] = r²θ³/(12n²) + O(n⁻⁴) → 二阶收敛。
func TestConvergenceOrder(t *testing.T) {
	theta := 2.0
	ns := []int{250, 500, 1000, 2000, 4000, 8000}
	for _, r := range []float64{0.5, 500} {
		sec, _ := SectorArea(r, theta)
		seg, _ := SegmentArea(r, theta)
		var lx, ly, lxx, lxy float64 // log2 空间最小二乘拟合斜率
		k := 0
		t.Logf("r=%v θ=%v: n, 扇形误差, 弓形误差", r, theta)
		for _, n := range ns {
			es := math.Abs(polySectorArea(Point{}, r, theta, n) - sec)
			eg := math.Abs(polySegmentArea(Point{}, r, theta, n) - seg)
			t.Logf("  n=%6d 扇形误差=%.6e 弓形误差=%.6e", n, es, eg)
			x, y := math.Log2(float64(n)), math.Log2(es)
			lx += x
			ly += y
			lxx += x * x
			lxy += x * y
			k++
		}
		fk := float64(k)
		slope := (fk*lxy - lx*ly) / (fk*lxx - lx*lx)
		t.Logf("r=%v 实测阶数 = %.4f (推导阶数 2)", r, -slope)
		if math.Abs(slope+2) > 0.05 {
			t.Fatalf("r=%v 实测阶数 %.4f 偏离推导阶数 2", r, -slope)
		}
		// 十万点: 相对误差低于阈值
		const n = 100000
		es := math.Abs(polySectorArea(Point{}, r, theta, n)-sec) / sec
		eg := math.Abs(polySegmentArea(Point{}, r, theta, n)-seg) / seg
		t.Logf("r=%v n=%d 扇形相对误差=%.3e 弓形相对误差=%.3e", r, n, es, eg)
		if es > 1e-8 || eg > 1e-8 {
			t.Fatalf("r=%v 十万点相对误差 扇形=%.3e 弓形=%.3e 超阈值 1e-8", r, es, eg)
		}
	}
}
