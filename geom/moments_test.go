package geom

import (
	"errors"
	"math"
	"testing"
	"time"
)

// rot 旋转 deg 度后平移 (dx,dy)。
func rot(pts []Point, deg, dx, dy float64) []Point {
	r := deg * math.Pi / 180
	c, s := math.Cos(r), math.Sin(r)
	out := make([]Point, len(pts))
	for i, p := range pts {
		out[i] = Point{c*p.X - s*p.Y + dx, s*p.X + c*p.Y + dy}
	}
	return out
}

func rectPts(w, h float64) []Point {
	return []Point{{-w / 2, -h / 2}, {w / 2, -h / 2}, {w / 2, h / 2}, {-w / 2, h / 2}}
}

// angDiff 无向轴夹角，[0,90]。
func angDiff(a, b float64) float64 {
	d := math.Abs(math.Mod(a-b, 180))
	if d > 90 {
		d = 180 - d
	}
	return d
}

func closePt(a, b Point, tol float64) bool {
	return math.Hypot(a.X-b.X, a.Y-b.Y) <= tol
}

// 凸多边形：40x10 矩形旋转 30° 平移到 (5,7)。
func TestConvexRectMoments(t *testing.T) {
	poly := rot(rectPts(40, 10), 30, 5, 7)
	m, err := PolygonMoments(poly)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(m.Area-400) > 1e-9 {
		t.Fatalf("面积错误 %v", m.Area)
	}
	if !closePt(m.Centroid, Point{5, 7}, 1e-9) {
		t.Fatalf("质心错误 %v", m.Centroid)
	}
	if d := angDiff(m.AxisAngle, 30); d > 1e-9 {
		t.Fatalf("主轴角错误 %v 偏差 %v", m.AxisAngle, d)
	}
	c, long, short, ra, err := MinAreaRect(poly)
	if err != nil {
		t.Fatal(err)
	}
	if !closePt(c, Point{5, 7}, 1e-9) || math.Abs(long-40) > 1e-9 || math.Abs(short-10) > 1e-9 {
		t.Fatalf("最小矩形错误 center=%v long=%v short=%v", c, long, short)
	}
	if d := angDiff(m.AxisAngle, ra); d > 1 {
		t.Fatalf("主轴与长边方向偏差超一度: %v°", d)
	}
	t.Logf("凸矩形   %s 长边方向=%.4f°", m, ra)
}

// 凹多边形：关于 x 轴对称的飞镖形，质心-主轴变换不变量前后对照。
func TestTransformInvariance(t *testing.T) {
	dart := []Point{{0, 0}, {12, -2}, {9, 0}, {12, 2}} // 凹点 (9,0)
	m0, err := PolygonMoments(dart)
	if err != nil {
		t.Fatal(err)
	}
	moved := rot(dart, 55, 100, -37)
	m1, err := PolygonMoments(moved)
	if err != nil {
		t.Fatal(err)
	}
	r := 55 * math.Pi / 180
	wantC := Point{
		math.Cos(r)*m0.Centroid.X - math.Sin(r)*m0.Centroid.Y + 100,
		math.Sin(r)*m0.Centroid.X + math.Cos(r)*m0.Centroid.Y - 37,
	}
	if !closePt(m1.Centroid, wantC, 1e-9) {
		t.Fatalf("质心未随变换移动: got %v want %v", m1.Centroid, wantC)
	}
	wantAxis := math.Mod(m0.AxisAngle+55, 180)
	if d := angDiff(m1.AxisAngle, wantAxis); d > 1e-9 {
		t.Fatalf("主轴角未随旋转移动: got %v want %v", m1.AxisAngle, wantAxis)
	}
	if _, _, _, ra, err := MinAreaRect(dart); err != nil || angDiff(m0.AxisAngle, ra) > 1 {
		t.Fatalf("凹多边形主轴与长边方向不符: 主轴=%v 长边=%v err=%v", m0.AxisAngle, ra, err)
	}
	t.Logf("前后对照 前: %s", m0)
	t.Logf("前后对照 后: %s 期望质心=(%.6f,%.6f) 期望主轴角=%.4f°", m1, wantC.X, wantC.Y, wantAxis)
}

// 带洞：40x10 外环逆时针 + 8x4 洞顺时针，同中心同旋转。
func TestRegionWithHole(t *testing.T) {
	outer := rot(rectPts(40, 10), 30, 5, 7)
	hole := rot([]Point{{-4, -2}, {-4, 2}, {4, 2}, {4, -2}}, 30, 5, 7) // 顺时针
	m, err := RegionMoments(outer, hole)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(m.Area-368) > 1e-9 {
		t.Fatalf("带洞面积错误 %v", m.Area)
	}
	if !closePt(m.Centroid, Point{5, 7}, 1e-9) {
		t.Fatalf("带洞质心错误 %v", m.Centroid)
	}
	if d := angDiff(m.AxisAngle, 30); d > 1e-9 {
		t.Fatalf("带洞主轴角错误 %v", m.AxisAngle)
	}
	if _, _, _, ra, err := MinAreaRect(outer); err != nil || angDiff(m.AxisAngle, ra) > 1 {
		t.Fatalf("带洞主轴与长边方向不符: 主轴=%v 长边=%v err=%v", m.AxisAngle, ra, err)
	}
	// 洞方向反了应报 ErrBadWinding
	if _, err := RegionMoments(outer, rot(rectPts(8, 4), 30, 5, 7)); !errors.Is(err, ErrBadWinding) {
		t.Fatalf("洞方向错误应报 ErrBadWinding, got %v", err)
	}
	t.Logf("带洞     %s", m)
}

// 退化与错误输入：负面积、顶点不足、退化共线、零面积，各报点号。
func TestErrorCases(t *testing.T) {
	var ve *VertexError

	_, err := PolygonMoments([]Point{{0, 0}, {1, 1}})
	if !errors.Is(err, ErrTooFewPts) || !errors.As(err, &ve) || ve.Index != 2 {
		t.Fatalf("顶点不足: %v", err)
	}
	t.Logf("顶点不足   -> %v", err)

	cw := []Point{{0, 0}, {0, 10}, {10, 10}, {10, 0}} // 顺时针
	_, err = PolygonMoments(cw)
	if !errors.Is(err, ErrNegativeArea) || !errors.As(err, &ve) {
		t.Fatalf("负面积: %v", err)
	}
	t.Logf("面积为负   -> %v", err)

	col := []Point{{0, 0}, {1, 1}, {2, 2}, {3, 3}} // 共线
	_, err = PolygonMoments(col)
	if !errors.Is(err, ErrDegenerate) || !errors.As(err, &ve) {
		t.Fatalf("退化共线: %v", err)
	}
	t.Logf("退化共线   -> %v", err)

	bow := []Point{{0, 0}, {2, 2}, {2, 0}, {0, 2}} // 自交，有向面积 0
	_, err = PolygonMoments(bow)
	if !errors.Is(err, ErrDegenerate) || !errors.As(err, &ve) {
		t.Fatalf("零面积: %v", err)
	}
	t.Logf("零面积     -> %v", err)
}

// 两档顶点数：旋转 25° 的椭圆离散多边形，各报质心、主轴角与耗时。
func TestVertexCountTiers(t *testing.T) {
	for _, n := range []int{360, 7200} {
		pts := make([]Point, n)
		for i := range pts {
			th := 2 * math.Pi * float64(i) / float64(n)
			pts[i] = Point{50 * math.Cos(th), 10 * math.Sin(th)}
		}
		pts = rot(pts, 25, 3, -2)

		start := time.Now()
		m, err := PolygonMoments(pts)
		elMom := time.Since(start)
		if err != nil {
			t.Fatal(err)
		}
		start = time.Now()
		_, _, _, ra, err := MinAreaRect(pts)
		elRect := time.Since(start)
		if err != nil {
			t.Fatal(err)
		}
		if !closePt(m.Centroid, Point{3, -2}, 1e-6) {
			t.Fatalf("n=%d 质心错误 %v", n, m.Centroid)
		}
		if d := angDiff(m.AxisAngle, 25); d > 1 {
			t.Fatalf("n=%d 主轴角偏差超一度: %v°", n, d)
		}
		if d := angDiff(m.AxisAngle, ra); d > 1 {
			t.Fatalf("n=%d 主轴与长边方向偏差超一度: %v°", n, d)
		}
		t.Logf("顶点数=%d 质心=(%.6f,%.6f) 主轴角=%.4f° 长边方向=%.4f° 矩耗时=%v 矩形耗时=%v",
			n, m.Centroid.X, m.Centroid.Y, m.AxisAngle, ra, elMom, elRect)
	}
}
