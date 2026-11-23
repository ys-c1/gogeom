// cut_test.go 按比切分的验收测试：
// 凸块、凹块、细长、带洞、近退化五种形状，以及比值越界、自交、共线三类
// 非法输入（校验 EdgeError 边号）；面积守恒与简单性两条独立断言。
package geom

import (
	"errors"
	"math"
	"testing"
)

const ratioTol = 1e-7

// shapes 五类形状（逆时针外轮廓；带洞者孔洞顺时针）。
func shapes() map[string]*Polygon {
	return map[string]*Polygon{
		// 凸块：10x10 正方形。
		"convex": {Outer: []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}},
		// 凹块：L 形。
		"concave": {Outer: []Point{{0, 0}, {10, 0}, {10, 4}, {4, 4}, {4, 10}, {0, 10}}},
		// 细长：100:1 矩形。
		"slender": {Outer: []Point{{0, 0}, {100, 0}, {100, 1}, {0, 1}}},
		// 带洞：10x10 方框挖内嵌大矩形洞 x∈[3,7],y∈[2,8]（面积 24，净 76）。
		// 洞在切深方向足够宽，使 0.4/0.5/0.6 三个比例的切线都贯穿孔洞，
		// 两块均可缝合成无洞的简单多边形。
		"holed": {
			Outer: []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}},
			Holes: [][]Point{{{3, 2}, {3, 8}, {7, 8}, {7, 2}}},
		},
		// 近退化：极扁四边形（高 1e-6）。
		"nearDegenerate": {Outer: []Point{{0, 0}, {10, 0}, {10, 1e-6}, {0, 0.5e-6}}},
	}
}

// ratiosFor 每个形状各自合法的测试比例。带洞形状只能取“切线贯穿孔洞”
// 的比例；其它比例会因单一直线无法把含洞块切成简单环而被正确拒绝
// （见 TestHoleNotPiercedRejected）。
func ratiosFor(name string) []float64 {
	if name == "holed" {
		return []float64{0.4, 0.5, 0.6}
	}
	return []float64{0.25, 0.5, 0.73}
}

func TestCutAllShapes(t *testing.T) {
	for name, p := range shapes() {
		for _, r := range ratiosFor(name) {
			res, err := CutPolygon(p, r, ratioTol)
			if err != nil {
				t.Fatalf("[%s ratio=%v] 切分失败: %v", name, r, err)
			}
			// 验收 1：面积比落入容差。
			if res.AbsErr > ratioTol {
				t.Fatalf("[%s ratio=%v] 面积比误差 %g 超容差", name, r, res.AbsErr)
			}
			// 验收 2（独立）：两块面积之和等于原面积。
			if !res.AreaSumOK {
				t.Fatalf("[%s ratio=%v] 面积不守恒: %g+%g != %g",
					name, r, res.AreaNeg, res.AreaPos, res.AreaAll)
			}
			// 验收 3（独立）：两块都是简单多边形。
			if !res.PiecesOK || !IsSimple(res.Neg) || !IsSimple(res.Pos) {
				t.Fatalf("[%s ratio=%v] 切块存在非简单多边形", name, r)
			}
			// 每块面积必须为正。
			if res.AreaNeg <= 0 || res.AreaPos <= 0 {
				t.Fatalf("[%s ratio=%v] 出现零面积切块", name, r)
			}
			// 法向与方向单位正交。
			n, u := res.Line.N, res.Line.Dir
			if math.Abs(n.X*n.X+n.Y*n.Y-1) > 1e-12 ||
				math.Abs(u.X*u.X+u.Y*u.Y-1) > 1e-12 ||
				math.Abs(n.X*u.X+n.Y*u.Y) > 1e-12 {
				t.Fatalf("[%s] 切线方向未单位正交: %+v", name, res.Line)
			}
			t.Logf("[%s ratio=%.2f] %s  面积 neg=%.6f pos=%.6f 原=%.6f 实际比=%.8f 误差=%.2e",
				name, r, res.Line, res.AreaNeg, res.AreaPos, res.AreaAll, res.Ratio, res.AbsErr)
		}
	}
}

// TestCutLineEquationPoints 切线上的拼装点必须精确落在直线上。
func TestCutLineEquationPoints(t *testing.T) {
	for name, p := range shapes() {
		r := ratiosFor(name)[0] // 取各形状的合法比例
		res, err := CutPolygon(p, r, ratioTol)
		if err != nil {
			t.Fatalf("[%s] %v", name, err)
		}
		onLine := func(q Point) bool {
			return math.Abs(res.Line.N.X*q.X+res.Line.N.Y*q.Y-res.Line.D) <= 1e-7
		}
		// 两块各自应至少有两个切线上的点（共享切点）。
		var nNeg, nPos int
		for _, q := range res.Neg {
			if onLine(q) {
				nNeg++
			}
		}
		for _, q := range res.Pos {
			if onLine(q) {
				nPos++
			}
		}
		if nNeg < 2 || nPos < 2 {
			t.Fatalf("[%s] 切块上的切点数不足 neg=%d pos=%d", name, nNeg, nPos)
		}
	}
}

// TestRatioOutOfRange 比值越界拒绝（无边可指，返回哨兵错误）。
func TestRatioOutOfRange(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	for _, r := range []float64{0, 1, -0.5, 1.2, math.NaN()} {
		if _, err := CutRing(sq, r, ratioTol); !errors.Is(err, ErrRatioOutOfRange) {
			t.Fatalf("ratio=%v 应拒绝为 ErrRatioOutOfRange，得到 %v", r, err)
		}
	}
}

// TestSelfIntersectRejected 蝴蝶结自交：报 ErrSelfIntersect 与边号。
func TestSelfIntersectRejected(t *testing.T) {
	bowtie := []Point{{0, 0}, {10, 10}, {10, 0}, {0, 10}}
	_, err := CutRing(bowtie, 0.5, ratioTol)
	var ee *EdgeError
	if !errors.As(err, &ee) || !errors.Is(ee.Kind, ErrSelfIntersect) {
		t.Fatalf("自交多边形应报 EdgeError/ErrSelfIntersect，得到 %v", err)
	}
	if ee.Edge != 0 || ee.Ring != 0 {
		t.Fatalf("自交边号应为 0，得到 ring=%d edge=%d", ee.Ring, ee.Edge)
	}
	t.Logf("自交拒绝信息: %v", err)
}

// TestCollinearRejected 连续三顶点共线：报 ErrCollinear 与边号。
func TestCollinearRejected(t *testing.T) {
	// 在正方形上边中点插一个共线点 (5,10)，它是顶点 4，出射边为边 4。
	coll := []Point{{0, 0}, {10, 0}, {10, 10}, {5, 10}, {0, 10}}
	_, err := CutRing(coll, 0.5, ratioTol)
	var ee *EdgeError
	if !errors.As(err, &ee) || !errors.Is(ee.Kind, ErrCollinear) {
		t.Fatalf("共线输入应报 EdgeError/ErrCollinear，得到 %v", err)
	}
	if ee.Edge != 3 {
		t.Fatalf("共线边号应为 3（中间顶点 (5,10) 的入射/出射共线边），得到 %d", ee.Edge)
	}
	t.Logf("共线拒绝信息: %v", err)

	// 重复点（零长边）同样按共线拒绝，边号即零长边。
	dup := []Point{{0, 0}, {10, 0}, {10, 0}, {10, 10}, {0, 10}}
	_, err = CutRing(dup, 0.5, ratioTol)
	if !errors.As(err, &ee) || !errors.Is(ee.Kind, ErrCollinear) || ee.Edge != 1 {
		t.Fatalf("重复点应报边 1 共线，得到 %v", err)
	}
}

// TestHoleErrors 孔洞越界与孔洞交叠均拒绝并报孔洞环号。
func TestHoleErrors(t *testing.T) {
	out := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	bad := &Polygon{Outer: out, Holes: [][]Point{{{20, 20}, {21, 20}, {21, 21}}}}
	_, err := CutPolygon(bad, 0.5, ratioTol)
	var ee *EdgeError
	if !errors.As(err, &ee) || !errors.Is(ee.Kind, ErrInvalidHole) || ee.Ring != 1 {
		t.Fatalf("越界孔洞应报环 1，得到 %v", err)
	}
}

// TestHoleNotPiercedRejected 当目标比例对应的切线够不到孔洞时，含洞那一块
// 无法成为无洞简单环，必须拒绝（几何上单一直线不可切），并指明孔洞环号。
func TestHoleNotPiercedRejected(t *testing.T) {
	p := shapes()["holed"]
	// 0.1 的切线落在洞（x∈[3,7]）左侧，整块洞都在另一侧。
	_, err := CutPolygon(p, 0.1, ratioTol)
	if !errors.Is(err, ErrUncuttable) {
		t.Fatalf("切线未穿洞应报 ErrUncuttable，得到 %v", err)
	}
	t.Logf("未穿洞拒绝信息: %v", err)
	// 对照组：穿洞比例应成功。
	if _, err := CutPolygon(p, 0.5, ratioTol); err != nil {
		t.Fatalf("穿洞比例 0.5 应成功，得到 %v", err)
	}
}

// TestExistingAPIUnchanged 原有接口行为保持不变。
func TestExistingAPIUnchanged(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	a, err := Area(sq)
	if err != nil || a != 100 {
		t.Fatalf("Area 行为变化: %v %v", a, err)
	}
	if !Inside(sq, 5, 5) || Inside(sq, 15, 5) {
		t.Fatal("Inside 行为变化")
	}
	if _, err := Area(sq[:2]); err == nil {
		t.Fatal("Area 短输入行为变化")
	}
	pts, err := Flatten(Point{}, Point{Y: 10}, Point{X: 10, Y: 10}, Point{X: 10}, 4)
	if err != nil || len(pts) != 5 {
		t.Fatalf("Flatten 行为变化: %v", err)
	}
}
