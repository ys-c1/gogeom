package geom

import (
	"errors"
	"math"
	"strings"
	"testing"
)

// arcCubic 构造圆心 (cx,cy)、半径 r、从角 a0 到 a1 的圆弧三次贝塞尔逼近。
func arcCubic(cx, cy, r, a0, a1 float64) Cubic {
	k := 4.0 / 3 * math.Tan((a1-a0)/4)
	p := func(a float64) Point { return Point{cx + r*math.Cos(a), cy + r*math.Sin(a)} }
	d := func(a float64) Point { return Point{-math.Sin(a), math.Cos(a)} }
	return Cubic{
		P0: p(a0),
		P1: padd(p(a0), pscale(d(a0), k*r)),
		P2: psub(p(a1), pscale(d(a1), k*r)),
		P3: p(a1),
	}
}

// shift 平移曲线使起点落在 to。
func shift(c Cubic, to Point) Cubic {
	d := psub(to, c.P0)
	return Cubic{padd(c.P0, d), padd(c.P1, d), padd(c.P2, d), padd(c.P3, d)}
}

// 直线接弧：直线末端水平，弧的起点切向抬高 30°，接缝原本有折角。
func lineArcCase() []Cubic {
	line := Cubic{Point{0, 0}, Point{10, 0}, Point{20, 0}, Point{30, 0}}
	arc := arcCubic(0, 0, 20, -math.Pi/2+math.Pi/6, 0) // 起点切向与水平成 30°
	return []Cubic{line, shift(arc, Point{30, 0})}
}

// 弧接弧：两段半径不同的弧，第二段绕接缝旋转 25° 制造折角。
func arcArcCase() []Cubic {
	a := arcCubic(0, 0, 30, -math.Pi/2, 0) // 终点 (30,0)，切向竖直向上
	b := arcCubic(0, 0, 12, -math.Pi/2+deg(25), deg(60))
	return []Cubic{a, shift(b, Point{30, 0})}
}

// sCase：S 形，两段弯曲方向相反，接缝两侧曲率异号。
func sCase() []Cubic {
	up := Cubic{Point{0, 0}, Point{10, 12}, Point{20, 12}, Point{30, 0}}
	down := Cubic{Point{30, 0}, Point{38, -14}, Point{52, -14}, Point{60, 0}}
	return []Cubic{up, down}
}

// foldCase：回折，第二段起点切向朝后（接近 180° 折返）。
func foldCase() []Cubic {
	a := Cubic{Point{0, 0}, Point{8, 8}, Point{22, 8}, Point{30, 0}}
	b := Cubic{Point{30, 0}, Point{22, -6}, Point{44, -10}, Point{55, -2}}
	return []Cubic{a, b}
}

func deg(r float64) float64 { return r * math.Pi / 180 }

// checkJoined 校验三档拼接后的两个量：切向夹角须为零，曲率差须在上界内。
func checkJoined(t *testing.T, name string, segs []Cubic) *JoinReport {
	t.Helper()
	out, rep, err := Join(segs, JoinCurvatureIter)
	if err != nil {
		t.Fatalf("%s: Join 报错 %v", name, err)
	}
	if len(rep.Levels) != 3 {
		t.Fatalf("%s: 应有三档统计, 得 %d", name, len(rep.Levels))
	}
	for _, s := range rep.Seams {
		if s.TangentAngle > 1e-12 {
			t.Errorf("%s 接缝%d: 切向夹角非零 %g", name, s.Seam, s.TangentAngle)
		}
		if s.CurvatureDiff > s.CurvatureBound {
			t.Errorf("%s 接缝%d: 曲率差 %g 超出上界 %g", name, s.Seam, s.CurvatureDiff, s.CurvatureBound)
		}
	}
	// 接缝处 C0 连续
	for i := 0; i+1 < len(out); i++ {
		if out[i].P3 != out[i+1].P0 {
			t.Errorf("%s 接缝%d: 端点未吸附", name, i)
		}
	}
	// 三档应收敛：复核轮不再产生修正
	if rep.Levels[2].Corrections != 0 {
		t.Errorf("%s: 三档未收敛, 复核轮仍修正 %d 次", name, rep.Levels[2].Corrections)
	}
	// 控制点不得漂移出发散：应保持在输入包围盒的合理扩界内
	min, max := bboxOf(segs)
	size := math.Max(max.X-min.X, max.Y-min.Y)
	for i, c := range out {
		for _, p := range []Point{c.P0, c.P1, c.P2, c.P3} {
			if math.Abs(p.X-min.X) > 10*size+1 || math.Abs(p.X-max.X) > 10*size+1 ||
				math.Abs(p.Y-min.Y) > 10*size+1 || math.Abs(p.Y-max.Y) > 10*size+1 {
				t.Errorf("%s 段%d: 控制点 %v 漂移出界", name, i, p)
			}
		}
	}
	return rep
}

func bboxOf(segs []Cubic) (min, max Point) {
	min = Point{math.Inf(1), math.Inf(1)}
	max = Point{math.Inf(-1), math.Inf(-1)}
	for _, c := range segs {
		for _, p := range []Point{c.P0, c.P1, c.P2, c.P3} {
			min.X = math.Min(min.X, p.X)
			min.Y = math.Min(min.Y, p.Y)
			max.X = math.Max(max.X, p.X)
			max.Y = math.Max(max.Y, p.Y)
		}
	}
	return
}

func TestJoinLineArc(t *testing.T) { checkJoined(t, "直线接弧", lineArcCase()) }
func TestJoinArcArc(t *testing.T)  { checkJoined(t, "弧接弧", arcArcCase()) }
func TestJoinS(t *testing.T)       { checkJoined(t, "S形", sCase()) }
func TestJoinFold(t *testing.T)    { checkJoined(t, "回折", foldCase()) }

// 多段链（三段）同样收敛。
func TestJoinChain(t *testing.T) {
	segs := append(sCase(), shift(arcCubic(0, 0, 15, deg(200), deg(320)), Point{60, 0}))
	rep := checkJoined(t, "三段链", segs)
	if len(rep.Seams) != 2 {
		t.Fatalf("应有 2 条接缝, 得 %d", len(rep.Seams))
	}
}

// 三档深度各自上报修正次数与耗时。
func TestJoinDepthLevels(t *testing.T) {
	for _, d := range []JoinDepth{JoinTangent, JoinCurvature, JoinCurvatureIter} {
		_, rep, err := Join(arcArcCase(), d)
		if err != nil {
			t.Fatal(err)
		}
		if len(rep.Levels) != int(d) {
			t.Fatalf("深度%d: 应有 %d 档统计, 得 %d", d, d, len(rep.Levels))
		}
		for _, l := range rep.Levels {
			if l.Elapsed < 0 {
				t.Fatalf("深度%d: 耗时异常", l.Depth)
			}
		}
	}
	// 一档只切向：夹角为零，但曲率差一般不为零
	_, rep, _ := Join(arcArcCase(), JoinTangent)
	if rep.Seams[0].TangentAngle > 1e-12 {
		t.Fatal("一档后切向夹角应为零")
	}
	if rep.Seams[0].CurvatureDiff == 0 {
		t.Fatal("一档不做曲率修正, 曲率差应非零")
	}
}

// 退化输入：段数不足、导数为零、控制点重合，各报段号。
func TestJoinDegenerate(t *testing.T) {
	good := Cubic{Point{0, 0}, Point{10, 0}, Point{20, 0}, Point{30, 0}}

	if _, _, err := Join([]Cubic{good}, JoinCurvatureIter); !errors.Is(err, ErrTooFewSegments) {
		t.Fatalf("段数不足应报 ErrTooFewSegments, 得 %v", err)
	}

	zeroEnd := Cubic{Point{0, 0}, Point{5, 5}, Point{10, 10}, Point{10, 10}} // P2==P3
	_, _, err := Join([]Cubic{good, zeroEnd}, JoinCurvatureIter)
	if !errors.Is(err, ErrZeroDerivative) || !strings.Contains(err.Error(), "第 1 段") {
		t.Fatalf("导数为零应报段号 1, 得 %v", err)
	}

	zeroStart := Cubic{Point{3, 3}, Point{3, 3}, Point{8, 8}, Point{12, 8}} // P0==P1
	_, _, err = Join([]Cubic{zeroStart, good, good}, JoinCurvatureIter)
	if !errors.Is(err, ErrZeroDerivative) || !strings.Contains(err.Error(), "第 0 段") {
		t.Fatalf("导数为零应报段号 0, 得 %v", err)
	}

	collapsed := Cubic{Point{5, 5}, Point{5, 5}, Point{5, 5}, Point{5, 5}}
	_, _, err = Join([]Cubic{good, good, collapsed}, JoinCurvatureIter)
	if !errors.Is(err, ErrCoincidentControls) || !strings.Contains(err.Error(), "第 2 段") {
		t.Fatalf("控制点重合应报段号 2, 得 %v", err)
	}
}

// 拼接不改动入参。
func TestJoinDoesNotMutateInput(t *testing.T) {
	segs := foldCase()
	orig := append([]Cubic(nil), segs...)
	if _, _, err := Join(segs, JoinCurvatureIter); err != nil {
		t.Fatal(err)
	}
	for i := range segs {
		if segs[i] != orig[i] {
			t.Fatalf("入参第 %d 段被修改", i)
		}
	}
}

// TestJoinTable 打印各用例的两个量与修正次数表（go test -v 查看）。
func TestJoinTable(t *testing.T) {
	cases := []struct {
		name string
		segs []Cubic
	}{
		{"直线接弧", lineArcCase()},
		{"弧接弧", arcArcCase()},
		{"S形", sCase()},
		{"回折", foldCase()},
	}
	for _, c := range cases {
		_, rep, err := Join(c.segs, JoinCurvatureIter)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		t.Logf("\n[%s]\n%s", c.name, rep.Table())
	}
}
