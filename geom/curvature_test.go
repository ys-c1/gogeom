package geom

import (
	"errors"
	"math"
	"testing"
)

// 五个代表性曲线（控制点）。
var (
	// 圆弧：单位圆四分之一弧的贝塞尔逼近，关于 t=0.5 对称。
	arcH = 4 * (math.Sqrt2 - 1) / 3
	arc  = [4]Point{{1, 0}, {1, arcH}, {arcH, 1}, {0, 1}}
	// 尖角：B'(t) = 3·((1−t)(1−2t), t(2t−1))，t=0.5 处速度为零、
	// 二阶导数非零，为真尖点，|κ|→∞。
	cusp = [4]Point{{0, 0}, {1, 0}, {0.5, -0.5}, {0.5, 0.5}}
	// S 形：中心对称，κ(t) = −κ(1−t)，t=0.5 为拐点。
	ess = [4]Point{{0, 0}, {1, 2}, {2, -2}, {3, 0}}
	// 回折：B'(t) = B'(1−t)、B''(t) = −B''(1−t)，故 κ(t) = −κ(1−t)；
	// t≈0.26、0.74 处切线竖直，急转弯。
	loop = [4]Point{{0, 0}, {2.5, 0}, {-1.5, 2}, {1, 2}}
	// 近直线：x'≡3、y''≡−6e-3 常数，κ 仅在 t=0.5 有一个极小，
	// 手算 κ(0.5) = −1.8e-2 / 9^{3/2} = −6.667e-4。
	flat = [4]Point{{0, 0}, {1, 1e-3}, {2, 1e-3}, {3, 0}}
)

func extremaOf(t *testing.T, c [4]Point, n int) []Extremum {
	t.Helper()
	ex, err := CurvatureExtrema(c[0], c[1], c[2], c[3], n)
	if err != nil {
		t.Fatalf("CurvatureExtrema 报错: %v", err)
	}
	if len(ex) == 0 {
		t.Fatal("未找到任何曲率极值")
	}
	return ex
}

// checkMonotone 核验极值点两侧曲率单调：极大值左升右降，极小值反之。
func checkMonotone(t *testing.T, c [4]Point, e Extremum) {
	t.Helper()
	if e.Cusp {
		// 尖点：|κ| 向尖点位置单调增大
		var prev float64
		for k := 3; k >= 1; k-- {
			d := float64(k) * 0.01
			kl := math.Abs(Curvature(c[0], c[1], c[2], c[3], e.T-d))
			kr := math.Abs(Curvature(c[0], c[1], c[2], c[3], e.T+d))
			if k < 3 && (kl <= prev || kr <= prev) {
				t.Fatalf("尖点两侧 |κ| 未单调增大: t=%v", e.T)
			}
			prev = math.Min(kl, kr)
		}
		return
	}
	for _, d := range []float64{0.005, 0.01, 0.02} {
		kl := Curvature(c[0], c[1], c[2], c[3], e.T-d)
		kr := Curvature(c[0], c[1], c[2], c[3], e.T+d)
		if e.Max && (kl >= e.K || kr >= e.K) {
			t.Fatalf("极大值两侧不单调: t=%v K=%v 左=%v 右=%v", e.T, e.K, kl, kr)
		}
		if !e.Max && (kl <= e.K || kr <= e.K) {
			t.Fatalf("极小值两侧不单调: t=%v K=%v 左=%v 右=%v", e.T, e.K, kl, kr)
		}
	}
}

// checkConverge 两档采样密度下极值位置必须收敛到同一处。
func checkConverge(t *testing.T, c [4]Point) {
	t.Helper()
	a := extremaOf(t, c, 64)
	b := extremaOf(t, c, 256)
	if len(a) != len(b) {
		t.Fatalf("两档密度极值个数不同: %d vs %d", len(a), len(b))
	}
	for i := range a {
		if math.Abs(a[i].T-b[i].T) > 1e-6 {
			t.Fatalf("极值位置不收敛: n=64 得 %v, n=256 得 %v", a[i].T, b[i].T)
		}
		if a[i].Max != b[i].Max || a[i].Cusp != b[i].Cusp {
			t.Fatalf("极值类型不一致: %+v vs %+v", a[i], b[i])
		}
	}
}

func TestArcExtremum(t *testing.T) {
	ex := extremaOf(t, arc, 64)
	// 手算：对称性 ⇒ 极值在 t=0.5。解析推导得
	//   x'(0.5) = −y'(0.5) = 0.75(h−2)，S = 1.125(2−h)²
	//   x''(0.5) = y''(0.5) = −3h，N = x'y''−y'x'' = 4.5h(2−h)
	//   κ(0.5) = N / S^{3/2} = 4.5h / (1.125^{3/2}·(2−h)²)
	want := 4.5 * arcH / (math.Pow(1.125, 1.5) * (2 - arcH) * (2 - arcH))
	var mid *Extremum
	for i := range ex {
		if math.Abs(ex[i].T-0.5) < 1e-9 {
			mid = &ex[i]
		}
		checkMonotone(t, arc, ex[i])
	}
	// 贝塞尔圆弧逼近的曲率在中点取局部极小、两侧各有一个极大
	if mid == nil || mid.Max || mid.Cusp {
		t.Fatalf("圆弧中点应为极小值: %+v", ex)
	}
	// 两条各自核验：手算解析值 vs 库解析式
	if math.Abs(mid.K-want) > 1e-9 {
		t.Fatalf("极值曲率与手算不符: 得 %v 手算 %v", mid.K, want)
	}
	if got := Curvature(arc[0], arc[1], arc[2], arc[3], 0.5); math.Abs(got-want) > 1e-9 {
		t.Fatalf("Curvature 解析式与手算不符: 得 %v 手算 %v", got, want)
	}
	// 逼近单位圆，曲率应接近 1
	if math.Abs(mid.K-1) > 0.01 {
		t.Fatalf("圆弧曲率应接近 1/R=1: %v", mid.K)
	}
	checkConverge(t, arc)
}

func TestCuspExtremum(t *testing.T) {
	ex := extremaOf(t, cusp, 64)
	var cp *Extremum
	for i := range ex {
		if ex[i].Cusp {
			cp = &ex[i]
		}
		checkMonotone(t, cusp, ex[i])
	}
	if cp == nil {
		t.Fatalf("应检出尖点: %+v", ex)
	}
	// 手算：B'(0.5) = 3·(0.5·0, 0.5·0) = 0 ⇒ 尖点在 t=0.5
	if math.Abs(cp.T-0.5) > 1e-9 {
		t.Fatalf("尖点位置与手算不符: %v", cp.T)
	}
	if !math.IsInf(cp.K, 1) {
		t.Fatalf("尖点曲率应为 +Inf: %v", cp.K)
	}
	// 尖点处曲率不存在
	if got := Curvature(cusp[0], cusp[1], cusp[2], cusp[3], 0.5); !math.IsNaN(got) {
		t.Fatalf("尖点处 Curvature 应为 NaN: %v", got)
	}
	checkConverge(t, cusp)
}

func TestSShapeExtrema(t *testing.T) {
	ex := extremaOf(t, ess, 64)
	var max, min *Extremum
	for i := range ex {
		checkMonotone(t, ess, ex[i])
		if ex[i].Max && (max == nil || ex[i].K > max.K) {
			max = &ex[i]
		}
		if !ex[i].Max && !ex[i].Cusp && (min == nil || ex[i].K < min.K) {
			min = &ex[i]
		}
	}
	if max == nil || min == nil {
		t.Fatalf("S 形应各有一个极大、极小: %+v", ex)
	}
	// 中心对称 ⇒ κ(t) = −κ(1−t)：极值位置互补、曲率相反
	if math.Abs(max.T+min.T-1) > 1e-6 {
		t.Fatalf("极值位置不互补: %v, %v", max.T, min.T)
	}
	if math.Abs(max.K+min.K) > 1e-6 {
		t.Fatalf("极值曲率不相反: %v, %v", max.K, min.K)
	}
	// 拐点在中点
	if k := Curvature(ess[0], ess[1], ess[2], ess[3], 0.5); math.Abs(k) > 1e-9 {
		t.Fatalf("S 形中点曲率应为零: %v", k)
	}
	checkConverge(t, ess)
}

func TestLoopBackExtrema(t *testing.T) {
	ex := extremaOf(t, loop, 64)
	if len(ex) < 2 {
		t.Fatalf("回折应有多个极值: %+v", ex)
	}
	for _, e := range ex {
		checkMonotone(t, loop, e)
	}
	// κ(t) = −κ(1−t)：极值关于 t=0.5 成对、位置互补、曲率相反
	for _, e := range ex {
		found := false
		for _, g := range ex {
			if math.Abs(e.T+g.T-1) < 1e-6 && math.Abs(e.K+g.K) < 1e-6 {
				found = true
			}
		}
		if !found {
			t.Fatalf("回折极值不关于中点反对称: %+v", e)
		}
	}
	checkConverge(t, loop)
}

func TestNearLineExtrema(t *testing.T) {
	ex := extremaOf(t, flat, 64)
	for _, e := range ex {
		checkMonotone(t, flat, e)
		if math.Abs(e.K) > 0.01 {
			t.Fatalf("近直线曲率应接近零: %+v", e)
		}
	}
	// 手算：κ(0.5) = x'y''/|x'|³ = 3·(−6e-3)/27 = −6.667e-4（极小）
	var mid *Extremum
	for i := range ex {
		if math.Abs(ex[i].T-0.5) < 1e-9 {
			mid = &ex[i]
		}
	}
	if mid == nil || mid.Max {
		t.Fatalf("近直线中点应为极小值: %+v", ex)
	}
	if want := -1.8e-2 / 27; math.Abs(mid.K-want) > 1e-9 {
		t.Fatalf("近直线极值曲率与手算不符: 得 %v 手算 %v", mid.K, want)
	}
	checkConverge(t, flat)
}

func TestExtremaRejectsBadInput(t *testing.T) {
	// 段数不足
	_, err := CurvatureExtrema(arc[0], arc[1], arc[2], arc[3], 4)
	var ce *CurveError
	if !errors.As(err, &ce) || !errors.Is(err, ErrBadSegments) || ce.Segment != 4 {
		t.Fatalf("段数不足应报 ErrBadSegments 且段号为 4: %v", err)
	}
	// 导数恒为零（四点全重合）
	p := Point{1, 1}
	_, err = CurvatureExtrema(p, p, p, p, 64)
	if !errors.As(err, &ce) || !errors.Is(err, ErrZeroDeriv) {
		t.Fatalf("四点重合应报 ErrZeroDeriv: %v", err)
	}
	// 控制点重合（相邻两点），段号为重合控制点段序号
	_, err = CurvatureExtrema(Point{0, 0}, Point{1, 1}, Point{1, 1}, Point{2, 0}, 64)
	if !errors.As(err, &ce) || !errors.Is(err, ErrCoincident) || ce.Segment != 1 {
		t.Fatalf("p1==p2 应报 ErrCoincident 且段号为 1: %v", err)
	}
}

// TestCurvatureTable 收尾输出：各曲线的极值位置与曲率表。
func TestCurvatureTable(t *testing.T) {
	cases := []struct {
		name string
		c    [4]Point
	}{
		{"圆弧", arc}, {"尖角", cusp}, {"S形", ess}, {"回折", loop}, {"近直线", flat},
	}
	for _, tc := range cases {
		ex := extremaOf(t, tc.c, 256)
		t.Logf("%s: %d 个极值", tc.name, len(ex))
		for _, e := range ex {
			kind := "极小"
			if e.Max {
				kind = "极大"
			}
			if e.Cusp {
				kind = "尖点"
			}
			t.Logf("  t=%.6f  κ=%+.6f  %s", e.T, e.K, kind)
		}
	}
}
