package geom

import (
	"errors"
	"math"
	"strings"
	"testing"
)

// 生成一段圆弧（不含重复的终点）。
func arcPoints(cx, cy, r float64, n int, a0, a1 float64) []Point {
	out := make([]Point, n+1)
	for i := 0; i <= n; i++ {
		a := a0 + (a1-a0)*float64(i)/float64(n)
		out[i] = Point{cx + r*math.Cos(a), cy + r*math.Sin(a)}
	}
	return out
}

func sawtooth(n int) []Point {
	out := make([]Point, n)
	for i := range out {
		x := float64(i)
		y := 0.0
		if i%2 == 1 {
			y = 1
		}
		out[i] = Point{x, y}
	}
	return out
}

func bowtie() []Point {
	// 自交：两条对角线交叉，最后回到起点附近形成闭合环。
	return []Point{{0, 0}, {10, 10}, {10, 0}, {0, 10}, {0, 0}}
}

func assertPointError(t *testing.T, err error, want error, wantIdx int) {
	t.Helper()
	if err == nil {
		t.Fatalf("期望错误 %v，得到 nil", want)
	}
	var pe *PointError
	if !errors.As(err, &pe) {
		t.Fatalf("错误应带点号，得到 %T: %v", err, err)
	}
	if !errors.Is(err, want) {
		t.Fatalf("期望 %v，得到 %v", want, err)
	}
	if pe.Index != wantIdx {
		t.Fatalf("点号期望 %d，得到 %d（%v）", wantIdx, pe.Index, err)
	}
}

func TestSimplifyErrorsWithPointIndex(t *testing.T) {
	tri := []Point{{0, 0}, {1, 0}, {0, 1}}
	if _, err := Simplify(tri, -1); err == nil {
		t.Fatal("负容差应报错")
	} else {
		assertPointError(t, err, ErrNegativeTolerance, 0)
	}
	if _, err := Simplify(nil, 1); err == nil {
		t.Fatal("无点应报错")
	} else {
		assertPointError(t, err, ErrTooFewPts, 0)
	}
	two := []Point{{0, 0}, {1, 1}}
	if _, err := Simplify(two, 1); err == nil {
		t.Fatal("两点应报错")
	} else {
		assertPointError(t, err, ErrTooFewPts, 2)
	}
	line := []Point{{0, 0}, {1, 1}, {2, 2}, {3, 3}}
	if _, err := Simplify(line, 0.5); err == nil {
		t.Fatal("全共线应报错")
	} else {
		assertPointError(t, err, ErrCollinear, 0)
	}
	same := make([]Point, 5)
	if _, err := Simplify(same, 0.5); err == nil {
		t.Fatal("全部重合应报错")
	} else {
		assertPointError(t, err, ErrCollinear, 0)
	}
}

func TestDeletedRangesContiguousAndQueryable(t *testing.T) {
	got := DeletedRanges(10, []int{0, 3, 7, 9})
	want := []IndexRange{{1, 2}, {4, 6}, {8, 8}}
	if len(got) != len(want) {
		t.Fatalf("区间数期望 %d，得到 %d: %v", len(want), len(got), got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("区间 %d 期望 %v，得到 %v", i, want[i], got[i])
		}
		if got[i].Len() != want[i].End-want[i].Start+1 {
			t.Fatalf("区间 %d 长度错误: %v", i, got[i])
		}
	}
	if r := DeletedRanges(4, nil); len(r) != 1 || r[0] != (IndexRange{0, 3}) {
		t.Fatalf("无保留点时应删整段: %v", r)
	}
	if r := DeletedRanges(0, nil); len(r) != 0 {
		t.Fatalf("空路径不应有区间: %v", r)
	}
	// 全部保留、全部删除的边界情况。
	if r := DeletedRanges(3, []int{0, 1, 2}); len(r) != 0 {
		t.Fatalf("全保留不应有被删区间: %v", r)
	}
	if r := DeletedRanges(3, []int{0, 2}); len(r) != 1 || r[0] != (IndexRange{1, 1}) {
		t.Fatalf("中间删除区间错误: %v", r)
	}
}

func TestPointErrorString(t *testing.T) {
	pe := &PointError{Err: ErrTooFewPts, Index: 7}
	if got := pe.Error(); !strings.Contains(got, "7") || !strings.Contains(got, ErrTooFewPts.Error()) {
		t.Fatalf("错误串应同时含原因与点号: %q", got)
	}
	if !errors.Is(pe, ErrTooFewPts) {
		t.Fatal("应能 errors.Is 到原始错误")
	}
}

func TestMaxDeviationDegenerate(t *testing.T) {
	pts := []Point{{0, 0}, {1, 1}, {2, 0}}
	if d, _ := MaxDeviation(pts, nil); !math.IsInf(d, 1) {
		t.Fatalf("空简化路径应返回 +Inf，得到 %v", d)
	}
	// 单点简化路径按点距计算。
	if d, idx := MaxDeviation(pts, []Point{{0, 0}}); d != 2 || idx != 2 {
		t.Fatalf("单点对照错误 d=%v idx=%d", d, idx)
	}
}

func TestSimplifyDegenerateClosed(t *testing.T) {
	// 首尾相同但只有两点不同 → 全部落在同一直线，应报共线。
	tri := []Point{{0, 0}, {1, 1}, {0, 0}}
	_, err := Simplify(tri, 0.1)
	assertPointError(t, err, ErrCollinear, 0)

	// 合法的最小闭环 {A,B,C,A}，首尾差必须为零。
	ring := []Point{{0, 0}, {10, 0}, {5, 8}, {0, 0}}
	res, err := Simplify(ring, 0.1)
	if err != nil {
		t.Fatal(err)
	}
	if !res.Closed || res.ClosureGap != 0 {
		t.Fatalf("最小闭环处理错误 closed=%v gap=%v", res.Closed, res.ClosureGap)
	}
	if res.RetentionRatio <= 0 || res.RetentionRatio > 1 {
		t.Fatalf("保留率超出范围: %v", res.RetentionRatio)
	}
}

func TestDeletedRangesMatchResult(t *testing.T) {
	pts := arcPoints(0, 0, 100, 200, 0, math.Pi/2)
	res, err := Simplify(pts, 1)
	if err != nil {
		t.Fatal(err)
	}
	// Deleted 必须与 Kept 完全互补，且每个被删点号恰好落在一个区间里。
	inDeleted := 0
	for i := 0; i < len(pts); i++ {
		if res.DeletedContains(i) {
			inDeleted++
		}
	}
	if inDeleted+res.RetainedPoints != len(pts) {
		t.Fatalf("被删 %d + 保留 %d != 总点数 %d", inDeleted, res.RetainedPoints, len(pts))
	}
	for _, k := range res.Kept {
		if res.DeletedContains(k) {
			t.Fatalf("保留点 %d 不应出现在被删区间里", k)
		}
	}
	// 区间按点号升序、互不重叠，且相邻区间之间恰好夹着保留点。
	prevEnd := -1
	for _, rg := range res.Deleted {
		if rg.End < rg.Start || rg.Start <= prevEnd {
			t.Fatalf("区间不连续或重叠: %v (前一个结束于 %d)", res.Deleted, prevEnd)
		}
		prevEnd = rg.End
	}
}

func TestSimplifyStraightLineIsCollinear(t *testing.T) {
	pts := make([]Point, 50)
	for i := range pts {
		pts[i] = Point{float64(i), 2 * float64(i)}
	}
	_, err := Simplify(pts, 0.1)
	assertPointError(t, err, ErrCollinear, 0)
}

func TestSimplifyArc(t *testing.T) {
	pts := arcPoints(0, 0, 100, 400, 0, math.Pi/2)
	orig := append([]Point(nil), pts...)
	res, err := Simplify(pts, 2)
	if err != nil {
		t.Fatal(err)
	}
	// 原始路径不能被改动。
	for i := range orig {
		if pts[i] != orig[i] {
			t.Fatalf("原始路径第 %d 点被修改", i)
		}
	}
	if res.MaxDeviation > res.Tolerance {
		t.Fatalf("最大偏差 %v 超过容差 %v", res.MaxDeviation, res.Tolerance)
	}
	if res.RetainedPoints >= len(pts) {
		t.Fatalf("圆弧应被抽稀，保留 %d/%d", res.RetainedPoints, len(pts))
	}
	if res.RetentionRatio != float64(res.RetainedPoints)/float64(len(pts)) {
		t.Fatalf("保留率不一致: %v vs %d/%d", res.RetentionRatio, res.RetainedPoints, len(pts))
	}
	if res.Kept[0] != 0 || res.Kept[len(res.Kept)-1] != len(pts)-1 {
		t.Fatalf("首尾点必须保留: %v", res.Kept)
	}
	// 容差越大，保留点越少（单调不增）。
	prev := res.RetainedPoints
	for _, tol := range []float64{5, 10, 20, 40} {
		r2, err := Simplify(pts, tol)
		if err != nil {
			t.Fatal(err)
		}
		if r2.RetainedPoints > prev {
			t.Fatalf("容差 %v 保留点 %d 多于更小容差时的 %d", tol, r2.RetainedPoints, prev)
		}
		if r2.MaxDeviation > r2.Tolerance {
			t.Fatalf("容差 %v 时偏差 %v 超限", tol, r2.MaxDeviation)
		}
		prev = r2.RetainedPoints
	}
}

func TestSimplifySawtooth(t *testing.T) {
	pts := sawtooth(101)
	res, err := Simplify(pts, 0.4)
	if err != nil {
		t.Fatal(err)
	}
	// 锯齿幅值为 1，容差 0.4 必须保留所有尖峰。
	for i := 1; i < len(pts)-1; i += 2 {
		if res.DeletedContains(i) {
			t.Fatalf("尖峰 %d 被误删", i)
		}
	}
	if res.MaxDeviation > 0.4 {
		t.Fatalf("偏差 %v 超限", res.MaxDeviation)
	}
	// 零容差时所有偏离直线的点都应保留。
	zero, err := Simplify(pts, 0)
	if err != nil {
		t.Fatal(err)
	}
	if zero.RetainedPoints != len(pts) {
		t.Fatalf("零容差应保留全部点，得到 %d/%d", zero.RetainedPoints, len(pts))
	}
}

func TestSimplifyClosedLoop(t *testing.T) {
	ring := arcPoints(0, 0, 50, 360, 0, 2*math.Pi)
	ring[len(ring)-1] = ring[0] // 消除浮点舍入，确保首尾逐位相等
	if ring[0] != ring[len(ring)-1] {
		t.Fatal("构造的闭环首尾不相等")
	}
	orig := append([]Point(nil), ring...)
	res, err := Simplify(ring, 1.5)
	if err != nil {
		t.Fatal(err)
	}
	if !res.Closed {
		t.Fatal("应识别为闭合环")
	}
	if res.ClosureGap != 0 {
		t.Fatalf("闭环首尾差值必须为零，得到 %v", res.ClosureGap)
	}
	if res.Points[0] != res.Points[len(res.Points)-1] {
		t.Fatalf("闭环简化后首尾点不相等: %v vs %v", res.Points[0], res.Points[len(res.Points)-1])
	}
	if res.MaxDeviation > res.Tolerance {
		t.Fatalf("闭环偏差 %v 超过容差 %v", res.MaxDeviation, res.Tolerance)
	}
	if res.RetentionRatio <= 0 || res.RetentionRatio >= 1 {
		t.Fatalf("保留率应在 (0,1) 内，得到 %v", res.RetentionRatio)
	}
	for i := range orig {
		if ring[i] != orig[i] {
			t.Fatalf("原始闭环第 %d 点被修改", i)
		}
	}
}

func TestSimplifySelfIntersecting(t *testing.T) {
	pts := bowtie()
	res, err := Simplify(pts, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	if res.MaxDeviation > res.Tolerance {
		t.Fatalf("自交路径偏差 %v 超限", res.MaxDeviation)
	}
	if res.Points[0] != res.Points[len(res.Points)-1] {
		t.Fatal("自交闭环首尾应相等")
	}
	// 大容差把拐角也压掉时，仍需通过逐点核对。
	res2, err := Simplify(pts, 100)
	if err != nil {
		t.Fatal(err)
	}
	if res2.MaxDeviation > res2.Tolerance {
		t.Fatalf("大容差偏差 %v 超限", res2.MaxDeviation)
	}
}

func TestVerifyCatchesOvershoot(t *testing.T) {
	pts := arcPoints(0, 0, 100, 64, 0, math.Pi/2)
	// 只用两端点的弦，弧顶偏差约为 100*(1-cos(pi/4)) ≈ 29.3。
	if _, err := Verify(pts, []Point{pts[0], pts[len(pts)-1]}, 1); err == nil {
		t.Fatal("应检出超限")
	} else {
		var pe *PointError
		if !errors.As(err, &pe) || !errors.Is(err, ErrToleranceExceeded) {
			t.Fatalf("应返回带点号的超限错误，得到 %v", err)
		}
		if pe.Index <= 0 || pe.Index >= len(pts)-1 {
			t.Fatalf("超限点号应位于弧中，得到 %d", pe.Index)
		}
	}
	if d, err := Verify(pts, []Point{pts[0], pts[len(pts)-1]}, 30); err != nil {
		t.Fatalf("上界 30 应通过，偏差 %v: %v", d, err)
	}
	if _, err := Verify(pts, pts, -1); err == nil {
		t.Fatal("负上界应报错")
	} else {
		assertPointError(t, err, ErrNegativeTolerance, 0)
	}
}

func TestSimplifyTwoDensitiesReport(t *testing.T) {
	// 两档点数：分别报告保留点数与耗时。
	for _, n := range []int{1000, 10000} {
		pts := arcPoints(0, 0, 500, n, 0, 2*math.Pi*0.75)
		res, err := Simplify(pts, 0.05)
		if err != nil {
			t.Fatal(err)
		}
		if res.MaxDeviation > res.Tolerance {
			t.Fatalf("n=%d 偏差 %v 超限", n, res.MaxDeviation)
		}
		t.Logf("n=%d 保留=%d (%.2f%%) 抽稀耗时=%v 核对耗时=%v 最大偏差=%.6f",
			n, res.RetainedPoints, 100*res.RetentionRatio, res.Elapsed, res.VerifyElapsed, res.MaxDeviation)
	}
}

func BenchmarkSimplify10k(b *testing.B) {
	pts := arcPoints(0, 0, 500, 10000, 0, 2*math.Pi*0.75)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := Simplify(pts, 0.05); err != nil {
			b.Fatal(err)
		}
	}
}
