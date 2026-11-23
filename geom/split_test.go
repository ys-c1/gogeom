package geom

import (
	"errors"
	"math"
	"testing"
)

const tol = 1e-9

func checkSplit(t *testing.T, p Polygon, n int, wantTotal float64) *SplitResult {
	t.Helper()
	res, err := SplitEqualArea(p, n)
	if err != nil {
		t.Fatalf("切分失败: %v", err)
	}
	if math.Abs(res.Total-wantTotal) > tol {
		t.Fatalf("总面积错误: got %v want %v", res.Total, wantTotal)
	}
	if err := res.Verify(1e-6); err != nil {
		t.Fatalf("核验失败: %v", err)
	}
	target := wantTotal / float64(n)
	for i, pc := range res.Pieces {
		if math.Abs(pc.Area-target) > 1e-6 {
			t.Errorf("第 %d 块面积 %v，目标 %v", i, pc.Area, target)
		}
	}
	return res
}

func TestSplitConvex(t *testing.T) {
	sq := Polygon{Outer: []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}}
	res := checkSplit(t, sq, 4, 100)
	want := []float64{2.5, 5, 7.5}
	for i, c := range res.Cuts {
		if math.Abs(c-want[i]) > 1e-9 {
			t.Errorf("切线 %d 位置 %v，期望 %v", i, c, want[i])
		}
	}
}

func TestSplitConcave(t *testing.T) {
	// L 形：4x4 去掉右上 2x2，面积 12
	l := Polygon{Outer: []Point{{0, 0}, {4, 0}, {4, 2}, {2, 2}, {2, 4}, {0, 4}}}
	res := checkSplit(t, l, 4, 12)
	want := []float64{0.75, 1.5, 2.5}
	for i, c := range res.Cuts {
		if math.Abs(c-want[i]) > 1e-9 {
			t.Errorf("切线 %d 位置 %v，期望 %v", i, c, want[i])
		}
	}
}

func TestSplitWithHole(t *testing.T) {
	box := Polygon{
		Outer: []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}},
		Holes: [][]Point{{{4, 4}, {4, 6}, {6, 6}, {6, 4}}}, // 顺时针，面积 4
	}
	res := checkSplit(t, box, 3, 96)
	for i, pc := range res.Pieces {
		if math.Abs(pc.Area-32) > 1e-6 {
			t.Errorf("第 %d 块面积 %v，目标 32", i, pc.Area)
		}
	}
}

func TestSplitThin(t *testing.T) {
	thin := Polygon{Outer: []Point{{0, 0}, {1000, 0}, {1000, 1}, {0, 1}}}
	res := checkSplit(t, thin, 5, 1000)
	for i, c := range res.Cuts {
		want := float64(i+1) * 200
		if math.Abs(c-want) > 1e-9 {
			t.Errorf("切线 %d 位置 %v，期望 %v", i, c, want)
		}
	}
}

func TestSplitErrors(t *testing.T) {
	sq := Polygon{Outer: []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}}
	if _, err := SplitEqualArea(sq, 0); !errors.Is(err, ErrZeroParts) {
		t.Errorf("块数为零应报 ErrZeroParts，got %v", err)
	}
	// 退化：三点共线，面积为零
	degen := Polygon{Outer: []Point{{0, 0}, {1, 0}, {2, 0}}}
	if _, err := SplitEqualArea(degen, 2); !errors.Is(err, ErrZeroArea) {
		t.Errorf("面积为零应报 ErrZeroArea，got %v", err)
	}
	// 自交：蝴蝶结，边 0 与边 2 相交
	bowtie := Polygon{Outer: []Point{{0, 0}, {2, 2}, {2, 0}, {0, 2}}}
	_, err := SplitEqualArea(bowtie, 2)
	var si *SelfIntersectionError
	if !errors.As(err, &si) {
		t.Fatalf("自交应报 SelfIntersectionError，got %v", err)
	}
	if si.RingA != 0 || si.EdgeA != 0 || si.EdgeB != 2 {
		t.Errorf("自交边号错误: %+v", si)
	}
}

func TestCutsNonIntersecting(t *testing.T) {
	// 切线为互相平行的竖直线且位置严格递增，故两两不相交
	arrow := Polygon{Outer: []Point{{0, 0}, {6, 0}, {6, 2}, {8, 2}, {4, 5}, {0, 2}, {2, 2}, {2, 0}}}
	res := checkSplit(t, arrow, 6, 20)
	for i := 1; i < len(res.Cuts); i++ {
		if !(res.Cuts[i] > res.Cuts[i-1]) {
			t.Fatalf("切线位置未严格递增: %v", res.Cuts)
		}
	}
}

func TestSplitSumConservation(t *testing.T) {
	// 带孔 + 凹 + 多块的守恒与极差核验
	p := Polygon{
		Outer: []Point{{0, 0}, {12, 0}, {12, 8}, {8, 8}, {8, 4}, {4, 4}, {4, 8}, {0, 8}},
		Holes: [][]Point{{{1, 1}, {1, 3}, {3, 3}, {3, 1}}},
	}
	res, err := SplitEqualArea(p, 7)
	if err != nil {
		t.Fatalf("切分失败: %v", err)
	}
	if err := res.Verify(1e-6); err != nil {
		t.Fatalf("核验失败: %v", err)
	}
	var sum float64
	for _, pc := range res.Pieces {
		sum += pc.Area
	}
	if d := math.Abs(sum - res.Total); d > 1e-9 {
		t.Fatalf("面积不守恒: 和 %v 原面积 %v 差 %v", sum, res.Total, d)
	}
}
