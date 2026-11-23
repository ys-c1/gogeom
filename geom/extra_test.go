package geom

import (
	"errors"
	"strings"
	"testing"
)

func TestRasterSetGetAndBounds(t *testing.T) {
	r := &Raster{MinX: 0, MinY: 0, Pixel: 1, W: 4, H: 3, Cover: make([]bool, 12)}
	// 越界 Get 返回 false，Set 不 panic
	if r.Get(-1, 0) || r.Get(4, 0) || r.Get(0, 3) {
		t.Fatal("越界 Get 应为 false")
	}
	r.Set(-5, 0)
	r.Set(0, 0)
	r.Set(3.9, 2.9)
	if !r.Get(0.5, 0.5) || !r.Get(3.5, 2.5) {
		t.Fatal("Set 后 Get 应为 true")
	}
	if r.Count() != 2 {
		t.Fatalf("Count=%d want 2", r.Count())
	}
	// 连通分量：两个对角块 + 一个孤立点
	r2 := &Raster{MinX: 0, MinY: 0, Pixel: 1, W: 5, H: 5, Cover: make([]bool, 25)}
	for _, idx := range []int{0, 1, 5, 6} { // 2x2 块
		r2.Cover[idx] = true
	}
	r2.Cover[24] = true // 孤立点
	if c := r2.Components(); c != 2 {
		t.Fatalf("Components=%d want 2", c)
	}
}

func TestInsidePolyWrapper(t *testing.T) {
	sq := []Point{{X: 0, Y: 0}, {X: 10, Y: 0}, {X: 10, Y: 10}, {X: 0, Y: 10}}
	if !InsidePoly(sq, 5, 5) || InsidePoly(sq, 15, 5) {
		t.Fatal("InsidePoly 包装行为错误")
	}
}

func TestMoreErrorBranches(t *testing.T) {
	// 顶点不足
	if _, err := StrokePolygons([]Point{{X: 0, Y: 0}}, StrokeParams{Width: 10}); !errors.Is(err, ErrTooFewPts) {
		t.Fatalf("顶点不足错误不对: %v", err)
	}
	// 非正线宽
	path := []Point{{X: 0, Y: 0}, {X: 10, Y: 0}, {X: 10, Y: 10}}
	if _, err := StrokePolygons(path, StrokeParams{Width: 0}); !errors.Is(err, ErrDegenerate) {
		t.Fatalf("非正线宽应报退化: %v", err)
	}
	if _, err := StrokePolygons(path, StrokeParams{Width: -2}); !errors.Is(err, ErrDegenerate) {
		t.Fatalf("负线宽应报退化: %v", err)
	}
	// 非法样式（负值）
	if _, err := StrokePolygons(path, StrokeParams{Width: 10, Style: -1}); !errors.Is(err, ErrBadJoinStyle) {
		t.Fatalf("负样式值应报错: %v", err)
	}
	// 零长边
	zeroEdge := []Point{{X: 0, Y: 0}, {X: 0, Y: 0}, {X: 10, Y: 10}}
	if _, err := StrokePolygons(zeroEdge, StrokeParams{Width: 10, Style: JoinRound, Radius: []float64{5}}); !errors.Is(err, ErrDegenerate) {
		t.Fatalf("零长边应报退化: %v", err)
	}
	// JoinGapPixels 的错误分支
	v := Point{X: 10, Y: 0}
	t0 := Point{X: 1, Y: 0}
	if _, _, _, err := JoinGapPixels(v, t0, t0, 5, 5, JoinRound, 64, 0.2); !errors.Is(err, ErrCollinearJoin) {
		t.Fatalf("共线 gap 应报错: %v", err)
	}
	if _, _, _, err := JoinGapPixels(v, t0, Point{X: 0, Y: 1}, 5, 0, JoinRound, 64, 0.2); !errors.Is(err, ErrZeroRadius) {
		t.Fatalf("零半径 gap 应报错: %v", err)
	}
	if _, _, _, err := JoinGapPixels(v, t0, Point{X: 0, Y: 1}, 5, 5, JoinStyle(11), 64, 0.2); !errors.Is(err, ErrBadJoinStyle) {
		t.Fatalf("非法样式 gap 应报错: %v", err)
	}
	// Stroke 端到端错误透传
	if _, err := Stroke(path, StrokeParams{Width: 10, Style: JoinRound, Radius: []float64{0}}, 0.2); !errors.Is(err, ErrZeroRadius) {
		t.Fatalf("Stroke 未透传零半径错误: %v", err)
	}
	// 半径超 h 被钳制：不报错且面积等于 r=h
	_, _, _, err := JoinGapPixels(v, Point{X: 1, Y: 0}, Point{X: 0, Y: 1}, 5, 99, JoinRound, 64, 0.2)
	if err != nil {
		t.Fatalf("半径钳制不应报错: %v", err)
	}
	// RadiusTurns 短序列
	if turns := RadiusTurns([]float64{1, 2}); turns != nil {
		t.Fatalf("短序列应返回 nil, got %v", turns)
	}
	// 错误信息包含拐点号的字符串检查
	_, err = StrokePolygons(path, StrokeParams{Width: 10, Style: JoinBevel, Radius: []float64{0}})
	if err == nil || !strings.Contains(err.Error(), "拐点 1") {
		t.Fatalf("错误信息缺拐点号: %v", err)
	}
}
