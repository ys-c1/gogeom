package geom

import (
	"errors"
	"math"
	"testing"
)

func sumRingArea(t *testing.T, rings [][]Point) float64 {
	t.Helper()
	var s float64
	for _, r := range rings {
		a, err := Area(r)
		if err != nil {
			t.Fatalf("环面积异常 %v", err)
		}
		s += math.Abs(a)
	}
	return s
}

func perimeter(p []Point) float64 {
	var s float64
	for i := range p {
		a, b := p[i], p[(i+1)%len(p)]
		s += math.Hypot(b.X-a.X, b.Y-a.Y)
	}
	return s
}

// rasterRef 逐像素数格子的参照：返回交集近似面积、连通分量数和格宽。
// 采样点在格心，测试用例的边界都避开格心，不影响计数。
func rasterRef(a, b []Point, n int) (area float64, comps int, cell float64) {
	minX, minY := math.Inf(1), math.Inf(1)
	maxX, maxY := math.Inf(-1), math.Inf(-1)
	for _, ps := range [][]Point{a, b} {
		for _, p := range ps {
			minX = math.Min(minX, p.X)
			minY = math.Min(minY, p.Y)
			maxX = math.Max(maxX, p.X)
			maxY = math.Max(maxY, p.Y)
		}
	}
	cell = math.Max(maxX-minX, maxY-minY) / float64(n)
	inside := make([][]bool, n)
	var cnt int
	for iy := 0; iy < n; iy++ {
		inside[iy] = make([]bool, n)
		for ix := 0; ix < n; ix++ {
			x := minX + (float64(ix)+0.5)*cell
			y := minY + (float64(iy)+0.5)*cell
			if Inside(a, x, y) && Inside(b, x, y) {
				inside[iy][ix] = true
				cnt++
			}
		}
	}
	// 4-连通 flood fill 数分量
	seen := make([][]bool, n)
	for i := range seen {
		seen[i] = make([]bool, n)
	}
	for iy := 0; iy < n; iy++ {
		for ix := 0; ix < n; ix++ {
			if !inside[iy][ix] || seen[iy][ix] {
				continue
			}
			comps++
			stack := [][2]int{{iy, ix}}
			seen[iy][ix] = true
			for len(stack) > 0 {
				c := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				for _, d := range [][2]int{{1, 0}, {-1, 0}, {0, 1}, {0, -1}} {
					ny, nx := c[0]+d[0], c[1]+d[1]
					if ny >= 0 && ny < n && nx >= 0 && nx < n && inside[ny][nx] && !seen[ny][nx] {
						seen[ny][nx] = true
						stack = append(stack, [2]int{ny, nx})
					}
				}
			}
		}
	}
	return float64(cnt) * cell * cell, comps, cell
}

// checkRaster 两条独立核验：环数等于参照连通分量数，面积与数格子一致。
func checkRaster(t *testing.T, a, b []Point, rings [][]Point) {
	t.Helper()
	refArea, comps, cell := rasterRef(a, b, 400)
	if len(rings) != comps {
		t.Errorf("环数 %d != 参照连通分量数 %d", len(rings), comps)
	}
	got := sumRingArea(t, rings)
	tol := 1.5 * (perimeter(a) + perimeter(b)) * cell
	if math.Abs(got-refArea) > tol {
		t.Errorf("交区面积 %v 与数格子参照 %v 偏差超过 %v", got, refArea, tol)
	}
	t.Logf("环数=%d 面积=%.6f 参照面积=%.6f 参照分量=%d", len(rings), got, refArea, comps)
}

func TestIntersectCases(t *testing.T) {
	sq := func(x0, y0, x1, y1 float64) []Point {
		return []Point{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}
	}
	cases := []struct {
		name     string
		a, b     []Point
		wantRing int
		wantArea float64
	}{
		{"相离", sq(0, 0, 1, 1), sq(3, 3, 4, 4), 0, 0},
		{"相交", sq(0, 0, 2, 2), sq(1, 1, 3, 3), 1, 1},
		{"包含", sq(0, 0, 4, 4), sq(1, 1, 2, 2), 1, 1},
		{"被包含", sq(1, 1, 2, 2), sq(0, 0, 4, 4), 1, 1},
		{"共线贴边", sq(0, 0, 1, 1), sq(1, 0, 2, 1), 0, 0},
		{"完全重合", sq(0, 0, 2, 2), sq(0, 0, 2, 2), 1, 4},
		{"两分量", // 横条跨过 U 形两臂，交集为两个矩形
			[]Point{{-1, 2.5}, {7, 2.5}, {7, 3.5}, {-1, 3.5}},
			[]Point{{0, 0}, {6, 0}, {6, 4}, {4, 4}, {4, 2}, {2, 2}, {2, 4}, {0, 4}},
			2, 4},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			rings, err := Intersect(c.a, c.b)
			if err != nil {
				t.Fatalf("求交失败 %v", err)
			}
			if len(rings) != c.wantRing {
				t.Fatalf("环数 %d，期望 %d", len(rings), c.wantRing)
			}
			if got := sumRingArea(t, rings); math.Abs(got-c.wantArea) > 1e-9 {
				t.Fatalf("面积 %v，期望 %v", got, c.wantArea)
			}
			checkRaster(t, c.a, c.b, rings)
		})
	}
}

// 绕向自推：同一组点顺时针/逆时针输入结果必须一致。
func TestIntersectWinding(t *testing.T) {
	ccwSq := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	cwSq := []Point{{0, 0}, {0, 2}, {2, 2}, {2, 0}}
	other := []Point{{1, 1}, {3, 1}, {3, 3}, {1, 3}}
	r1, err1 := Intersect(ccwSq, other)
	r2, err2 := Intersect(cwSq, other)
	if err1 != nil || err2 != nil {
		t.Fatalf("求交失败 %v %v", err1, err2)
	}
	if len(r1) != 1 || len(r2) != 1 {
		t.Fatalf("环数异常 %d %d", len(r1), len(r2))
	}
	a1, _ := Area(r1[0])
	a2, _ := Area(r2[0])
	if a1 != 1 || a2 != 1 {
		t.Fatalf("面积异常 %v %v", a1, a2)
	}
}

func TestIntersectReject(t *testing.T) {
	ok := []Point{{0, 0}, {4, 0}, {0, 4}}
	cases := []struct {
		name string
		poly []Point
		want error
		idx  int
	}{
		{"顶点不足", []Point{{0, 0}, {1, 1}}, ErrTooFewPts, 2},
		{"顶点重合", []Point{{0, 0}, {2, 0}, {2, 0}, {0, 2}}, ErrDupVertex, 2},
		{"自交未拆分", []Point{{0, 0}, {2, 2}, {2, 0}, {0, 2}}, ErrSelfIntersect, 2},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := Intersect(c.poly, ok)
			if err == nil {
				t.Fatal("应拒绝未拒绝")
			}
			if !errors.Is(err, c.want) {
				t.Fatalf("错误类型 %v，期望 %v", err, c.want)
			}
			var ve *VertexError
			if !errors.As(err, &ve) || ve.Index != c.idx {
				t.Fatalf("序号异常 %+v，期望 %d", err, c.idx)
			}
			t.Logf("拒绝：%v", err)
		})
	}
}
