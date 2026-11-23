package geom

import (
	"errors"
	"math"
	"testing"
)

// 容差：面积相对误差阈值
const relTol = 1e-9

func ringN(n int, fn func(i int) Point) []Point {
	out := make([]Point, n)
	for i := 0; i < n; i++ {
		out[i] = fn(i)
	}
	return out
}

// assertReport 逐块核对：三角形数 = 轮数；面积守恒；无重叠。
func assertReport(t *testing.T, name string, rep *TriReport, wantTris int, wantArea float64) {
	t.Helper()
	if len(rep.Triangles) != wantTris {
		t.Errorf("%s: 三角形数 = %d, 期望 %d", name, len(rep.Triangles), wantTris)
	}
	if rep.Rounds != len(rep.Triangles) {
		t.Errorf("%s: 耳切轮数 %d != 三角形数 %d（每轮应产出一个三角形）",
			name, rep.Rounds, len(rep.Triangles))
	}
	if rep.OverlapPairs != 0 {
		t.Errorf("%s: 存在 %d 对三角形内部重叠", name, rep.OverlapPairs)
	}
	// 逐块重算面积，确保报告值与独立计算一致
	var sum float64
	for _, tr := range rep.Triangles {
		sum += absAreaTri(tr[0], tr[1], tr[2])
	}
	if math.Abs(sum-rep.TotalTriArea) > relTol*math.Max(1, rep.TotalTriArea) {
		t.Errorf("%s: 逐块面积和 %.9f 与报告 %.9f 不符", name, sum, rep.TotalTriArea)
	}
	if math.Abs(rep.TotalTriArea-wantArea) > relTol*math.Max(1, wantArea) {
		t.Errorf("%s: 三角化面积 %.9f 与原多边形面积 %.9f 不一致，差 %g",
			name, rep.TotalTriArea, wantArea, rep.AreaDiff)
	}
	if rep.AreaRelDiff > relTol {
		t.Errorf("%s: 面积相对差 %g 超过阈值 %g", name, rep.AreaRelDiff, relTol)
	}
}

func TestTriangulateConvex(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	rep, err := TriangulateReport(sq)
	if err != nil {
		t.Fatalf("方块三角化失败: %v", err)
	}
	assertReport(t, "square", rep, 2, 100)
	if rep.EstRounds != 2 || rep.Rounds != 2 {
		t.Errorf("方块轮数估算/实际 = %d/%d", rep.EstRounds, rep.Rounds)
	}

	// 凸六边形
	hex := ringN(6, func(i int) Point {
		a := float64(i) * math.Pi / 3
		return Point{10 * math.Cos(a), 10 * math.Sin(a)}
	})
	rep2, err := TriangulateReport(hex)
	if err != nil {
		t.Fatalf("六边形三角化失败: %v", err)
	}
	want := ringAbsArea(hex)
	assertReport(t, "hexagon", rep2, 4, want)
}

func TestTriangulateConcave(t *testing.T) {
	// L 形
	L := []Point{{0, 0}, {10, 0}, {10, 4}, {4, 4}, {4, 10}, {0, 10}}
	rep, err := TriangulateReport(L)
	if err != nil {
		t.Fatalf("L 形三角化失败: %v", err)
	}
	assertReport(t, "L-shape", rep, 4, ringAbsArea(L))
	if rep.SplitFaces != 1 {
		t.Errorf("简单凹多边形 SplitFaces 应为 1，得到 %d", rep.SplitFaces)
	}

	// 星形（简单凹多边形，16 顶点）
	star := ringN(16, func(i int) Point {
		a := float64(i) * math.Pi / 8
		r := 10.0
		if i%2 == 1 {
			r = 4
		}
		return Point{r * math.Cos(a), r * math.Sin(a)}
	})
	rep2, err := TriangulateReport(star)
	if err != nil {
		t.Fatalf("星形三角化失败: %v", err)
	}
	assertReport(t, "star", rep2, 14, ringAbsArea(star))
}

func TestTriangulateWithHoles(t *testing.T) {
	outer := []Point{{0, 0}, {20, 0}, {20, 20}, {0, 20}}
	hole := []Point{{6, 6}, {6, 14}, {14, 14}, {14, 6}}
	rep, err := TriangulateRingsReport(outer, [][]Point{hole})
	if err != nil {
		t.Fatalf("单洞三角化失败: %v", err)
	}
	want := ringAbsArea(outer) - ringAbsArea(hole)
	assertReport(t, "one-hole", rep, len(rep.Triangles), want)
	if len(rep.Bridges) != 1 {
		t.Errorf("单洞应有 1 条桥，得到 %d", len(rep.Bridges))
	}

	// 双洞
	h1 := []Point{{2, 6}, {2, 14}, {7, 14}, {7, 6}}
	h2 := []Point{{13, 6}, {13, 14}, {18, 14}, {18, 6}}
	rep2, err := TriangulateRingsReport(outer, [][]Point{h1, h2})
	if err != nil {
		t.Fatalf("双洞三角化失败: %v", err)
	}
	want2 := ringAbsArea(outer) - ringAbsArea(h1) - ringAbsArea(h2)
	assertReport(t, "two-hole", rep2, len(rep2.Triangles), want2)
	if len(rep2.Bridges) != 2 {
		t.Errorf("双洞应有 2 条桥，得到 %d", len(rep2.Bridges))
	}
}

func TestTriangulateDegenerate(t *testing.T) {
	// 顶点不足三
	if _, err := Triangulate([]Point{{0, 0}, {1, 1}}); err == nil {
		t.Fatal("顶点不足三应报错")
	} else {
		var ve *VertexError
		if !errors.As(err, &ve) || ve.Code != "TooFew" {
			t.Fatalf("应为 TooFew VertexError，得到 %v", err)
		}
		if len(ve.Vertices) != 2 {
			t.Errorf("TooFew 应报出全部顶点号，得到 %v", ve.Vertices)
		}
	}

	// 全共线
	line := []Point{{0, 0}, {5, 0}, {10, 0}, {15, 0}}
	_, err := Triangulate(line)
	var ve *VertexError
	if !errors.As(err, &ve) || ve.Code != "AllCollinear" {
		t.Fatalf("全共线应报 AllCollinear，得到 %v", err)
	}
	if len(ve.Vertices) != 4 {
		t.Errorf("全共线应报出全部顶点号，得到 %v", ve.Vertices)
	}

	// 近似共线：一个顶点仅微小凸起（凸起量小于坐标尺度对应的共线容差），
	// 应被清理而非报错
	near := []Point{{0, 0}, {10, 0}, {10, 10}, {5, 10 + 1e-11}, {0, 10}}
	rep, err := TriangulateReport(near)
	if err != nil {
		t.Fatalf("近似共线不应报错: %v", err)
	}
	if len(rep.CollinearRemoved) == 0 {
		t.Error("近似共线顶点应被清理并记录顶点号")
	}
	assertReport(t, "near-collinear", rep, len(rep.Triangles), ringAbsArea(near))
}

func TestTriangulateSelfIntersect(t *testing.T) {
	// 蝴蝶结：两条对角线规范相交，拆成两个三角形面
	bow := []Point{{0, 0}, {10, 10}, {10, 0}, {0, 10}}
	rep, err := TriangulateReport(bow)
	if err != nil {
		t.Fatalf("蝴蝶结应能拆分: %v", err)
	}
	if rep.SplitFaces != 2 {
		t.Errorf("蝴蝶结应拆成 2 个面，得到 %d", rep.SplitFaces)
	}
	// 奇偶规则：两个三角形，每个面积 25，共 50
	assertReport(t, "bowtie", rep, 2, 50)

	// 五角星（5 顶点自交），奇偶填充，中心五边形留空
	pent := ringN(5, func(i int) Point {
		a := float64(i)*4*math.Pi/5 - math.Pi/2
		return Point{10 * math.Cos(a), 10 * math.Sin(a)}
	})
	rep2, err := TriangulateReport(pent)
	if err != nil {
		t.Fatalf("五角星应能拆分: %v", err)
	}
	if rep2.SplitFaces != 5 {
		t.Errorf("五角星应拆成 5 个面，得到 %d", rep2.SplitFaces)
	}
	if rep2.OverlapPairs != 0 {
		t.Errorf("五角星三角网格存在重叠: %d 对", rep2.OverlapPairs)
	}
	// 奇偶规则下五角星面积 = 5 个尖角三角形之和
	assertReport(t, "pentagram", rep2, 5, rep2.OriginalArea)
}

func TestTriangulateSelfIntersectUnsplit(t *testing.T) {
	// 共线区间重叠边：边 (0,0)->(10,0) 与非相邻边 (7,0)->(3,0) 在 [3,7] 上
	// 实质重叠，平面半边无法把它拆成合法面，必须带顶点号报 SelfIntersect
	bad := []Point{{0, 0}, {10, 0}, {10, 10}, {7, 0}, {3, 0}, {0, 10}}
	_, err := Triangulate(bad)
	var ve *VertexError
	if !errors.As(err, &ve) || ve.Code != "SelfIntersect" {
		t.Fatalf("共线重叠自交应报 SelfIntersect，得到 %v", err)
	}
	if len(ve.Vertices) == 0 {
		t.Error("自交未拆分应报出相关顶点号")
	}
	// 报出的顶点号必须都是合法下标
	for _, idx := range ve.Vertices {
		if idx < 0 || idx >= len(bad) {
			t.Errorf("报出的顶点号 %d 越界（输入共 %d 个顶点）", idx, len(bad))
		}
	}
}

func TestRoundsEstimation(t *testing.T) {
	// 简单 n 边形：轮数恒为 n-2
	for _, n := range []int{3, 4, 5, 8} {
		poly := ringN(n, func(i int) Point {
			a := float64(i) * 2 * math.Pi / float64(n)
			return Point{math.Cos(a), math.Sin(a)}
		})
		rep, err := TriangulateReport(poly)
		if err != nil {
			t.Fatalf("n=%d 三角化失败: %v", n, err)
		}
		if rep.EstRounds != n-2 || rep.Rounds != n-2 {
			t.Errorf("n=%d 轮数估算/实际 = %d/%d，期望均为 %d",
				n, rep.EstRounds, rep.Rounds, n-2)
		}
	}
}
