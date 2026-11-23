package geom

import (
	"math"
	"math/rand"
	"sort"
	"testing"
)

// expectDecompose 对一个合法用例执行完整的独立验证。
func expectDecompose(t *testing.T, name string, poly []Point, wantReflex int) {
	t.Helper()

	got, err := ReflexVertices(poly)
	if err != nil {
		t.Fatalf("%s: 反射顶点检测失败: %v", name, err)
	}
	if len(got) != wantReflex {
		t.Fatalf("%s: 反射顶点数 = %d (%v), 期望 %d", name, len(got), got, wantReflex)
	}

	pieces, err := ConvexDecompose(poly)
	if err != nil {
		t.Fatalf("%s: 凸分解失败: %v", name, err)
	}

	// 每块必须是凸的。
	for i, p := range pieces {
		if len(p) < 3 {
			t.Fatalf("%s: 第 %d 块顶点不足", name, i)
		}
		if !IsConvex(p) {
			t.Fatalf("%s: 第 %d 块非凸: %v", name, i, p)
		}
	}

	// 守恒条件一（独立）：面积和等于原面积。
	origArea, sumArea, err := CheckAreaConservation(poly, pieces)
	if err != nil {
		t.Fatalf("%s: %v (原面积 %v, 块面积和 %v)", name, err, origArea, sumArea)
	}

	// 守恒条件二（独立）：块间不重叠。
	if err := CheckNoOverlap(pieces); err != nil {
		t.Fatalf("%s: %v", name, err)
	}

	// 越界复核：所有顶点必须在原多边形内部或边界上。
	if err := CheckPiecesWithin(poly, pieces); err != nil {
		t.Fatalf("%s: %v", name, err)
	}

	// 报告：块数与各块面积。
	t.Logf("%s: 反射顶点=%v 块数=%d 原面积=%.4f 块面积和=%.4f",
		name, got, len(pieces), origArea, sumArea)
	for i, p := range pieces {
		a, _ := Area(p)
		t.Logf("  块 %d: %d 顶点 面积=%.4f", i, len(p), math.Abs(a))
	}
}

// 单凹口：方盘上开一个矩形缺口，2 个反射顶点。
func TestDecomposeSingleNotch(t *testing.T) {
	poly := []Point{
		{0, 0}, {10, 0}, {10, 10},
		{7, 10}, {7, 4}, {3, 4}, {3, 10}, {0, 10},
	}
	expectDecompose(t, "单凹口", poly, 2)
}

// 多凹口：梳子形，4 个反射顶点。
func TestDecomposeMultiNotch(t *testing.T) {
	poly := []Point{
		{0, 0}, {12, 0}, {12, 8},
		{9, 8}, {9, 3}, {7, 3}, {7, 8},
		{5, 8}, {5, 3}, {3, 3}, {3, 8}, {0, 8},
	}
	expectDecompose(t, "多凹口", poly, 4)
}

// 螺旋：由 1 格宽的螺旋走廊单元并集取边界得到，保证简单多边形。
func TestDecomposeSpiral(t *testing.T) {
	cells := []struct{ x, y int }{
		{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0},
		{4, 1}, {4, 2}, {3, 2}, {2, 2}, {1, 2},
		{1, 3}, {1, 4}, {2, 4}, {3, 4},
	}
	poly := cellsToPolygon(cells)
	if err := ValidatePolygon(poly); err != nil {
		t.Fatalf("螺旋夹具自身不合法: %v", err)
	}
	expectDecompose(t, "螺旋", poly, len(reflexIndices(poly)))
}

// 细颈：两个大块通过 0.2 宽的颈相连。
func TestDecomposeThinNeck(t *testing.T) {
	poly := []Point{
		{0, 0}, {4, 0}, {4, 1.9}, {8, 1.9},
		{8, 0}, {12, 0}, {12, 4},
		{8, 4}, {8, 2.1}, {4, 2.1}, {4, 4}, {0, 4},
	}
	expectDecompose(t, "细颈", poly, 4)
}

// 近凸：浅 L 形（凹进仅 1 个单位），边上带共线顶点，仅 1 个反射顶点。
func TestDecomposeNearConvex(t *testing.T) {
	poly := []Point{
		{0, 0}, {5, 0}, {10, 0}, {10, 9}, {9, 9},
		{9, 10}, {5, 10}, {0, 10}, {0, 5},
	}
	expectDecompose(t, "近凸", poly, 1)
}

// 纯凸多边形：不应切分，原样返回一块。
func TestDecomposeConvex(t *testing.T) {
	poly := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	pieces, err := ConvexDecompose(poly)
	if err != nil {
		t.Fatalf("凸多边形分解失败: %v", err)
	}
	if len(pieces) != 1 {
		t.Fatalf("凸多边形块数 = %d, 期望 1", len(pieces))
	}
	orig, sum, err := CheckAreaConservation(poly, pieces)
	if err != nil || orig != 100 || sum != 100 {
		t.Fatalf("凸多边形面积校验失败: %v %v %v", orig, sum, err)
	}
}

// 拒绝一：顶点重复，报出两个重复顶点号。
func TestRejectDuplicate(t *testing.T) {
	poly := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}, {10, 0}}
	_, err := ConvexDecompose(poly)
	var ie *IndexError
	if !asIndexError(err, &ie) || ie.Op != "duplicate" {
		t.Fatalf("应报顶点重复，实际: %v", err)
	}
	if !sameInts(ie.Index, []int{1, 4}) {
		t.Fatalf("重复顶点号错误: %v", ie.Index)
	}
	t.Logf("顶点重复正确拒绝: %v", err)
}

// 拒绝二：自交（蝴蝶结），报出相交的两条边的顶点号。
func TestRejectSelfIntersect(t *testing.T) {
	bowtie := []Point{{0, 0}, {10, 10}, {10, 0}, {0, 10}}
	_, err := ConvexDecompose(bowtie)
	var ie *IndexError
	if !asIndexError(err, &ie) || ie.Op != "self-intersect" {
		t.Fatalf("应报自交，实际: %v", err)
	}
	t.Logf("自交正确拒绝: %v", err)
}

// 拒绝三：对角线相交，报出两条对角线的四个顶点号。
func TestRejectCrossingDiagonals(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	err := CheckDiagonals(sq, [][2]int{{0, 2}, {1, 3}})
	var ie *IndexError
	if !asIndexError(err, &ie) || ie.Op != "diagonal-cross" {
		t.Fatalf("应报对角线相交，实际: %v", err)
	}
	if len(ie.Index) != 4 {
		t.Fatalf("应报四个顶点号，实际: %v", ie.Index)
	}
	t.Logf("对角线相交正确拒绝: %v", err)
}

// 越界即失败：从顶点引向多边形外部的对角线必须被拒绝。
func TestRejectDiagonalOutside(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	if err := CheckDiagonals(sq, [][2]int{{0, 99}}); err == nil {
		t.Fatal("越界顶点号应报错")
	}
	// 近凸用例上，所有切分对角线都必须落在内部（由分解过程内部保证并复核）。
	poly := []Point{
		{0, 0}, {5, 0}, {10, 0}, {10, 9}, {9, 9},
		{9, 10}, {5, 10}, {0, 10}, {0, 5},
	}
	pieces, err := ConvexDecompose(poly)
	if err != nil {
		t.Fatalf("分解失败: %v", err)
	}
	if err := CheckPiecesWithin(poly, pieces); err != nil {
		t.Fatalf("存在越界对角线: %v", err)
	}
}

// 分解结果用到的对角线必须两两不相交。
func TestDecomposeDiagonalsDoNotCross(t *testing.T) {
	cases := map[string][]Point{
		"单凹口": {{0, 0}, {10, 0}, {10, 10}, {7, 10}, {7, 4}, {3, 4}, {3, 10}, {0, 10}},
		"细颈":  {{0, 0}, {4, 0}, {4, 1.9}, {8, 1.9}, {8, 0}, {12, 0}, {12, 4}, {8, 4}, {8, 2.1}, {4, 2.1}, {4, 4}, {0, 4}},
	}
	for name, poly := range cases {
		pieces, err := ConvexDecompose(poly)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		diags := diagonalsOf(poly, pieces)
		if err := CheckDiagonals(poly, diags); err != nil {
			t.Fatalf("%s: 分解对角线不合法: %v", name, err)
		}
		if err := CheckNoOverlap(pieces); err != nil {
			t.Fatalf("%s: %v", name, err)
		}
	}
}

// 顺时针多边形必须与逆时针得到同样正确的结果。
func TestDecomposeClockwise(t *testing.T) {
	fwd := []Point{{0, 0}, {10, 0}, {10, 10}, {7, 10}, {7, 4}, {3, 4}, {3, 10}, {0, 10}}
	rev := reversePoints(fwd)
	p1, err := ConvexDecompose(fwd)
	if err != nil {
		t.Fatal(err)
	}
	p2, err := ConvexDecompose(rev)
	if err != nil {
		t.Fatalf("顺时针分解失败: %v", err)
	}
	if len(p1) != len(p2) {
		t.Fatalf("顺/逆时针块数不一致: %d vs %d", len(p1), len(p2))
	}
	if _, _, err := CheckAreaConservation(rev, p2); err != nil {
		t.Fatal(err)
	}
	if err := CheckNoOverlap(p2); err != nil {
		t.Fatal(err)
	}
}

// T 型接触：末顶点 (10,5) 压在非相邻边 (10,0)->(10,10) 上，必须按自交拒绝并报顶点号。
func TestRejectVertexOnEdge(t *testing.T) {
	poly := []Point{
		{0, 0}, {10, 0}, {10, 10},
		{0, 10}, {0, 5}, {10, 5},
	}
	err := ValidatePolygon(poly)
	var ie *IndexError
	if !asIndexError(err, &ie) || ie.Op != "self-intersect" {
		t.Fatalf("应按自交拒绝 T 型接触，实际: %v", err)
	}
	t.Logf("T 型接触正确拒绝: %v", err)
}

// 反射顶点数为零的输入必须一块不切。
func TestReflexFreeInput(t *testing.T) {
	tri := []Point{{0, 0}, {6, 0}, {3, 8}}
	r, err := ReflexVertices(tri)
	if err != nil || len(r) != 0 {
		t.Fatalf("三角形反射顶点应为空: %v %v", r, err)
	}
	pieces, err := ConvexDecompose(tri)
	if err != nil || len(pieces) != 1 {
		t.Fatalf("三角形应原样返回一块: %d %v", len(pieces), err)
	}
}

// 模糊测试：随机简单多边形（极坐标扰动的星形），每次分解都要通过全部守恒检查。
func FuzzDecompose(f *testing.F) {
	f.Add(int64(1), 8)
	f.Fuzz(func(t *testing.T, seed int64, n int) {
		if n < 3 || n > 40 {
			return
		}
		rng := rand.New(rand.NewSource(seed))
		poly := randomSimplePolygon(rng, n)
		if err := ValidatePolygon(poly); err != nil {
			t.Skip("生成器产出非法多边形")
		}
		reflex := reflexIndices(poly)
		pieces, err := ConvexDecompose(poly)
		if err != nil {
			t.Fatalf("分解失败: %v", err)
		}
		// 仅允许顶点间对角线（无 Steiner 点）时，凸分块数的经典上界为 2r+1。
		if len(pieces) > 2*len(reflex)+1 {
			t.Fatalf("块数 %d 超过无 Steiner 上界 2r+1 = %d", len(pieces), 2*len(reflex)+1)
		}
		// 另一条结构性不变量：凸分块数不可能超过三角化的三角形数 n-2。
		if len(pieces) > len(poly)-2 {
			t.Fatalf("块数 %d 超过三角化上界 n-2 = %d", len(pieces), len(poly)-2)
		}
		for i, p := range pieces {
			if !IsConvex(p) {
				t.Fatalf("第 %d 块非凸", i)
			}
		}
		if _, _, err := CheckAreaConservation(poly, pieces); err != nil {
			t.Fatal(err)
		}
		if err := CheckNoOverlap(pieces); err != nil {
			t.Fatal(err)
		}
		if err := CheckPiecesWithin(poly, pieces); err != nil {
			t.Fatal(err)
		}
	})
}

func reversePoints(in []Point) []Point {
	out := make([]Point, len(in))
	for i := range in {
		out[i] = in[len(in)-1-i]
	}
	return out
}

// randomSimplePolygon 生成以原点为中心、按极角排序的扰动星形简单多边形。
func randomSimplePolygon(rng *rand.Rand, n int) []Point {
	pts := make([]Point, n)
	step := 2 * math.Pi / float64(n)
	for i := 0; i < n; i++ {
		theta := float64(i)*step + rng.Float64()*step*0.5
		r := 5 + rng.Float64()*5
		pts[i] = Point{r * math.Cos(theta), r * math.Sin(theta)}
	}
	// 按极角排序保证不自交（星形多边形族）。
	sort.Slice(pts, func(a, b int) bool {
		return math.Atan2(pts[a].Y, pts[a].X) < math.Atan2(pts[b].Y, pts[b].X)
	})
	return pts
}

// ---- 辅助 ----

func asIndexError(err error, target **IndexError) bool {
	for err != nil {
		if ie, ok := err.(*IndexError); ok {
			*target = ie
			return true
		}
		type wrapper interface{ Unwrap() error }
		w, ok := err.(wrapper)
		if !ok {
			return false
		}
		err = w.Unwrap()
	}
	return false
}

func sameInts(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	x := append([]int(nil), a...)
	y := append([]int(nil), b...)
	sort.Ints(x)
	sort.Ints(y)
	for i := range x {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}

// diagonalsOf 从分块的共享边反推出分解所用对角线。
func diagonalsOf(poly []Point, pieces [][]Point) [][2]int {
	type edge struct{ a, b int }
	idx := map[Point]int{}
	for i, p := range poly {
		idx[p] = i
	}
	count := map[edge]int{}
	add := func(i, j int) {
		if i > j {
			i, j = j, i
		}
		count[edge{i, j}]++
	}
	for _, piece := range pieces {
		for k := range piece {
			a, b := piece[k], piece[(k+1)%len(piece)]
			add(idx[a], idx[b])
		}
	}
	var out [][2]int
	for e, c := range count {
		if c >= 2 {
			out = append(out, [2]int{e.a, e.b}) // 共享边即对角线
		}
	}
	return out
}

// cellsToPolygon 把轴对齐单位方格并集转成 CCW 边界多边形（夹具生成器）。
func cellsToPolygon(cells []struct{ x, y int }) []Point {
	type cell struct{ x, y int }
	set := map[cell]bool{}
	for _, c := range cells {
		set[cell{c.x, c.y}] = true
	}
	type edge struct{ ax, ay, bx, by int }
	type ukey struct{ ax, ay, bx, by int } // 无向规范键
	canonical := func(e edge) ukey {
		if e.ax > e.bx || (e.ax == e.bx && e.ay > e.by) {
			return ukey{e.bx, e.by, e.ax, e.ay}
		}
		return ukey{e.ax, e.ay, e.bx, e.by}
	}
	count := map[ukey]int{}
	var directed []edge
	toggle := func(ax, ay, bx, by int) {
		e := edge{ax, ay, bx, by}
		count[canonical(e)]++
		directed = append(directed, e)
	}
	for c := range set {
		x, y := c.x, c.y
		// 逆时针方格：底、右、顶、左
		toggle(x, y, x+1, y)
		toggle(x+1, y, x+1, y+1)
		toggle(x+1, y+1, x, y+1)
		toggle(x, y+1, x, y)
	}
	// 仅保留出现一次的边（内部边被两个方格各引用一次，共两次）。
	var boundary []edge
	for _, e := range directed {
		if count[canonical(e)] == 1 {
			boundary = append(boundary, e)
		}
	}
	// 把有向边按首尾相接串成环。
	start := boundary[0]
	used := make([]bool, len(boundary))
	used[0] = true
	cur := start
	poly := []Point{{float64(cur.ax), float64(cur.ay)}}
	for {
		poly = append(poly, Point{float64(cur.bx), float64(cur.by)})
		if cur.bx == start.ax && cur.by == start.ay {
			break // 回到起点，环闭合
		}
		nextIdx, found := -1, false
		for k, e := range boundary {
			if used[k] {
				continue
			}
			if e.ax == cur.bx && e.ay == cur.by {
				nextIdx, found = k, true
				break
			}
		}
		if !found {
			break
		}
		used[nextIdx] = true
		cur = boundary[nextIdx]
	}
	// 闭合时末点与首点重复，去掉重复的闭合点。
	if len(poly) > 1 && poly[len(poly)-1] == poly[0] {
		poly = poly[:len(poly)-1]
	}
	return poly
}
