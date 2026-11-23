package geom

import (
	"errors"
	"fmt"
	"math"
	"strings"
	"testing"
)

// turnPath 构造两段折线：v 为拐点，t0=(1,0)，出边与入边夹角 deg（左转）。
func turnPath(deg float64, length float64) []Point {
	rad := deg * math.Pi / 180
	v := Point{X: 20, Y: 0}
	t1 := Point{X: math.Cos(rad), Y: math.Sin(rad)}
	return []Point{
		{X: v.X - length, Y: 0},
		v,
		{X: v.X + length*t1.X, Y: v.Y + length*t1.Y},
	}
}

// isConvex 判定简单多边形是否凸（允许共线边）。
// 容差按边长按尺度自适应，避免近 180° 巨型多边形上的浮点噪声误判。
func isConvex(poly []Point) bool {
	n := len(poly)
	scale := 0.0
	for i := 0; i < n; i++ {
		a := poly[i]
		b := poly[(i+1)%n]
		scale = math.Max(scale, math.Hypot(b.X-a.X, b.Y-a.Y))
	}
	eps := scale * scale * 1e-9
	sign := 0
	for i := 0; i < n; i++ {
		a := poly[i]
		b := poly[(i+1)%n]
		c := poly[(i+2)%n]
		cross := (b.X-a.X)*(c.Y-b.Y) - (b.Y-a.Y)*(c.X-b.X)
		s := 0
		if cross > eps {
			s = 1
		} else if cross < -eps {
			s = -1
		}
		if s != 0 && sign != 0 && s != sign {
			return false
		}
		if s != 0 {
			sign = s
		}
	}
	return true
}

func polyLen(poly []Point) float64 {
	var l float64
	for i := range poly {
		a, b := poly[i], poly[(i+1)%len(poly)]
		l += math.Hypot(b.X-a.X, b.Y-a.Y)
	}
	return l
}

// TestAngleCategories 覆盖锐角/直角/钝角/近180度：三种接头均生成凸多边形，
// 且面积与手算解析公式一致。
func TestAngleCategories(t *testing.T) {
	h := 5.0
	v := Point{X: 20, Y: 0}
	t0 := Point{X: 1, Y: 0}
	cases := []struct {
		name string
		deg  float64
	}{
		{"锐角30", 30}, {"锐角45", 45}, {"锐角60", 60},
		{"直角90", 90},
		{"钝角120", 120}, {"钝角150", 150},
		{"近180度179", 179},
	}
	for _, c := range cases {
		rad := c.deg * math.Pi / 180
		t1 := Point{X: math.Cos(rad), Y: math.Sin(rad)}
		mp := miterJoin(v, t0, t1, h)
		rp := roundJoin(v, t0, t1, h, h, 256)
		bp := bevelJoin(v, t0, t1, h, h)
		for name, p := range map[string][]Point{"miter": mp, "round": rp, "bevel": bp} {
			if !isConvex(p) {
				t.Fatalf("%s/%s 多边形非凸: %v", c.name, name, p)
			}
		}
		wantMiter := h * h * math.Tan(rad/2)
		wantSector := 0.5 * h * h * rad
		if got := absArea(mp); math.Abs(got-wantMiter) > 1e-8*wantMiter {
			t.Fatalf("%s miter 面积 %v want %v", c.name, got, wantMiter)
		}
		if got := absArea(bp); math.Abs(got-wantMiter) > 1e-8*wantMiter {
			t.Fatalf("%s bevel(r=h) 面积 %v want %v", c.name, got, wantMiter)
		}
		// 圆弧离散误差上界：n 段弦切割的总面积 ≤ r²α³/(6n²)
		n := 256
		discBound := h * h * rad * rad * rad / (6 * float64(n) * float64(n))
		if got := absArea(rp); wantSector-got > discBound || got > wantSector+1e-9 {
			t.Fatalf("%s round 面积 %v 超出离散上界 [%v, %v]", c.name, got, wantSector-discBound, wantSector)
		}
		// 右转镜像面积一致
		t1r := Point{X: math.Cos(rad), Y: -math.Sin(rad)}
		if got := absArea(miterJoin(v, t0, t1r, h)); math.Abs(got-wantMiter) > 1e-8*wantMiter {
			t.Fatalf("%s 右转 miter 不对称", c.name)
		}
	}
}

// TestCollinearAndNear180 共线（180度同向与折返）报错且带拐点号；
// 近180度（179）正常处理。
func TestCollinearAndNear180(t *testing.T) {
	h := 5.0
	// 路径含两个内部拐点：i=1 正常 90 度，i=2 共线
	path := []Point{{0, 0}, {10, 0}, {10, 10}, {20, 10}, {30, 10}}
	for _, style := range []JoinStyle{JoinRound, JoinBevel, JoinMiter} {
		params := StrokeParams{Width: 2 * h, Style: style, Radius: []float64{h, h, h}, ArcSegs: 64}
		_, err := StrokePolygons(path, params)
		if err == nil {
			t.Fatalf("样式 %d 共线拐点应报错", style)
		}
		if !errors.Is(err, ErrCollinearJoin) || !strings.Contains(err.Error(), "拐点 3") {
			t.Fatalf("共线错误信息缺拐点号: %v", err)
		}
	}
	// 折返（180度反向）同样报共线
	fold := []Point{{0, 0}, {10, 0}, {0, 0}}
	_, err := StrokePolygons(fold, StrokeParams{Width: 2 * h, Style: JoinRound, Radius: []float64{h}})
	if !errors.Is(err, ErrCollinearJoin) {
		t.Fatalf("折返应报共线，got %v", err)
	}
	// 179 度近180：不报错，且直折尖角足够长
	near := turnPath(179, 30)
	polys, err := StrokePolygons(near, StrokeParams{Width: 2 * h, Style: JoinMiter})
	if err != nil {
		t.Fatalf("近180度不应报错: %v", err)
	}
	if len(polys) != 3 { // 2 条边带 + 1 个接头
		t.Fatalf("近180度多边形数 %d", len(polys))
	}
}

// TestJoinErrors 半径为零、样式非法均返回带拐点号的错误。
func TestJoinErrors(t *testing.T) {
	path := turnPath(90, 20)
	// 半径为零（圆角/斜切），拐点号为 1
	for _, style := range []JoinStyle{JoinRound, JoinBevel} {
		_, err := StrokePolygons(path, StrokeParams{Width: 10, Style: style, Radius: []float64{0}})
		if !errors.Is(err, ErrZeroRadius) || !strings.Contains(err.Error(), "拐点 1") {
			t.Fatalf("零半径错误不对: %v", err)
		}
		// 负半径同效
		_, err = StrokePolygons(path, StrokeParams{Width: 10, Style: style, Radius: []float64{-2}})
		if !errors.Is(err, ErrZeroRadius) {
			t.Fatalf("负半径应报零半径错误, got %v", err)
		}
	}
	// 非法样式
	_, err := StrokePolygons(path, StrokeParams{Width: 10, Style: JoinStyle(99)})
	if !errors.Is(err, ErrBadJoinStyle) {
		t.Fatalf("非法样式应报错, got %v", err)
	}
	// 半径数与拐点数不一致
	_, err = StrokePolygons(path, StrokeParams{Width: 10, Style: JoinRound, Radius: []float64{1, 2}})
	if err == nil {
		t.Fatal("半径数不匹配应报错")
	}
}

// TestDefaultMiterUnchanged 直折为默认分支：零值参数即直折，
// 且输出与显式 JoinMiter 完全一致（确定性快照）。
func TestDefaultMiterUnchanged(t *testing.T) {
	path := []Point{{0, 0}, {20, 0}, {20, 20}, {40, 20}}
	got, err := StrokePolygons(path, StrokeParams{Width: 10})
	if err != nil {
		t.Fatal(err)
	}
	want, err := StrokePolygons(path, StrokeParams{Width: 10, Style: JoinMiter})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != len(want) {
		t.Fatalf("默认分支多边形数 %d != %d", len(got), len(want))
	}
	var sb strings.Builder
	for k := range got {
		if len(got[k]) != len(want[k]) {
			t.Fatalf("多边形 %d 顶点数不一致", k)
		}
		for i := range got[k] {
			if got[k][i] != want[k][i] {
				t.Fatalf("默认分支输出变化: poly %d pt %d %v != %v", k, i, got[k][i], want[k][i])
			}
			fmt.Fprintf(&sb, "%d:%d:(%.4f,%.4f)\n", k, i, got[k][i].X, got[k][i].Y)
		}
	}
	// 金快照：路径 (0,0)->(20,0)->(20,20)->(40,20)、线宽 10 的直折输出。
	// 3 条边带（CCW）+ 2 个直折接头（90° 外正方形），共 5 个多边形。
	snapshot := strings.Join([]string{
		"0:0:(0.0000,5.0000)", "0:1:(20.0000,5.0000)", "0:2:(20.0000,-5.0000)", "0:3:(0.0000,-5.0000)",
		"1:0:(15.0000,0.0000)", "1:1:(15.0000,20.0000)", "1:2:(25.0000,20.0000)", "1:3:(25.0000,0.0000)",
		"2:0:(20.0000,25.0000)", "2:1:(40.0000,25.0000)", "2:2:(40.0000,15.0000)", "2:3:(20.0000,15.0000)",
		// 接头1（v=(20,0)，外正方形 q0-tip-q1-v）
		"3:0:(20.0000,5.0000)", "3:1:(15.0000,5.0000)", "3:2:(15.0000,0.0000)", "3:3:(20.0000,0.0000)",
		// 接头2（v=(20,20)，外正方形 q0-tip-q1-v）
		"4:0:(25.0000,20.0000)", "4:1:(25.0000,15.0000)", "4:2:(20.0000,15.0000)", "4:3:(20.0000,20.0000)",
	}, "\n") + "\n"
	if sb.String() != snapshot {
		t.Fatalf("直折快照漂移:\n%s", sb.String())
	}
	// 栅格化结果同样一致
	r1, _ := Stroke(path, StrokeParams{Width: 10}, 0.5)
	r2, _ := Stroke(path, StrokeParams{Width: 10, Style: JoinMiter}, 0.5)
	if r1.Count() != r2.Count() {
		t.Fatalf("默认/显式直折像素数 %d != %d", r1.Count(), r2.Count())
	}
}

// TestRasterAreaBound 接头覆盖像素与手算面积对照，误差给严格上界：
// 中心采样栅格化的对称差被边界管道包含，|Np²-A| ≤ L·2p + πp²。
func TestRasterAreaBound(t *testing.T) {
	h := 5.0
	v := Point{X: 30, Y: 30}
	t0 := Point{X: 1, Y: 0}
	pixel := 0.1
	cases := []struct {
		deg  float64
		make func(rad float64) []Point
		area func(rad float64) float64
	}{
		{90, func(a float64) []Point { return miterJoin(v, t0, Point{math.Cos(a), math.Sin(a)}, h) },
			func(a float64) float64 { return h * h * math.Tan(a/2) }},
		{60, func(a float64) []Point { return bevelJoin(v, t0, Point{math.Cos(a), math.Sin(a)}, h, h/2) },
			func(a float64) float64 { return h * (h / 2) * math.Tan(a/2) }},
	}
	for _, c := range cases {
		rad := c.deg * math.Pi / 180
		poly := c.make(rad)
		analytic := c.area(rad)
		if got := absArea(poly); math.Abs(got-analytic) > 1e-7*analytic {
			t.Fatalf("deg %v 解析面积不一致 %v vs %v", c.deg, got, analytic)
		}
		r := RasterizePolys([][]Point{poly}, pixel)
		areaPixels := float64(r.Count()) * pixel * pixel
		bound := polyLen(poly)*2*pixel + math.Pi*pixel*pixel
		if math.Abs(areaPixels-analytic) > bound {
			t.Fatalf("deg %v 像素面积 %v 与解析 %v 差 %v 超上界 %v",
				c.deg, areaPixels, analytic, math.Abs(areaPixels-analytic), bound)
		}
	}
	// 圆角扇形：像素面积对照离散多边形面积（边界管道上界），
	// 多边形面积再对照理想扇形（离散上界）。
	rad := 90.0 * math.Pi / 180
	t1 := Point{X: math.Cos(rad), Y: math.Sin(rad)}
	for _, n := range []int{32, 128, 512} {
		poly := roundJoin(v, t0, t1, h, h, n)
		r := RasterizePolys([][]Point{poly}, pixel)
		areaPixels := float64(r.Count()) * pixel * pixel
		polyArea := absArea(poly)
		bound := polyLen(poly)*2*pixel + math.Pi*pixel*pixel
		if math.Abs(areaPixels-polyArea) > bound {
			t.Fatalf("n=%d 像素面积 %v 与多边形 %v 差超上界 %v", n, areaPixels, polyArea, bound)
		}
		ideal := 0.5 * h * h * rad
		discBound := h * h * rad * rad * rad / (6 * float64(n) * float64(n))
		if ideal-polyArea > discBound {
			t.Fatalf("n=%d 离散误差超上界: %v > %v", n, ideal-polyArea, discBound)
		}
	}
}

// TestGapPixelsMonotone 缺口像素随半径单调不增：半径越大，接头越接近直折信封。
func TestGapPixelsMonotone(t *testing.T) {
	h := 5.0
	v := Point{X: 30, Y: 0}
	t0 := Point{X: 1, Y: 0}
	rad := 40.0 * math.Pi / 180
	t1 := Point{X: math.Cos(rad), Y: math.Sin(rad)}
	prev := math.MaxInt
	for k := 1; k <= 10; k++ {
		r := h * float64(k) / 10
		gap, _, _, err := JoinGapPixels(v, t0, t1, h, r, JoinRound, 128, 0.1)
		if err != nil {
			t.Fatal(err)
		}
		if gap > prev {
			t.Fatalf("半径 %v 缺口 %d 大于上一档 %d", r, gap, prev)
		}
		prev = gap
	}
	// 直折样式缺口为零
	gap, _, _, err := JoinGapPixels(v, t0, t1, h, h, JoinMiter, 128, 0.1)
	if err != nil || gap != 0 {
		t.Fatalf("直折缺口应为 0, got %d err %v", gap, err)
	}
}

// TestCurvedRadiusNoGap 半径沿路径连成曲线（先增后减，含拐点）时，
// 整条描边栅格化后保持单连通，无像素缝隙；曲线拐点被标出。
func TestCurvedRadiusNoGap(t *testing.T) {
	// 6 个内部拐点的折线
	path := []Point{{0, 0}, {10, 0}, {20, 5}, {30, 0}, {40, 5}, {50, 0}, {60, 5}, {70, 0}}
	h := 5.0
	// 正弦曲线半径 profile：r = h·(0.55 + 0.45·sin)，二阶差分在中间变号
	radii := make([]float64, len(path)-2)
	for i := range radii {
		radii[i] = h * (0.55 + 0.45*math.Sin(float64(i)*math.Pi/3))
	}
	turns := RadiusTurns(radii)
	if len(turns) == 0 {
		t.Fatal("曲线半径 profile 应至少有一个拐点")
	}
	params := StrokeParams{Width: 2 * h, Style: JoinRound, Radius: radii, ArcSegs: 64}
	r, err := Stroke(path, params, 0.25)
	if err != nil {
		t.Fatal(err)
	}
	if c := r.Components(); c != 1 {
		t.Fatalf("曲线半径描边出现缝隙，连通分量=%d（应为 1）", c)
	}
	// 斜切同样无缝
	params.Style = JoinBevel
	r2, err := Stroke(path, params, 0.25)
	if err != nil {
		t.Fatal(err)
	}
	if c := r2.Components(); c != 1 {
		t.Fatalf("斜切曲线半径描边连通分量=%d", c)
	}
}

// TestRadiusTurns 拐点检测：单调序列无拐点，S 形序列拐点位置正确。
func TestRadiusTurns(t *testing.T) {
	if turns := RadiusTurns([]float64{1, 2, 3, 4, 5}); len(turns) != 0 {
		t.Fatalf("单调序列不应有拐点, got %v", turns)
	}
	// 凸序列（二阶差分恒正）无拐点
	if turns := RadiusTurns([]float64{1, 4, 9, 16, 25}); len(turns) != 0 {
		t.Fatalf("凸序列不应有拐点, got %v", turns)
	}
	// S 形：先凸后凹，拐点在 i=4（二阶差分由正变号处）
	s := []float64{0, 1, 3, 6, 9, 11, 12, 12.5, 12.75}
	turns := RadiusTurns(s)
	found := false
	for _, idx := range turns {
		if idx == 4 {
			found = true
		}
	}
	if !found {
		t.Fatalf("S 形序列拐点 4 未标出, got %v", turns)
	}
}
