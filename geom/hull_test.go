package geom

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"testing"
)

// ---------- 测试数据生成 ----------

type genFunc func(n int, seed int64) []Point

func uniformSquare(n int, seed int64) []Point {
	r := rand.New(rand.NewSource(seed))
	out := make([]Point, n)
	for i := range out {
		out[i] = Point{r.Float64() * 100, r.Float64() * 100}
	}
	return out
}

// clustered 生成三团高斯聚集点。
func clustered(n int, seed int64) []Point {
	r := rand.New(rand.NewSource(seed))
	centers := []Point{{20, 20}, {78, 24}, {52, 82}}
	out := make([]Point, n)
	for i := range out {
		c := centers[i%len(centers)]
		out[i] = Point{c.X + r.NormFloat64()*6, c.Y + r.NormFloat64()*6}
	}
	return out
}

func twoClusters(n int, seed int64) []Point {
	r := rand.New(rand.NewSource(seed))
	out := make([]Point, n)
	for i := range out {
		if i%2 == 0 {
			out[i] = Point{15 + r.NormFloat64()*4, 50 + r.NormFloat64()*4}
		} else {
			out[i] = Point{85 + r.NormFloat64()*4, 50 + r.NormFloat64()*4}
		}
	}
	return out
}

// collinear 生成严格落在同一条直线上的点（整数坐标保证无浮点舍入）。
func collinear(n int, seed int64) []Point {
	r := rand.New(rand.NewSource(seed))
	out := make([]Point, n)
	for i := range out {
		x := float64(r.Intn(1000))
		out[i] = Point{x, 3*x + 7}
	}
	return out
}

// rings 在 rings 个同心圆上均分生成点；环数固定，点数只影响每环密度。
func rings(n, ringCount int, seed int64) []Point {
	_ = seed
	if ringCount < 2 {
		ringCount = 2
	}
	out := make([]Point, 0, n)
	for i := 0; i < ringCount; i++ {
		rad := 10 + 18*float64(i)
		cnt := n / ringCount
		if i == ringCount-1 {
			cnt += n % ringCount
		}
		off := 0.37 * float64(i)
		for j := 0; j < cnt; j++ {
			a := 2*math.Pi*float64(j)/float64(cnt) + off
			out = append(out, Point{rad * math.Cos(a), rad * math.Sin(a)})
		}
	}
	return out
}

func countUnique(pts []Point) int {
	m := map[Point]struct{}{}
	for _, p := range pts {
		m[p] = struct{}{}
	}
	return len(m)
}

func ringGen(ringCount int) genFunc {
	return func(n int, seed int64) []Point { return rings(n, ringCount, seed) }
}

// ---------- 不变量校验 ----------

func checkLayerInvariants(t *testing.T, name string, pts []Point) *LayerSet {
	t.Helper()
	set, err := Peel(pts)
	if err != nil {
		t.Fatalf("%s: Peel 失败 %v", name, err)
	}
	if err := ValidateLayers(set.Layers); err != nil {
		t.Fatalf("%s: 层校验失败 %v", name, err)
	}
	if set.StoppedAt != len(set.Layers) {
		t.Fatalf("%s: StoppedAt=%d 与层数 %d 不一致", name, set.StoppedAt, len(set.Layers))
	}
	if set.Remaining >= 3 {
		t.Fatalf("%s: 停止时剩余点 %d 应不足三个", name, set.Remaining)
	}
	// 不变量一：各层面积随层号单调不增。
	for i := 1; i < len(set.Layers); i++ {
		prev, cur := set.Layers[i-1].Area, set.Layers[i].Area
		if cur > prev+1e-9*math.Max(1, math.Abs(prev)) {
			t.Fatalf("%s: 第 %d 层面积 %.12g 大于上一层 %.12g", name, i, cur, prev)
		}
	}
	// 点数守恒：各层剥去点数 + 剩余点数 == 互异点数。
	consumed := set.Remaining
	for _, l := range set.Layers {
		consumed += l.Size
		if len(l.Hull) < 3 || l.Area <= 0 {
			t.Fatalf("%s: 第 %d 层退化 hull=%d area=%g", name, l.Index, len(l.Hull), l.Area)
		}
		if l.Size < len(l.Hull) {
			t.Fatalf("%s: 第 %d 层剥去点数 %d 小于 hull 顶点数 %d", name, l.Index, l.Size, len(l.Hull))
		}
	}
	if want := countUnique(pts); consumed != want {
		t.Fatalf("%s: 点数不守恒，消耗 %d，互异 %d", name, consumed, want)
	}
	// 每层 hull 都应为逆时针凸多边形（拐角严格左转）。
	for _, l := range set.Layers {
		for i := range l.Hull {
			a, b, c := l.Hull[i], l.Hull[(i+1)%len(l.Hull)], l.Hull[(i+2)%len(l.Hull)]
			if cross(a, b, c) <= 0 {
				t.Fatalf("%s: 第 %d 层 hull 非严格凸 / 非逆时针", name, l.Index)
			}
		}
	}
	return set
}

// TestPeelInvariantsTwoTiers 在两档点数下验证面积单调不增。
func TestPeelInvariantsTwoTiers(t *testing.T) {
	tiers := []int{300, 2400}
	gens := []struct {
		name string
		gen  genFunc
	}{
		{"随机", uniformSquare},
		{"聚集", clustered},
		{"双簇", twoClusters},
		{"同心环", ringGen(6)},
	}
	for _, tier := range tiers {
		for _, g := range gens {
			g, tier := g, tier
			t.Run(fmt.Sprintf("%s/n=%d", g.name, tier), func(t *testing.T) {
				pts := g.gen(tier, 42)
				set := checkLayerInvariants(t, g.name, pts)
				t.Logf("%s n=%d: 层数=%d 停止层=%d 剩余=%d", g.name, tier, len(set.Layers), set.StoppedAt, set.Remaining)
			})
		}
	}
}

// TestLayerCountNotLinear 验证层数不随点数线性增长：
// 点数放大 8 倍，层数应亚线性增长（拟合幂次明显小于 1）。
func TestLayerCountNotLinear(t *testing.T) {
	const small, large = 300, 2400
	const growth = float64(large) / float64(small)
	gens := []struct {
		name string
		gen  genFunc
	}{
		{"随机", uniformSquare},
		{"聚集", clustered},
		{"双簇", twoClusters},
		{"同心环", ringGen(6)},
	}
	for _, g := range gens {
		g := g
		t.Run(g.name, func(t *testing.T) {
			ls := checkLayerInvariants(t, g.name, g.gen(small, 7))
			ll := checkLayerInvariants(t, g.name, g.gen(large, 7))
			if len(ls.Layers) == 0 {
				t.Fatalf("%s: 小规模层数为零", g.name)
			}
			ratio := float64(len(ll.Layers)) / float64(len(ls.Layers))
			// 层数随点数的幂次；线性增长对应 1.0。
			exponent := math.Log(ratio) / math.Log(growth)
			if exponent >= 0.85 {
				t.Fatalf("%s: 点数 %.0f 倍时层数由 %d 增至 %d（%.2f 倍，幂次 %.2f），接近线性增长",
					g.name, growth, len(ls.Layers), len(ll.Layers), ratio, exponent)
			}
			// 同心环层数应恒等于环数，不随每环点数增长。
			if g.name == "同心环" && (len(ls.Layers) != 6 || len(ll.Layers) != 6) {
				t.Fatalf("同心环层数应为 6，实得 %d / %d", len(ls.Layers), len(ll.Layers))
			}
			if len(ll.Layers) >= large/8 {
				t.Fatalf("%s: 层数 %d 过大，接近线性", g.name, len(ll.Layers))
			}
			t.Logf("%s: n=%d 层数=%d, n=%d 层数=%d, 比值=%.2f, 增长幂次=%.2f",
				g.name, small, len(ls.Layers), large, len(ll.Layers), ratio, exponent)
		})
	}
}

// ---------- 单层凸包行为 ----------

func TestConvexHullUnchanged(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	pts := append([]Point(nil), sq...)
	// 加入大量内部点与边上的共线点，单层凸包只应保留四个拐角。
	pts = append(pts, Point{5, 5}, Point{2, 2}, Point{8, 8}, Point{3, 7}, Point{0, 5}, Point{5, 0})
	hull, err := ConvexHull(pts)
	if err != nil {
		t.Fatalf("ConvexHull 失败 %v", err)
	}
	want := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	if len(hull) != len(want) {
		t.Fatalf("单层凸包应只有 %d 个拐角，实得 %d: %v", len(want), len(hull), hull)
	}
	for i := range want {
		if hull[i] != want[i] {
			t.Fatalf("单层凸包第 %d 点错误: 实得 %v 期望 %v", i, hull[i], want[i])
		}
	}
	// 单层凸包必须与分层剥离的第 0 层一致。
	layers, err := Peel(pts)
	if err != nil {
		t.Fatalf("Peel 失败 %v", err)
	}
	if len(layers.Layers[0].Hull) != len(hull) {
		t.Fatalf("第 0 层与单层凸包不一致: %v vs %v", layers.Layers[0].Hull, hull)
	}
	for i := range hull {
		if layers.Layers[0].Hull[i] != hull[i] {
			t.Fatalf("第 0 层第 %d 点与单层凸包不一致", i)
		}
	}
}

func TestConvexHullErrors(t *testing.T) {
	cases := []struct {
		name  string
		pts   []Point
		want  error
		layer int
	}{
		{"顶点不足", []Point{{0, 0}, {1, 1}}, ErrTooFewPts, 0},
		{"重复退化", []Point{{1, 1}, {1, 1}, {1, 1}}, ErrDegenerate, 0},
		{"全部共线", []Point{{0, 0}, {1, 1}, {2, 2}, {3, 3}}, ErrAllCollinear, 0},
	}
	for _, c := range cases {
		_, err := ConvexHull(c.pts)
		var le *LayerError
		if !errors.As(err, &le) || !errors.Is(err, c.want) || le.Layer != c.layer {
			t.Fatalf("%s: 期望第 %d 层 %v，实得 %v", c.name, c.layer, c.want, err)
		}
	}
}

// ---------- 分层错误均报出层号 ----------

func TestLayerErrorsReportLayer(t *testing.T) {
	// 外圈方形剥掉后只剩三个共线点：共线错误应报到第 1 层。
	squareAndLine := []Point{
		{0, 0}, {10, 0}, {10, 10}, {0, 10},
		{2, 5}, {5, 5}, {8, 5},
	}
	_, err := Peel(squareAndLine)
	var le *LayerError
	if !errors.As(err, &le) || !errors.Is(err, ErrAllCollinear) || le.Layer != 1 {
		t.Fatalf("共线应报第 1 层，实得 %v", err)
	}

	// 剥到第 1 层时剩余点仅 1 个，请求该层应报出层号 1。
	five := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}, {5, 5}}
	set, err := Peel(five)
	if err != nil {
		t.Fatalf("Peel 失败 %v", err)
	}
	if set.Remaining != 1 || set.StoppedAt != 1 {
		t.Fatalf("应剥 1 层后剩 1 点，实得层数 %d 剩余 %d", set.StoppedAt, set.Remaining)
	}
	if _, err := TrimLayers(set.Layers, 2); !errors.As(err, &le) || !errors.Is(err, ErrTooFewRemaining) || le.Layer != 1 {
		t.Fatalf("剩点不足应报第 1 层，实得 %v", err)
	}
	if _, err := LayerAt(set.Layers, 1); !errors.As(err, &le) || !errors.Is(err, ErrTooFewRemaining) || le.Layer != 1 {
		t.Fatalf("剩点不足取第 1 层应报第 1 层，实得 %v", err)
	}
	if _, err := CoverageAt(set.Layers, 2); !errors.As(err, &le) || !errors.Is(err, ErrTooFewRemaining) || le.Layer != 1 {
		t.Fatalf("覆盖层数不足应报首个不可用层 1，实得 %v", err)
	}
	if _, err := Peel([]Point{{0, 0}, {1, 1}}); !errors.As(err, &le) || !errors.Is(err, ErrTooFewPts) || le.Layer != 0 {
		t.Fatalf("初始点数不足应报第 0 层，实得 %v", err)
	}

	// 层号跳号：索引 0 后直接给 2，应在缺号处报出层号 1。
	gapped := []HullLayer{
		{Index: 0, Hull: []Point{{0, 0}, {1, 0}, {0, 1}}, Area: 0.5, Size: 3},
		{Index: 2, Hull: []Point{{0, 0}, {1, 0}, {0, 1}}, Area: 0.5, Size: 3},
	}
	if err := ValidateLayers(gapped); !errors.As(err, &le) || !errors.Is(err, ErrLayerGap) || le.Layer != 1 {
		t.Fatalf("跳号应报第 1 层，实得 %v", err)
	}
	if _, err := LayerAt(gapped, 0); !errors.As(err, &le) || !errors.Is(err, ErrLayerGap) || le.Layer != 1 {
		t.Fatalf("LayerAt 遇跳号应报第 1 层，实得 %v", err)
	}
}

// ---------- 面积覆盖反解与裁剪 ----------

func TestCoverageSolveAndTrim(t *testing.T) {
	pts := uniformSquare(800, 3)
	set, err := Peel(pts)
	if err != nil {
		t.Fatalf("Peel 失败 %v", err)
	}
	if len(set.Layers) < 3 {
		t.Fatalf("样本层数过少: %d", len(set.Layers))
	}
	outer := set.Layers[0].Area
	// 覆盖按嵌套环带计算：保留 keep 层覆盖 (A[0]-A[keep])/A[0]，到最内层为 1。
	// 另用 Band 独立复算，二者应一致。
	coverage := func(keep int) float64 {
		if keep == len(set.Layers) {
			return 1
		}
		return (outer - set.Layers[keep].Area) / outer
	}
	bandTotal, bandCum := 0.0, 0.0
	for _, l := range set.Layers {
		bandTotal += l.Band
	}
	for i, l := range set.Layers {
		bandCum += l.Band
		if math.Abs(bandCum/bandTotal-coverage(i+1)) > 1e-12 {
			t.Fatalf("环带累加口径与面积差口径不一致于 keep=%d", i+1)
		}
		if i+1 < len(set.Layers) && math.Abs(l.Band-(l.Area-set.Layers[i+1].Area)) > 1e-9*math.Max(1, l.Area) {
			t.Fatalf("第 %d 层环带面积不符: %g vs %g", i, l.Band, l.Area-set.Layers[i+1].Area)
		}
	}
	if math.Abs(bandTotal-outer) > 1e-9*outer {
		t.Fatalf("环带总和 %g 应等于最外层面积 %g", bandTotal, outer)
	}
	if coverage(0) != 0 {
		t.Fatalf("保留 0 层覆盖应为 0，实得 %v", coverage(0))
	}
	if coverage(len(set.Layers)) != 1 {
		t.Fatalf("保留全部层覆盖应为 1")
	}
	prevKeep := 0
	for _, target := range []float64{0.10, 0.35, 0.60, 0.85, 0.95, 1.0} {
		plan, err := SolveCoverage(set.Layers, target)
		if err != nil {
			t.Fatalf("目标 %v 反解失败 %v", target, err)
		}
		if plan.Keep < 1 || plan.Keep <= prevKeep {
			t.Fatalf("目标 %v 反解层数 %d 应随目标单调增（上一档 %d）", target, plan.Keep, prevKeep)
		}
		prevKeep = plan.Keep
		if plan.LastLayer != plan.Keep-1 {
			t.Fatalf("目标 %v: LastLayer=%d 与 Keep=%d 不符", target, plan.LastLayer, plan.Keep)
		}
		if plan.Peeled != len(set.Layers)-plan.Keep {
			t.Fatalf("目标 %v: Peeled=%d 与 Keep=%d 不符", target, plan.Peeled, plan.Keep)
		}
		if plan.Coverage < target-1e-12 || plan.Coverage > 1+1e-12 {
			t.Fatalf("目标 %v: 实际覆盖 %.15g 越界", target, plan.Coverage)
		}
		// 独立复算并与 CoverageAt 对齐。
		if math.Abs(coverage(plan.Keep)-plan.Coverage) > 1e-12 {
			t.Fatalf("目标 %v: 覆盖复算 %.15g 与反解 %.15g 不符", target, coverage(plan.Keep), plan.Coverage)
		}
		if got, err := CoverageAt(set.Layers, plan.Keep); err != nil || math.Abs(got-plan.Coverage) > 1e-12 {
			t.Fatalf("目标 %v: CoverageAt 与反解不符 %v %v", target, got, err)
		}
		// 少保留一层则达不到目标（Keep 必须最小）。
		if plan.Keep > 1 && coverage(plan.Keep-1) >= target-1e-12 {
			t.Fatalf("目标 %v: 反解层数 %d 不是最小（少一层已覆盖 %.15g）",
				target, plan.Keep, coverage(plan.Keep-1))
		}
	}
	// 覆盖目标非法。
	for _, bad := range []float64{0, -0.1, 1.01, math.NaN()} {
		if _, err := SolveCoverage(set.Layers, bad); !errors.Is(err, ErrBadCoverage) {
			t.Fatalf("目标 %v 应报 ErrBadCoverage，实得 %v", bad, err)
		}
	}
	// 按反解层号裁剪输出。
	trimmed, plan, err := PeelToCoverage(pts, 0.8)
	if err != nil {
		t.Fatalf("PeelToCoverage 失败 %v", err)
	}
	if len(trimmed) != plan.Keep {
		t.Fatalf("裁剪层数 %d 与反解 %d 不符", len(trimmed), plan.Keep)
	}
	if err := ValidateLayers(trimmed); err != nil {
		t.Fatalf("裁剪结果不合法 %v", err)
	}
	if _, err := TrimLayers(set.Layers, len(set.Layers)+1); !errors.Is(err, ErrTooFewRemaining) {
		t.Fatalf("超额裁剪应报错，实得 %v", err)
	}
}

// ---------- 收尾报告：反解层号与面积表 ----------

func TestCoverageReport(t *testing.T) {
	datasets := []struct {
		name string
		pts  []Point
	}{
		{"随机", uniformSquare(900, 11)},
		{"聚集", clustered(900, 11)},
		{"双簇", twoClusters(900, 11)},
		{"同心环", rings(900, 6, 11)},
	}
	for _, ds := range datasets {
		set, err := Peel(ds.pts)
		if err != nil {
			t.Fatalf("%s: Peel 失败 %v", ds.name, err)
		}
		t.Logf("== %s n=%d 层数=%d 剩余=%d ==", ds.name, len(ds.pts), len(set.Layers), set.Remaining)
		outer := set.Layers[0].Area
		for _, l := range set.Layers {
			t.Logf("  层 %d: 面积=%10.4f 占外层面积=%6.2f%% 剥去点数=%d 拐角=%d",
				l.Index, l.Area, 100*l.Area/outer, l.Size, len(l.Hull))
		}
		for _, target := range []float64{0.5, 0.8, 0.95} {
			plan, err := SolveCoverage(set.Layers, target)
			if err != nil {
				t.Fatalf("%s: 反解 %v 失败 %v", ds.name, target, err)
			}
			t.Logf("  目标 %.2f -> 保留 %d 层（最深第 %d 层），实际覆盖 %.4f",
				target, plan.Keep, plan.LastLayer, plan.Coverage)
		}
	}
}

func TestCollinearReportsLayer(t *testing.T) {
	for _, n := range []int{300, 2400} {
		pts := collinear(n, 5)
		_, err := Peel(pts)
		var le *LayerError
		if !errors.As(err, &le) || !errors.Is(err, ErrAllCollinear) || le.Layer != 0 {
			t.Fatalf("全共线 n=%d 应报第 0 层共线错误，实得 %v", n, err)
		}
	}
	// 外圈非退化、内层共线时，层号应指向内层。
	pts := []Point{
		{0, 0}, {10, 0}, {10, 10}, {0, 10},
		{2, 5}, {5, 5}, {8, 5},
	}
	_, err := Peel(pts)
	var le *LayerError
	if !errors.As(err, &le) || !errors.Is(err, ErrAllCollinear) || le.Layer != 1 {
		t.Fatalf("内层共线应报第 1 层，实得 %v", err)
	}
}

// ---------- 与独立实现交叉校验 ----------

// jarvisHull 用礼物包装法独立求凸包拐角（逆时针），用于交叉校验分层结果。
// 共线时取最远点入包，因此只保留严格拐角。
func jarvisHull(pts []Point) []Point {
	n := len(pts)
	if n < 3 {
		return nil
	}
	start := 0
	for i := 1; i < n; i++ {
		if pts[i].X < pts[start].X || (pts[i].X == pts[start].X && pts[i].Y < pts[start].Y) {
			start = i
		}
	}
	hull := []Point{pts[start]}
	cur := start
	for {
		next := -1
		for i := 0; i < n; i++ {
			if i == cur {
				continue
			}
			if next == -1 {
				next = i
				continue
			}
			c := cross(pts[cur], pts[next], pts[i])
			if c < 0 {
				next = i
			} else if c == 0 {
				// 共线取更远者。
				d1 := math.Hypot(pts[next].X-pts[cur].X, pts[next].Y-pts[cur].Y)
				d2 := math.Hypot(pts[i].X-pts[cur].X, pts[i].Y-pts[cur].Y)
				if d2 > d1 {
					next = i
				}
			}
		}
		if next == start {
			break
		}
		hull = append(hull, pts[next])
		cur = next
		if len(hull) > n {
			return nil
		}
	}
	if len(hull) < 3 {
		return nil
	}
	return hull
}

// strictlyInside 判断点是否严格在凸多边形内部（不在任一边上）。
func strictlyInside(hull []Point, p Point) bool {
	if len(hull) < 3 {
		return false
	}
	for i := range hull {
		a, b := hull[i], hull[(i+1)%len(hull)]
		if cross(a, b, p) <= 0 {
			return false
		}
	}
	return true
}

func TestPeelMatchesIndependentHull(t *testing.T) {
	datasets := []struct {
		name string
		pts  []Point
	}{
		{"随机", uniformSquare(240, 21)},
		{"聚集", clustered(240, 21)},
		{"双簇", twoClusters(240, 21)},
		{"同心环", rings(240, 6, 21)},
	}
	for _, ds := range datasets {
		set, err := Peel(ds.pts)
		if err != nil {
			t.Fatalf("%s: Peel 失败 %v", ds.name, err)
		}
		rem := append([]Point(nil), ds.pts...)
		// 去重。
		{
			m := map[Point]bool{}
			u := rem[:0]
			for _, p := range rem {
				if !m[p] {
					m[p] = true
					u = append(u, p)
				}
			}
			rem = u
		}
		idx := 0
		for len(rem) >= 3 {
			h := jarvisHull(rem)
			if len(h) < 3 {
				break
			}
			if idx >= len(set.Layers) {
				t.Fatalf("%s: 独立实现多出第 %d 层", ds.name, idx)
			}
			l := set.Layers[idx]
			if len(l.Hull) != len(h) {
				t.Fatalf("%s 第 %d 层顶点数不符: peel=%d jarvis=%d", ds.name, idx, len(l.Hull), len(h))
			}
			// 同为逆时针，起点可能不同：逐个按顺序比较（两者均从最左点起）。
			for i := range h {
				if l.Hull[i] != h[i] {
					t.Fatalf("%s 第 %d 层第 %d 点不符: peel=%v jarvis=%v", ds.name, idx, i, l.Hull[i], h[i])
				}
			}
			if got := math.Abs(polygonArea(h)); math.Abs(got-l.Area) > 1e-9*math.Max(1, got) {
				t.Fatalf("%s 第 %d 层面积不符: peel=%.12g jarvis=%.12g", ds.name, idx, l.Area, got)
			}
			var rest []Point
			for _, p := range rem {
				// 剥去落在 hull 边界上的点（含拐角与边上共线点），只留严格内部点。
				if strictlyInside(h, p) {
					rest = append(rest, p)
				}
			}
			rem = rest
			idx++
		}
		if idx != len(set.Layers) {
			t.Fatalf("%s: 层数不符: peel=%d jarvis=%d", ds.name, len(set.Layers), idx)
		}
		if len(rem) != set.Remaining {
			t.Fatalf("%s: 剩余点数不符: peel=%d jarvis=%d", ds.name, set.Remaining, len(rem))
		}
	}
}

// TestCoverageMonotone 验证覆盖比例随保留层数单调不减，且边界值为 0 / 1。
func TestCoverageMonotone(t *testing.T) {
	datasets := [][]Point{
		uniformSquare(500, 33),
		clustered(500, 33),
		twoClusters(500, 33),
		rings(500, 6, 33),
	}
	for di, pts := range datasets {
		set, err := Peel(pts)
		if err != nil {
			t.Fatalf("数据集 %d: Peel 失败 %v", di, err)
		}
		prev := -1.0
		for keep := 0; keep <= len(set.Layers); keep++ {
			cov, err := CoverageAt(set.Layers, keep)
			if err != nil {
				t.Fatalf("数据集 %d keep=%d: CoverageAt 失败 %v", di, keep, err)
			}
			if keep == 0 && cov != 0 {
				t.Fatalf("数据集 %d: keep=0 覆盖应为 0，实得 %v", di, cov)
			}
			if keep == len(set.Layers) && math.Abs(cov-1) > 1e-15 {
				t.Fatalf("数据集 %d: keep=全部 覆盖应为 1，实得 %v", di, cov)
			}
			if cov < prev-1e-12 {
				t.Fatalf("数据集 %d: keep=%d 覆盖 %.15g 小于上一层 %.15g", di, keep, cov, prev)
			}
			prev = cov
		}
	}
}

// ---------- 边界输入 ----------

func TestPeelEdgeInputs(t *testing.T) {
	// 空输入与 nil。
	for _, pts := range [][]Point{nil, {}} {
		if _, err := Peel(pts); !errors.Is(err, ErrTooFewPts) {
			t.Fatalf("空输入应报 ErrTooFewPts，实得 %v", err)
		}
		if _, err := ConvexHull(pts); !errors.Is(err, ErrTooFewPts) {
			t.Fatalf("空输入 ConvexHull 应报 ErrTooFewPts，实得 %v", err)
		}
	}
	// 恰为三角形：一层，无剩余。
	tri := []Point{{0, 0}, {4, 0}, {0, 4}}
	set, err := Peel(tri)
	if err != nil || len(set.Layers) != 1 || set.Remaining != 0 || set.StoppedAt != 1 {
		t.Fatalf("三角形应恰一层且无剩余: %+v %v", set, err)
	}
	// 重复点被合并为一个互异点：(5,5) 落在方形内部，剥掉外圈后剩这一点，无法再成层。
	dup := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}, {5, 5}, {5, 5}, {5, 5}}
	set, err = Peel(dup)
	if err != nil {
		t.Fatalf("含重复点 Peel 失败 %v", err)
	}
	if len(set.Layers) != 1 || set.Layers[0].Size != 4 || set.Remaining != 1 {
		t.Fatalf("重复点应合并: 层数=%d size=%d 剩余=%d", len(set.Layers), set.Layers[0].Size, set.Remaining)
	}
	// 输入切片不得被修改。
	orig := append([]Point(nil), dup...)
	if _, err := Peel(dup); err != nil {
		t.Fatal(err)
	}
	for i := range orig {
		if dup[i] != orig[i] {
			t.Fatalf("Peel 修改了调用方切片")
		}
	}
	// TrimLayers 返回副本，不与原层共享底层数组。
	layers := []HullLayer{
		{Index: 0, Hull: []Point{{0, 0}, {2, 0}, {0, 2}}, Area: 2, Band: 1, Size: 3},
		{Index: 1, Hull: []Point{{0, 0}, {1, 0}, {0, 1}}, Area: 0.5, Band: 0.5, Size: 3},
	}
	got, err := TrimLayers(layers, 1)
	if err != nil {
		t.Fatal(err)
	}
	got[0].Area = 99
	if layers[0].Area != 2 {
		t.Fatalf("TrimLayers 未返回副本")
	}
}
