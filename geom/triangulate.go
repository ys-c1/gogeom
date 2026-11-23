package geom

import (
	"math"
	"sort"
)

// Triangle 是三角化输出的基本单元，字段为三个顶点（坐标取自输入）。
type Triangle [3]Point

// TriReport 是三角化的完整报告：
// Triangles 为输出三角网格；Rounds 为实际耳切轮数（每轮切除一只耳朵、
// 产出一个三角形）；EstRounds 为按顶点数估算的轮数；RoundsNote 说明
// 实际轮次与估算不一致的原因；CollinearRemoved 记录被近似共线清理掉的
// 输入顶点号；Merged 记录带洞桥接/自交拆分过程中使用的顶点对；
// SplitFaces 为自交拆分后参与填充的面数（简单多边形为 1）。
type TriReport struct {
	Triangles        []Triangle
	Rounds           int
	EstRounds        int
	RoundsNote       string
	CollinearRemoved []int
	Bridges          [][2]int
	SplitFaces       int
	TotalTriArea     float64
	OriginalArea     float64
	AreaDiff         float64
	AreaRelDiff      float64
	OverlapPairs     int
}

const (
	// epsRel 为相对坐标尺度的容差，用于近似共线判定与浮点比较
	epsRel = 1e-10
)

func epsFor(pts []Point) float64 {
	var minx, miny, maxx, maxy float64
	for i, p := range pts {
		if i == 0 {
			minx, maxx, miny, maxy = p.X, p.X, p.Y, p.Y
			continue
		}
		if p.X < minx {
			minx = p.X
		}
		if p.X > maxx {
			maxx = p.X
		}
		if p.Y < miny {
			miny = p.Y
		}
		if p.Y > maxy {
			maxy = p.Y
		}
	}
	span := math.Max(maxx-minx, maxy-miny)
	if span == 0 {
		return epsRel
	}
	return epsRel * span
}

// cross 返回向量 (b-a)×(c-a)
func cross(a, b, c Point) float64 {
	return (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
}

func orient(a, b, c Point, eps float64) int {
	v := cross(a, b, c)
	// 共线判据：c 到直线 ab 的垂直距离 |cross|/|ab| <= eps。
	// tol 只按单条边长度缩放，任意坐标尺度下阈值统一，不会把
	// 正常转角（如正方形对角的交叉）误吞为共线。
	lab := math.Hypot(b.X-a.X, b.Y-a.Y)
	tol := eps * math.Max(1, lab)
	if v > tol {
		return 1
	}
	if v < -tol {
		return -1
	}
	return 0
}

func absAreaTri(a, b, c Point) float64 {
	return math.Abs(cross(a, b, c)) / 2
}

// pointInTriContain 判定 p 是否严格落在三角形内部（边界上不算严格内部）。
// 三条边的定向必须同号且均非零。
func pointInTriContain(p, a, b, c Point, eps float64) bool {
	o1 := orient(a, b, p, eps)
	o2 := orient(b, c, p, eps)
	o3 := orient(c, a, p, eps)
	return (o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0)
}

// dedupPoints 合并坐标完全相同的相邻顶点，返回新切片与被删除的原始下标。
// 首点与末点相同（闭合环）时也会去掉末点。
func dedupPoints(pts []Point, eps float64) ([]Point, []int) {
	var out []Point
	var removed []int
	n := len(pts)
	for i, p := range pts {
		if len(out) > 0 {
			last := out[len(out)-1]
			if math.Abs(last.X-p.X) <= eps && math.Abs(last.Y-p.Y) <= eps {
				removed = append(removed, i)
				continue
			}
		}
		out = append(out, p)
	}
	if len(out) > 1 {
		first, last := out[0], out[len(out)-1]
		if math.Abs(first.X-last.X) <= eps && math.Abs(first.Y-last.Y) <= eps {
			// 末点是首点的重复，需要映射回原始下标：找最后一个输出点的原始号
			removed = append(removed, n-1)
			out = out[:len(out)-1]
		}
	}
	return out, removed
}

// removeCollinear 清理近似共线（含重复）顶点，返回清理后的环、被删顶点号。
// 若清理后不足三个非共线点，ok 为 false（即全共线）。
func removeCollinear(pts []Point) (out []Point, removed []int, ok bool) {
	if len(pts) < 3 {
		return pts, nil, false
	}
	eps := epsFor(pts)
	work, dupRemoved := dedupPoints(pts, eps)
	removed = append(removed, dupRemoved...)
	if len(work) < 3 {
		return work, removed, false
	}

	// 反复扫描，直到一轮中没有点可删（删除会让相邻点变成新的共线三元组）
	for {
		any := false
		m := len(work)
		var next []Point
		var nextRemoved []int
		// 用“原始点 -> 原下标”的映射：work 中的点仍能在 pts 中按坐标查找，
		// 为稳妥起见这里维护索引链。
		for i := 0; i < m; i++ {
			a := work[(i-1+m)%m]
			b := work[i]
			c := work[(i+1)%m]
			if orient(a, b, c, eps) == 0 {
				// 共线：仅当 b 落在 a-c 线段上才删除（环退化时可能反向）
				if onSegment(a, b, c, eps) {
					idx := indexOfPoint(pts, b, eps)
					nextRemoved = append(nextRemoved, idx)
					any = true
					continue
				}
			}
			next = append(next, b)
		}
		work = next
		removed = append(removed, nextRemoved...)
		if !any || len(work) < 3 {
			break
		}
	}
	if len(work) < 3 {
		return work, removed, false
	}
	// 全共线判定：所有相邻三元组定向均为零（不能用总面积==0，
	// 自交多边形如蝴蝶结的有向面积也可能恰好为 0）
	if allCollinearOriented(work, epsFor(work)) {
		return work, removed, false
	}
	return work, removed, true
}

// allCollinearOriented 检查环上每个顶点与其前后邻居是否全部共线
func allCollinearOriented(pts []Point, eps float64) bool {
	n := len(pts)
	for i := 0; i < n; i++ {
		if orient(pts[(i-1+n)%n], pts[i], pts[(i+1)%n], eps) != 0 {
			return false
		}
	}
	return true
}

func indexOfPoint(pts []Point, q Point, eps float64) int {
	for i, p := range pts {
		if math.Abs(p.X-q.X) <= eps && math.Abs(p.Y-q.Y) <= eps {
			return i
		}
	}
	return -1
}

// onSegment 判定 b 是否落在 a-c 线段上（假定三者共线）
func onSegment(a, b, c Point, eps float64) bool {
	return math.Min(a.X, c.X)-eps <= b.X && b.X <= math.Max(a.X, c.X)+eps &&
		math.Min(a.Y, c.Y)-eps <= b.Y && b.Y <= math.Max(a.Y, c.Y)+eps
}

// signedArea2 返回环的有向面积两倍（CCW 为正）
func signedArea2(ring []Point) float64 {
	var s float64
	n := len(ring)
	for i := 0; i < n; i++ {
		a, b := ring[i], ring[(i+1)%n]
		s += a.X*b.Y - b.X*a.Y
	}
	return s
}

// ensureCCW 把环调整为逆时针
func ensureCCW(ring []Point) []Point {
	if signedArea2(ring) < 0 {
		out := make([]Point, len(ring))
		for i := range ring {
			out[i] = ring[len(ring)-1-i]
		}
		return out
	}
	return ring
}

// earInFill 在耳三角形内部做网格采样，确认所有采样点都在填充区内；
// 采样点避开边界，只要有一点落在洞区即判定该耳跨洞、不可切除。
func earInFill(a, b, c Point, fillTest func(Point) bool) bool {
	const n = 8
	for i := 1; i < n; i++ {
		for j := 1; j < n-i; j++ {
			u := float64(i) / n
			v := float64(j) / n
			w := 1 - u - v
			if w <= 0 {
				continue
			}
			p := Point{
				X: a.X*w + b.X*u + c.X*v,
				Y: a.Y*w + b.Y*u + c.Y*v,
			}
			if !fillTest(p) {
				return false
			}
		}
	}
	return true
}

// earClip 对一条简单、无自交、非共线的环执行耳切。
// fillTest 为可选的填充区谓词：非 nil 时，每只候选耳的质心必须落在
// 填充区内，用于带洞合并环排除跨越洞区的假耳；普通简单环传 nil。
// 若因数值原因找不到合格耳朵，进入“最瘦耳朵”兜底并在 rescueNote 中说明。
func earClip(ring []Point, origIdx []int, fillTest func(Point) bool) (tris []Triangle, rounds int, rescueNote string) {
	n := len(ring)
	eps := epsFor(ring)
	// 活动链：prev/next 下标对应当前切片位置，删除采用标记法
	alive := make([]bool, n)
	prev := make([]int, n)
	next := make([]int, n)
	for i := 0; i < n; i++ {
		alive[i] = true
		prev[i] = (i - 1 + n) % n
		next[i] = (i + 1) % n
	}

	isEar := func(i int) bool {
		p, pr, nx := i, prev[i], next[i]
		a, b, c := ring[pr], ring[p], ring[nx]
		if orient(a, b, c, eps) <= 0 {
			return false // 凸耳要求严格左转（CCW 环）
		}
		// 带洞合并环：耳三角形内多点采样，任一点落到洞区即说明该耳跨洞
		if fillTest != nil && !earInFill(a, b, c, fillTest) {
			return false
		}
		// 耳内不得包含任何其它活动顶点
		for k := 0; k < n; k++ {
			if !alive[k] || k == p || k == pr || k == nx {
				continue
			}
			if pointInTriContain(ring[k], a, b, c, eps) {
				return false
			}
		}
		return true
	}

	aliveCount := n
	totalRounds := 0
	for aliveCount > 3 {
		found := -1
		for i := 0; i < n; i++ {
			if alive[i] && isEar(i) {
				found = i
				break
			}
		}
		if found < 0 {
			// 兜底：在所有凸顶点里取耳面积最小者强制切除。
			// 正常简单多边形必有至少两只耳朵，走到这里只可能是数值共线，
			// 属于实际轮次估算之外的异常路径，需在报告中说明。
			best := -1
			bestArea := math.Inf(1)
			for i := 0; i < n; i++ {
				if !alive[i] {
					continue
				}
				a, b, c := ring[prev[i]], ring[i], ring[next[i]]
				if orient(a, b, c, eps) <= 0 {
					continue
				}
				ar := absAreaTri(a, b, c)
				if ar < bestArea {
					bestArea, best = ar, i
				}
			}
			if best < 0 {
				// 连凸顶点都没有：取面积最小的三元组，保证流程可终止
				for i := 0; i < n; i++ {
					if !alive[i] {
						continue
					}
					ar := absAreaTri(ring[prev[i]], ring[i], ring[next[i]])
					if ar < bestArea {
						bestArea, best = ar, i
					}
				}
			}
			found = best
			rescueNote = "存在数值近似共线顶点，标准耳判定失效，采用最小面积耳兜底切除"
		}
		pr, nx := prev[found], next[found]
		tris = append(tris, Triangle{ring[pr], ring[found], ring[nx]})
		totalRounds++
		alive[found] = false
		next[pr] = nx
		prev[nx] = pr
		aliveCount--
	}
	// 剩余三个活动点构成最后一个三角形，作为最后一轮切除，
	// 因此简单多边形的耳切轮数恒为 n-2，与按顶点数的估算一致。
	var tail []int
	for i := 0; i < n; i++ {
		if alive[i] {
			tail = append(tail, i)
		}
	}
	if len(tail) == 3 {
		tris = append(tris, Triangle{ring[tail[0]], ring[tail[1]], ring[tail[2]]})
		totalRounds++
	}
	return tris, totalRounds, rescueNote
}

// ---------------------------------------------------------------------------
// 自交检测与拆分（平面子图面遍历 + 奇偶填充规则）
// ---------------------------------------------------------------------------

type dEdge struct {
	from, to int
	rev      int // 反向半边下标
}

type segHit struct {
	x, y     float64
	t1, t2   float64
	edgeA    int
	edgeB    int
	endpoint bool // 交点恰为某条边端点
}

// segmentIntersect 求两条线段的交点。interior 为 true 表示交点落在两条
// 线段内部（规范相交，t,u 均严格介于 0,1 之间）；endpoint=true 表示交点
// 只在某条线段的端点上；collinear 表示两线段共线（含重叠）。
func segmentIntersect(p1, p2, q1, q2 Point, eps float64) (hit segHit, interior, endpoint, collinear bool) {
	rx, ry := p2.X-p1.X, p2.Y-p1.Y
	sx, sy := q2.X-q1.X, q2.Y-q1.Y
	lr := math.Hypot(rx, ry)
	ls := math.Hypot(sx, sy)
	rcs := rx*sy - ry*sx
	qpX, qpY := q1.X-p1.X, q1.Y-p1.Y
	// 共线判据与 orient 一致：交点到两直线的垂直距离均不超过 eps，
	// 即 |rcs|/|r| 与 |rcs|/|s| 都 <= eps，取较短边做保守阈值。
	lmin := math.Max(1, math.Min(lr, ls))
	if math.Abs(rcs) <= eps*lmin {
		if orient(p1, p2, q1, eps) == 0 && orient(p1, p2, q2, eps) == 0 {
			// 同一直线：仅当投影区间实质重叠才算共线重叠，
			// 否则（如同一直线上互不接触的两段）按不相交处理
			if overlap1D(p1, p2, q1, q2, eps) {
				return segHit{}, false, false, true
			}
		}
		return segHit{}, false, false, false
	}
	t := (qpX*sy - qpY*sx) / rcs
	u := (qpX*ry - qpY*rx) / rcs
	if t < -epsRel || t > 1+epsRel || u < -epsRel || u > 1+epsRel {
		return segHit{}, false, false, false
	}
	t = math.Max(0, math.Min(1, t))
	u = math.Max(0, math.Min(1, u))
	hit = segHit{x: p1.X + t*rx, y: p1.Y + t*ry, t1: t, t2: u}
	// 参数偏离 0/1 的程度必须超过按线段长度缩放的容差，才算真正进入内部
	inSeg := func(v, len float64) bool {
		return v > eps*math.Max(1, len) && v < 1-eps*math.Max(1, len)
	}
	inT, inU := inSeg(t, lr), inSeg(u, ls)
	interior = inT && inU
	endpoint = !interior
	return hit, interior, endpoint, false
}

// splitSelfIntersecting 把自交环拆成若干简单面，按奇偶规则保留填充面。
// 失败时返回带顶点号的 VertexError。
func splitSelfIntersecting(ring []Point) (faces [][]Point, err error) {
	n := len(ring)
	eps := epsFor(ring)

	type node struct {
		pt Point
	}
	var nodes []node
	nodeAt := make(map[[2]int64]int)
	coordKey := func(p Point) [2]int64 {
		s := math.Max(1, 1/eps)
		return [2]int64{int64(math.Round(p.X * s)), int64(math.Round(p.Y * s))}
	}
	getNode := func(p Point) int {
		k := coordKey(p)
		if id, ok := nodeAt[k]; ok {
			return id
		}
		id := len(nodes)
		nodeAt[k] = id
		nodes = append(nodes, node{pt: p})
		return id
	}

	// 1) 收集所有边的交点
	type edgeSplit struct {
		t  float64
		pt Point
	}
	splits := make([][]edgeSplit, n)
	for i := 0; i < n; i++ {
		splits[i] = []edgeSplit{{t: 0, pt: ring[i]}, {t: 1, pt: ring[(i+1)%n]}}
	}
	var badVerts []int
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			// 相邻边只在公共端点处相接，跳过
			if j == i+1 || (i == 0 && j == n-1) {
				continue
			}
			a, b := ring[i], ring[(i+1)%n]
			c, d := ring[j], ring[(j+1)%n]
			hit, isInterior, isEndpoint, col := segmentIntersect(a, b, c, d, eps)
			if col {
				badVerts = append(badVerts, i, (i+1)%n, j, (j+1)%n)
				continue
			}
			if !isInterior && !isEndpoint {
				continue
			}
			if isEndpoint {
				// 顶点落在非相邻边内部（T 接）：仍然需要在该边处切开
			}
			splits[i] = append(splits[i], edgeSplit{t: hit.t1, pt: Point{hit.x, hit.y}})
			splits[j] = append(splits[j], edgeSplit{t: hit.t2, pt: Point{hit.x, hit.y}})
		}
	}
	if len(badVerts) > 0 {
		return nil, errSelfIntersect(uniqSorted(badVerts), "存在共线重叠边，自交无法拆分")
	}

	// 2) 组装半边：每个平面段同时加入正反两条半边
	var edges []dEdge
	outHalf := make(map[int][]int) // node -> 以该点为起点的半边
	addEdge := func(u, v int) {
		if u == v {
			return
		}
		hf, hr := len(edges), len(edges)+1
		edges = append(edges,
			dEdge{from: u, to: v, rev: hr},
			dEdge{from: v, to: u, rev: hf},
		)
		outHalf[u] = append(outHalf[u], hf)
		outHalf[v] = append(outHalf[v], hr)
	}
	for i := 0; i < n; i++ {
		ss := splits[i]
		sort.Slice(ss, func(a, b int) bool { return ss[a].t < ss[b].t })
		var ids []int
		for _, sp := range ss {
			ids = append(ids, getNode(sp.pt))
		}
		// 去重相邻同点
		for k := 0; k+1 < len(ids); k++ {
			if ids[k] != ids[k+1] {
				addEdge(ids[k], ids[k+1])
			}
		}
	}
	// 自洽性检查：每条半边的反向边必须存在
	for h := range edges {
		r := edges[h].rev
		if r >= len(edges) || edges[r].from != edges[h].to || edges[r].to != edges[h].from {
			return nil, errSelfIntersect(nil, "自交拆分得到非闭合平面边")
		}
	}

	// 3) 每个节点的出边按极角排序（CCW）
	ang := func(u, v int) float64 {
		a, b := nodes[u].pt, nodes[v].pt
		return math.Atan2(b.Y-a.Y, b.X-a.X)
	}
	for u, hs := range outHalf {
		sort.Slice(hs, func(a, b int) bool {
			return ang(u, edges[hs[a]].to) < ang(u, edges[hs[b]].to)
		})
		outHalf[u] = hs
	}
	// nextInFace: 到达 v（半边 h 的终点）后，取 v 处反向边在 CCW 序中的前一条
	nextInFace := make([]int, len(edges))
	for _, hs := range outHalf {
		for _, h := range hs {
			r := edges[h].rev
			// r 是以 v 为起点的出边之一
			v := edges[h].to
			list := outHalf[v]
			pos := indexOfInt(list, r)
			if pos < 0 {
				return nil, errSelfIntersect(nil, "自交拆半边遍历失败")
			}
			nextInFace[h] = list[(pos-1+len(list))%len(list)]
		}
	}

	// 4) 遍历所有面
	visited := make([]bool, len(edges))
	var rawFaces [][]int
	for h := 0; h < len(edges); h++ {
		if visited[h] {
			continue
		}
		var face []int
		cur := h
		for !visited[cur] {
			visited[cur] = true
			face = append(face, cur)
			cur = nextInFace[cur]
		}
		rawFaces = append(rawFaces, face)
	}

	// 5) 面定向与奇偶判定
	for _, face := range rawFaces {
		var s float64
		for _, h := range face {
			a, b := nodes[edges[h].from].pt, nodes[edges[h].to].pt
			s += a.X*b.Y - b.X*a.Y
		}
		if s <= 0 {
			continue // 外面（CW）或退化面
		}
		loop := make([]Point, 0, len(face))
		for _, h := range face {
			loop = append(loop, nodes[edges[h].from].pt)
		}
		// 取面内一点做奇偶测试
		c := centroid(loop)
		if parityInside(ring, c, eps) {
			faces = append(faces, loop)
		}
	}
	if len(faces) == 0 {
		return nil, errSelfIntersect(allIndices(n), "自交拆分后没有可填充面")
	}
	return faces, nil
}

func centroid(loop []Point) Point {
	var c Point
	for _, p := range loop {
		c.X += p.X
		c.Y += p.Y
	}
	n := float64(len(loop))
	c.X /= n
	c.Y /= n
	return c
}

// parityInside 对原始（可能自交的）环做射线奇偶判定
func parityInside(ring []Point, p Point, eps float64) bool {
	inside := false
	n := len(ring)
	for i := 0; i < n; i++ {
		a, b := ring[i], ring[(i+1)%n]
		if (a.Y > p.Y) == (b.Y > p.Y) {
			continue
		}
		t := (p.Y - a.Y) / (b.Y - a.Y)
		x := a.X + (b.X-a.X)*t
		if x > p.X {
			inside = !inside
		}
	}
	return inside
}

func allIndices(n int) []int {
	out := make([]int, n)
	for i := range out {
		out[i] = i
	}
	return out
}

func uniqSorted(in []int) []int {
	seen := make(map[int]bool)
	var out []int
	for _, v := range in {
		if v >= 0 && !seen[v] {
			seen[v] = true
			out = append(out, v)
		}
	}
	sort.Ints(out)
	return out
}

func indexOfInt(list []int, v int) int {
	for i, x := range list {
		if x == v {
			return i
		}
	}
	return -1
}

// hasProperIntersection 检测环是否存在非相邻边的规范相交（交点落在两边
// 内部）或共线实质重叠。端点接触（含 T 接、桥缝共享顶点）默认也算自交，
// 因为输入多边形出现这种情况需要走自交拆分；但带洞桥接产生的合并环合法
// 地共享桥接点，此时传 ignoreEndpoint=true 只认规范相交与共线重叠。
func hasProperIntersection(ring []Point) (bool, []int) {
	return intersectCheck(ring, false)
}

func intersectCheck(ring []Point, ignoreEndpoint bool) (bool, []int) {
	n := len(ring)
	eps := epsFor(ring)
	var verts []int
	proper := false
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if j == i+1 || (i == 0 && j == n-1) {
				continue
			}
			a, b := ring[i], ring[(i+1)%n]
			c, d := ring[j], ring[(j+1)%n]
			// 桥接零宽缝：两条边端点完全相同、方向相反（a==d 且 b==c），
			// 这是带洞合并环的合法结构，不算自交
			if ptsEq(a, d, eps) && ptsEq(b, c, eps) {
				continue
			}
			_, isInterior, isEndpoint, col := segmentIntersect(a, b, c, d, eps)
			if col {
				// 共线分两种：重叠才是自交；仅延长线共线、无重叠则忽略
				if overlap1D(a, b, c, d, eps) {
					proper = true
					verts = append(verts, i, (i+1)%n, j, (j+1)%n)
				}
				continue
			}
			if isInterior {
				proper = true
				verts = append(verts, i, (i+1)%n, j, (j+1)%n)
				continue
			}
			if isEndpoint && !ignoreEndpoint {
				// 非相邻边仅端点相接：重复顶点接触（如蝴蝶结尖点）算自交，
				// 需拆分；T 接（边内部点落在另一边上）同理。
				proper = true
				if ptsEq(a, c, eps) || ptsEq(a, d, eps) || ptsEq(b, c, eps) || ptsEq(b, d, eps) {
					verts = append(verts, i, j)
				} else {
					verts = append(verts, i, (i+1)%n, j, (j+1)%n)
				}
			}
		}
	}
	return proper, uniqSorted(verts)
}

func ptsEq(a, b Point, eps float64) bool {
	return math.Abs(a.X-b.X) <= eps && math.Abs(a.Y-b.Y) <= eps
}

// ---------------------------------------------------------------------------
// 带洞多边形：桥接合并
// ---------------------------------------------------------------------------

// mergeHoles 将 outer（CCW）与若干洞（内部统一为 CW）通过互相可见的
// 桥接点逐洞合并为一条简单环。每个候选桥都实际拼接并验证（无自交且
// 有向面积守恒），选取第一个通过验证的桥。桥接点对的顶点号映射回
// 输入的 outer/hole 切片。
func mergeHoles(outer []Point, holes [][]Point) ([]Point, [][2]int, error) {
	cur := ensureCCW(outer)
	eps := epsFor(cur)
	var bridges [][2]int

	for _, holeIn := range holes {
		hole := holeIn
		if signedArea2(hole) > 0 {
			rev := make([]Point, len(hole))
			for i := range hole {
				rev[i] = hole[len(hole)-1-i]
			}
			hole = rev
		}
		// 合并前的目标有向面积（桥为零宽缝，不改变有向面积）
		wantArea := signedArea2(cur) + signedArea2(hole)

		merged, ci, hj, ok := spliceWithValidBridge(cur, hole, wantArea, eps)
		if !ok {
			return nil, nil, errSelfIntersect(nil, "洞与外环找不到可见桥接点")
		}
		// 映射回原始 outer 的顶点号（按坐标查找）
		bridges = append(bridges, [2]int{indexOfPoint(outer, cur[ci], eps), hj})
		cur = merged
	}
	return cur, bridges, nil
}

// spliceWithValidBridge 枚举 (外边界点, 洞点) 候选，按距离升序逐个拼接，
// 返回第一个“无自交且有向面积守恒”的合并环及所用桥接点下标。
func spliceWithValidBridge(cur, hole []Point, wantArea, eps float64) ([]Point, int, int, bool) {
	type cand struct {
		ci, hj int
		d      float64
	}
	var cands []cand
	for hj, hp := range hole {
		for ci, cp := range cur {
			cands = append(cands, cand{ci, hj, math.Hypot(cp.X-hp.X, cp.Y-hp.Y)})
		}
	}
	sort.Slice(cands, func(a, b int) bool { return cands[a].d < cands[b].d })

	for _, c := range cands {
		merged := spliceRing(cur, hole, c.ci, c.hj)
		// 面积守恒（桥是零宽缝，面积应严格等于外环+洞的有向面积）
		if math.Abs(signedArea2(merged)-wantArea) > eps*math.Max(1, math.Abs(wantArea))*1e3 {
			continue
		}
		// 合并环不得自交（桥缝共享顶点属合法结构，忽略端点接触）
		if inter, _ := intersectCheck(merged, true); inter {
			continue
		}
		return merged, c.ci, c.hj, true
	}
	return nil, 0, 0, false
}

// spliceRing 按桥接点 (cur[ci], hole[hj]) 拼接。hole 为顺时针（CW）洞，
// 其数组相邻顺序即 CW 走向（每个顶点处右转）。从外环桥接点 q=cur[ci]
// 经桥进入洞点 p=hole[hj]，沿数组顺序（CW）绕洞一周回到 p，再经桥
// 返回 q，继续外环。两条桥段 q→p 与 p→q 方向相反、共线，面积相互
// 抵消，故合并环有向面积 = 外环 + 洞：
//
//	cur[..ci] + p + hole[hj+1..] + hole[..hj](=p) + q + cur[ci+1..]
//
// 桥缝两端的 p、q 各出现两次，属于带洞合并环合法的弱简单结构。
func spliceRing(cur, hole []Point, ci, hj int) []Point {
	hn := len(hole)
	var out []Point
	out = append(out, cur[:ci+1]...) // 外环走到 q（含 q）
	out = append(out, hole[hj])      // 桥 q → p
	for k := 1; k <= hn; k++ {       // CW 绕洞：hj+1, ..., 回到 hj(=p)
		out = append(out, hole[(hj+k)%hn])
	}
	out = append(out, cur[ci])       // 桥 p → q
	out = append(out, cur[ci+1:]...) // 继续外环
	return out
}

// overlap1D 判定两条已确认共线的线段是否有实质性重叠（仅端点接触不算）。
// 与线段方向无关：把四个端点都投影到 a→b 方向轴上，比较两条线段各自
// 投影区间的交叠长度。
func overlap1D(a, b, c, d Point, eps float64) bool {
	dx, dy := b.X-a.X, b.Y-a.Y
	lab := math.Hypot(dx, dy)
	if lab <= eps {
		return false
	}
	ux, uy := dx/lab, dy/lab
	proj := func(p Point) float64 { return (p.X-a.X)*ux + (p.Y-a.Y)*uy }
	p1, p2 := proj(a), proj(b)
	p3, p4 := proj(c), proj(d)
	lo := math.Max(math.Min(p1, p2), math.Min(p3, p4))
	hi := math.Min(math.Max(p1, p2), math.Max(p3, p4))
	return hi-lo > eps
}

// ---------------------------------------------------------------------------
// 校验：面积逐块核对、重叠逐对核对
// ---------------------------------------------------------------------------

// verifyTriangulation 逐三角形核对面积求和，并逐对核对内部不重叠。
func verifyTriangulation(tris []Triangle, originalAbsArea float64) (total float64, overlapPairs int, diff, rel float64) {
	for _, t := range tris {
		total += absAreaTri(t[0], t[1], t[2])
	}
	for i := 0; i < len(tris); i++ {
		for j := i + 1; j < len(tris); j++ {
			if trianglesOverlap(tris[i], tris[j]) {
				overlapPairs++
			}
		}
	}
	diff = math.Abs(total - originalAbsArea)
	if originalAbsArea > 0 {
		rel = diff / originalAbsArea
	}
	return total, overlapPairs, diff, rel
}

// trianglesOverlap 判定两个三角形是否存在二维内部重叠（仅共享边/顶点不算）。
// 做法：取所有顶点，检查“一个三角形的顶点严格落在另一个内部”，以及
// 边与边的规范相交。
func trianglesOverlap(t1, t2 Triangle) bool {
	eps := epsFor(append(t1[:], t2[:]...))
	a := t1[:]
	b := t2[:]
	for _, p := range a {
		if pointInTriContain(p, b[0], b[1], b[2], eps) {
			return true
		}
	}
	for _, p := range b {
		if pointInTriContain(p, a[0], a[1], a[2], eps) {
			return true
		}
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			_, isInterior, _, col := segmentIntersect(a[i], a[(i+1)%3], b[j], b[(j+1)%3], eps)
			if col {
				// 共线边（含共享边、同直线上的重合段）只有一维交叠，
				// 不产生二维面积重叠，忽略
				continue
			}
			if isInterior {
				return true
			}
		}
	}
	return false
}

// ringAbsArea 取有向面积绝对值
func ringAbsArea(ring []Point) float64 {
	return math.Abs(signedArea2(ring)) / 2
}
