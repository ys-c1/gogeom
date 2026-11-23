package geom

import (
	"errors"
	"math"
)

// 输入与分解过程中的错误。带顶点号的错误统一用 IndexError 包装。
var (
	ErrDuplicateVertex   = errors.New("顶点重复")
	ErrSelfIntersect     = errors.New("多边形自交")
	ErrDegeneratePolygon = errors.New("退化多边形")
	ErrDiagonalOutside   = errors.New("对角线穿出多边形")
	ErrDiagonalCross     = errors.New("对角线相交")
)

// IndexError 携带出错的顶点号（沿用输入多边形的 0 基下标）。
type IndexError struct {
	Op    string // duplicate | self-intersect | degenerate | diagonal-cross | diagonal-outside
	Index []int  // 相关顶点号
	Err   error  // 底层错误
}

func (e *IndexError) Error() string {
	return e.Op + " 顶点号 " + joinInts(e.Index) + ": " + e.Err.Error()
}

func (e *IndexError) Unwrap() error { return e.Err }

func joinInts(xs []int) string {
	if len(xs) == 0 {
		return "[]"
	}
	s := "[" + itoa(xs[0])
	for _, x := range xs[1:] {
		s += " " + itoa(x)
	}
	return s + "]"
}

func itoa(x int) string {
	if x == 0 {
		return "0"
	}
	neg := x < 0
	if neg {
		x = -x
	}
	var b [24]byte
	i := len(b)
	for x > 0 {
		i--
		b[i] = byte('0' + x%10)
		x /= 10
	}
	if neg {
		i--
		b[i] = '-'
	}
	return string(b[i:])
}

// epsScale 给出与坐标尺度匹配的容差，保证在小坐标和大坐标下都有效。
func epsScale(a, b Point) float64 {
	m := math.Max(math.Abs(a.X), math.Abs(a.Y))
	m = math.Max(m, math.Max(math.Abs(b.X), math.Abs(b.Y)))
	return 1e-12 * math.Max(1, m)
}

// cross2 计算 (b-a)×(c-a)。
func cross2(a, b, c Point) float64 {
	return (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
}

func ptEq(a, b Point, eps float64) bool {
	return math.Abs(a.X-b.X) <= eps && math.Abs(a.Y-b.Y) <= eps
}

// turn 返回 a->b->c 的转向：>0 左转（CCW），<0 右转，0 共线。
func turn(a, b, c Point) float64 {
	return cross2(a, b, c)
}

// properIntersect 判断两条线段是否在内部真正相交（不含端点接触）。
func properIntersect(a, b, c, d Point) bool {
	p := cross2(a, b, c)
	q := cross2(a, b, d)
	r := cross2(c, d, a)
	s := cross2(c, d, b)
	return ((p > 0 && q < 0) || (p < 0 && q > 0)) &&
		((r > 0 && s < 0) || (r < 0 && s > 0))
}

// onSeg 判断 p 是否落在闭线段 ab 上（含端点）。
func onSeg(a, b, p Point) bool {
	if math.Abs(cross2(a, b, p)) > epsScale(a, b) {
		return false
	}
	return math.Min(a.X, b.X)-epsScale(a, b) <= p.X &&
		p.X <= math.Max(a.X, b.X)+epsScale(a, b) &&
		math.Min(a.Y, b.Y)-epsScale(a, b) <= p.Y &&
		p.Y <= math.Max(a.Y, b.Y)+epsScale(a, b)
}

// ValidatePolygon 校验简单多边形：顶点数、重复顶点、自交（含顶点压边）。
// 出错时返回带顶点号的 IndexError。
func ValidatePolygon(poly []Point) error {
	n := len(poly)
	if n < 3 {
		return ErrTooFewPts
	}
	// 顶点重复（按坐标）。
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if ptEq(poly[i], poly[j], epsScale(poly[i], poly[j])) {
				return &IndexError{Op: "duplicate", Index: []int{i, j}, Err: ErrDuplicateVertex}
			}
		}
	}
	// 非相邻边求交；相邻边仅在端点相接，端点重复已在上面拒绝。
	// 自交检查先于零面积检查：蝴蝶结等自交多边形的有符号面积恰好为零。
	for i := 0; i < n; i++ {
		a, b := poly[i], poly[(i+1)%n]
		for j := i + 1; j < n; j++ {
			c, d := poly[j], poly[(j+1)%n]
			adjacent := j == i+1 || (i == 0 && j == n-1)
			if adjacent {
				continue
			}
			if properIntersect(a, b, c, d) {
				return &IndexError{Op: "self-intersect", Index: []int{i, (i + 1) % n, j, (j + 1) % n}, Err: ErrSelfIntersect}
			}
			// 顶点压在非关联边上（T 型接触，破坏简单多边形）。
			if onSeg(c, d, a) {
				return &IndexError{Op: "self-intersect", Index: []int{i, j, (j + 1) % n}, Err: ErrSelfIntersect}
			}
			if onSeg(c, d, b) {
				return &IndexError{Op: "self-intersect", Index: []int{(i + 1) % n, j, (j + 1) % n}, Err: ErrSelfIntersect}
			}
			if onSeg(a, b, c) {
				return &IndexError{Op: "self-intersect", Index: []int{j, i, (i + 1) % n}, Err: ErrSelfIntersect}
			}
			if onSeg(a, b, d) {
				return &IndexError{Op: "self-intersect", Index: []int{(j + 1) % n, i, (i + 1) % n}, Err: ErrSelfIntersect}
			}
		}
	}
	// 退化：共线点或零长边已被重复顶点覆盖一部分，这里再查连续三点共线导致的零面积。
	area2 := 0.0
	maxCoord := 0.0
	for i := 0; i < n; i++ {
		area2 += poly[i].X*poly[(i+1)%n].Y - poly[(i+1)%n].X*poly[i].Y
		maxCoord = math.Max(maxCoord, math.Max(math.Abs(poly[i].X), math.Abs(poly[i].Y)))
	}
	if math.Abs(area2) <= 1e-12*maxCoord*maxCoord {
		return &IndexError{Op: "degenerate", Index: []int{0, 1, 2}, Err: ErrDegeneratePolygon}
	}
	return nil
}

// polygonArea2 返回有符号二倍面积（正为 CCW）。
func polygonArea2(poly []Point) float64 {
	var s float64
	for i := range poly {
		a, b := poly[i], poly[(i+1)%len(poly)]
		s += a.X*b.Y - b.X*a.Y
	}
	return s
}

// reflexIndices 返回反射（凹）顶点在当前 poly 中的下标。poly 须为已校验的简单多边形。
func reflexIndices(poly []Point) []int {
	n := len(poly)
	ccw := polygonArea2(poly) > 0
	var out []int
	for i := 0; i < n; i++ {
		a, b, c := poly[(i-1+n)%n], poly[i], poly[(i+1)%n]
		t := turn(a, b, c)
		if ccw && t < -epsScale(a, b) {
			out = append(out, i)
		} else if !ccw && t > epsScale(a, b) {
			out = append(out, i)
		}
	}
	return out
}

// ReflexVertices 校验后返回所有反射顶点的原始下标。
func ReflexVertices(poly []Point) ([]int, error) {
	if err := ValidatePolygon(poly); err != nil {
		return nil, err
	}
	return reflexIndices(poly), nil
}

// IsConvex 判断已校验多边形是否为凸（允许共线边）。
func IsConvex(poly []Point) bool {
	if len(poly) < 3 {
		return false
	}
	sign := 0.0
	for i := range poly {
		t := turn(poly[i], poly[(i+1)%len(poly)], poly[(i+2)%len(poly)])
		if math.Abs(t) <= epsScale(poly[i], poly[(i+1)%len(poly)]) {
			continue
		}
		if sign == 0 {
			sign = t
		} else if (sign > 0) != (t > 0) {
			return false
		}
	}
	return true
}

// pointOnBoundary 判断 p 是否落在多边形边界上。
func pointOnBoundary(poly []Point, p Point) bool {
	for i := range poly {
		if onSeg(poly[i], poly[(i+1)%len(poly)], p) {
			return true
		}
	}
	return false
}

// internalDiagonal 判断 poly[i] 到 poly[j] 的线段是否为合法内部对角线：
// 不穿边界、不被其他顶点压、中点严格位于多边形内部。
func internalDiagonal(poly []Point, i, j int) bool {
	n := len(poly)
	if i == j {
		return false
	}
	a, b := poly[i], poly[j]
	// 不能被任何非端点顶点压住（避免对角线沿边穿点）。
	for k := 0; k < n; k++ {
		if k == i || k == j {
			continue
		}
		if onSeg(a, b, poly[k]) {
			return false
		}
	}
	// 不得与任何非相邻边真正相交。
	for k := 0; k < n; k++ {
		e1, e2 := poly[k], poly[(k+1)%n]
		if k == i || k == j || (k+1)%n == i || (k+1)%n == j {
			continue
		}
		if properIntersect(a, b, e1, e2) {
			return false
		}
	}
	// 中点必须严格在多边形内部（在边界上即视为穿出）。
	mid := Point{(a.X + b.X) / 2, (a.Y + b.Y) / 2}
	if pointOnBoundary(poly, mid) {
		return false
	}
	if !Inside(poly, mid.X, mid.Y) {
		return false
	}
	return true
}

// ConvexDecompose 对简单多边形做凸分解。
//
// 流程：先用耳切法三角化（对角线均落在多边形内部），再沿共享对角线贪心
// 合并相邻凸块——只要合并后的多边形仍为凸就合并。最终每块都是凸的，
// 所用对角线全部是三角化对角线的子集，因此两两不相交且不穿出边界。
//
// 校验失败（顶点重复、自交）或分解失败时返回带顶点号的错误。
func ConvexDecompose(poly []Point) ([][]Point, error) {
	if err := ValidatePolygon(poly); err != nil {
		return nil, err
	}
	tris, err := triangulate(poly)
	if err != nil {
		return nil, err
	}
	// 独立复核：三角化用到的每条对角线都必须落在多边形内部且两两不相交。
	triDiags := triangulationDiagonals(poly, tris)
	if err := checkDiagList(poly, triDiags); err != nil {
		return nil, err
	}
	for _, d := range triDiags {
		if !internalDiagonal(poly, d[0], d[1]) {
			return nil, &IndexError{Op: "diagonal-outside", Index: d[:], Err: ErrDiagonalOutside}
		}
	}
	idxPieces := mergePieces(poly, tris)
	// 分块之间的共享边即实际使用的对角线，做一次两两不相交复核。
	diags := sharedDiagonals(idxBoundaryEdges(poly), idxPieces)
	if err := checkDiagList(poly, diags); err != nil {
		return nil, err
	}
	pieces := make([][]Point, len(idxPieces))
	for i, idx := range idxPieces {
		p := make([]Point, len(idx))
		for k, vi := range idx {
			p[k] = poly[vi]
		}
		pieces[i] = p
	}
	return pieces, nil
}

// triangulate 耳切法三角化，返回 CCW 的原多边形下标三元组序列。
func triangulate(poly []Point) ([][3]int, error) {
	n := len(poly)
	ring := make([]int, n)
	for i := range ring {
		ring[i] = i
	}
	// 统一成 CCW，耳判定只处理左转。
	if polygonArea2(poly) < 0 {
		for i, j := 0, n-1; i < j; i, j = i+1, j-1 {
			ring[i], ring[j] = ring[j], ring[i]
		}
	}
	var tris [][3]int
	for len(ring) > 3 {
		found := -1
		for k := range ring {
			if isEar(poly, ring, k) {
				found = k
				break
			}
		}
		if found < 0 {
			return nil, &IndexError{Op: "diagonal-outside", Index: append([]int(nil), ring...), Err: ErrDiagonalOutside}
		}
		m := len(ring)
		a := ring[(found-1+m)%m]
		b := ring[found]
		c := ring[(found+1)%m]
		tris = append(tris, [3]int{a, b, c})
		ring = append(ring[:found], ring[found+1:]...)
	}
	tris = append(tris, [3]int{ring[0], ring[1], ring[2]})
	return tris, nil
}

// isEar 判断 ring[k] 是否为耳：凸顶点、(prev,next) 是内部对角线、
// 耳三角形内不含其他顶点。
func isEar(poly []Point, ring []int, k int) bool {
	m := len(ring)
	a, b, c := ring[(k-1+m)%m], ring[k], ring[(k+1)%m]
	pa, pb, pc := poly[a], poly[b], poly[c]
	if cross2(pa, pb, pc) <= epsScale(pa, pb) {
		return false // 共线或反射顶点不可能是耳
	}
	if !internalDiagonal(poly, a, c) {
		return false
	}
	eps := epsScale(pa, pc)
	for _, vi := range ring {
		if vi == a || vi == b || vi == c {
			continue
		}
		p := poly[vi]
		// CCW 三角形的严格内部（含边，共线点也算被耳包住，跳过该耳）。
		if cross2(pa, pb, p) >= -eps &&
			cross2(pb, pc, p) >= -eps &&
			cross2(pc, pa, p) >= -eps {
			return false
		}
	}
	return true
}

// mergePieces 从三角形集合出发，反复合并共享一条边且并集仍为凸的两块。
func mergePieces(poly []Point, tris [][3]int) [][]int {
	pieces := make([][]int, 0, len(tris))
	for _, t := range tris {
		pieces = append(pieces, []int{t[0], t[1], t[2]})
	}
	for {
		merged := false
		for i := 0; i < len(pieces) && !merged; i++ {
			for j := i + 1; j < len(pieces); j++ {
				u, v, ok := sharedEdgeIdx(pieces[i], pieces[j])
				if !ok {
					continue
				}
				union := mergeTwo(pieces[i], pieces[j], u, v)
				if isConvexIdx(poly, union) {
					pieces[i] = union
					pieces = append(pieces[:j], pieces[j+1:]...)
					merged = true
					break
				}
			}
		}
		if !merged {
			break
		}
	}
	return pieces
}

// sharedEdgeIdx 返回两块共享边的两个原多边形下标（无共享边时 ok=false）。
func sharedEdgeIdx(a, b []int) (u, v int, ok bool) {
	for i := 0; i < len(a); i++ {
		x, y := a[i], a[(i+1)%len(a)]
		for j := 0; j < len(b); j++ {
			p, q := b[j], b[(j+1)%len(b)]
			if (x == p && y == q) || (x == q && y == p) {
				return x, y, true
			}
		}
	}
	return 0, 0, false
}

// mergeTwo 沿共享边 (u,v) 拼接两个 CCW 凸多边形，返回并集的 CCW 下标环。
func mergeTwo(a, b []int, u, v int) []int {
	// 找到 a 中 u->v 的方向，并集 = a 上 v 到 u 的链 + b 上 u 到 v 的链。
	au := indexOfInt(a, u)
	if a[(au+1)%len(a)] != v {
		a, b = b, a // 保证 a 含 u->v
		au = indexOfInt(a, u)
	}
	var out []int
	k := (au + 1) % len(a) // 从 v 出发
	for {
		out = append(out, a[k])
		if a[k] == u {
			break
		}
		k = (k + 1) % len(a)
	}
	bu := indexOfInt(b, u)
	kb := bu
	for {
		kb = (kb + 1) % len(b)
		if b[kb] == v {
			break
		}
		out = append(out, b[kb])
	}
	return out
}

func indexOfInt(xs []int, x int) int {
	for i, v := range xs {
		if v == x {
			return i
		}
	}
	return -1
}

// isConvexIdx 用原多边形坐标判断下标环是否为凸（允许共线）。
func isConvexIdx(poly []Point, idx []int) bool {
	n := len(idx)
	sign := 0.0
	for i := 0; i < n; i++ {
		t := cross2(poly[idx[i]], poly[idx[(i+1)%n]], poly[idx[(i+2)%n]])
		if math.Abs(t) <= epsScale(poly[idx[i]], poly[idx[(i+1)%n]]) {
			continue
		}
		if sign == 0 {
			sign = t
		} else if (sign > 0) != (t > 0) {
			return false
		}
	}
	return true
}

// triangulationDiagonals 提取三角形集合中不属于原多边形边界的边（即对角线）。
func triangulationDiagonals(poly []Point, tris [][3]int) [][2]int {
	boundary := idxBoundaryEdges(poly)
	count := map[[2]int]int{}
	var keys [][2]int
	note := func(u, v int) {
		if u > v {
			u, v = v, u
		}
		e := [2]int{u, v}
		if _, seen := count[e]; !seen {
			keys = append(keys, e)
		}
		count[e]++
	}
	for _, t := range tris {
		note(t[0], t[1])
		note(t[1], t[2])
		note(t[2], t[0])
	}
	var out [][2]int
	for _, e := range keys {
		if count[e] == 1 && !boundary[e] {
			out = append(out, e)
		}
	}
	return out
}

// idxBoundaryEdges 给出原多边形边界边的规范键集合。
func idxBoundaryEdges(poly []Point) map[[2]int]bool {
	set := map[[2]int]bool{}
	for i := range poly {
		u, v := i, (i+1)%len(poly)
		if u > v {
			u, v = v, u
		}
		set[[2]int{u, v}] = true
	}
	return set
}

// sharedDiagonals 返回分块之间被两块共享、且不属于原边界的边（即对角线）。
func sharedDiagonals(boundary map[[2]int]bool, pieces [][]int) [][2]int {
	count := map[[2]int]int{}
	for _, p := range pieces {
		for i := range p {
			u, v := p[i], p[(i+1)%len(p)]
			if u > v {
				u, v = v, u
			}
			count[[2]int{u, v}]++
		}
	}
	var out [][2]int
	for e, c := range count {
		if c == 2 && !boundary[e] {
			out = append(out, [2]int{e[0], e[1]})
		}
	}
	return out
}

// checkDiagList 检查对角线两两是否真正相交（端点相接合法）。
func checkDiagList(poly []Point, diags [][2]int) error {
	for i := 0; i < len(diags); i++ {
		for j := i + 1; j < len(diags); j++ {
			a, b := poly[diags[i][0]], poly[diags[i][1]]
			c, d := poly[diags[j][0]], poly[diags[j][1]]
			if properIntersect(a, b, c, d) {
				return &IndexError{
					Op:    "diagonal-cross",
					Index: []int{diags[i][0], diags[i][1], diags[j][0], diags[j][1]},
					Err:   ErrDiagonalCross,
				}
			}
		}
	}
	return nil
}

// CheckDiagonals 校验外部给出的对角线集合：每条必须是合法内部对角线，
// 且两两不得相交（共享端点除外）。违反时返回带四个顶点号的 IndexError。
func CheckDiagonals(poly []Point, diags [][2]int) error {
	if err := ValidatePolygon(poly); err != nil {
		return err
	}
	n := len(poly)
	for _, d := range diags {
		if d[0] < 0 || d[0] >= n || d[1] < 0 || d[1] >= n || d[0] == d[1] {
			return &IndexError{Op: "diagonal-outside", Index: []int{d[0], d[1]}, Err: ErrDiagonalOutside}
		}
		if !internalDiagonal(poly, d[0], d[1]) {
			return &IndexError{Op: "diagonal-outside", Index: []int{d[0], d[1]}, Err: ErrDiagonalOutside}
		}
	}
	return checkDiagList(poly, diags)
}

// CheckAreaConservation 分开验证守恒条件一：各块面积之和等于原多边形面积。
func CheckAreaConservation(original []Point, pieces [][]Point) (origArea, sumArea float64, err error) {
	o, err := Area(original)
	if err != nil {
		return 0, 0, err
	}
	origArea = math.Abs(o)
	for _, p := range pieces {
		a, err := Area(p)
		if err != nil {
			return origArea, 0, err
		}
		sumArea += math.Abs(a)
	}
	tol := 1e-9 * math.Max(1, origArea)
	if math.Abs(origArea-sumArea) > tol {
		return origArea, sumArea, &IndexError{
			Op:  "area-mismatch",
			Err: errors.New("各块面积之和不等于原面积"),
		}
	}
	return origArea, sumArea, nil
}

// pointStrictInside 判断 p 是否严格位于多边形内部（边界上不算内部）。
func pointStrictInside(poly []Point, p Point) bool {
	if pointOnBoundary(poly, p) {
		return false
	}
	return Inside(poly, p.X, p.Y)
}

// CheckNoOverlap 分开验证守恒条件二：块与块之间不得重叠。
// 判定：任一顶点不得严格落在另一块内部；任两对边不得真正相交。
// 共享边、共线贴合等仅边界接触的情况合法。
func CheckNoOverlap(pieces [][]Point) error {
	for i := 0; i < len(pieces); i++ {
		for j := i + 1; j < len(pieces); j++ {
			a, b := pieces[i], pieces[j]
			for _, v := range a {
				if pointStrictInside(b, v) {
					return &IndexError{Op: "overlap", Index: []int{i, j}, Err: errors.New("分块重叠")}
				}
			}
			for _, v := range b {
				if pointStrictInside(a, v) {
					return &IndexError{Op: "overlap", Index: []int{i, j}, Err: errors.New("分块重叠")}
				}
			}
			for x := 0; x < len(a); x++ {
				e1, e2 := a[x], a[(x+1)%len(a)]
				for y := 0; y < len(b); y++ {
					f1, f2 := b[y], b[(y+1)%len(b)]
					if properIntersect(e1, e2, f1, f2) {
						return &IndexError{Op: "overlap", Index: []int{i, j}, Err: errors.New("分块边相交")}
					}
				}
			}
		}
	}
	return nil
}

// CheckPiecesWithin 验证所有分块都位于原多边形内部（对角线未越界的独立复核）。
func CheckPiecesWithin(original []Point, pieces [][]Point) error {
	for pi, p := range pieces {
		for _, v := range p {
			if !Inside(original, v.X, v.Y) && !pointOnBoundary(original, v) {
				return &IndexError{Op: "diagonal-outside", Index: []int{pi}, Err: ErrDiagonalOutside}
			}
		}
	}
	return nil
}
