// cut.go 按面积比切分多边形：切线方向与位置均由方程解出，无任何二分/迭代试凑。
//
// 切分直线 n·x = d，其中 n 为单位法向，u 为切线方向（n 逆时针转 90°）。
//   - 方向 n：对多边形（含孔洞，孔洞以有向面积自然带负权）求面积加权的
//     二阶矩协方差矩阵，其特征向量由 2×2 对称矩阵的闭式公式直接解出；
//     取最大特征值方向（主轴）为法向，短轴为切线方向。
//   - 位置 d：各顶点在 n 上的投影把轴切成若干板；板内截线穿过多边形的
//     截面长度 L(d) 是 d 的一次函数，故累计面积 F(d) 是 d 的二次函数。
//     线性扫描定位目标板后，直接解二次方程得到 d，无二分。
package geom

import (
	"errors"
	"fmt"
	"math"
)

// 输入校验相关错误。
var (
	ErrRatioOutOfRange = errors.New("面积比必须落在 (0,1) 开区间")
	ErrSelfIntersect   = errors.New("多边形自交")
	ErrCollinear       = errors.New("存在连续三顶点共线")
	ErrInvalidHole     = errors.New("孔洞非法（越界或相互交叠）")
	ErrUncuttable      = errors.New("该方向下无法按比例切分")
	ErrAssembly        = errors.New("切块拼装失败")
)

const (
	epsAbs = 1e-10 // 绝对容差（与常规坐标量级配套）
	epsRel = 1e-9  // 相对容差
)

// Polygon 表示带洞多边形：Outer 为外轮廓环，Holes 为孔洞环。
// 各环顶点顺序任意，构造时会自动规整为外逆内顺。
type Polygon struct {
	Outer []Point
	Holes [][]Point
}

// CutLine 是切分直线 n·x = d，Dir 为直线方向（单位向量，Dir ⟂ N）。
type CutLine struct {
	N   Point   // 单位法向
	Dir Point   // 切线方向（单位向量，Dir = N 逆时针转 90°）
	D   float64 // 有向偏移：N·X = D
}

// String 输出切线的一般式方程。
func (c CutLine) String() string {
	return fmt.Sprintf("%.10g*x + %.10g*y = %.10g", c.N.X, c.N.Y, c.D)
}

// CutResult 是一次切分的全部产物。
type CutResult struct {
	Line    CutLine // 切线方程
	Neg     []Point // N·x <= D 一侧的简单多边形
	Pos     []Point // N·x >= D 一侧的简单多边形
	AreaNeg float64 // Neg 面积（正数）
	AreaPos float64 // Pos 面积（正数）
	AreaAll float64 // 原多边形面积（正数，孔洞已扣除）
	Ratio   float64 // 实际面积比 AreaNeg / AreaAll
	AbsErr  float64 // |实际比 - 目标比|
	// 两条独立验收结论（面积守恒、简单性分开验）：
	AreaSumOK bool // 两块面积之和 == 原面积
	PiecesOK  bool // 两块都是简单多边形
}

// EdgeError 定位到具体环与边的输入错误。边号从 0 起；Ring 0 为外轮廓，
// Ring>=1 为第 Ring 个孔洞。
type EdgeError struct {
	Ring int
	Edge int // 0 基边号，边 i 连接顶点 i 与顶点 (i+1)%n
	Kind error
}

func (e *EdgeError) Error() string {
	if e.Ring == 0 {
		return fmt.Sprintf("%s：边 %d（顶点 %d→%d）", e.Kind, e.Edge, e.Edge, (e.Edge + 1))
	}
	return fmt.Sprintf("%s：环 %d 边 %d（顶点 %d→%d）", e.Kind, e.Ring, e.Edge, e.Edge, (e.Edge + 1))
}

func (e *EdgeError) Unwrap() error { return e.Kind }

// ---------------------------------------------------------------------------
// 基础几何工具
// ---------------------------------------------------------------------------

func cross2(ax, ay, bx, by float64) float64 { return ax*by - ay*bx }

// orient 三顶点叉积 (b-a)×(c-a)。
func orient(a, b, c Point) float64 {
	return cross2(b.X-a.X, b.Y-a.Y, c.X-a.X, c.Y-a.Y)
}

// tol 按参与量的量级给出共线/重合判定的绝对阈值。
func tol(vals ...float64) float64 {
	m := epsAbs
	for _, v := range vals {
		if s := math.Abs(v) * epsRel; s > m {
			m = s
		}
	}
	return m
}

func ptsEq(a, b Point) bool {
	return math.Hypot(a.X-b.X, a.Y-b.Y) <= tol(a.X, a.Y, b.X, b.Y)
}

func reverse(p []Point) {
	for i, j := 0, len(p)-1; i < j; i, j = i+1, j-1 {
		p[i], p[j] = p[j], p[i]
	}
}

func signedArea(ring []Point) float64 {
	var s float64
	for i := range ring {
		a, b := ring[i], ring[(i+1)%len(ring)]
		s += a.X*b.Y - b.X*a.Y
	}
	return s / 2
}

// ---------------------------------------------------------------------------
// 输入校验：比值越界 / 自交 / 连续共线，均报出边号
// ---------------------------------------------------------------------------

// validateRing 校验单环：最少 3 点、零长边、连续共线、自交。边号 0 基。
func validateRing(ring []Point, ringIdx int) error {
	n := len(ring)
	if n < 3 {
		return &EdgeError{Ring: ringIdx, Edge: 0, Kind: ErrTooFewPts}
	}
	// 零长边（重复点）按退化共线处理，报该边号。
	for i := 0; i < n; i++ {
		a, b := ring[i], ring[(i+1)%n]
		if math.Hypot(b.X-a.X, b.Y-a.Y) <= tol(a.X, a.Y, b.X, b.Y) {
			return &EdgeError{Ring: ringIdx, Edge: i, Kind: ErrCollinear}
		}
	}
	// 连续三顶点共线：报中间顶点出射边的边号。
	for i := 0; i < n; i++ {
		a, b, c := ring[i], ring[(i+1)%n], ring[(i+2)%n]
		if math.Abs(orient(a, b, c)) <= tol(a.X, a.Y, b.X, b.Y, c.X, c.Y) {
			return &EdgeError{Ring: ringIdx, Edge: (i + 1) % n, Kind: ErrCollinear}
		}
	}
	// 自交：任意两条不相邻边不得规范相交，也不得在对方内部相触。
	for i := 0; i < n; i++ {
		a, b := ring[i], ring[(i+1)%n]
		for j := i + 2; j < n; j++ {
			if i == 0 && j == n-1 {
				continue // 首尾相邻
			}
			c, d := ring[j], ring[(j+1)%n]
			if properIntersect(a, b, c, d) ||
				onInterior(a, b, c) || onInterior(a, b, d) ||
				onInterior(c, d, a) || onInterior(c, d, b) {
				return &EdgeError{Ring: ringIdx, Edge: i, Kind: ErrSelfIntersect}
			}
		}
	}
	return nil
}

// properIntersect 两线段是否规范相交（端点落在对方内部由 onInterior 另判）。
func properIntersect(a, b, c, d Point) bool {
	o1, o2 := orient(a, b, c), orient(a, b, d)
	o3, o4 := orient(c, d, a), orient(c, d, b)
	t := tol(a.X, a.Y, b.X, b.Y, c.X, c.Y, d.X, d.Y)
	return ((o1 > t && o2 < -t) || (o1 < -t && o2 > t)) &&
		((o3 > t && o4 < -t) || (o3 < -t && o4 > t))
}

// onInterior 点 p 是否严格落在开线段 (a,b) 内部。
func onInterior(a, b, p Point) bool {
	t1 := tol(a.X, a.Y, b.X, b.Y, p.X, p.Y)
	if math.Abs(orient(a, b, p)) > t1 {
		return false
	}
	return (p.X > math.Min(a.X, b.X)+t1 || p.Y > math.Min(a.Y, b.Y)+t1) &&
		(p.X < math.Max(a.X, b.X)-t1 || p.Y < math.Max(a.Y, b.Y)-t1)
}

// validatePolygon 校验整多边形并把各环规整为外逆内顺，返回环集与净面积。
func validatePolygon(p *Polygon) ([][]Point, float64, error) {
	if err := validateRing(p.Outer, 0); err != nil {
		return nil, 0, err
	}
	rings := make([][]Point, 0, 1+len(p.Holes))
	outer := append([]Point(nil), p.Outer...)
	if signedArea(outer) < 0 {
		reverse(outer)
	}
	rings = append(rings, outer)

	total := math.Abs(signedArea(outer))
	for h, hole := range p.Holes {
		if err := validateRing(hole, h+1); err != nil {
			return nil, 0, err
		}
		hr := append([]Point(nil), hole...)
		if signedArea(hr) > 0 {
			reverse(hr)
		}
		if !Inside(outer, hr[0].X, hr[0].Y) {
			return nil, 0, &EdgeError{Ring: h + 1, Edge: 0, Kind: ErrInvalidHole}
		}
		rings = append(rings, hr)
		total -= math.Abs(signedArea(hr))
	}
	// 孔洞两两不得交叠。
	for i := 1; i < len(rings); i++ {
		for j := i + 1; j < len(rings); j++ {
			if ringsIntersect(rings[i], rings[j]) || Inside(rings[i], rings[j][0].X, rings[j][0].Y) {
				return nil, 0, &EdgeError{Ring: j, Edge: 0, Kind: ErrInvalidHole}
			}
		}
	}
	if total <= 0 {
		return nil, 0, ErrDegenerate
	}
	return rings, total, nil
}

// ringInside b 的首点是否严格在 a 内（简单环条件下足以代表包含关系）。
func ringsIntersect(a, b []Point) bool {
	for i := 0; i < len(a); i++ {
		a1, a2 := a[i], a[(i+1)%len(a)]
		for j := 0; j < len(b); j++ {
			b1, b2 := b[j], b[(j+1)%len(b)]
			if properIntersect(a1, a2, b1, b2) {
				return true
			}
		}
	}
	return false
}

// IsSimple 判定一个环是否为简单多边形：>=3 个互异顶点、无自交、
// 无非相邻边相触、无连续共线。供验收与外部调用。
func IsSimple(ring []Point) bool {
	return validateRing(ring, 0) == nil
}

// ---------------------------------------------------------------------------
// 方向：面积加权二阶矩的特征向量（闭式解，无迭代）
// ---------------------------------------------------------------------------

// principalAxis 求法向 n（主轴，最大特征值方向）与切线方向 u（短轴）。
// 三角形扇分解求二阶矩；孔洞环有向面积为负，自动带负权。
func principalAxis(rings [][]Point) (n, u Point) {
	var wSum, wxSum, wySum float64
	for _, r := range rings {
		for i := 0; i < len(r); i++ {
			a, b := r[i], r[(i+1)%len(r)]
			w := a.X*b.Y - b.X*a.Y // 2 倍有向三角形面积
			wSum += w
			wxSum += (a.X + b.X) * w
			wySum += (a.Y + b.Y) * w
		}
	}
	cx, cy := wxSum/(3*wSum), wySum/(3*wSum)

	// 关于形心的协方差矩阵（三角扇闭式，公共常数因子不影响特征向量方向）。
	var sxx, sxy, syy float64
	for _, r := range rings {
		for i := 0; i < len(r); i++ {
			a, b := r[i], r[(i+1)%len(r)]
			ax, ay := a.X-cx, a.Y-cy
			bx, by := b.X-cx, b.Y-cy
			w := ax*by - bx*ay
			sxx += (ax*ax + ax*bx + bx*bx) * w
			syy += (ay*ay + ay*by + by*by) * w
			sxy += (2*ax*ay + ax*by + bx*ay + 2*bx*by) * w
		}
	}
	sxx /= 12
	syy /= 12
	sxy /= 24

	// 2×2 对称矩阵特征向量闭式解：atan2 给出主轴（大特征值）方向角。
	theta := 0.5 * math.Atan2(2*sxy, sxx-syy)
	n = Point{math.Cos(theta), math.Sin(theta)}
	u = Point{-n.Y, n.X}
	return n, u
}

// ---------------------------------------------------------------------------
// 位置：板截面一次函数 + 板内二次面积方程（闭式解，无二分）
// ---------------------------------------------------------------------------

// crossSection 求直线 n·x = d（方向 u）与多边形的全部交点在 u 坐标上的
// 有序投影；排序后相邻两点配对即材料区间。
// 切线恰好穿过顶点或与边重合时，按“每条边贡献其落在切线上的端点”收集，
// 再按坐标容差去重，保证偶数个交点且可正确配对。
func crossSection(rings [][]Point, n, u Point, d float64) []float64 {
	proj := func(p Point) float64 { return p.X*n.X + p.Y*n.Y }
	upos := func(p Point) float64 { return p.X*u.X + p.Y*u.Y }
	var hits []float64
	add := func(p Point) {
		s := upos(p)
		for _, h := range hits {
			if math.Abs(h-s) <= tol(s, h) {
				return
			}
		}
		hits = append(hits, s)
	}
	for _, r := range rings {
		m := len(r)
		for i := 0; i < m; i++ {
			a, b := r[i], r[(i+1)%m]
			pa, pb := proj(a), proj(b)
			tv := tol(pa, pb, d)
			sa, sb := 0, 0
			if pa > d+tv {
				sa = 1
			} else if pa < d-tv {
				sa = -1
			}
			if pb > d+tv {
				sb = 1
			} else if pb < d-tv {
				sb = -1
			}
			switch {
			case sa == -sb && sa != 0: // 严格跨线：插值交点
				t := (d - pa) / (pb - pa)
				add(Point{a.X + (b.X-a.X)*t, a.Y + (b.Y-a.Y)*t})
			case sa == 0 && sb == 0: // 整条边在线上：两端都计入
				add(a)
				add(b)
			case sa == 0:
				add(a)
			case sb == 0:
				add(b)
			}
		}
	}
	for i := 1; i < len(hits); i++ {
		for j := i; j > 0 && hits[j-1] > hits[j]; j-- {
			hits[j-1], hits[j] = hits[j], hits[j-1]
		}
	}
	return hits
}

func sectionLength(hits []float64) float64 {
	var L float64
	for i := 0; i+1 < len(hits); i += 2 {
		L += hits[i+1] - hits[i]
	}
	return L
}

// sweepSolve 在固定法向下解 d，使负侧面积 = target。
// 板 [d0,d1] 内 L(d)=αd+β，累计面积
//
//	F(d)=F(d0)+ α/2*(d²-d0²) + β*(d-d0)
//
// 目标板用线性扫描定位（非二分），板内直接解二次方程。
func sweepSolve(rings [][]Point, n, u Point, target float64) (float64, error) {
	projSet := map[float64]struct{}{}
	for _, r := range rings {
		for _, p := range r {
			projSet[p.X*n.X+p.Y*n.Y] = struct{}{}
		}
	}
	ds := make([]float64, 0, len(projSet))
	for d := range projSet {
		ds = append(ds, d)
	}
	for i := 1; i < len(ds); i++ {
		for j := i; j > 0 && ds[j-1] > ds[j]; j-- {
			ds[j-1], ds[j] = ds[j], ds[j-1]
		}
	}
	uniq := ds[:0] // 合并数值上重合的板界
	for _, d := range ds {
		if len(uniq) == 0 || math.Abs(d-uniq[len(uniq)-1]) > tol(d) {
			uniq = append(uniq, d)
		}
	}
	ds = uniq

	var F float64 // 进入当前板时的累计负侧面积
	for k := 0; k+1 < len(ds); k++ {
		d0, d1 := ds[k], ds[k+1]
		if d1-d0 <= tol(d0, d1) {
			continue
		}
		// 板内 1/4、3/4 处取两个截面，精确求一次系数 α、β。
		mid := 0.5 * (d0 + d1)
		half := 0.25 * (d1 - d0)
		LA := sectionLength(crossSection(rings, n, u, mid-half))
		LB := sectionLength(crossSection(rings, n, u, mid+half))
		alpha := (LB - LA) / (2 * half)
		beta := LA - alpha*(mid-half)

		slabArea := alpha/2*(d1*d1-d0*d0) + beta*(d1-d0)
		if slabArea < 0 {
			slabArea = 0
		}
		if F+slabArea+tol(target) < target {
			F += slabArea
			continue
		}
		// 目标落在本板：令 h=d-d0，解 α/2*h² + β*h - want = 0。
		want := target - F
		A := alpha / 2
		B := beta
		var h float64
		if math.Abs(A) <= tol(alpha) {
			if math.Abs(B) <= tol(beta) {
				if math.Abs(want) > tol(target) {
					return 0, ErrUncuttable
				}
			} else {
				h = want / B // 退化为一次方程
			}
		} else {
			disc := B*B + 4*A*want // C=-want，故判别式 = B²-4AC = B²+4A·want
			if disc < -tol(target) {
				return 0, ErrUncuttable
			}
			if disc < 0 {
				disc = 0
			}
			sq := math.Sqrt(disc)
			h = pickRoot((-B-sq)/(2*A), (-B+sq)/(2*A), d1-d0)
		}
		if h < -tol(d0, d1) || h > (d1-d0)+tol(d0, d1) {
			return 0, ErrUncuttable
		}
		h = math.Max(0, math.Min(h, d1-d0))
		return d0 + h, nil
	}
	return 0, ErrUncuttable
}

// pickRoot 在两个二次方程根里选落在板宽 [0,w] 内的那个。
func pickRoot(r1, r2, w float64) float64 {
	in := func(r float64) bool { return r >= -epsRel && r <= w+epsRel }
	switch {
	case in(r1):
		return r1
	case in(r2):
		return r2
	default:
		if math.Abs(r1) < math.Abs(r1-w) {
			return 0
		}
		return w
	}
}

// ---------------------------------------------------------------------------
// 切块拼装：Sutherland–Hodgman 半平面裁剪 + 切边双零宽桥缝合
// ---------------------------------------------------------------------------

// sideVal 点到切线的有向距离 n·p - d。
func sideVal(p, n Point, d float64) float64 { return p.X*n.X + p.Y*n.Y - d }

// clipRing 用半平面 n·x <= d（negSide=true）或 n·x >= d（negSide=false）
// 对单个环做 Sutherland–Hodgman 裁剪。保持输入环的顶点朝向，交点由线性
// 插值得到。返回首尾不重复的闭合顶点序列；环完全落在另一侧时返回 nil。
func clipRing(ring []Point, n Point, d float64, negSide bool) []Point {
	inside := func(p Point) bool {
		v := sideVal(p, n, d)
		if negSide {
			return v <= epsAbs
		}
		return v >= -epsAbs
	}
	m := len(ring)
	var out []Point
	for i := 0; i < m; i++ {
		cur := ring[i]
		nxt := ring[(i+1)%m]
		curIn, nxtIn := inside(cur), inside(nxt)
		intersect := func() Point {
			vc, vn := sideVal(cur, n, d), sideVal(nxt, n, d)
			t := vc / (vc - vn)
			return Point{cur.X + (nxt.X-cur.X)*t, cur.Y + (nxt.Y-cur.Y)*t}
		}
		switch {
		case curIn && nxtIn:
			out = append(out, nxt)
		case curIn && !nxtIn:
			out = append(out, intersect())
		case !curIn && nxtIn:
			out = append(out, intersect(), nxt)
		}
	}
	cleaned := out[:0] // 去除连续重复点
	for _, p := range out {
		if len(cleaned) == 0 || !ptsEq(cleaned[len(cleaned)-1], p) {
			cleaned = append(cleaned, p)
		}
	}
	if len(cleaned) >= 2 && ptsEq(cleaned[0], cleaned[len(cleaned)-1]) {
		cleaned = cleaned[:len(cleaned)-1]
	}
	if len(cleaned) < 3 {
		return nil
	}
	return cleaned
}

// indexOfPoint 返回环中与 q 重合（容差）的顶点下标，找不到返回 -1。
func indexOfPoint(ring []Point, q Point) int {
	for i, p := range ring {
		if ptsEq(p, q) {
			return i
		}
	}
	return -1
}

// ringPath 沿环向从下标 from 走到下标 to，返回含两端的顶点序列。
func ringPath(ring []Point, from, to int) []Point {
	var path []Point
	for k := from; ; k = (k + 1) % len(ring) {
		path = append(path, ring[k])
		if k == to {
			return path
		}
	}
}

// dedupMerge 拼接多段路径时去掉相邻重复点。
func dedupMerge(segs ...[]Point) []Point {
	var out []Point
	push := func(seg []Point, skipFirst bool) {
		start := 0
		if skipFirst {
			start = 1
		}
		for k := start; k < len(seg); k++ {
			if len(out) == 0 || !ptsEq(out[len(out)-1], seg[k]) {
				out = append(out, seg[k])
			}
		}
	}
	for i, s := range segs {
		push(s, i > 0)
	}
	if len(out) >= 2 && ptsEq(out[0], out[len(out)-1]) {
		out = out[:len(out)-1]
	}
	return out
}

// onLinePts 返回环上落在切线上的顶点，按它们在环中出现的先后顺序排列。
func onLinePts(ring []Point, n, u Point, d float64) []Point {
	var pts []Point
	for _, p := range ring {
		if math.Abs(sideVal(p, n, d)) <= tol(p.X, p.Y, d) {
			pts = append(pts, p)
		}
	}
	return pts
}

// stitchOne 用两条沿切线的桥把一个洞环缝进外轮廓，返回缝合后的新外轮廓。
//
// 裁剪后外轮廓为逆时针（CCW）、洞环为顺时针（CW），每个环在切线上有一对
// 切点，按 u 坐标分低端 lo、高端 hi。外轮廓切点与洞切点的 u 坐标不同，
// 二者之间是切线上的真实线段（桥），不能省略。缝合环由四段构成：
//
//	外实体段 ＋ 桥(外lo→洞lo) ＋ 洞实体段 ＋ 桥(洞hi→外hi)
//
// 外、洞实体段各有两种环向。用面积不变量
//
//	signedArea(缝合环) == signedArea(外轮廓) + signedArea(洞环)
//
// 在四种组合里选出正确者；该判据与切线朝向、负/正侧均无关，天然保证
// 面积守恒且产出简单环。
func stitchOne(outer, hole []Point, n, u Point, d float64) []Point {
	oOn := onLinePts(outer, n, u, d)
	hOn := onLinePts(hole, n, u, d)
	if len(oOn) < 2 || len(hOn) < 2 {
		return outer
	}
	ukey := func(p Point) float64 { return p.X*u.X + p.Y*u.Y }
	loHi := func(a, b Point) (lo, hi Point) {
		if ukey(a) <= ukey(b) {
			return a, b
		}
		return b, a
	}
	oLo, oHi := loHi(oOn[0], oOn[1])
	hLo, hHi := loHi(hOn[0], hOn[1])
	iLo, iHi := indexOfPoint(outer, oLo), indexOfPoint(outer, oHi)
	jLo, jHi := indexOfPoint(hole, hLo), indexOfPoint(hole, hHi)
	if iLo < 0 || iHi < 0 || jLo < 0 || jHi < 0 {
		return outer
	}
	want := signedArea(outer) + signedArea(hole)
	oSegs := [][]Point{ringPath(outer, iHi, iLo), ringPath(outer, iLo, iHi)}
	hSegs := [][]Point{ringPath(hole, jLo, jHi), ringPath(hole, jHi, jLo)}
	best := outer
	bestDiff := math.Inf(1)
	for _, oSeg := range oSegs {
		// 外段两端必为 oHi、oLo；据此决定桥接顺序。
		oStartHi := ptsEq(oSeg[0], oHi)
		for _, hSeg := range hSegs {
			var merged []Point
			merged = append(merged, oSeg...)
			if oStartHi {
				// 外段 hi→lo：桥到洞 lo，洞段 lo→hi，桥回外 hi。
				if !ptsEq(hSeg[0], hLo) {
					continue
				}
				merged = append(merged, hLo)
				merged = append(merged, hSeg[1:]...)
				merged = append(merged, oHi)
			} else {
				// 外段 lo→hi：桥到洞 hi，洞段 hi→lo，桥回外 lo。
				if !ptsEq(hSeg[0], hHi) {
					continue
				}
				merged = append(merged, hHi)
				merged = append(merged, hSeg[1:]...)
				merged = append(merged, oLo)
			}
			merged = dedupMerge(merged, []Point{merged[0]})
			if diff := math.Abs(signedArea(merged) - want); diff < bestDiff {
				bestDiff, best = diff, merged
			}
		}
	}
	return best
}

// assemble 用半平面裁剪产出负、正两个简单多边形。
// rings[0] 为逆时针外轮廓，其余为顺时针孔洞；前提是每条孔洞都被切线贯穿。
func assemble(rings [][]Point, n, u Point, d float64) (neg, pos []Point, err error) {
	negOuter := clipRing(rings[0], n, d, true)
	posOuter := clipRing(rings[0], n, d, false)
	if negOuter == nil || posOuter == nil {
		return nil, nil, ErrAssembly
	}
	for hi := 1; hi < len(rings); hi++ {
		if h := clipRing(rings[hi], n, d, true); h != nil {
			negOuter = stitchOne(negOuter, h, n, u, d)
		}
		if h := clipRing(rings[hi], n, d, false); h != nil {
			posOuter = stitchOne(posOuter, h, n, u, d)
		}
	}
	if len(negOuter) < 3 || len(posOuter) < 3 {
		return nil, nil, ErrAssembly
	}
	return negOuter, posOuter, nil
}

// ---------------------------------------------------------------------------
// 对外入口
// ---------------------------------------------------------------------------

// CutPolygon 按目标面积比 ratio 切分带洞多边形，tolRatio 为面积比绝对容差。
// 返回的 CutResult 中切线方程与两块面积一并给出。
// 非法输入（比值越界 / 自交 / 连续共线）以 *EdgeError 报出边号。
func CutPolygon(p *Polygon, ratio, tolRatio float64) (*CutResult, error) {
	if p == nil {
		return nil, ErrDegenerate
	}
	if math.IsNaN(ratio) || ratio <= 0 || ratio >= 1 {
		return nil, ErrRatioOutOfRange
	}
	if tolRatio <= 0 || math.IsNaN(tolRatio) {
		tolRatio = epsRel
	}
	rings, areaAll, err := validatePolygon(p)
	if err != nil {
		return nil, err
	}
	n, u := principalAxis(rings)
	d, err := sweepSolve(rings, n, u, ratio*areaAll)
	if err != nil {
		return nil, err
	}
	line := CutLine{N: n, Dir: u, D: d}

	// 孔洞必须被切线穿过，否则含洞那一块不可能是简单多边形（单环）。
	for hi := 1; hi < len(rings); hi++ {
		if len(crossSection(rings[hi:hi+1], n, u, d)) == 0 {
			return nil, fmt.Errorf("%w：环 %d 的孔洞未被切线穿过，单一直线无法把它切成两个简单多边形", ErrUncuttable, hi)
		}
	}

	neg, pos, err := assemble(rings, n, u, d)
	if err != nil {
		return nil, fmt.Errorf("%w：%v", ErrAssembly, err)
	}
	// 统一输出为逆时针。
	if signedArea(neg) < 0 {
		reverse(neg)
	}
	if signedArea(pos) < 0 {
		reverse(pos)
	}

	aNeg := math.Abs(signedArea(neg))
	aPos := math.Abs(signedArea(pos))
	got := aNeg / areaAll
	res := &CutResult{
		Line:    line,
		Neg:     neg,
		Pos:     pos,
		AreaNeg: aNeg,
		AreaPos: aPos,
		AreaAll: areaAll,
		Ratio:   got,
		AbsErr:  math.Abs(got - ratio),
	}
	// 两条独立验收，分开记录、分开判定。
	res.AreaSumOK = math.Abs((aNeg+aPos)-areaAll) <= tol(areaAll)*10
	res.PiecesOK = IsSimple(neg) && IsSimple(pos)
	if res.AbsErr > tolRatio {
		return res, fmt.Errorf("%w：实际比误差 %g 超过容差 %g", ErrUncuttable, res.AbsErr, tolRatio)
	}
	if !res.AreaSumOK || !res.PiecesOK {
		return res, fmt.Errorf("%w：面积守恒=%v 简单性=%v", ErrAssembly, res.AreaSumOK, res.PiecesOK)
	}
	return res, nil
}

// CutRing 是无洞多边形的便捷入口，与原库 []Point 风格一致。
func CutRing(ring []Point, ratio, tolRatio float64) (*CutResult, error) {
	return CutPolygon(&Polygon{Outer: ring}, ratio, tolRatio)
}
