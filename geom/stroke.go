package geom

import (
	"fmt"
	"math"
)

// JoinStyle 描边拐角接头样式。
type JoinStyle int

const (
	JoinMiter JoinStyle = iota // 直折（尖角），默认分支
	JoinRound                  // 圆角
	JoinBevel                  // 斜切
)

// StrokeParams 控制描边：线宽、接头样式、每顶点接头半径 r、圆角离散段数。
// Radius 长度必须为 0 或 len(path)-2（每个内部拐点一个半径）；为空时取半线宽 h。
// 半径为 0 时圆角/斜切接头退化为无接头，直接报错。
type StrokeParams struct {
	Width   float64
	Style   JoinStyle
	Radius  []float64
	ArcSegs int
}

// Raster 以像素网格中心点采样得到的覆盖结果。
type Raster struct {
	MinX, MinY float64
	Pixel      float64
	W, H       int
	Cover      []bool
}

// Count 返回覆盖像素总数。
func (r *Raster) Count() int {
	n := 0
	for _, v := range r.Cover {
		if v {
			n++
		}
	}
	return n
}

func (r *Raster) index(x, y float64) (int, int, bool) {
	ix := int(math.Floor((x - r.MinX) / r.Pixel))
	iy := int(math.Floor((y - r.MinY) / r.Pixel))
	if ix < 0 || iy < 0 || ix >= r.W || iy >= r.H {
		return 0, 0, false
	}
	return ix, iy, true
}

// Set 标记世界坐标 (x,y) 所在像素为覆盖。
func (r *Raster) Set(x, y float64) {
	if ix, iy, ok := r.index(x, y); ok {
		r.Cover[iy*r.W+ix] = true
	}
}

// Get 查询世界坐标 (x,y) 所在像素是否覆盖。
func (r *Raster) Get(x, y float64) bool {
	if ix, iy, ok := r.index(x, y); ok {
		return r.Cover[iy*r.W+ix]
	}
	return false
}

// Components 对覆盖像素做 4 邻域洪水填充，返回连通分量个数。
func (r *Raster) Components() int {
	seen := make([]bool, len(r.Cover))
	comps := 0
	var stack []int
	for i, on := range r.Cover {
		if !on || seen[i] {
			continue
		}
		comps++
		stack = append(stack, i)
		seen[i] = true
		for len(stack) > 0 {
			cur := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			x, y := cur%r.W, cur/r.W
			if x > 0 {
				stack = pushNeighbor(r, seen, stack, cur-1)
			}
			if x+1 < r.W {
				stack = pushNeighbor(r, seen, stack, cur+1)
			}
			if y > 0 {
				stack = pushNeighbor(r, seen, stack, cur-r.W)
			}
			if y+1 < r.H {
				stack = pushNeighbor(r, seen, stack, cur+r.W)
			}
		}
	}
	return comps
}

func pushNeighbor(r *Raster, seen []bool, stack []int, q int) []int {
	if r.Cover[q] && !seen[q] {
		seen[q] = true
		return append(stack, q)
	}
	return stack
}

func hypot(p Point) float64 { return math.Hypot(p.X, p.Y) }

func leftNormal(t Point) Point  { return Point{-t.Y, t.X} }
func rightNormal(t Point) Point { return Point{t.Y, -t.X} }

func rotate(v Point, ang float64) Point {
	c, s := math.Cos(ang), math.Sin(ang)
	return Point{v.X*c - v.Y*s, v.X*s + v.Y*c}
}

// edgeTangent 返回边 a→b 的单位切向与长度。
func edgeTangent(a, b Point) (t Point, length float64) {
	d := Point{b.X - a.X, b.Y - a.Y}
	l := hypot(d)
	return Point{d.X / l, d.Y / l}, l
}

// joinAngle 计算内部拐点处的接头圆心角 α∈(0,π]，
// 即两条边切向夹角（= 两条外法线夹角）；共线返回 π，折返返回 0。
func joinAngle(t0, t1 Point) float64 {
	c := t0.X*t1.X + t0.Y*t1.Y
	c = math.Max(-1, math.Min(1, c))
	return math.Acos(c)
}

// isCollinear 判定两切向是否共线（同向 α≈π，或反向 α≈0）。
func isCollinear(alpha float64) bool {
	return math.Abs(alpha-math.Pi) < 1e-12 || alpha < 1e-12
}

// bandQuad 一条边的带状四边形（含两个平头端盖）。
func bandQuad(a, b Point, h float64) []Point {
	t, _ := edgeTangent(a, b)
	lp, rp := leftNormal(t), rightNormal(t)
	return []Point{
		{a.X + h*lp.X, a.Y + h*lp.Y},
		{b.X + h*lp.X, b.Y + h*lp.Y},
		{b.X + h*rp.X, b.Y + h*rp.Y},
		{a.X + h*rp.X, a.Y + h*rp.Y},
	}
}

// joinOuterNormal 取拐点的外法线（外楔形角平分线），
// 方向由 (入边切向)×(出边切向) 的符号决定。
func joinOuterNormal(t0, t1 Point) Point {
	cross := t0.X*t1.Y - t0.Y*t1.X
	m := Point{leftNormal(t0).X + leftNormal(t1).X, leftNormal(t0).Y + leftNormal(t1).Y}
	ml := hypot(m)
	if ml == 0 {
		return leftNormal(t0) // 严格折返的退化兜底
	}
	m = Point{m.X / ml, m.Y / ml}
	if cross < 0 {
		m = Point{-m.X, -m.Y}
	}
	return m
}

// joinSide 返回该拐点外侧所用的边法线函数（左法线或右法线）与旋转符号。
func joinSide(t0, t1 Point) (func(Point) Point, float64) {
	if t0.X*t1.Y-t0.Y*t1.X >= 0 {
		return leftNormal, 1
	}
	return rightNormal, -1
}

// miterJoin 直折接头多边形（外尖角四边形），半径不参与。
func miterJoin(v, t0, t1 Point, h float64) []Point {
	alpha := joinAngle(t0, t1)
	nOut := joinOuterNormal(t0, t1)
	side, _ := joinSide(t0, t1)
	n0, n1 := side(t0), side(t1)
	d := h / math.Cos(alpha/2) // 尖角沿外法线外推距离
	tip := Point{v.X + nOut.X*d, v.Y + nOut.Y*d}
	q0 := Point{v.X + n0.X*h, v.Y + n0.Y*h}
	q1 := Point{v.X + n1.X*h, v.Y + n1.Y*h}
	return []Point{q0, tip, q1, v}
}

// bevelJoin 斜切接头多边形：斜切弦 c0-c1 与 v 平行外移，半径 r 为弦到 v 的
// 法向距离。c0 = v+r·n0-d·t0，c1 = v+r·n1+d·t1，d = r·tan(α/2)
// （两点到 v 距离均为 r/cos(α/2)）。r=0 退化为 v（无接头），r=h 时弦即 q0-q1。
// 凸多边形顶点顺序：q0 → c0 → c1 → q1 → v。
func bevelJoin(v, t0, t1 Point, h, r float64) []Point {
	alpha := joinAngle(t0, t1)
	side, _ := joinSide(t0, t1)
	n0, n1 := side(t0), side(t1)
	d := r * math.Tan(alpha/2)
	c0 := Point{v.X + r*n0.X - d*t0.X, v.Y + r*n0.Y - d*t0.Y}
	c1 := Point{v.X + r*n1.X + d*t1.X, v.Y + r*n1.Y + d*t1.Y}
	q0 := Point{v.X + h*n0.X, v.Y + h*n0.Y}
	q1 := Point{v.X + h*n1.X, v.Y + h*n1.Y}
	return []Point{q0, c0, c1, q1, v}
}

// roundJoin 圆角接头多边形：v 为圆心、半径 r 的外楔形扇形。
// arcSegs 控制圆弧离散段数（偶数，关于角平分线对称）。
func roundJoin(v, t0, t1 Point, h, r float64, arcSegs int) []Point {
	alpha := joinAngle(t0, t1)
	side, sign := joinSide(t0, t1)
	n0, n1 := side(t0), side(t1)
	a := Point{v.X + n0.X*h, v.Y + n0.Y*h}
	b := Point{v.X + n1.X*h, v.Y + n1.Y*h}
	p0 := Point{v.X + n0.X*r, v.Y + n0.Y*r}
	half := arcSegs / 2
	if half < 1 {
		half = 1
	}
	poly := []Point{a, p0}
	// 前半弧：n0 → 角平分线
	for k := 1; k <= half; k++ {
		dir := rotate(n0, sign*alpha*float64(k)/float64(2*half))
		poly = append(poly, Point{v.X + dir.X*r, v.Y + dir.Y*r})
	}
	// 后半弧：角平分线 → n1
	for k := half - 1; k >= 0; k-- {
		dir := rotate(n1, -sign*alpha*float64(k)/float64(2*half))
		poly = append(poly, Point{v.X + dir.X*r, v.Y + dir.Y*r})
	}
	poly = append(poly, b, v)
	return poly
}

// joinPolys 生成路径上所有内部拐点的接头多边形，并返回每个拐点的圆心角。
func joinPolys(path []Point, p StrokeParams) ([][]Point, []float64, error) {
	arcSegs := p.ArcSegs
	if arcSegs <= 0 {
		arcSegs = 64
	}
	n := len(path)
	tan := make([]Point, n-1)
	for i := 0; i+1 < n; i++ {
		t, l := edgeTangent(path[i], path[i+1])
		if l == 0 {
			return nil, nil, fmt.Errorf("%w: 顶点 %d 处零长边", ErrDegenerate, i)
		}
		tan[i] = t
	}
	polys := make([][]Point, 0, n-2)
	angles := make([]float64, 0, n-2)
	for i := 1; i < n-1; i++ {
		t0, t1 := tan[i-1], tan[i]
		alpha := joinAngle(t0, t1)
		// 共线（同向 180° 或反向折返）无法形成接头，三种样式统一报错并带拐点号
		if isCollinear(alpha) {
			return nil, nil, fmt.Errorf("%w: 拐点 %d", ErrCollinearJoin, i)
		}
		r := p.Width / 2
		if len(p.Radius) > 0 {
			r = p.Radius[i-1]
		}
		switch p.Style {
		case JoinMiter:
			polys = append(polys, miterJoin(path[i], t0, t1, p.Width/2))
		case JoinRound:
			if r <= 0 {
				return nil, nil, fmt.Errorf("%w: 拐点 %d", ErrZeroRadius, i)
			}
			if r > p.Width/2 {
				r = p.Width / 2
			}
			polys = append(polys, roundJoin(path[i], t0, t1, p.Width/2, r, arcSegs))
		case JoinBevel:
			if r <= 0 {
				return nil, nil, fmt.Errorf("%w: 拐点 %d", ErrZeroRadius, i)
			}
			if r > p.Width/2 {
				r = p.Width / 2
			}
			polys = append(polys, bevelJoin(path[i], t0, t1, p.Width/2, r))
		default:
			return nil, nil, fmt.Errorf("%w: 拐点 %d 的样式值 %d", ErrBadJoinStyle, i, p.Style)
		}
		angles = append(angles, alpha)
	}
	return polys, angles, nil
}

// StrokePolygons 返回描边覆盖区域的全部多边形（边带 + 接头），可直接用于填充。
// 默认（Style == JoinMiter）不读取 Radius，输出保持直折分支不变。
func StrokePolygons(path []Point, p StrokeParams) ([][]Point, error) {
	if len(path) < 2 {
		return nil, ErrTooFewPts
	}
	if p.Width <= 0 {
		return nil, fmt.Errorf("%w: 线宽必须为正", ErrDegenerate)
	}
	if p.Style < JoinMiter || p.Style > JoinBevel {
		return nil, fmt.Errorf("%w: 样式值 %d", ErrBadJoinStyle, p.Style)
	}
	if len(p.Radius) != 0 && len(p.Radius) != len(path)-2 {
		return nil, fmt.Errorf("%w: 半径数 %d 与拐点数 %d 不一致", ErrDegenerate, len(p.Radius), len(path)-2)
	}
	h := p.Width / 2
	out := make([][]Point, 0, 2*len(path)-1)
	for i := 0; i+1 < len(path); i++ {
		out = append(out, bandQuad(path[i], path[i+1], h))
	}
	polys, _, err := joinPolys(path, p)
	if err != nil {
		return nil, err
	}
	out = append(out, polys...)
	return out, nil
}

// rasterBounds 计算覆盖所有多边形的像素网格边界（向外扩 1 像素）。
func rasterBounds(polys [][]Point, pixel float64) (minX, minY float64, w, h int) {
	minX, minY = math.Inf(1), math.Inf(1)
	maxX, maxY := math.Inf(-1), math.Inf(-1)
	for _, poly := range polys {
		for _, q := range poly {
			minX = math.Min(minX, q.X)
			minY = math.Min(minY, q.Y)
			maxX = math.Max(maxX, q.X)
			maxY = math.Max(maxY, q.Y)
		}
	}
	minX -= pixel
	minY -= pixel
	w = int(math.Ceil((maxX-minX)/pixel)) + 1
	h = int(math.Ceil((maxY-minY)/pixel)) + 1
	return
}

// RasterizePolys 以 pixel 为像素边长，用网格中心采样把多边形栅格化（并集）。
func RasterizePolys(polys [][]Point, pixel float64) *Raster {
	minX, minY, w, h := rasterBounds(polys, pixel)
	r := &Raster{MinX: minX, MinY: minY, Pixel: pixel, W: w, H: h, Cover: make([]bool, w*h)}
	for j := 0; j < h; j++ {
		y := minY + (float64(j)+0.5)*pixel
		for i := 0; i < w; i++ {
			x := minX + (float64(i)+0.5)*pixel
			for _, poly := range polys {
				if Inside(poly, x, y) {
					r.Cover[j*w+i] = true
					break
				}
			}
		}
	}
	return r
}

// Stroke 一步完成多边形构造与栅格化，返回覆盖网格。
func Stroke(path []Point, p StrokeParams, pixel float64) (*Raster, error) {
	polys, err := StrokePolygons(path, p)
	if err != nil {
		return nil, err
	}
	return RasterizePolys(polys, pixel), nil
}

// JoinGapPixels 计算单个接头相对直折尖角信封的缺口像素数：
// 直折多边形覆盖、但当前接头（半径 r）未覆盖的网格中心数。
// 同一像素网格下结果可直接比较。
func JoinGapPixels(v, t0, t1 Point, h, r float64, style JoinStyle, arcSegs int, pixel float64) (gap int, miterPoly, joinPoly []Point, err error) {
	alpha := joinAngle(t0, t1)
	if isCollinear(alpha) {
		return 0, nil, nil, fmt.Errorf("%w: 共线接头", ErrCollinearJoin)
	}
	miterPoly = miterJoin(v, t0, t1, h)
	switch style {
	case JoinMiter:
		joinPoly = miterPoly
	case JoinRound:
		if r <= 0 {
			return 0, nil, nil, fmt.Errorf("%w", ErrZeroRadius)
		}
		if r > h {
			r = h
		}
		joinPoly = roundJoin(v, t0, t1, h, r, arcSegs)
	case JoinBevel:
		if r <= 0 {
			return 0, nil, nil, fmt.Errorf("%w", ErrZeroRadius)
		}
		if r > h {
			r = h
		}
		joinPoly = bevelJoin(v, t0, t1, h, r)
	default:
		return 0, nil, nil, fmt.Errorf("%w: 样式值 %d", ErrBadJoinStyle, style)
	}
	minX, minY, w, hg := rasterBounds([][]Point{miterPoly}, pixel)
	for j := 0; j < hg; j++ {
		y := minY + (float64(j)+0.5)*pixel
		for i := 0; i < w; i++ {
			x := minX + (float64(i)+0.5)*pixel
			if Inside(miterPoly, x, y) && !Inside(joinPoly, x, y) {
				gap++
			}
		}
	}
	return gap, miterPoly, joinPoly, nil
}

// RadiusTurns 标记半径序列的曲线拐点（二阶差分符号变化处）。
// 入参为每个内部拐点的半径，返回序列下标（对应路径顶点 i+1）；单调序列返回空。
func RadiusTurns(radii []float64) []int {
	if len(radii) < 3 {
		return nil
	}
	scale := 0.0
	for _, r := range radii {
		scale = math.Max(scale, math.Abs(r))
	}
	eps := scale * 1e-9
	var out []int
	prev := 0
	for i := 1; i+1 < len(radii); i++ {
		d2 := radii[i+1] - 2*radii[i] + radii[i-1]
		sign := 0
		if d2 > eps {
			sign = 1
		} else if d2 < -eps {
			sign = -1
		}
		if sign != 0 && prev != 0 && sign != prev {
			out = append(out, i)
		}
		if sign != 0 {
			prev = sign
		}
	}
	return out
}

// 以下导出包装供命令行工具与外部调用使用，内部逻辑与同名未导出函数一致。

// MiterPoly 返回直折接头多边形。
func MiterPoly(v, t0, t1 Point, h float64) []Point { return miterJoin(v, t0, t1, h) }

// RoundPoly 返回圆角接头多边形。
func RoundPoly(v, t0, t1 Point, h, r float64, arcSegs int) []Point {
	return roundJoin(v, t0, t1, h, r, arcSegs)
}

// BevelPoly 返回斜切接头多边形。
func BevelPoly(v, t0, t1 Point, h, r float64) []Point {
	return bevelJoin(v, t0, t1, h, r)
}

// InsidePoly 判定点是否在多边形内（射线法）。
func InsidePoly(poly []Point, x, y float64) bool { return Inside(poly, x, y) }
