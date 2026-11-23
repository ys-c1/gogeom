// 渲染器：把多边形的“整块填充（before）”与“三角网格填充（after）”
// 渲染成 360x360 P5 PGM，灰度约定与参考产物一致：
// 背景 245、整块填充 90、三角分片按序循环 70/110/150/190、轮廓 30。
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"

	"gogeom/geom"
)

const (
	W, H    = 360, 360
	bg      = 245
	fill    = 90
	outline = 30
)

var triShades = []byte{70, 110, 150, 190}

type canvas struct {
	buf [W * H]byte
}

func newCanvas() *canvas {
	c := &canvas{}
	for i := range c.buf {
		c.buf[i] = bg
	}
	return c
}

func (c *canvas) set(x, y int, v byte) {
	if x >= 0 && x < W && y >= 0 && y < H {
		c.buf[y*W+x] = v
	}
}

func (c *canvas) get(x, y int) byte {
	if x < 0 || x >= W || y < 0 || y >= H {
		return bg
	}
	return c.buf[y*W+x]
}

// fillPoly 用射线奇偶规则光栅化填充多边形（整块）
func (c *canvas) fillPoly(ring []geom.Point, v byte) {
	n := len(ring)
	for py := 0; py < H; py++ {
		y := float64(py) + 0.5
		var xs []float64
		for i := 0; i < n; i++ {
			a, b := ring[i], ring[(i+1)%n]
			if (a.Y > y) == (b.Y > y) {
				continue
			}
			t := (y - a.Y) / (b.Y - a.Y)
			xs = append(xs, a.X+(b.X-a.X)*t)
		}
		// 简单插入排序
		for i := 1; i < len(xs); i++ {
			for j := i; j > 0 && xs[j-1] > xs[j]; j-- {
				xs[j-1], xs[j] = xs[j], xs[j-1]
			}
		}
		for k := 0; k+1 < len(xs); k += 2 {
			lo := int(xs[k] + 0.5)
			hi := int(xs[k+1] + 0.5)
			for x := lo; x < hi; x++ {
				c.set(x, py, v)
			}
		}
	}
}

// fillTri 重心坐标光栅化三角形
func (c *canvas) fillTri(t geom.Triangle, v byte) {
	minX := maxi(0, mini3(t[0].X, t[1].X, t[2].X))
	maxX := mini(W-1, maxi3(t[0].X, t[1].X, t[2].X))
	minY := maxi(0, mini3(t[0].Y, t[1].Y, t[2].Y))
	maxY := mini(H-1, maxi3(t[0].Y, t[1].Y, t[2].Y))
	den := (t[1].Y-t[2].Y)*(t[0].X-t[2].X) + (t[2].X-t[1].X)*(t[0].Y-t[2].Y)
	if den == 0 {
		return
	}
	for py := minY; py <= maxY; py++ {
		for px := minX; px <= maxX; px++ {
			x, y := float64(px)+0.5, float64(py)+0.5
			a := ((t[1].Y-t[2].Y)*(x-t[2].X) + (t[2].X-t[1].X)*(y-t[2].Y)) / den
			b := ((t[2].Y-t[0].Y)*(x-t[2].X) + (t[0].X-t[2].X)*(y-t[2].Y)) / den
			cc := 1 - a - b
			if a >= 0 && b >= 0 && cc >= 0 {
				c.set(px, py, v)
			}
		}
	}
}

// drawLine 整数 DDA 画轮廓
func (c *canvas) drawLine(p0, p1 geom.Point, v byte) {
	dx := p1.X - p0.X
	dy := p1.Y - p0.Y
	steps := int(maxf(mathAbs(dx), mathAbs(dy))) + 1
	for i := 0; i <= steps; i++ {
		t := float64(i) / float64(steps)
		c.set(int(p0.X+dx*t+0.5), int(p0.Y+dy*t+0.5), v)
	}
}

func (c *canvas) drawRing(ring []geom.Point) {
	for i := 0; i < len(ring); i++ {
		c.drawLine(ring[i], ring[(i+1)%len(ring)], outline)
	}
}

func (c *canvas) writePGM(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	fmt.Fprintf(w, "P5\n%d %d\n255\n", W, H)
	_, err = w.Write(c.buf[:])
	return err
}

func mini3(a, b, c float64) int { return int(minf(minf(a, b), c)) }
func maxi3(a, b, c float64) int { return int(maxf(maxf(a, b), c) + 0.999) }
func mini(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func maxi(a, b int) int {
	if a > b {
		return a
	}
	return b
}
func minf(a, b float64) float64 {
	if a < b {
		return a
	}
	return b
}
func maxf(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}
func mathAbs(a float64) float64 {
	if a < 0 {
		return -a
	}
	return a
}

// renderCase 渲染单个用例的 before/after
func renderCase(dir, tag string, rings [][]geom.Point, holes [][]geom.Point) error {
	var allRings [][]geom.Point
	allRings = append(allRings, rings...)
	allRings = append(allRings, holes...)

	// before：整块奇偶填充
	before := newCanvas()
	for _, r := range allRings {
		before.fillPoly(r, fill)
	}
	// 洞在奇偶规则下二次扫过会变回背景，单独清洞
	for _, h := range holes {
		before.fillPoly(h, bg)
	}
	for _, r := range allRings {
		before.drawRing(r)
	}
	if err := before.writePGM(filepath.Join(dir, tag+"_before.pgm")); err != nil {
		return err
	}

	// after：三角网格
	var tris []geom.Triangle
	if len(holes) > 0 {
		t, err := geom.TriangulateRings(rings[0], holes)
		if err != nil {
			return err
		}
		tris = t
	} else {
		t, err := geom.Triangulate(rings[0])
		if err != nil {
			return err
		}
		tris = t
	}
	after := newCanvas()
	for i, t := range tris {
		after.fillTri(t, triShades[i%len(triShades)])
	}
	for _, r := range allRings {
		after.drawRing(r)
	}
	return after.writePGM(filepath.Join(dir, tag+"_after.pgm"))
}

func main() {
	dir := os.Args[1]
	if err := os.MkdirAll(dir, 0o755); err != nil {
		panic(err)
	}
	type tc struct {
		tag   string
		rings [][]geom.Point
		holes [][]geom.Point
	}
	cases := []tc{
		{"tri_00", [][]geom.Point{{{X: 40, Y: 40}, {X: 320, Y: 40}, {X: 320, Y: 320}, {X: 40, Y: 320}}}, nil},
		{"tri_01", [][]geom.Point{octagon()}, nil},
		{"tri_02", [][]geom.Point{{{X: 40, Y: 320}, {X: 40, Y: 40}, {X: 150, Y: 40}, {X: 150, Y: 180}, {X: 320, Y: 180}, {X: 320, Y: 320}}}, nil},
		{"tri_03", [][]geom.Point{star8()}, nil},
		{"tri_04", [][]geom.Point{{{X: 40, Y: 40}, {X: 320, Y: 40}, {X: 320, Y: 320}, {X: 40, Y: 320}}},
			[][]geom.Point{{{X: 130, Y: 130}, {X: 130, Y: 230}, {X: 230, Y: 230}, {X: 230, Y: 130}}}},
		{"tri_05", [][]geom.Point{{{X: 20, Y: 110}, {X: 340, Y: 110}, {X: 340, Y: 250}, {X: 20, Y: 250}}},
			[][]geom.Point{
				{{X: 50, Y: 140}, {X: 50, Y: 220}, {X: 120, Y: 220}, {X: 120, Y: 140}},
				{{X: 200, Y: 140}, {X: 200, Y: 220}, {X: 300, Y: 220}, {X: 300, Y: 140}},
			}},
		{"tri_06", [][]geom.Point{circleN(180, 180, 150, 48)}, [][]geom.Point{circleN(180, 180, 70, 32)}},
		{"tri_07", [][]geom.Point{{{X: 40, Y: 40}, {X: 320, Y: 320}, {X: 320, Y: 40}, {X: 40, Y: 320}}}, nil},
		{"tri_08", [][]geom.Point{pentagram()}, nil},
		{"tri_09", [][]geom.Point{{{X: 40, Y: 40}, {X: 320, Y: 40}, {X: 320, Y: 320}, {X: 180, Y: 320 + 1e-6}, {X: 40, Y: 320}}}, nil},
	}
	for _, c := range cases {
		if err := renderCase(dir, c.tag, c.rings, c.holes); err != nil {
			fmt.Printf("%s 渲染失败: %v\n", c.tag, err)
		} else {
			fmt.Printf("%s 完成\n", c.tag)
		}
	}
}

func octagon() []geom.Point {
	out := make([]geom.Point, 8)
	for i := range out {
		a := float64(i) * math.Pi / 4
		out[i] = geom.Point{X: 180 + 140*math.Cos(a), Y: 180 + 140*math.Sin(a)}
	}
	return out
}

func star8() []geom.Point {
	out := make([]geom.Point, 16)
	for i := range out {
		a := float64(i) * math.Pi / 8
		r := 140.0
		if i%2 == 1 {
			r = 60
		}
		out[i] = geom.Point{X: 180 + r*math.Cos(a), Y: 180 + r*math.Sin(a)}
	}
	return out
}

func circleN(cx, cy, r float64, n int) []geom.Point {
	out := make([]geom.Point, n)
	for i := range out {
		a := float64(i) * 2 * math.Pi / float64(n)
		out[i] = geom.Point{X: cx + r*math.Cos(a), Y: cy + r*math.Sin(a)}
	}
	return out
}

func pentagram() []geom.Point {
	out := make([]geom.Point, 5)
	for i := range out {
		a := float64(i)*4*math.Pi/5 - math.Pi/2
		out[i] = geom.Point{X: 180 + 150*math.Cos(a), Y: 180 + 150*math.Sin(a)}
	}
	return out
}
