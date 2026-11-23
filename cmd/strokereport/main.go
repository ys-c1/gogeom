// strokereport 生成描边接头的验证报告：
//  1. 两档线宽的覆盖像素与耗时（三种样式）
//  2. 单接头覆盖像素 vs 手算面积，误差与严格上界
//  3. 十档半径扫描的缺口像素（圆角/斜切）
//  4. 半径连成曲线时的无缝（单连通）验证与曲线拐点
//  5. 可视化 PNG 输出到 -out 目录
package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
	"path/filepath"
	"time"

	"gogeom/geom"
)

type gapRow struct {
	k         int
	r         float64
	roundGap  int
	bevelGap  int
	roundPoly []geom.Point
	bevelPoly []geom.Point
}

func main() {
	outDir := flag.String("out", "artifacts", "报告与图像输出目录")
	pixel := flag.Float64("pixel", 0.2, "栅格化像素边长")
	flag.Parse()
	if err := os.MkdirAll(*outDir, 0o755); err != nil {
		panic(err)
	}

	const arcSegs = 128
	widths := []float64{10, 30} // 两档线宽

	// ---- 1. 两档线宽：单接头像素与耗时 ----
	fmt.Println("================ 1. 两档线宽：单接头覆盖像素与耗时 ================")
	fmt.Printf("像素边长 p=%.2f，孤立 40° 拐角；圆角半径 r=半线宽，斜切 r=半线宽/2，圆弧 %d 段\n", *pixel, arcSegs)
	fmt.Println("线宽  半宽h  样式   覆盖像素  栅格化耗时(ms)  构造耗时(µs)")
	jv := geom.Point{X: 40, Y: 40}
	jt0 := geom.Point{X: 1, Y: 0}
	jRad := 40.0 * math.Pi / 180
	jt1 := geom.Point{X: math.Cos(jRad), Y: math.Sin(jRad)}
	for _, w := range widths {
		h := w / 2
		for _, st := range []struct {
			name string
			s    geom.JoinStyle
		}{
			{"直折", geom.JoinMiter}, {"圆角", geom.JoinRound}, {"斜切", geom.JoinBevel},
		} {
			t0 := time.Now()
			var poly []geom.Point
			switch st.s {
			case geom.JoinMiter:
				poly = geom.MiterPoly(jv, jt0, jt1, h)
			case geom.JoinRound:
				poly = geom.RoundPoly(jv, jt0, jt1, h, h, arcSegs)
			case geom.JoinBevel:
				poly = geom.BevelPoly(jv, jt0, jt1, h, h/2)
			}
			build := time.Since(t0)
			var best time.Duration
			var r *geom.Raster
			for k := 0; k < 30; k++ {
				t1 := time.Now()
				r = geom.RasterizePolys([][]geom.Point{poly}, *pixel)
				if d := time.Since(t1); k == 0 || d < best {
					best = d
				}
			}
			fmt.Printf("%4.0f  %5.1f  %-4s  %8d  %12.3f  %12.3f\n",
				w, h, st.name, r.Count(), float64(best.Microseconds())/1000, float64(build.Microseconds()))
			writePNG(filepath.Join(*outDir, fmt.Sprintf("join_w%02.0f_%s.png", w, st.name)), r)
		}
	}

	// ---- 2. 接头覆盖像素 vs 手算面积 ----
	fmt.Println()
	fmt.Println("================ 2. 单接头：像素面积 vs 手算面积 ================")
	fmt.Println("（p=0.1；对称差严格上界 = 周长·2p + πp²，离散扇形另给 O(r²α³/n²) 上界）")
	fmt.Println("角度  样式   半径   手算面积    像素面积     |差|       严格上界     上界内?")
	h := 5.0
	v := geom.Point{X: 30, Y: 30}
	t0v := geom.Point{X: 1, Y: 0}
	const pFine = 0.1
	for _, deg := range []float64{30, 60, 90, 120, 150} {
		rad := deg * math.Pi / 180
		t1 := geom.Point{X: math.Cos(rad), Y: math.Sin(rad)}
		type entry struct {
			name     string
			poly     []geom.Point
			analytic float64
		}
		entries := []entry{
			{"直折", geom.MiterPoly(v, t0v, t1, h), h * h * math.Tan(rad/2)},
			{"斜切", geom.BevelPoly(v, t0v, t1, h, h/2), h * (h / 2) * math.Tan(rad/2)},
		}
		// 圆角：理想扇形面积（离散多边形另算上界）
		rp := geom.RoundPoly(v, t0v, t1, h, h, arcSegs)
		entries = append(entries, entry{"圆角", rp, 0.5 * h * h * rad})
		for _, e := range entries {
			perim := perimeter(e.poly)
			r := geom.RasterizePolys([][]geom.Point{e.poly}, pFine)
			pixArea := float64(r.Count()) * pFine * pFine
			bound := perim*2*pFine + math.Pi*pFine*pFine
			diff := math.Abs(pixArea - e.analytic)
			ok := diff <= bound
			mark := "是"
			if !ok {
				mark = "否"
			}
			fmt.Printf("%3.0f°  %-4s  %4.1f  %10.4f  %10.4f  %9.2e  %10.2e   %s\n",
				deg, e.name, h, e.analytic, pixArea, diff, bound, mark)
		}
		// 圆角离散误差上界单列
		n := arcSegs
		discBound := h * h * rad * rad * rad / (6 * float64(n) * float64(n))
		polyA := polygonArea(rp)
		ideal := 0.5 * h * h * rad
		fmt.Printf("       └ 圆角离散: 多边形面积=%.4f 理想扇形=%.4f |差|=%.2e 离散上界=%.2e n=%d\n",
			polyA, ideal, ideal-polyA, discBound, n)
	}

	// ---- 3. 十档半径扫描：缺口像素 ----
	fmt.Println()
	fmt.Println("================ 3. 十档半径扫描：相对直折信封的缺口像素 ================")
	gapV := geom.Point{X: 40, Y: 0}
	gapT0 := geom.Point{X: 1, Y: 0}
	gapRad := 40.0 * math.Pi / 180
	gapT1 := geom.Point{X: math.Cos(gapRad), Y: math.Sin(gapRad)}
	fmt.Printf("拐角 α=40°，半线宽 h=%.1f，p=%.2f；缺口=直折覆盖而当前接头未覆盖的网格中心数\n", h, *pixel)
	fmt.Println("档位  半径r    圆角缺口  斜切缺口")
	var gaps []gapRow
	prevRound, prevBevel := math.MaxInt, math.MaxInt
	for k := 1; k <= 10; k++ {
		r := h * float64(k) / 10
		rg, _, rpoly, err := geom.JoinGapPixels(gapV, gapT0, gapT1, h, r, geom.JoinRound, arcSegs, *pixel)
		if err != nil {
			panic(err)
		}
		bg, _, bpoly, err := geom.JoinGapPixels(gapV, gapT0, gapT1, h, r, geom.JoinBevel, arcSegs, *pixel)
		if err != nil {
			panic(err)
		}
		if rg > prevRound || bg > prevBevel {
			panic(fmt.Sprintf("第 %d 档缺口非单调: round %d→%d bevel %d→%d", k, prevRound, rg, prevBevel, bg))
		}
		prevRound, prevBevel = rg, bg
		gaps = append(gaps, gapRow{k, r, rg, bg, rpoly, bpoly})
		fmt.Printf("%3d  %6.3f  %8d  %8d\n", k, r, rg, bg)
	}
	// 十档半径可视化拼图（圆角）
	montagePath := filepath.Join(*outDir, "radius_scan_round.png")
	writeMontage(montagePath, gaps, h)

	// ---- 4. 半径连成曲线：无缝验证 + 拐点 ----
	fmt.Println()
	fmt.Println("================ 4. 曲线半径 profile：无缝与拐点 ================")
	curvePath := zigzagPath(35, 80) // 7 个内部拐点
	radii := make([]float64, len(curvePath)-2)
	for i := range radii {
		radii[i] = h * (0.55 + 0.45*math.Sin(float64(i)*math.Pi/3.5))
	}
	turns := geom.RadiusTurns(radii)
	fmt.Printf("半径 profile（%d 个拐点）: ", len(radii))
	for i, r := range radii {
		fmt.Printf("%.2f", r)
		if i < len(radii)-1 {
			fmt.Print(" → ")
		}
	}
	fmt.Println()
	fmt.Printf("曲线拐点（半径序列下标，对应路径顶点 i+1）: %v\n", turns)
	for _, st := range []struct {
		name string
		s    geom.JoinStyle
	}{
		{"圆角", geom.JoinRound}, {"斜切", geom.JoinBevel},
	} {
		params := geom.StrokeParams{Width: 2 * h, Style: st.s, Radius: radii, ArcSegs: arcSegs}
		r, err := geom.Stroke(curvePath, params, *pixel)
		if err != nil {
			panic(err)
		}
		comps := r.Components()
		fmt.Printf("%s 曲线描边：覆盖像素=%d，4 邻域连通分量=%d（1 = 无缝隙）\n", st.name, r.Count(), comps)
		if comps != 1 {
			panic(fmt.Sprintf("%s 曲线半径描边出现缝隙", st.name))
		}
		writePNG(filepath.Join(*outDir, "curve_radius_"+st.name+".png"), r)
	}

	// ---- 5. 错误分支抽查 ----
	fmt.Println()
	fmt.Println("================ 5. 错误分支（均带拐点号） ================")
	errPath := []geom.Point{{X: 0, Y: 0}, {X: 20, Y: 0}, {X: 20, Y: 20}}
	check := func(name string, p geom.StrokeParams) {
		_, err := geom.StrokePolygons(errPath, p)
		fmt.Printf("%-12s → %v\n", name, err)
	}
	check("半径为零", geom.StrokeParams{Width: 10, Style: geom.JoinRound, Radius: []float64{0}})
	check("样式非法", geom.StrokeParams{Width: 10, Style: geom.JoinStyle(7)})
	collinear := []geom.Point{{X: 0, Y: 0}, {X: 10, Y: 0}, {X: 20, Y: 0}}
	_, err := geom.StrokePolygons(collinear, geom.StrokeParams{Width: 10, Style: geom.JoinBevel, Radius: []float64{5}})
	fmt.Printf("%-12s → %v\n", "两边共线", err)

	fmt.Println()
	fmt.Println("报告完成，图像输出于", *outDir)
}

// zigzagPath 生成 n 段长度交替上下的等长折线，拐角角度由 deg 控制。
func zigzagPath(deg, length float64) []geom.Point {
	rad := deg * math.Pi / 180
	pts := []geom.Point{{X: 0, Y: 0}}
	dir := geom.Point{X: 1, Y: 0}
	up := true
	// 每次把方向旋转 ±(180-deg)，使相邻边夹角恰为 deg
	turn := math.Pi - rad
	for len(pts) < 8 {
		if up {
			dir = rotateVec(dir, turn)
		} else {
			dir = rotateVec(dir, -turn)
		}
		last := pts[len(pts)-1]
		pts = append(pts, geom.Point{X: last.X + length*dir.X, Y: last.Y + length*dir.Y})
		up = !up
	}
	return pts
}

func rotateVec(v geom.Point, a float64) geom.Point {
	c, s := math.Cos(a), math.Sin(a)
	return geom.Point{X: v.X*c - v.Y*s, Y: v.X*s + v.Y*c}
}

func perimeter(p []geom.Point) float64 {
	var l float64
	for i := range p {
		a, b := p[i], p[(i+1)%len(p)]
		l += math.Hypot(b.X-a.X, b.Y-a.Y)
	}
	return l
}

func polygonArea(p []geom.Point) float64 {
	var s float64
	for i := range p {
		a, b := p[i], p[(i+1)%len(p)]
		s += a.X*b.Y - b.X*a.Y
	}
	return math.Abs(s / 2)
}

// writePNG 把覆盖网格渲染为 PNG（深色背景 + 亮色像素）。
func writePNG(path string, r *geom.Raster) {
	img := image.NewGray(image.Rect(0, 0, r.W, r.H))
	bg := uint8(20)
	for i := range img.Pix {
		img.Pix[i] = bg
	}
	for j := 0; j < r.H; j++ {
		for i := 0; i < r.W; i++ {
			if r.Cover[j*r.W+i] {
				img.Pix[j*r.W+i] = 235
			}
		}
	}
	f, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	if err := png.Encode(f, img); err != nil {
		panic(err)
	}
}

// writeMontage 把十档半径的圆角接头横向拼成一张 PNG。
func writeMontage(path string, gaps []gapRow, h float64) {
	const cell = 90
	n := len(gaps)
	img := image.NewGray(image.Rect(0, 0, cell*n, cell))
	for i := range img.Pix {
		img.Pix[i] = 20
	}
	white := color.Gray{Y: 235}
	for idx, g := range gaps {
		// 以拐点 (40,0) 为中心，取半径 h+余量 的窗口
		const cx, cy, half = 40.0, 0.0, 7.0
		for py := 0; py < cell; py++ {
			for px := 0; px < cell; px++ {
				wx := cx - half + 2*half*float64(px)/cell
				wy := cy - half + 2*half*float64(py)/cell
				if geom.InsidePoly(g.roundPoly, wx, wy) {
					img.SetGray(idx*cell+px, py, white)
				}
			}
		}
	}
	f, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	if err := png.Encode(f, img); err != nil {
		panic(err)
	}
}
