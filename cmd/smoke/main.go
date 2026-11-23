package main

import (
	"fmt"
	"math"

	"gogeom/geom"
)

func ring(n int, fn func(i int) geom.Point) []geom.Point {
	out := make([]geom.Point, n)
	for i := 0; i < n; i++ {
		out[i] = fn(i)
	}
	return out
}

func show(name string, pts []geom.Point, holes ...[]geom.Point) {
	var rep *geom.TriReport
	var err error
	if len(holes) > 0 {
		rep, err = geom.TriangulateRingsReport(pts, holes)
	} else {
		rep, err = geom.TriangulateReport(pts)
	}
	if err != nil {
		fmt.Printf("%-18s ERROR %v\n", name, err)
		if ve, ok := err.(*geom.VertexError); ok {
			fmt.Printf("%-18s   verts=%v\n", name, ve.Vertices)
		}
		return
	}
	fmt.Printf("%-18s tris=%-3d rounds=%-3d est=%-3d faces=%d overlap=%d areaOrig=%.6f areaTri=%.6f diff=%.2e rel=%.2e\n",
		name, len(rep.Triangles), rep.Rounds, rep.EstRounds, rep.SplitFaces, rep.OverlapPairs,
		rep.OriginalArea, rep.TotalTriArea, rep.AreaDiff, rep.AreaRelDiff)
	fmt.Printf("%-18s note=%s removed=%v bridges=%v\n", "", rep.RoundsNote, rep.CollinearRemoved, rep.Bridges)
}

func main() {
	// 0 方块（凸）
	sq := []geom.Point{{X: 40, Y: 40}, {X: 320, Y: 40}, {X: 320, Y: 320}, {X: 40, Y: 320}}
	show("square", sq)

	// 1 凸八边形
	oct := ring(8, func(i int) geom.Point {
		a := float64(i) * math.Pi / 4
		return geom.Point{X: 180 + 140*math.Cos(a), Y: 180 + 140*math.Sin(a)}
	})
	show("convex-octagon", oct)

	// 2 L 形（凹）
	L := []geom.Point{{X: 40, Y: 320}, {X: 40, Y: 40}, {X: 150, Y: 40}, {X: 150, Y: 180}, {X: 320, Y: 180}, {X: 320, Y: 320}}
	show("concave-L", L)

	// 3 八角星（简单凹多边形，非自交）
	star8 := ring(16, func(i int) geom.Point {
		a := float64(i) * math.Pi / 8
		r := 140.0
		if i%2 == 1 {
			r = 60
		}
		return geom.Point{X: 180 + r*math.Cos(a), Y: 180 + r*math.Sin(a)}
	})
	show("star8-simple", star8)

	// 4 单洞方块
	outer := []geom.Point{{X: 40, Y: 40}, {X: 320, Y: 40}, {X: 320, Y: 320}, {X: 40, Y: 320}}
	hole := []geom.Point{{X: 130, Y: 130}, {X: 130, Y: 230}, {X: 230, Y: 230}, {X: 230, Y: 130}}
	show("one-hole", outer, hole)

	// 5 双洞
	o5 := []geom.Point{{X: 20, Y: 110}, {X: 340, Y: 110}, {X: 340, Y: 250}, {X: 20, Y: 250}}
	h1 := []geom.Point{{X: 50, Y: 140}, {X: 50, Y: 220}, {X: 120, Y: 220}, {X: 120, Y: 140}}
	h2 := []geom.Point{{X: 200, Y: 140}, {X: 200, Y: 220}, {X: 300, Y: 220}, {X: 300, Y: 140}}
	show("two-hole", o5, h1, h2)

	// 6 圆环（48 边外环 + 32 边洞）
	ringOut := ring(48, func(i int) geom.Point {
		a := float64(i) * 2 * math.Pi / 48
		return geom.Point{X: 180 + 150*math.Cos(a), Y: 180 + 150*math.Sin(a)}
	})
	ringIn := ring(32, func(i int) geom.Point {
		a := float64(i) * 2 * math.Pi / 32
		return geom.Point{X: 180 + 70*math.Cos(a), Y: 180 + 70*math.Sin(a)}
	})
	show("annulus", ringOut, ringIn)

	// 7 蝴蝶结（自交）
	bow := []geom.Point{{X: 40, Y: 40}, {X: 320, Y: 320}, {X: 320, Y: 40}, {X: 40, Y: 320}}
	show("bowtie", bow)

	// 8 五角星（自交）
	pent := ring(10, func(i int) geom.Point {
		a := float64(i)*math.Pi/5 - math.Pi/2
		r := 140.0
		if i%2 == 1 {
			r = 55
		}
		return geom.Point{X: 180 + r*math.Cos(a), Y: 180 + r*math.Sin(a)}
	})
	show("pentagram", pent)

	// 9 近似共线方块（一个顶点略微凸起）
	near := []geom.Point{{X: 40, Y: 40}, {X: 320, Y: 40}, {X: 320, Y: 320}, {X: 180, Y: 320 + 1e-7}, {X: 40, Y: 320}}
	show("near-collinear", near)

	// 10 全共线（退化）
	line := []geom.Point{{X: 40, Y: 180}, {X: 120, Y: 180}, {X: 200, Y: 180}, {X: 300, Y: 180}}
	show("all-collinear", line)

	// 错误：顶点不足
	show("too-few", []geom.Point{{X: 0, Y: 0}, {X: 1, Y: 1}})

	// 共线重叠边自交
	bad := []geom.Point{{X: 0, Y: 0}, {X: 10, Y: 0}, {X: 10, Y: 10}, {X: 5, Y: 0}, {X: 0, Y: 10}}
	show("overlap-edge", bad)
}
