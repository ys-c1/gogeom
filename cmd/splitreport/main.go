// splitreport 对等面积切分跑全部用例并打印对照表：
// 各块面积、极差、核验结果，以及每个用例一行
// “各块面积之和与原面积的差”。
package main

import (
	"fmt"
	"strings"

	"gogeom/geom"
)

func pts(xy ...float64) []geom.Point {
	out := make([]geom.Point, len(xy)/2)
	for i := range out {
		out[i] = geom.Point{X: xy[2*i], Y: xy[2*i+1]}
	}
	return out
}

func main() {
	type kase struct {
		name string
		poly geom.Polygon
		n    int
	}
	cases := []kase{
		{"凸多边形(正方形10x10)", geom.Polygon{Outer: pts(0, 0, 10, 0, 10, 10, 0, 10)}, 4},
		{"凹多边形(L形)", geom.Polygon{Outer: pts(0, 0, 4, 0, 4, 2, 2, 2, 2, 4, 0, 4)}, 4},
		{"带洞(方框含2x2孔)", geom.Polygon{
			Outer: pts(0, 0, 10, 0, 10, 10, 0, 10),
			Holes: [][]geom.Point{pts(4, 4, 4, 6, 6, 6, 6, 4)},
		}, 3},
		{"细长(1000x1)", geom.Polygon{Outer: pts(0, 0, 1000, 0, 1000, 1, 0, 1)}, 5},
		{"退化(三点共线)", geom.Polygon{Outer: pts(0, 0, 1, 0, 2, 0)}, 3},
		{"块数为零", geom.Polygon{Outer: pts(0, 0, 10, 0, 10, 10, 0, 10)}, 0},
		{"自交(蝴蝶结)", geom.Polygon{Outer: pts(0, 0, 2, 2, 2, 0, 0, 2)}, 2},
	}

	fmt.Println("| 用例 | 块数 | 原面积 | 目标块面积 | 各块面积 | 极差 | 核验 |")
	fmt.Println("|---|---|---|---|---|---|---|")
	for _, c := range cases {
		res, err := geom.SplitEqualArea(c.poly, c.n)
		if err != nil {
			fmt.Printf("| %s | %d | - | - | 报错：%v | - | - |\n", c.name, c.n, err)
			continue
		}
		areas := make([]string, len(res.Pieces))
		for i, pc := range res.Pieces {
			areas[i] = fmt.Sprintf("%.6f", pc.Area)
		}
		verdict := "通过"
		if verr := res.Verify(1e-6); verr != nil {
			verdict = "失败: " + verr.Error()
		}
		fmt.Printf("| %s | %d | %.6f | %.6f | %s | %.3e | %s |\n",
			c.name, c.n, res.Total, res.Total/float64(c.n),
			strings.Join(areas, ", "), res.Spread, verdict)
		fmt.Printf("[%s] 各块面积之和与原面积的差 = %+.3e\n", c.name, res.SumDiff)
	}
}
