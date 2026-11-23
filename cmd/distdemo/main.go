// distdemo 计算五组样本曲线的最近距离，打印最近点对与收敛步数。
// 运行：go run ./cmd/distdemo
package main

import (
	"fmt"

	"gogeom/geom"
)

func main() {
	fmt.Printf("%-6s %-12s %-24s %-24s %s\n", "样本", "最近距离", "曲线A最近点", "曲线B最近点", "收敛步数")
	for _, sp := range geom.DistanceSamples() {
		r, err := geom.CurveDistance(sp.A, sp.B)
		if err != nil {
			fmt.Printf("%-6s 错误: %v\n", sp.Name, err)
			continue
		}
		fmt.Printf("%-6s %-12.9f (% .9f, % .9f)   (% .9f, % .9f)   %d\n",
			sp.Name, r.Dist, r.PA.X, r.PA.Y, r.PB.X, r.PB.Y, r.Steps)
	}
}
