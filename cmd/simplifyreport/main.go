// Command simplifyreport 生成路径简化在十档点数下的数据表。
//
// 用法：go run ./cmd/simplifyreport [-out docs/simplify_report.md]
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"

	"gogeom/geom"
)

var tiers = []int{100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000}

type row struct {
	kind string
	tol  float64
	n    int
	res  *geom.Result
	err  error
}

// wavyPath 是开放的波浪路径：弯曲处需要保留点，缓弯处可大量抽稀。
func wavyPath(n int) []geom.Point {
	pts := make([]geom.Point, n)
	for i := range pts {
		x := 1000 * float64(i) / float64(n-1)
		y := 20*math.Sin(x/40) + 5*math.Sin(x/7)
		pts[i] = geom.Point{X: x, Y: y}
	}
	return pts
}

// circlePath 是闭合环，共 n 个点且首尾逐位相同。
func circlePath(n int) []geom.Point {
	if n < 3 {
		n = 3
	}
	pts := make([]geom.Point, n)
	for i := 0; i < n-1; i++ {
		a := 2 * math.Pi * float64(i) / float64(n-1)
		pts[i] = geom.Point{X: 500 + 400*math.Cos(a), Y: 500 + 400*math.Sin(a)}
	}
	pts[n-1] = pts[0]
	return pts
}

func run(kind string, tol float64, makePts func(int) []geom.Point) []row {
	rows := make([]row, 0, len(tiers))
	for _, n := range tiers {
		res, err := geom.Simplify(makePts(n), tol)
		rows = append(rows, row{kind: kind, tol: tol, n: n, res: res, err: err})
	}
	return rows
}

func main() {
	out := flag.String("out", "docs/simplify_report.md", "输出文件；- 表示只打印到标准输出")
	flag.Parse()

	const (
		openTol   = 0.10
		closedTol = 0.10
	)
	open := run("开放波浪", openTol, wavyPath)
	closed := run("闭合圆环", closedTol, circlePath)

	var b strings.Builder
	b.WriteString("# 路径简化十档数据表\n\n")
	b.WriteString("由 `go run ./cmd/simplifyreport` 生成。每次运行耗时会随机器负载浮动。\n\n")
	b.WriteString("抽稀算法：Douglas–Peucker，按最大偏差控制；逐点核对原始路径到简化路径的最小距离。\n\n")
	b.WriteString("## 开放波浪路径（容差 ")
	b.WriteString(fmt.Sprintf("%.2f", openTol))
	b.WriteString("）\n\n")
	writeTable(&b, open)
	b.WriteString("\n## 闭合圆环（容差 ")
	b.WriteString(fmt.Sprintf("%.2f", closedTol))
	b.WriteString("）\n\n")
	writeTable(&b, closed)
	b.WriteString("\n")
	writeTwoDensities(&b, open, closed)
	writeSummary(&b, open, closed)

	text := b.String()
	if *out == "-" {
		fmt.Print(text)
		return
	}
	if err := os.MkdirAll(filepath.Dir(*out), 0o755); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if err := os.WriteFile(*out, []byte(text), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	fmt.Println(text)
}

func writeTable(b *strings.Builder, rows []row) {
	b.WriteString("| 档 | 原始点数 | 保留点数 | 保留率 | 最大偏差 | 闭环首尾差 | 抽稀耗时 | 核对耗时 |\n")
	b.WriteString("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
	for i, r := range rows {
		if r.err != nil {
			b.WriteString(fmt.Sprintf("| %d | %d | — | — | — | — | — | 错误：%v |\n", i+1, r.n, r.err))
			continue
		}
		gap := "—"
		if r.res.Closed {
			gap = fmt.Sprintf("%.3e", r.res.ClosureGap)
		}
		b.WriteString(fmt.Sprintf("| %d | %d | %d | %.2f%% | %.4f | %s | %v | %v |\n",
			i+1, r.n, r.res.RetainedPoints, 100*r.res.RetentionRatio,
			r.res.MaxDeviation, gap, r.res.Elapsed, r.res.VerifyElapsed))
	}
}

// writeTwoDensities 单独报告两档点数（1 千 / 1 万）的保留点数与耗时。
func writeTwoDensities(b *strings.Builder, groups ...[]row) {
	const (
		lo = 1000
		hi = 10000
	)
	b.WriteString("## 两档点数对照（1 千 / 1 万）\n\n")
	b.WriteString("| 数据集 | 原始点数 | 保留点数 | 保留率 | 最大偏差 | 抽稀耗时 | 核对耗时 |\n")
	b.WriteString("|---|---:|---:|---:|---:|---:|---:|\n")
	for _, rows := range groups {
		name := rows[0].kind
		for _, r := range rows {
			if r.n != lo && r.n != hi {
				continue
			}
			if r.err != nil {
				b.WriteString(fmt.Sprintf("| %s | %d | — | — | — | — | 错误：%v |\n", name, r.n, r.err))
				continue
			}
			b.WriteString(fmt.Sprintf("| %s | %d | %d | %.2f%% | %.4f | %v | %v |\n",
				name, r.n, r.res.RetainedPoints, 100*r.res.RetentionRatio,
				r.res.MaxDeviation, r.res.Elapsed, r.res.VerifyElapsed))
		}
	}
	b.WriteString("\n")
}

func writeSummary(b *strings.Builder, open, closed []row) {
	b.WriteString("## 收尾登记\n\n")
	b.WriteString("| 数据集 | 首档保留率 | 末档保留率 | 首档最大偏差 | 末档最大偏差 | 十档最大偏差 | 十档最大偏差是否 ≤ 容差 |\n")
	b.WriteString("|---|---:|---:|---:|---:|---:|:---:|\n")
	register := func(name string, rows []row) {
		first, last := rows[0], rows[len(rows)-1]
		if first.err != nil || last.err != nil {
			b.WriteString(fmt.Sprintf("| %s | — | — | — | — | — | 错误 |\n", name))
			return
		}
		worst := 0.0
		for _, r := range rows {
			if r.err == nil && r.res.MaxDeviation > worst {
				worst = r.res.MaxDeviation
			}
		}
		b.WriteString(fmt.Sprintf("| %s | %.2f%% | %.2f%% | %.4f | %.4f | %.4f | %v |\n",
			name, 100*first.res.RetentionRatio, 100*last.res.RetentionRatio,
			first.res.MaxDeviation, last.res.MaxDeviation, worst,
			worst <= rows[0].tol))
	}
	register("开放波浪", open)
	register("闭合圆环", closed)
}
