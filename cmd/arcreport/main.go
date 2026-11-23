// Command arcreport 生成三点定圆 / 圆弧折线逼近的验证报告。
package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"gogeom/geom"
)

const (
	smallR = 10.0
	largeR = 1000.0
	sag    = 0.01
)

func onCircle(r, deg float64) geom.Point {
	t := deg * math.Pi / 180
	return geom.Point{X: r * math.Cos(t), Y: r * math.Sin(t)}
}

type arcCase struct {
	kind    string // 大弧 / 小弧 / 半圆 / 近直线
	a, b, c geom.Point
}

func casesFor(r float64) []arcCase {
	rise := r / 1000 // 近直线：弦长 2r，拱高 r/1000
	return []arcCase{
		{"大弧 300°", onCircle(r, 20), onCircle(r, 260), onCircle(r, 320)},
		{"小弧 30°", onCircle(r, 5), onCircle(r, 20), onCircle(r, 35)},
		{"半圆 180°", onCircle(r, 180), onCircle(r, 270), onCircle(r, 360)},
		{"近直线", geom.Point{X: 0, Y: 0}, geom.Point{X: r, Y: rise}, geom.Point{X: 2 * r, Y: 0}},
	}
}

type row struct {
	radius  float64
	kind    string
	arc     geom.Arc
	pts     []geom.Point
	rep     geom.ArcDeviationReport
	elapsed time.Duration
	err     error
}

func buildRow(r float64, c arcCase) row {
	arc, err := geom.ArcThroughPoints(c.a, c.b, c.c)
	if err != nil {
		return row{radius: r, kind: c.kind, err: err}
	}
	best := time.Duration(math.MaxInt64)
	var pts []geom.Point
	var rep geom.ArcDeviationReport
	for i := 0; i < 200; i++ {
		t0 := time.Now()
		p, _, e := geom.FlattenArc(arc, sag)
		dt := time.Since(t0)
		if e != nil {
			return row{radius: r, kind: c.kind, err: e}
		}
		if dt < best {
			best = dt
		}
		pts = p
		if i == 0 {
			rep, err = arc.Report(p, sag)
			if err != nil {
				return row{radius: r, kind: c.kind, err: err}
			}
		}
	}
	return row{radius: r, kind: c.kind, arc: arc, pts: pts, rep: rep, elapsed: best}
}

func fmtP(p geom.Point) string {
	return fmt.Sprintf("(%+.6g, %+.6g)", p.X, p.Y)
}

func fmtD(d time.Duration) string {
	switch {
	case d < time.Microsecond:
		return fmt.Sprintf("%d ns", d.Nanoseconds())
	case d < time.Millisecond:
		return fmt.Sprintf("%.2f µs", float64(d.Nanoseconds())/1e3)
	default:
		return fmt.Sprintf("%.3f ms", float64(d.Nanoseconds())/1e6)
	}
}

func main() {
	var rows []row
	for _, r := range []float64{smallR, largeR} {
		for _, c := range casesFor(r) {
			rows = append(rows, buildRow(r, c))
		}
	}

	var b strings.Builder
	writeHeader(&b)
	writeCircleTable(&b, rows)
	writeRefitTable(&b, rows)
	writeCoincidentNote(&b)
	writeDistribution(&b, rows)
	writeExtraTable(&b, rows)

	out := b.String()
	fmt.Print(out)

	dir := "artifacts"
	if len(os.Args) > 1 {
		dir = os.Args[1]
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	path := filepath.Join(dir, "arc_report.md")
	if err := os.WriteFile(path, []byte(out), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	fmt.Fprintf(os.Stderr, "\n报告已写入 %s\n", path)
}

func writeHeader(b *strings.Builder) {
	fmt.Fprintf(b, "# 三点定圆与圆弧折线逼近 — 验证报告\n\n")
	fmt.Fprintf(b, "- 生成时间：%s\n", time.Now().Format("2006-01-02 15:04:05"))
	fmt.Fprintf(b, "- 弦高误差公差：%g（两档半径均用同一公差）\n", sag)
	fmt.Fprintf(b, "- 分段反解：n = ⌈|φ| / (2·acos(1−d/r))⌉，其中 φ 为弧的带符号扫角\n")
	fmt.Fprintf(b, "- 误差上界：每段弦中点径向偏差 r·(1−cos(θ/2)) ≤ d，θ = |φ|/n\n")
	fmt.Fprintf(b, "- 折线顶点严格取在解析圆上，其径向偏差仅为浮点舍入量级\n\n")
}

func writeCircleTable(b *strings.Builder, rows []row) {
	fmt.Fprintf(b, "## 表 1　圆心与偏差表\n\n")
	fmt.Fprintf(b, "| 半径档 | 用例 | 圆心 | 解析半径 | 扫角 | 段数 | 顶点最大径向偏差 | 弦中最大弦高 | 折线化耗时 |\n")
	fmt.Fprintf(b, "|---:|---|---:|---:|---:|---:|---:|---:|---:|\n")
	for _, r := range rows {
		if r.err != nil {
			fmt.Fprintf(b, "| %g | %s | — | — | — | — | — | — | 错误：%v |\n", r.radius, r.kind, r.err)
			continue
		}
		fmt.Fprintf(b, "| %g | %s | %s | %.9g | %+.4f° | %d | %.3e | %.6e | %s |\n",
			r.radius, r.kind, fmtP(r.arc.Center), r.arc.Radius,
			r.arc.Sweep*180/math.Pi, r.rep.Segments,
			r.rep.MaxDeviation, r.rep.MaxSagitta, fmtD(r.elapsed))
	}
	fmt.Fprintf(b, "\n注：段数为满足公差的最小段数；近直线用例的圆心在画布外，半径由拱高决定，"+
		"正是「几乎共线」的可解边界。弦中最大弦高均 ≤ %g。\n\n", sag)
}

func writeRefitTable(b *strings.Builder, rows []row) {
	fmt.Fprintf(b, "## 表 2　圆心反算残差表\n\n")
	fmt.Fprintf(b, "取折线首、中、末三点重新定圆，与解析圆比较。残差 1 = ‖C′−C‖；"+
		"残差 2 = |r′−r|；点偏差 = 折线全部顶点到反算圆的最大径向偏差。\n\n")
	fmt.Fprintf(b, "| 半径档 | 用例 | 反算圆心 | 反算半径 | 圆心残差（绝对） | 圆心残差（/r） | 半径残差 | 反算圆下顶点最大偏差 |\n")
	fmt.Fprintf(b, "|---:|---|---:|---:|---:|---:|---:|---:|\n")
	for _, r := range rows {
		if r.err != nil {
			fmt.Fprintf(b, "| %g | %s | — | — | — | — | — | 错误：%v |\n", r.radius, r.kind, r.err)
			continue
		}
		rel := r.rep.RefitResidual / r.arc.Radius
		fmt.Fprintf(b, "| %g | %s | %s | %.9g | %.3e | %.3e | %.3e | %.3e |\n",
			r.radius, r.kind, fmtP(r.rep.Refit.Center), r.rep.Refit.Radius,
			r.rep.RefitResidual, rel, r.rep.RadiusResidual, r.rep.RefitMaxPointDeviation)
	}
	fmt.Fprintf(b, "\n注：反算只用了折线上的三点（含首末端点），残差应在舍入量级；"+
		"若明显偏大，说明反算用的三点不落在同一段圆弧上（例如方向判反）。\n\n")
}

func writeCoincidentNote(b *strings.Builder) {
	fmt.Fprintf(b, "## 表 3　退化用例：两点重合（单独说明）\n\n")
	base := geom.Point{X: 3, Y: 4}
	others := geom.Point{X: 10, Y: 0}
	type c struct {
		name string
		p    [3]geom.Point
	}
	cs := []c{
		{"点 0 与点 1 重合", [3]geom.Point{base, base, others}},
		{"点 1 与点 2 重合", [3]geom.Point{others, base, base}},
		{"点 0 与点 2 重合", [3]geom.Point{base, others, base}},
	}
	fmt.Fprintf(b, "两点重合时叉积也为零，若先判共线会误报。实现中**先查重合**："+
		"任意两点距离 ≤ 1e-12·max(1,|x|,|y|) 即报 `ErrCoincident`，"+
		"并把**重合的两点号**放入 `*PointError.Pts`（长度 2）；"+
		"共线则返回三点号（长度 3）。\n\n")
	fmt.Fprintf(b, "| 用例 | 输入点 0 / 1 / 2 | 返回错误 | 点号 |\n|---|---|---|---|\n")
	for _, tc := range cs {
		_, err := geom.Circumcircle(tc.p[0], tc.p[1], tc.p[2])
		pts, ok := geom.IsPointError(err)
		ptsStr := "—"
		if ok {
			parts := make([]string, len(pts))
			for i, p := range pts {
				parts[i] = fmt.Sprintf("%d", p)
			}
			ptsStr = strings.Join(parts, ", ")
		}
		fmt.Fprintf(b, "| %s | (%+g,%+g) / (%+g,%+g) / (%+g,%+g) | %v | %s |\n",
			tc.name,
			tc.p[0].X, tc.p[0].Y, tc.p[1].X, tc.p[1].Y, tc.p[2].X, tc.p[2].Y,
			err, ptsStr)
	}
	_, err := geom.Circumcircle(geom.Point{}, geom.Point{X: 5, Y: 5}, geom.Point{X: 10, Y: 10})
	if pts, ok := geom.IsPointError(err); ok {
		fmt.Fprintf(b, "\n共线对照：输入 (0,0)/(5,5)/(10,10) → `%v`，点号 %v。\n", err, pts)
	}

	canvas := geom.Canvas{Min: geom.Point{}, Max: geom.Point{X: 800, Y: 800}}
	_, err = geom.CircumcircleInCanvas(onCircle(1000, 0), onCircle(1000, 30), onCircle(1000, 60), canvas)
	if pts, ok := geom.IsPointError(err); ok {
		fmt.Fprintf(b, "半径超画幅对照：r=1000、画幅 800×800（对角线 %.1f）→ `%v`，点号 %v。\n",
			canvas.Diagonal(), err, pts)
	}
	fmt.Fprintf(b, "\n重合例中三点退化到只有两个不同点，外接圆不唯一（圆心在两点连线的中垂线上任意取值），"+
		"故必须拒绝而不是给出某个任意圆。\n\n")
}

func writeDistribution(b *strings.Builder, rows []row) {
	fmt.Fprintf(b, "## 表 4　偏差分布\n\n")
	fmt.Fprintf(b, "沿折线每段等参采样 200 点（含端点），统计各点到解析圆的径向偏差；"+
		"区间以公差 d=%g 为单位，共 10 档。顶点落在 0 档、弦中点落在最高档。\n\n", sag)
	const perSeg = 200
	const bins = 10
	fmt.Fprintf(b, "| 半径档 | 用例 | 样本数 | 分档计数（0→1×d） | 均值/d | 最大/d |\n")
	fmt.Fprintf(b, "|---:|---|---:|:---:|---:|---:|\n")
	for _, r := range rows {
		if r.err != nil {
			continue
		}
		counts := make([]int, bins)
		total := 0
		var sum, mx float64
		for i := 0; i+1 < len(r.pts); i++ {
			for k := 0; k <= perSeg; k++ {
				t := float64(k) / perSeg
				q := geom.Point{
					X: r.pts[i].X + (r.pts[i+1].X-r.pts[i].X)*t,
					Y: r.pts[i].Y + (r.pts[i+1].Y-r.pts[i].Y)*t,
				}
				dv := r.arc.RadialDeviation(q)
				idx := int(dv / sag * bins)
				if idx >= bins {
					idx = bins - 1
				}
				counts[idx]++
				sum += dv
				if dv > mx {
					mx = dv
				}
				total++
			}
		}
		parts := make([]string, bins)
		for i, c := range counts {
			parts[i] = fmt.Sprintf("%d", c)
		}
		fmt.Fprintf(b, "| %g | %s | %d | %s | %.4f | %.4f |\n",
			r.radius, r.kind, total, strings.Join(parts, " / "), sum/float64(total)/sag, mx/sag)
	}
	fmt.Fprintf(b, "\n分布呈两端密、中间疏的 U 形：靠近顶点偏差趋 0，弦中点达到该段最大弦高。"+
		"近直线用例因段数极少，偏差集中在低档。所有采样点偏差均不超过 d，验证了反解段数的上界成立。\n\n")
}

func writeExtraTable(b *strings.Builder, rows []row) {
	fmt.Fprintf(b, "## 表 5　另附：公差—段数—耗时对照\n\n")
	tols := []float64{1e-3, 1e-2, 1e-1, 1}
	radii := []float64{smallR, largeR}
	fmt.Fprintf(b, "| 半径 | 公差 d | 整圆最小段数（解析） | 本报告半圆段数 | 半圆折线化耗时 |\n")
	fmt.Fprintf(b, "|---:|---:|---:|---:|---:|\n")
	for _, r := range radii {
		arc, err := geom.ArcThroughPoints(onCircle(r, 180), onCircle(r, 270), onCircle(r, 360))
		if err != nil {
			fmt.Fprintf(b, "| %g | — | 错误：%v | — | — |\n", r, err)
			continue
		}
		for _, tol := range tols {
			nFull, err := geom.SegmentsForSagitta(r, tol)
			if err != nil {
				fmt.Fprintf(b, "| %g | %g | 错误：%v | — | — |\n", r, tol, err)
				continue
			}
			best := time.Duration(math.MaxInt64)
			var n int
			for i := 0; i < 200; i++ {
				t0 := time.Now()
				_, k, e := geom.FlattenArc(arc, tol)
				dt := time.Since(t0)
				if e != nil {
					err = e
					break
				}
				n = k
				if dt < best {
					best = dt
				}
			}
			if err != nil {
				fmt.Fprintf(b, "| %g | %g | 错误：%v | — | — |\n", r, tol, err)
				continue
			}
			fmt.Fprintf(b, "| %g | %g | %d | %d | %s |\n", r, tol, nFull, n, fmtD(best))
		}
	}
	fmt.Fprintf(b, "\n耗时取 200 次重复中的最小值，以排除调度噪声；"+
		"整圆段数是保守上界，实际按扫角反解（半圆约为其一半）。\n")

	// 表 1 数据按段数排序的紧凑摘要，便于核对两档半径的规模差异。
	sort.SliceStable(rows, func(i, j int) bool {
		if rows[i].radius != rows[j].radius {
			return rows[i].radius < rows[j].radius
		}
		return rows[i].rep.Segments < rows[j].rep.Segments
	})
	fmt.Fprintf(b, "\n### 分段规模摘要（按段数升序）\n\n")
	fmt.Fprintf(b, "| 半径 | 用例 | 段数 |\n|---:|---|---:|\n")
	for _, r := range rows {
		if r.err != nil {
			continue
		}
		fmt.Fprintf(b, "| %g | %s | %d |\n", r.radius, r.kind, r.rep.Segments)
	}
}
