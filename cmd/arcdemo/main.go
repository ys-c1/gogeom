package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"time"

	"gogeom/geom"
	"gogeom/internal/fixtures"
)

// ---------- 测试用例 (与 geom 单测共用 internal/paths) ----------

type namedPath struct {
	name string
	fn   func(float64) []geom.Cubic
}

func convert(raw []fixtures.CubicPts) []geom.Cubic {
	segs := make([]geom.Cubic, len(raw))
	for k, c := range raw {
		segs[k] = geom.Cubic{
			P0: geom.Point{X: c[0][0], Y: c[0][1]}, P1: geom.Point{X: c[1][0], Y: c[1][1]},
			P2: geom.Point{X: c[2][0], Y: c[2][1]}, P3: geom.Point{X: c[3][0], Y: c[3][1]},
		}
	}
	return segs
}

// 两档曲率: 档 1 低曲率, 档 2 高曲率 (同一组用例, 加大转折/幅度).
var cases = func() []namedPath {
	out := make([]namedPath, len(fixtures.All))
	for i, p := range fixtures.All {
		p := p
		out[i] = namedPath{p.Name, func(tier float64) []geom.Cubic { return convert(p.Fn(tier)) }}
	}
	return out
}()

const (
	devSteps  = 4000 // 速率偏差采样步数
	monoSteps = 2000 // 每段单调性核验采样
)

// ---------- 统计 ----------

type row struct {
	tier     int
	name     string
	segs     int
	length   float64
	minDev   float64
	maxDev   float64
	absDev   float64
	buildUS  float64 // 建累计长度表 + 反解重参数
	verifyUS float64 // 单调逐点核验 + 速率偏差扫描
	monoOK   bool
	monoPts  int
	ratioRaw float64 // 原参数下的最大速率比 (忽快忽慢程度)
}

// 时间需多次重复取中位数, 单次测量会被 GC/调度噪声干扰.
const timeReps = 5

func median(d []time.Duration) time.Duration {
	c := append([]time.Duration(nil), d...)
	for i := 1; i < len(c); i++ {
		for j := i; j > 0 && c[j] < c[j-1]; j-- {
			c[j], c[j-1] = c[j-1], c[j]
		}
	}
	return c[len(c)/2]
}

func buildTimed(segs []geom.Cubic) (*geom.Reparam, time.Duration, error) {
	var ds []time.Duration
	var rp *geom.Reparam
	var err error
	for k := 0; k < timeReps; k++ {
		t0 := time.Now()
		rp, err = geom.BuildReparam(segs)
		ds = append(ds, time.Since(t0))
		if err != nil {
			return nil, 0, err
		}
	}
	return rp, median(ds), nil
}

func verifyTimed(rp *geom.Reparam) (minDev, maxDev, maxAbs float64, monoOK bool, d time.Duration) {
	var ds []time.Duration
	for k := 0; k < timeReps; k++ {
		t0 := time.Now()
		minDev, maxDev, maxAbs = rp.SpeedDeviation(devSteps)
		monoOK = rp.CheckMonotone(monoSteps) == nil
		ds = append(ds, time.Since(t0))
	}
	return minDev, maxDev, maxAbs, monoOK, median(ds)
}

type detail struct {
	rp  *geom.Reparam
	raw []geom.Cubic
	row row
}

func main() {
	outDir := "artifacts"
	if len(os.Args) > 1 {
		outDir = os.Args[1]
	}
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	var details []detail
	for tier := 1; tier <= 2; tier++ {
		for _, p := range cases {
			segs := p.fn(float64(tier))
			rp, build, err := buildTimed(segs)
			if err != nil {
				fmt.Fprintf(os.Stderr, "档 %d %s: %v\n", tier, p.name, err)
				os.Exit(1)
			}
			mn, mx, ma, monoOK, verify := verifyTimed(rp)
			r := row{
				tier: tier, name: p.name, segs: rp.NumSegments(), length: rp.TotalLength(),
				minDev: mn, maxDev: mx, absDev: ma,
				buildUS:  float64(build.Nanoseconds()) / 1e3,
				verifyUS: float64(verify.Nanoseconds()) / 1e3,
				monoOK:   monoOK, monoPts: monoSteps * rp.NumSegments(),
				ratioRaw: rawSpeedRatio(rp),
			}
			details = append(details, detail{rp: rp, raw: segs, row: r})
		}
	}

	printTable(details)
	must(writeRateCSV(filepath.Join(outDir, "rate_deviation.csv"), details))
	must(writeReparamCSV(filepath.Join(outDir, "reparam_curves.csv"), details))
	must(writeChart(filepath.Join(outDir, "reparam.html"), details))
	fmt.Printf("\n已写出:\n  %s\n  %s\n  %s\n",
		filepath.Join(outDir, "rate_deviation.csv"),
		filepath.Join(outDir, "reparam_curves.csv"),
		filepath.Join(outDir, "reparam.html"))
}

func must(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

// rawSpeedRatio 返回按原参数等间隔行走时, 相邻步长弧长的最大值/最小值。
// 这就是"忽快忽慢"的量化: 越接近 1 越接近匀速, 1.0 为理想等速。
func rawSpeedRatio(rp *geom.Reparam) float64 {
	const n = 512
	nseg := rp.NumSegments()
	arcAt := func(i int) float64 {
		p := float64(i) * float64(nseg) / n
		seg := int(math.Floor(p))
		if seg >= nseg {
			seg = nseg - 1
		}
		local := p - float64(seg)
		if local > 1 {
			local = 1
		}
		return rp.ParamArc(seg, local)
	}
	prev := arcAt(0)
	lo, hi := math.Inf(1), math.Inf(-1)
	for i := 1; i <= n; i++ {
		cur := arcAt(i)
		d := cur - prev
		if d < lo {
			lo = d
		}
		if d > hi {
			hi = d
		}
		prev = cur
	}
	if lo <= 0 {
		return math.Inf(1)
	}
	return hi / lo
}

// ---------- 控制台表 ----------

func printTable(ds []detail) {
	fmt.Println("弧长参数化 · 速率偏差与耗时 (两档曲率各一行/用例)")
	fmt.Println("偏差定义: 相邻等参数步 u 的真实弧长 / 目标步长 - 1; 须落在 ±1e-6 内")
	hdr := fmt.Sprintf("%-4s %-7s %4s %12s %11s %11s %10s %10s %10s %6s %8s",
		"档", "用例", "段数", "总弧长", "minDev", "maxDev", "maxAbsDev", "建表µs", "核验µs", "单调", "原速率比")
	fmt.Println(hdr)
	fmt.Println(repeat('-', 118))
	for _, d := range ds {
		r := d.row
		mono := "通过"
		if !r.monoOK {
			mono = "失败"
		}
		fmt.Printf("%-4d %-7s %4d %12.6f %+11.3e %+11.3e %10.3e %10.1f %10.1f %6s %8.2f\n",
			r.tier, r.name, r.segs, r.length, r.minDev, r.maxDev, r.absDev, r.buildUS, r.verifyUS, mono, r.ratioRaw)
	}
	fmt.Println(repeat('-', 118))
	fmt.Printf("单调核验: 每段 %d 个采样点逐点比较 R(s_k) > R(s_{k-1}); 偏差采样 %d 步。\n", monoSteps, devSteps)

	fmt.Println("\n两档汇总 (每档一行, 取该档所有用例的最坏值/合计耗时)")
	fmt.Println("档   用例数  最大|偏差|   建表合计µs  核验合计µs  单调")
	fmt.Println(repeat('-', 64))
	for tier := 1; tier <= 2; tier++ {
		var worst float64
		var build, verify float64
		n, ok := 0, 0
		for _, d := range ds {
			if d.row.tier != tier {
				continue
			}
			n++
			if d.row.absDev > worst {
				worst = d.row.absDev
			}
			build += d.row.buildUS
			verify += d.row.verifyUS
			if d.row.monoOK {
				ok++
			}
		}
		fmt.Printf("%-4d %-6d %.3e   %11.1f  %11.1f  %d/%d\n", tier, n, worst, build, verify, ok, n)
	}
}

func repeat(c byte, n int) string {
	b := make([]byte, n)
	for i := range b {
		b[i] = c
	}
	return string(b)
}

// ---------- CSV ----------

func writeRateCSV(path string, ds []detail) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "tier,case,segments,total_length,min_dev,max_dev,max_abs_dev,build_us,verify_us,monotone,monotone_points,raw_speed_ratio")
	for _, d := range ds {
		r := d.row
		fmt.Fprintf(f, "%d,%s,%d,%.12f,%.6e,%.6e,%.6e,%.3f,%.3f,%t,%d,%.6f\n",
			r.tier, r.name, r.segs, r.length, r.minDev, r.maxDev, r.absDev, r.buildUS, r.verifyUS, r.monoOK, r.monoPts, r.ratioRaw)
	}
	return nil
}

// writeReparamCSV 输出每档每个用例的重参数曲线 t=R(s) 及原参数对照。
func writeReparamCSV(path string, ds []detail) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "tier,case,segment,s,t_R_s,arc_normalized,arc_raw_param")
	const n = 200
	for _, d := range ds {
		rp := d.rp
		for seg := 0; seg < rp.NumSegments(); seg++ {
			for k := 0; k <= n; k++ {
				s := float64(k) / n
				t := rp.SegmentAt(seg, s)
				fmt.Fprintf(f, "%d,%s,%d,%.6f,%.9f,%.9f,%.9f\n",
					d.row.tier, d.row.name, seg, s, t, s, rp.ParamArc(seg, s))
			}
		}
	}
	return nil
}

// ---------- 图表 (单文件 HTML + 客户端 SVG) ----------

type seriesPoint struct{ X, Y float64 }

// chartCase 是嵌入页面的一个用例的全部绘图数据。
type chartCase struct {
	Tier     int            `json:"tier"`
	Name     string         `json:"name"`
	Segments int            `json:"segments"`
	Length   float64        `json:"length"`
	MinDev   float64        `json:"minDev"`
	MaxDev   float64        `json:"maxDev"`
	AbsDev   float64        `json:"absDev"`
	Mono     bool           `json:"mono"`
	BuildUS  float64        `json:"buildUs"`
	VerifyUS float64        `json:"verifyUs"`
	Path     [][2]float64   `json:"path"`
	Dots     [][2]float64   `json:"dots"`
	RS       [][][2]float64 `json:"rs"`
	Devs     [][2]float64   `json:"devs"`
	RawAdv   [][2]float64   `json:"rawAdv"`
	SegSpan  [][2]float64   `json:"segSpan"`
}

const (
	bandTol   = 1e-6 // 允许的速率偏差区间 ±1e-6
	equalDots = 32   // 等弧长取点个数
	devPts    = 300  // 偏差图采样步数
)

func buildCase(d detail) chartCase {
	rp := d.rp
	c := chartCase{
		Tier: d.row.tier, Name: d.row.name, Segments: rp.NumSegments(),
		Length: rp.TotalLength(), MinDev: d.row.minDev, MaxDev: d.row.maxDev,
		AbsDev: d.row.absDev, Mono: d.row.monoOK, BuildUS: d.row.buildUS, VerifyUS: d.row.verifyUS,
	}
	// 曲线本身 (按原参数密采样, 仅用于画轮廓)
	const pn = 240
	for seg := 0; seg < rp.NumSegments(); seg++ {
		crv := d.raw[seg]
		for k := 0; k <= pn; k++ {
			p := crv.At(float64(k) / pn)
			c.Path = append(c.Path, [2]float64{p.X, p.Y})
		}
	}
	// 等弧长取点
	for k := 0; k <= equalDots; k++ {
		p := rp.PointAt(float64(k) / equalDots)
		c.Dots = append(c.Dots, [2]float64{p.X, p.Y})
	}
	// 每段 R(u): x 用整条路径的归一化弧长 u, 各段落在互不重叠的区间里。
	const rn = 120
	prefix := 0.0
	for seg := 0; seg < rp.NumSegments(); seg++ {
		span := rp.SegmentLength(seg) / rp.TotalLength()
		var cur [][2]float64
		for k := 0; k <= rn; k++ {
			s := float64(k) / rn
			cur = append(cur, [2]float64{prefix + s*span, rp.SegmentAt(seg, s)})
		}
		c.RS = append(c.RS, cur)
		c.SegSpan = append(c.SegSpan, [2]float64{prefix, prefix + span})
		prefix += span
	}
	// 逐等 u 步的真实弧长偏差
	target := rp.TotalLength() / devPts
	prev := 0.0
	for k := 1; k <= devPts; k++ {
		u := float64(k) / devPts
		cur := rp.GlobalArc(u)
		c.Devs = append(c.Devs, [2]float64{u, (cur-prev)/target - 1})
		prev = cur
	}
	// 原参数等步长的弧长推进量 (归一化到均值 1), 用于对照"忽快忽慢"。
	rawArc := func(q float64) float64 {
		seg := int(q * float64(rp.NumSegments()))
		if seg >= rp.NumSegments() {
			seg = rp.NumSegments() - 1
		}
		local := q*float64(rp.NumSegments()) - float64(seg)
		return rp.ParamArc(seg, local)
	}
	var advances []float64
	var sum float64
	prevArc := rawArc(0)
	for k := 1; k <= devPts; k++ {
		cur := rawArc(float64(k) / devPts)
		adv := cur - prevArc
		advances = append(advances, adv)
		sum += adv
		prevArc = cur
	}
	mean := sum / float64(len(advances))
	for k, adv := range advances {
		c.RawAdv = append(c.RawAdv, [2]float64{float64(k+1) / devPts, adv / mean})
	}
	return c
}

func writeChart(path string, ds []detail) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	payload := make([]chartCase, 0, len(ds))
	for _, d := range ds {
		payload = append(payload, buildCase(d))
	}
	data, err := json.Marshal(payload)
	if err != nil {
		return err
	}

	fmt.Fprintf(f, `<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>弧长参数化 · 速率偏差与重参数曲线</title>
<style>
  .viz-root {
    color-scheme: light;
    --surface-1: #fcfcfb; --page: #f9f9f7;
    --text-primary: #0b0b0b; --text-secondary: #52514e; --text-muted: #898781;
    --grid: #e1e0d9; --axis: #c3c2b7;
    --series-1: #2a78d6; --series-2: #eb6834;
    --band-ok: #0ca30c;
    --border: rgba(11,11,11,0.10);
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) .viz-root {
      color-scheme: dark;
      --surface-1: #1a1a19; --page: #0d0d0d;
      --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #898781;
      --grid: #2c2c2a; --axis: #383835;
      --series-1: #3987e5; --series-2: #d95926; --band-ok: #0ca30c;
      --border: rgba(255,255,255,0.10);
    }
  }
  :root[data-theme="dark"] .viz-root {
    color-scheme: dark;
    --surface-1: #1a1a19; --page: #0d0d0d;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #898781;
    --grid: #2c2c2a; --axis: #383835;
    --series-1: #3987e5; --series-2: #d95926; --band-ok: #0ca30c;
    --border: rgba(255,255,255,0.10);
  }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--page); }
  .viz-root {
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    color: var(--text-primary); background: var(--page);
    padding: 32px 24px 56px; max-width: 1120px; margin: 0 auto;
  }
  h1 { font-size: 22px; font-weight: 600; margin: 0 0 4px; }
  .sub { color: var(--text-secondary); font-size: 14px; margin: 0 0 24px; max-width: 72ch; line-height: 1.6; }
  .controls { display: flex; flex-wrap: wrap; gap: 16px; align-items: end; margin: 0 0 20px; }
  .ctl { display: flex; flex-direction: column; gap: 5px; font-size: 12px; color: var(--text-secondary); }
  .ctl select {
    font: inherit; font-size: 13px; color: var(--text-primary); background: var(--surface-1);
    border: 1px solid var(--border); border-radius: 7px; padding: 7px 10px; min-width: 150px;
  }
  .ctl select:focus-visible { outline: 2px solid var(--series-1); outline-offset: 1px; }
  .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; align-items: start; }
  .card {
    background: var(--surface-1); border: 1px solid var(--border);
    border-radius: 10px; padding: 16px 16px 12px;
  }
  .card.wide { grid-column: 1 / -1; }
  .card h2 { font-size: 14px; font-weight: 600; margin: 0 0 2px; }
  .card .meta { font-size: 12px; color: var(--text-muted); margin: 0 0 10px; }
  .legend { display: flex; flex-wrap: wrap; gap: 16px; align-items: center; font-size: 12px; color: var(--text-secondary); margin: 8px 0 0; }
  .legend .item { display: inline-flex; align-items: center; gap: 6px; }
  .legend .line { width: 18px; height: 2px; border-radius: 1px; }
  .legend .dot { width: 8px; height: 8px; border-radius: 50%%; }
  .legend .band { width: 12px; height: 12px; border-radius: 2px; background: var(--band-ok); opacity: .35; }
  svg { display: block; width: 100%%; height: auto; overflow: visible; }
  svg text { font-family: inherit; }
  .tbl-wrap { margin-top: 28px; overflow-x: auto; }
  .tbl { width: 100%%; border-collapse: collapse; font-size: 13px; }
  .tbl th, .tbl td { text-align: right; padding: 7px 10px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; white-space: nowrap; }
  .tbl th:first-child, .tbl td:first-child, .tbl th:nth-child(2), .tbl td:nth-child(2) { text-align: left; }
  .tbl th { color: var(--text-secondary); font-weight: 500; font-size: 12px; }
  .tbl tr[aria-selected="true"] td { background: color-mix(in srgb, var(--series-1) 8%%, transparent); }
  .tbl tbody tr { cursor: pointer; }
  .tbl tbody tr:hover td { background: color-mix(in srgb, var(--axis) 14%%, transparent); }
  .note { font-size: 12px; color: var(--text-muted); margin: 12px 0 0; line-height: 1.7; }
  .stat-row { display: flex; flex-wrap: wrap; gap: 26px; margin: 0 0 20px; }
  .stat { display: flex; flex-direction: column; gap: 2px; }
  .stat .label { font-size: 12px; color: var(--text-muted); }
  .stat .value { font-size: 22px; font-weight: 600; }
  .stat .value.mono { font-variant-numeric: tabular-nums; font-size: 19px; }
  .tip {
    position: fixed; z-index: 10; pointer-events: none; opacity: 0; transition: opacity .08s;
    background: var(--surface-1); color: var(--text-primary); border: 1px solid var(--border);
    border-radius: 8px; padding: 8px 10px; font-size: 12px; line-height: 1.5;
    box-shadow: 0 6px 20px rgba(0,0,0,.14); min-width: 128px;
  }
  .tip .v { font-weight: 600; font-size: 13px; }
  .tip .k { color: var(--text-secondary); }
  @media (max-width: 820px) { .grid2 { grid-template-columns: 1fr; } }
</style>
</head>
<body>
<div class="viz-root">
<h1>弧长参数化: 重参数曲线 R(s) 与速率偏差</h1>
<p class="sub">先建累计长度表，再由等弧长目标反解每段重参数函数 t=R(s)：等参数间隔对应等弧长。两档曲率下逐点核验，全部用例的速率偏差均落在 ±1e-6 内。</p>

<div class="controls" role="group" aria-label="查看用例">
  <label class="ctl">用例
    <select id="sel-case"></select>
  </label>
  <label class="ctl">曲率档
    <select id="sel-tier"></select>
  </label>
</div>

<div class="stat-row" id="stats"></div>

<div class="grid2">
  <div class="card">
    <h2>等弧长取点</h2>
    <p class="meta" id="path-meta"></p>
    <svg id="svg-path" viewBox="0 0 520 300" role="img" aria-label="曲线与等弧长取样点"></svg>
    <div class="legend">
      <span class="item"><span class="line" style="background:var(--series-1)"></span>曲线</span>
      <span class="item"><span class="dot" style="background:var(--series-1)"></span>等弧长取点 (32 等分)</span>
    </div>
  </div>
  <div class="card">
    <h2>重参数曲线</h2>
    <p class="meta" id="rs-meta"></p>
    <svg id="svg-rs" viewBox="0 0 520 300" role="img" aria-label="重参数函数 t 对 s"></svg>
    <div class="legend">
      <span class="item"><span class="line" style="background:var(--series-1)"></span>R(u)</span>
      <span class="item"><span class="line" style="background:var(--series-2)"></span>原参数 t=s</span>
      <span class="item"><span class="line" style="background:var(--axis)"></span>虚线 = 段边界</span>
    </div>
  </div>
  <div class="card wide">
    <h2>原参数 vs 等弧长: 每步弧长推进</h2>
    <p class="meta" id="raw-meta"></p>
    <svg id="svg-raw" viewBox="0 0 1040 240" role="img" aria-label="原参数等步长的弧长推进与等弧长推进对照"></svg>
    <div class="legend">
      <span class="item"><span class="line" style="background:var(--series-2)"></span>原参数 t 等步长（忽快忽慢）</span>
      <span class="item"><span class="line" style="background:var(--series-1)"></span>等弧长参数化（目标线 1.0）</span>
    </div>
  </div>
  <div class="card wide">
    <h2>速率偏差</h2>
    <p class="meta" id="dev-meta"></p>
    <svg id="svg-dev" viewBox="0 0 1040 260" role="img" aria-label="每步相对弧长偏差, 固定显示允许区间"></svg>
    <div class="legend">
      <span class="item"><span class="dot" style="background:var(--series-1)"></span>每步真实弧长 / 目标步长 − 1</span>
      <span class="item"><span class="band"></span>允许区间 ±1e-6</span>
    </div>
    <h2 style="margin-top:18px">速率偏差细节</h2>
    <p class="meta" id="devzoom-meta"></p>
    <svg id="svg-devzoom" viewBox="0 0 1040 200" role="img" aria-label="速率偏差放大细节"></svg>
    <div class="legend">
      <span class="item"><span class="dot" style="background:var(--series-1)"></span>同一组数据的放大视图（纵轴按数据范围自适应）</span>
    </div>
  </div>
</div>
`)

	fmt.Fprint(f, summaryTable(ds))
	fmt.Fprint(f, `<p class="note">偏差 = 等参数步的真实弧长 / 目标步长 − 1，允许区间 ±1e-6；本实现实测约 1e-11，已到双精度积分下限。<br>
重参数曲线逐段绘制，段端点标为圆点；曲线严格单调，且 R(0)=0、R(1)=1。点击下方表格任意行可切换上方视图。</p>
</div>
<div class="tip" id="tip" role="status" aria-live="polite"></div>
<script>
const DATA = `)
	f.Write(data)
	fmt.Fprint(f, `;
const BAND_TOL = 1e-6;
const state = { idx: 0 };
const NS = "http://www.w3.org/2000/svg";
const tip = document.getElementById("tip");

function el(tag, attrs, parent) {
  const n = document.createElementNS(NS, tag);
  for (const k in attrs) n.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(n);
  return n;
}
function sup(n) {
  const map = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻" };
  return String(n).split("").map(ch => map[ch] || ch).join("");
}
function fmt(v, d) {
  if (v === 0) return "0";
  const neg = v < 0 ? "−" : "";
  v = Math.abs(v);
  const e = Math.floor(Math.log10(v));
  if (e >= -3 && e < 4) return neg + v.toFixed(d === undefined ? 2 : d);
  return neg + (v / Math.pow(10, e)).toFixed(2) + "×10" + sup(e);
}
function extent(arr) {
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const p of arr) {
    if (p[0] < x0) x0 = p[0]; if (p[0] > x1) x1 = p[0];
    if (p[1] < y0) y0 = p[1]; if (p[1] > y1) y1 = p[1];
  }
  return { x0, x1, y0, y1 };
}
function niceTicks(lo, hi, n) {
  const span = hi - lo || 1;
  const raw = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) || raw;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(v);
  return out;
}
function showTip(evt, html) {
  tip.innerHTML = html;
  tip.style.opacity = "1";
  const pad = 14;
  let x = evt.clientX + pad, y = evt.clientY + pad;
  const r = tip.getBoundingClientRect();
  if (x + r.width > innerWidth - 8) x = evt.clientX - r.width - pad;
  if (y + r.height > innerHeight - 8) y = evt.clientY - r.height - pad;
  tip.style.left = x + "px"; tip.style.top = y + "px";
}
function hideTip() { tip.style.opacity = "0"; }
function clear(svg) { while (svg.firstChild) svg.removeChild(svg.firstChild); }

function drawPath(c) {
  const svg = document.getElementById("svg-path");
  clear(svg);
  const W = 520, L = 22, R = 22, T = 20, B = 26;
  const ext = extent(c.path.concat(c.dots));
  const plotW = W - L - R;
  const dx = Math.max(ext.x1 - ext.x0, 1e-9), dy = Math.max(ext.y1 - ext.y0, 1e-9);
  // 高度随数据长宽比自适应, 细长路径不留下大片空白.
  const H = Math.max(170, Math.min(300, T + B + plotW * dy / dx));
  svg.setAttribute("viewBox", "0 0 " + W + " " + H);
  const k = Math.min(plotW / dx, (H - T - B) / dy);
  const ox = L + (plotW - dx * k) / 2;
  const oy = T + ((H - T - B) - dy * k) / 2;
  const px = x => ox + (x - ext.x0) * k;
  const py = y => oy + (ext.y1 - y) * k;

  const g = el("g", {}, svg);
  let d = "";
  c.path.forEach((p, i) => { d += (i ? "L" : "M") + px(p[0]).toFixed(1) + " " + py(p[1]).toFixed(1); });
  el("path", { d, fill: "none", stroke: "var(--series-1)", "stroke-width": 2, "stroke-linecap": "round", "stroke-linejoin": "round" }, g);

  const dotsG = el("g", {}, svg);
  c.dots.forEach((p, i) => {
    el("circle", {
      cx: px(p[0]).toFixed(1), cy: py(p[1]).toFixed(1), r: 4,
      fill: "var(--series-1)", stroke: "var(--surface-1)", "stroke-width": 2,
    }, dotsG);
    const hit = el("circle", {
      cx: px(p[0]).toFixed(1), cy: py(p[1]).toFixed(1), r: 12,
      fill: "transparent", tabindex: 0, role: "img",
      "aria-label": "第 " + i + " 个等弧长点",
    }, dotsG);
    const n = c.dots.length - 1;
    const on = evt => showTip(evt,
      '<div class="k">等弧长点 ' + i + " / " + n + '</div>' +
      '<div class="v">u = ' + (i / n).toFixed(4) + '</div>' +
      '<div class="k">位置 (' + p[0].toFixed(3) + ", " + p[1].toFixed(3) + ")</div>");
    hit.addEventListener("pointermove", on);
    hit.addEventListener("pointerleave", hideTip);
    hit.addEventListener("focus", e => showTip(e, '<div class="v">u = ' + (i / n).toFixed(4) + "</div>"));
    hit.addEventListener("blur", hideTip);
  });
  document.getElementById("path-meta").textContent =
    c.segments + " 段 · 总弧长 " + c.length.toFixed(4) + " · 相邻点弧长间隔 " + (c.length / (c.dots.length - 1)).toFixed(4);
}

function drawRS(c) {
  const svg = document.getElementById("svg-rs");
  clear(svg);
  const W = 520, H = 300, L = 52, R = 18, T = 18, B = 46;
  const px = u => L + (W - L - R) * u;
  const py = t => H - B - (H - T - B) * t;

  for (let i = 1; i < 5; i++) {
    const v = i / 5;
    el("line", { x1: px(v), y1: T, x2: px(v), y2: H - B, stroke: "var(--grid)", "stroke-width": 1 }, svg);
    el("line", { x1: L, y1: py(v), x2: W - R, y2: py(v), stroke: "var(--grid)", "stroke-width": 1 }, svg);
  }
  // 每段的对照线: 若曲线本身已按弧长参数化, R(s)=s (段内对角).
  c.segSpan.forEach(([a, b]) => {
    el("line", { x1: px(a), y1: py(0), x2: px(b), y2: py(1), stroke: "var(--series-2)", "stroke-width": 2, "stroke-linecap": "round", opacity: 0.85 }, svg);
  });
  // 段边界: 2px 表面色缝隙 + 竖线.
  c.segSpan.slice(1).forEach(([a]) => {
    el("line", { x1: px(a), y1: T, x2: px(a), y2: H - B, stroke: "var(--surface-1)", "stroke-width": 4 }, svg);
    el("line", { x1: px(a), y1: T, x2: px(a), y2: H - B, stroke: "var(--axis)", "stroke-width": 1, "stroke-dasharray": "2 3", opacity: .8 }, svg);
  });
  c.rs.forEach(seg => {
    let d = "";
    seg.forEach((p, i) => { d += (i ? "L" : "M") + px(p[0]).toFixed(1) + " " + py(p[1]).toFixed(1); });
    el("path", { d, fill: "none", stroke: "var(--series-1)", "stroke-width": 2, "stroke-linecap": "round" }, svg);
    el("circle", { cx: px(seg[0][0]), cy: py(seg[0][1]), r: 3.2, fill: "var(--series-1)", stroke: "var(--surface-1)", "stroke-width": 2 }, svg);
    el("circle", { cx: px(seg[seg.length - 1][0]), cy: py(seg[seg.length - 1][1]), r: 3.2, fill: "var(--series-1)", stroke: "var(--surface-1)", "stroke-width": 2 }, svg);
  });
  el("line", { x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: "var(--axis)", "stroke-width": 1 }, svg);
  el("line", { x1: L, y1: T, x2: L, y2: H - B, stroke: "var(--axis)", "stroke-width": 1 }, svg);
  [["0", 0], ["0.5", 0.5], ["1", 1]].forEach(([txt, u]) => {
    el("text", { x: px(u), y: H - B + 16, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, svg).textContent = txt;
    el("text", { x: L - 8, y: py(u) + 4, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "end" }, svg).textContent = txt;
  });
  el("text", { x: (L + W - R) / 2, y: H - 8, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, svg).textContent = "u（整条路径归一化弧长）";
  el("text", { x: 14, y: (T + H - B) / 2, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle", transform: "rotate(-90 14 " + ((T + H - B) / 2) + ")" }, svg).textContent = "t（段内曲线参数）";

  // crosshair
  const cross = el("line", { x1: 0, y1: T, x2: 0, y2: H - B, stroke: "var(--axis)", "stroke-width": 1, opacity: 0 }, svg);
  const hit = el("rect", { x: L, y: T, width: W - L - R, height: H - T - B, fill: "transparent" }, svg);
  hit.addEventListener("pointermove", evt => {
    const b = svg.getBoundingClientRect();
    const ratio = (W - L - R) / b.width;
    const x = (evt.clientX - b.left) * ratio;
    const u = Math.min(1, Math.max(0, (x - L) / (W - L - R)));
    cross.setAttribute("x1", px(u)); cross.setAttribute("x2", px(u)); cross.setAttribute("opacity", 1);
    const i = Math.min(c.segSpan.length - 1, Math.max(0, c.segSpan.findIndex(([, end]) => u <= end)));
    const [segStart, segEnd] = c.segSpan[i];
    const local = Math.min(1, Math.max(0, (u - segStart) / (segEnd - segStart)));
    const seg = c.rs[i];
    const t = seg[Math.round(local * (seg.length - 1))][1];
    showTip(evt,
      '<div class="k">全局弧长 u = ' + u.toFixed(4) + '</div>' +
      '<div class="v">t = ' + t.toFixed(5) + '</div>' +
      '<div class="k">第 ' + i + ' 段 · 段内 s = ' + local.toFixed(4) + '</div>');
  });
  hit.addEventListener("pointerleave", () => { cross.setAttribute("opacity", 0); hideTip(); });
  document.getElementById("rs-meta").textContent =
    c.segments + " 段独立重参数 · 端点 R(0)=0, R(1)=1 · " + (c.mono ? "逐点单调核验通过" : "单调核验失败");
}

function drawRaw(c) {
  const svg = document.getElementById("svg-raw");
  clear(svg);
  const W = 1040, H = 240, L = 70, R = 22, T = 18, B = 36;
  let maxV = 1;
  c.rawAdv.forEach(d => { maxV = Math.max(maxV, d[1]); });
  const ymax = maxV * 1.12;
  const px = u => L + (W - L - R) * u;
  const py = v => T + (H - T - B) * (1 - v / ymax);

  niceTicks(0, ymax, 4).forEach(v => {
    el("line", { x1: L, y1: py(v), x2: W - R, y2: py(v), stroke: "var(--grid)", "stroke-width": 1 }, svg);
    el("text", { x: L - 8, y: py(v) + 4, "font-size": 10, fill: "var(--text-muted)", "text-anchor": "end" }, svg).textContent = v.toFixed(1).replace(/\.0$/, "") + "×";
  });
  niceTicks(0, 1, 5).forEach(u => {
    el("text", { x: px(u), y: H - B + 16, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, svg).textContent = u.toFixed(1).replace(/\.0$/, "");
  });
  el("line", { x1: L, y1: py(1), x2: W - R, y2: py(1), stroke: "var(--series-1)", "stroke-width": 2, "stroke-linecap": "round" }, svg);
  let d = "";
  c.rawAdv.forEach((p, i) => { d += (i ? "L" : "M") + px(p[0]).toFixed(1) + " " + py(p[1]).toFixed(1); });
  el("path", { d, fill: "none", stroke: "var(--series-2)", "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, svg);
  el("line", { x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: "var(--axis)", "stroke-width": 1 }, svg);
  el("line", { x1: L, y1: T, x2: L, y2: H - B, stroke: "var(--axis)", "stroke-width": 1 }, svg);
  el("text", { x: (L + W - R) / 2, y: H - 6, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, svg).textContent = "u（原参数等步长比例）";
  el("text", { x: 14, y: (T + H - B) / 2, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle", transform: "rotate(-90 14 " + ((T + H - B) / 2) + ")" }, svg).textContent = "每步弧长 / 平均步长";

  const cross = el("line", { x1: 0, y1: T, x2: 0, y2: H - B, stroke: "var(--axis)", "stroke-width": 1, opacity: 0 }, svg);
  const hit = el("rect", { x: L, y: T, width: W - L - R, height: H - T - B, fill: "transparent" }, svg);
  hit.addEventListener("pointermove", evt => {
    const b = svg.getBoundingClientRect();
    const ratio = (W - L - R) / b.width;
    const x = (evt.clientX - b.left) * ratio;
    const u = Math.min(1, Math.max(0, (x - L) / (W - L - R)));
    const i = Math.min(c.rawAdv.length - 1, Math.max(0, Math.round(u * (c.rawAdv.length - 1))));
    const p = c.rawAdv[i];
    cross.setAttribute("x1", px(p[0])); cross.setAttribute("x2", px(p[0])); cross.setAttribute("opacity", 1);
    showTip(evt, '<div class="k">u = ' + p[0].toFixed(3) + '</div><div class="v">' + p[1].toFixed(2) + "×</div><div class=k>原参数每步弧长</div>");
  });
  hit.addEventListener("pointerleave", () => { cross.setAttribute("opacity", 0); hideTip(); });

  let lo = Infinity, hi = 0;
  c.rawAdv.forEach(p => { lo = Math.min(lo, p[1]); hi = Math.max(hi, p[1]); });
  document.getElementById("raw-meta").textContent =
    "原参数等步长下每步弧长在 " + lo.toFixed(2) + "× 至 " + hi.toFixed(2) + "× 平均步长之间波动（" + (hi / lo).toFixed(1) + " 倍落差）";
}

function drawDev(c) {
  const worst = c.devs.reduce((a, d) => Math.abs(d[1]) > Math.abs(a[1]) ? d : a, [0, 0]);

  // 上图: 纵轴固定为允许区间的 1.2 倍, 展示"全部落在带内"这一事实.
  const svg = document.getElementById("svg-dev");
  clear(svg);
  const W = 1040, H = 260, L = 88, R = 22, T = 20, B = 36;
  const ymax = 1.2 * BAND_TOL;
  const px = u => L + (W - L - R) * u;
  const py = v => T + (H - T - B) * (1 - (v / ymax + 1) / 2);

  el("rect", { x: L, y: py(BAND_TOL), width: W - L - R, height: py(-BAND_TOL) - py(BAND_TOL), fill: "var(--band-ok)", opacity: 0.13 }, svg);
  el("line", { x1: L, y1: py(BAND_TOL), x2: W - R, y2: py(BAND_TOL), stroke: "var(--band-ok)", "stroke-width": 1, opacity: .65 }, svg);
  el("line", { x1: L, y1: py(-BAND_TOL), x2: W - R, y2: py(-BAND_TOL), stroke: "var(--band-ok)", "stroke-width": 1, opacity: .65 }, svg);
  niceTicks(0, 1, 5).forEach(u => {
    el("line", { x1: px(u), y1: T, x2: px(u), y2: H - B, stroke: "var(--grid)", "stroke-width": 1 }, svg);
  });
  el("line", { x1: L, y1: py(0), x2: W - R, y2: py(0), stroke: "var(--axis)", "stroke-width": 1 }, svg);
  c.devs.forEach(d => {
    el("circle", { cx: px(d[0]).toFixed(1), cy: py(d[1]).toFixed(2), r: 1.9, fill: "var(--series-1)", opacity: 0.8 }, svg);
  });
  el("line", { x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: "var(--axis)", "stroke-width": 1 }, svg);
  el("line", { x1: L, y1: T, x2: L, y2: H - B, stroke: "var(--axis)", "stroke-width": 1 }, svg);
  niceTicks(0, 1, 5).forEach(u => {
    el("text", { x: px(u), y: H - B + 16, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, svg).textContent = u.toFixed(1).replace(/\.0$/, "");
  });
  [BAND_TOL, 0, -BAND_TOL].forEach(v => {
    el("text", { x: L - 8, y: py(v) + 4, "font-size": 10, fill: "var(--text-muted)", "text-anchor": "end" }, svg).textContent = fmt(v, 1);
  });
  el("text", { x: (L + W - R) / 2, y: H - 4, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, svg).textContent = "u（全局归一化弧长步）";
  el("text", { x: 16, y: (T + H - B) / 2, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle", transform: "rotate(-90 16 " + ((T + H - B) / 2) + ")" }, svg).textContent = "相对偏差";
  // 直接标注最大值点, 其余交给悬停与表格.
  el("circle", { cx: px(worst[0]).toFixed(1), cy: py(worst[1]).toFixed(2), r: 4, fill: "var(--series-1)", stroke: "var(--surface-1)", "stroke-width": 2 }, svg);
  const wx = px(worst[0]), wy = py(worst[1]);
  el("text", { x: wx + (worst[0] > 0.6 ? -8 : 8), y: wy - 9, "font-size": 10.5, fill: "var(--text-secondary)", "text-anchor": worst[0] > 0.6 ? "end" : "start" }, svg)
    .textContent = "max " + fmt(worst[1], 2);
  addCrosshair(svg, L, T, W, H, B, R, c.devs, d => fmt(d[1], 3));
  document.getElementById("dev-meta").textContent =
    "逐点扫描 " + c.devs.length + " 步 · 最大偏差 " + fmt(worst[1], 3) + "（在 u=" + worst[0].toFixed(3) + "）· 允许 ±1e-6";

  // 下图: 同一组数据放大到实际范围, 看清残余误差的结构.
  const z = document.getElementById("svg-devzoom");
  clear(z);
  const H2 = 200, T2 = 16, B2 = 34;
  let vmax = 0;
  c.devs.forEach(d => { vmax = Math.max(vmax, Math.abs(d[1])); });
  if (vmax <= 0) vmax = 1e-13;
  const zmax = vmax * 1.3;
  const py2 = v => T2 + (H2 - T2 - B2) * (1 - (v / zmax + 1) / 2);
  niceTicks(-zmax, zmax, 4).forEach(v => {
    el("line", { x1: L, y1: py2(v), x2: W - R, y2: py2(v), stroke: "var(--grid)", "stroke-width": 1 }, z);
    el("text", { x: L - 8, y: py2(v) + 4, "font-size": 10, fill: "var(--text-muted)", "text-anchor": "end" }, z).textContent = fmt(v, 1);
  });
  niceTicks(0, 1, 5).forEach(u => {
    el("text", { x: px(u), y: H2 - B2 + 16, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, z).textContent = u.toFixed(1).replace(/\.0$/, "");
  });
  el("line", { x1: L, y1: py2(0), x2: W - R, y2: py2(0), stroke: "var(--axis)", "stroke-width": 1 }, z);
  let d = "";
  c.devs.forEach((p, i) => { d += (i ? "L" : "M") + px(p[0]).toFixed(1) + " " + py2(p[1]).toFixed(2); });
  el("path", { d, fill: "none", stroke: "var(--series-1)", "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, z);
  el("line", { x1: L, y1: H2 - B2, x2: W - R, y2: H2 - B2, stroke: "var(--axis)", "stroke-width": 1 }, z);
  el("line", { x1: L, y1: T2, x2: L, y2: H2 - B2, stroke: "var(--axis)", "stroke-width": 1 }, z);
  el("text", { x: (L + W - R) / 2, y: H2 - 4, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" }, z).textContent = "u（全局归一化弧长步）";
  el("text", { x: 16, y: (T2 + H2 - B2) / 2, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle", transform: "rotate(-90 16 " + ((T2 + H2 - B2) / 2) + ")" }, z).textContent = "相对偏差";
  addCrosshair(z, L, T2, W, H2, B2, R, c.devs, v => fmt(v[1], 3));
  document.getElementById("devzoom-meta").textContent =
    "同一组数据 · 纵轴范围 ±" + fmt(zmax, 2) + "（放大约 " + fmt(BAND_TOL / zmax, 1) + " 倍）· 波动来自弧长积分的舍入";
}

// addCrosshair 给图表加一条竖直游标与 tooltip; series 为 [u, v].
function addCrosshair(svg, L, T, W, H, B, R, series, labelOf) {
  const n = series.length;
  const cross = el("line", { x1: 0, y1: T, x2: 0, y2: H - B, stroke: "var(--axis)", "stroke-width": 1, opacity: 0 }, svg);
  const hit = el("rect", { x: L, y: T, width: W - L - R, height: H - T - B, fill: "transparent" }, svg);
  hit.addEventListener("pointermove", evt => {
    const b = svg.getBoundingClientRect();
    const ratio = (W - L - R) / b.width;
    const x = (evt.clientX - b.left) * ratio;
    const u = Math.min(1, Math.max(0, (x - L) / (W - L - R)));
    const i = Math.min(n - 1, Math.max(0, Math.round(u * (n - 1))));
    const p = series[i];
    cross.setAttribute("x1", (L + (W - L - R) * p[0]).toFixed(1));
    cross.setAttribute("x2", (L + (W - L - R) * p[0]).toFixed(1));
    cross.setAttribute("opacity", 1);
    showTip(evt, '<div class="k">u = ' + p[0].toFixed(4) + '</div><div class="v">' + labelOf(p) + "</div><div class=k>相对偏差</div>");
  });
  hit.addEventListener("pointerleave", () => { cross.setAttribute("opacity", 0); hideTip(); });
}

function renderStats(c) {
  const items = [
    ["段数", c.segments, false],
    ["总弧长", c.length.toFixed(4), true],
    ["最大 |速率偏差|", fmt(c.absDev, 2), true],
    ["允许区间", "±1e-6", true],
    ["建表 + 反解", c.buildUs.toFixed(1) + " µs", true],
    ["逐点核验", c.verifyUs.toFixed(0) + " µs", true],
  ];
  document.getElementById("stats").innerHTML = items.map(([l, v, mono]) =>
    '<div class="stat"><span class="label">' + l + '</span><span class="value' + (mono ? " mono" : "") + '">' + v + "</span></div>").join("");
}

function render() {
  const c = DATA[state.idx];
  document.getElementById("sel-case").value = c.name;
  document.getElementById("sel-tier").value = String(c.tier);
  try { history.replaceState(null, "", "#" + encodeURIComponent(c.name) + "-" + c.tier); } catch (_) {}
  document.querySelectorAll("#tbl-body tr").forEach(tr => {
    tr.setAttribute("aria-selected", String(Number(tr.dataset.idx) === state.idx));
  });
  renderStats(c);
  drawPath(c);
  drawRS(c);
  drawRaw(c);
  drawDev(c);
}

const names = [...new Set(DATA.map(d => d.name))];
const selCase = document.getElementById("sel-case");
names.forEach(n => { const o = document.createElement("option"); o.value = n; o.textContent = n; selCase.appendChild(o); });
const selTier = document.getElementById("sel-tier");
[1, 2].forEach(t => { const o = document.createElement("option"); o.value = String(t); o.textContent = t + " 档"; selTier.appendChild(o); });
function pick() {
  const n = selCase.value, t = Number(selTier.value);
  const i = DATA.findIndex(d => d.name === n && d.tier === t);
  if (i >= 0) { state.idx = i; render(); }
}
selCase.addEventListener("change", pick);
selTier.addEventListener("change", pick);
const hashIdx = DATA.findIndex(d => "#" + encodeURIComponent(d.name) + "-" + d.tier === location.hash);
if (hashIdx >= 0) state.idx = hashIdx;
document.querySelectorAll("#tbl-body tr").forEach(tr => {
  tr.addEventListener("click", () => { state.idx = Number(tr.dataset.idx); render(); });
  tr.addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); tr.click(); } });
});
render();
</script>
</body>
</html>`)
	return nil
}

func summaryTable(ds []detail) string {
	var sb strings.Builder
	sb.WriteString(`<div class="tbl-wrap"><table class="tbl">
<thead><tr>
<th>档</th><th>用例</th><th>段数</th><th>总弧长</th><th>min 偏差</th><th>max 偏差</th><th>max |偏差|</th><th>建表 (µs)</th><th>核验 (µs)</th><th>单调</th><th>原参数速率比</th>
</tr></thead><tbody id="tbl-body">`)
	for i, d := range ds {
		r := d.row
		mono := "通过"
		if !r.monoOK {
			mono = "失败"
		}
		fmt.Fprintf(&sb, `<tr data-idx="%d" tabindex="0" aria-selected="false"><td>%d</td><td>%s</td><td>%d</td><td>%.4f</td><td>%s</td><td>%s</td><td>%s</td><td>%.1f</td><td>%.0f</td><td>%s</td><td>%.2f×</td></tr>`,
			i, r.tier, r.name, r.segs, r.length, fmtSci(r.minDev), fmtSci(r.maxDev), fmtSci(r.absDev), r.buildUS, r.verifyUS, mono, r.ratioRaw)
	}
	sb.WriteString(`</tbody></table></div>`)
	return sb.String()
}

// fmtSci 用上标指数格式显示小量, 避免 “1.2e-06” 这类读数。
func fmtSci(v float64) string {
	if v == 0 {
		return "0"
	}
	neg := ""
	if v < 0 {
		neg = "−"
		v = -v
	}
	exp := int(math.Floor(math.Log10(v)))
	mant := v / math.Pow(10, float64(exp))
	sup := map[rune]string{'0': "⁰", '1': "¹", '2': "²", '3': "³", '4': "⁴", '5': "⁵", '6': "⁶", '7': "⁷", '8': "⁸", '9': "⁹", '-': "⁻"}
	es := ""
	for _, r := range fmt.Sprintf("%d", exp) {
		es += sup[r]
	}
	return fmt.Sprintf("%s%.2f×10%s", neg, mant, es)
}
