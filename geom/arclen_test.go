package geom

import (
	"errors"
	"math"
	"strings"
	"testing"

	"gogeom/internal/fixtures"
)

// 测试用例来自 internal/paths, 见 fixtures.All。
type namedPath struct {
	name string
	fn   func(float64) []Cubic
}

func testPaths() []namedPath {
	out := make([]namedPath, len(fixtures.All))
	for i, p := range fixtures.All {
		p := p
		out[i] = namedPath{p.Name, func(tier float64) []Cubic {
			raw := p.Fn(tier)
			segs := make([]Cubic, len(raw))
			for k, c := range raw {
				segs[k] = Cubic{
					P0: Point{c[0][0], c[0][1]}, P1: Point{c[1][0], c[1][1]},
					P2: Point{c[2][0], c[2][1]}, P3: Point{c[3][0], c[3][1]},
				}
			}
			return segs
		}}
	}
	return out
}

const (
	deviationTol = 1e-6 // 给定区间: |速率偏差| <= 1e-6
	monoSteps    = 2000 // 逐点核验采样数
	devSteps     = 4000
)

func TestArcLengthReparamAllCases(t *testing.T) {
	for tier := 1.0; tier <= 2; tier++ {
		for _, tc := range testPaths() {
			tc := tc
			segs := tc.fn(tier)
			t.Run(tc.name, func(t *testing.T) {
				rp, err := BuildReparam(segs)
				if err != nil {
					t.Fatalf("档 %.0f 构建失败: %v", tier, err)
				}
				// 1) 重参数函数逐个采样点严格单调.
				if err := rp.CheckMonotone(monoSteps); err != nil {
					t.Fatal(err)
				}
				// 2) 等 u 间隔的真实弧长偏差落在区间 [-tol, tol].
				mn, mx, ma := rp.SpeedDeviation(devSteps)
				if ma > deviationTol {
					t.Fatalf("速率偏差 %.3e 超出区间 [%.0e, %.0e] (min=%g max=%g)", ma, -deviationTol, deviationTol, mn, mx)
				}
				// 3) 端点命中: u=0/1 对应整条路径首尾.
				p0, p1 := rp.PointAt(0), rp.PointAt(1)
				if p0 != segs[0].P0 || p1 != segs[len(segs)-1].P3 {
					t.Fatalf("端点不匹配: %v %v", p0, p1)
				}
				// 4) 等弧长取点: 每个点落在 i/steps 的目标弧长上 (弦长受曲率影响本就略短).
				steps := 64
				pts, err := rp.FlattenEqualArc(steps)
				if err != nil || len(pts) != steps+1 {
					t.Fatalf("等弧长细分异常 %v %d", err, len(pts))
				}
				for i := 0; i <= steps; i++ {
					u := float64(i) / float64(steps)
					if got := rp.GlobalArc(u); math.Abs(got-u*rp.TotalLength()) > 1e-9*rp.TotalLength() {
						t.Fatalf("第 %d 点弧长 %.12f 偏离目标 %.12f", i, got, u*rp.TotalLength())
					}
				}
				for i := 1; i < len(pts); i++ {
					if pts[i] == pts[i-1] {
						t.Fatalf("第 %d 与第 %d 个等弧长点重合", i-1, i)
					}
				}
			})
		}
	}
}

// TestReparamCurvesMonotonePointByPoint 直接检查每段 R(s) 的采样序列,
// 并核对 R(0)=0, R(1)=1, 以及 s 等间隔对应弧长等间隔.
func TestReparamCurvesMonotonePointByPoint(t *testing.T) {
	for _, tc := range testPaths() {
		segs := tc.fn(2)
		rp, err := BuildReparam(segs)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		for seg := 0; seg < rp.NumSegments(); seg++ {
			prev := rp.SegmentAt(seg, 0)
			for k := 1; k <= monoSteps; k++ {
				s := float64(k) / float64(monoSteps)
				cur := rp.SegmentAt(seg, s)
				if !(cur > prev) {
					t.Fatalf("%s 段 %d: R(%g)=%g 未严格大于 R(前点)=%g", tc.name, seg, s, cur, prev)
				}
				prev = cur
			}
			if math.Abs(rp.SegmentAt(seg, 1)-1) > 1e-9 {
				t.Fatalf("%s 段 %d: R(1) != 1", tc.name, seg)
			}
		}
	}
}

// TestErrorsReportSegmentIndex 三类错误都必须带段号.
func TestErrorsReportSegmentIndex(t *testing.T) {
	good := Cubic{P0: Point{0, 0}, P1: Point{1, 0}, P2: Point{2, 0}, P3: Point{3, 0}}
	zero := Cubic{P0: Point{5, 5}, P1: Point{5, 5}, P2: Point{5, 5}, P3: Point{5, 5}}
	degenerate := Cubic{P0: Point{0, 0}, P1: Point{0, 0}, P2: Point{2, 1}, P3: Point{3, 1}}

	cases := []struct {
		name string
		segs []Cubic
		want error
		idx  string
	}{
		{"长度为零", []Cubic{good, zero}, ErrZeroLength, "段 1"},
		{"控制点退化", []Cubic{good, good, degenerate}, ErrDegenerate, "段 2"},
	}
	for _, c := range cases {
		_, err := BuildReparam(c.segs)
		if !errors.Is(err, c.want) {
			t.Fatalf("%s: 期望 %v, 实际 %v", c.name, c.want, err)
		}
		if !strings.Contains(err.Error(), c.idx) {
			t.Fatalf("%s: 错误信息缺少段号 %q: %v", c.name, c.idx, err)
		}
	}

	// 曲线点数不足: 每段需 4 个控制点.
	if _, err := CubicsFromPoints([][]Point{{{0, 0}, {1, 0}, {2, 0}}}); !errors.Is(err, ErrTooFewPts) {
		t.Fatalf("点数不足应报 ErrTooFewPts, 实际 %v", err)
	} else if !strings.Contains(err.Error(), "段 0") {
		t.Fatalf("点数不足错误缺少段号: %v", err)
	}
	if _, err := BuildReparam(nil); !errors.Is(err, ErrTooFewPts) {
		t.Fatalf("空路径应报 ErrTooFewPts, 实际 %v", err)
	}
}

// TestOriginalParamAPIUnchanged 保证原有按参数取点接口行为不变.
func TestOriginalParamAPIUnchanged(t *testing.T) {
	p0, p1, p2, p3 := Point{0, 0}, Point{0, 10}, Point{10, 10}, Point{10, 0}
	for _, tval := range []float64{0, 0.25, 0.5, 0.75, 1} {
		if got, want := Bezier(p0, p1, p2, p3, tval), (Cubic{p0, p1, p2, p3}).At(tval); got != want {
			t.Fatalf("Bezier 与 Cubic.At 不一致 t=%g: %v %v", tval, got, want)
		}
	}
	pts, err := Flatten(p0, p1, p2, p3, 8)
	if err != nil || len(pts) != 9 || pts[0] != p0 || pts[8] != p3 {
		t.Fatalf("Flatten 行为改变: %v %d", err, len(pts))
	}
}
