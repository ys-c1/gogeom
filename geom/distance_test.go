package geom

import (
	"errors"
	"math"
	"strings"
	"testing"
)

func almost(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// 五组样本：平行、相交、交叠、分叉、远离。
func TestCurveDistanceSamples(t *testing.T) {
	samples := DistanceSamples()
	want := map[string]float64{
		"平行": 1,
		"相交": 0,
		"交叠": 0,
		"分叉": 0,
		"远离": math.Hypot(10-3, 4-0), // A 终点 (3,0) 到 B 起点 (10,4)
	}
	for _, sp := range samples {
		t.Run(sp.Name, func(t *testing.T) {
			r, err := CurveDistance(sp.A, sp.B)
			if err != nil {
				t.Fatalf("CurveDistance: %v", err)
			}
			if !almost(r.Dist, want[sp.Name], 1e-7) {
				t.Fatalf("距离错误: got %v want %v", r.Dist, want[sp.Name])
			}
			// 最近点对间距必须与报告的距离一致
			gap := math.Hypot(r.PA.X-r.PB.X, r.PA.Y-r.PB.Y)
			if !almost(gap, r.Dist, 1e-7) {
				t.Fatalf("最近点对不自洽: |PA-PB|=%v Dist=%v", gap, r.Dist)
			}
			if r.Steps < 1 {
				t.Fatalf("步数异常: %d", r.Steps)
			}
			t.Logf("Dist=%.9f PA=%v PB=%v Steps=%d", r.Dist, r.PA, r.PB, r.Steps)
		})
	}
}

// 验证一：最近点处两曲线的法向共线。
// 连接向量 w = PA-PB 同时垂直于两条曲线的切向（端点处只要求内侧那一半），
// 即两条法线都平行于 w，故共线。距离为零（相交/交叠/分叉）时 w 为零向量，
// 条件平凡成立，跳过。
func TestNormalCollinear(t *testing.T) {
	for _, sp := range DistanceSamples() {
		t.Run(sp.Name, func(t *testing.T) {
			r, err := CurveDistance(sp.A, sp.B)
			if err != nil {
				t.Fatalf("CurveDistance: %v", err)
			}
			if r.Dist < 1e-7 {
				t.Skip("距离为零，法向共线平凡成立")
			}
			w := Point{r.PA.X - r.PB.X, r.PA.Y - r.PB.Y}
			wn := math.Hypot(w.X, w.Y)
			// 参数在区间内部时，连接向量必须与该曲线切向正交
			if r.S > 1e-6 && r.S < 1-1e-6 {
				ta := sp.A[r.SegA].Deriv(r.S)
				dot := math.Abs(w.X*ta.X + w.Y*ta.Y)
				if dot > 1e-6*wn*math.Hypot(ta.X, ta.Y) {
					t.Fatalf("A 侧法向不共线: w·A'=%v", dot)
				}
			}
			if r.T > 1e-6 && r.T < 1-1e-6 {
				tb := sp.B[r.SegB].Deriv(r.T)
				dot := math.Abs(w.X*tb.X + w.Y*tb.Y)
				if dot > 1e-6*wn*math.Hypot(tb.X, tb.Y) {
					t.Fatalf("B 侧法向不共线: w·B'=%v", dot)
				}
			}
		})
	}
}

// 验证二：交替投影的距离历史单调不增。
func TestDistanceMonotone(t *testing.T) {
	for _, sp := range DistanceSamples() {
		t.Run(sp.Name, func(t *testing.T) {
			_, _, _, hist := segDist(sp.A[0], sp.B[0])
			if len(hist) == 0 {
				t.Fatal("无迭代历史")
			}
			for i := 1; i < len(hist); i++ {
				if hist[i] > hist[i-1]+1e-15 {
					t.Fatalf("第%d步距离上升: %v -> %v", i, hist[i-1], hist[i])
				}
			}
		})
	}
}

// 三种必须拦下的情形：曲线重合、控制点退化、段数为零，并指出第几段。
func TestCurveDistanceRejected(t *testing.T) {
	seg := Cubic{Point{0, 0}, Point{1, 0}, Point{2, 0}, Point{3, 0}}

	// 重合：两条完全相同的曲线
	if _, err := CurveDistance([]Cubic{seg}, []Cubic{seg}); !errors.Is(err, ErrCoincident) {
		t.Fatalf("重合应报 ErrCoincident, got %v", err)
	}

	// 控制点退化：曲线 A 第 1 段（下标 1）四个控制点缩成一点
	degen := Cubic{Point{5, 5}, Point{5, 5}, Point{5, 5}, Point{5, 5}}
	_, err := CurveDistance([]Cubic{seg, degen}, []Cubic{seg})
	if !errors.Is(err, ErrDegenerate) {
		t.Fatalf("退化应报 ErrDegenerate, got %v", err)
	}
	if !strings.Contains(err.Error(), "曲线A第1段") {
		t.Fatalf("错误应指出第几段: %v", err)
	}

	// 段数为零
	if _, err := CurveDistance(nil, []Cubic{seg}); !errors.Is(err, ErrBadSegments) {
		t.Fatalf("段数为零应报 ErrBadSegments, got %v", err)
	}
	if _, err := CurveDistance([]Cubic{seg}, nil); !errors.Is(err, ErrBadSegments) {
		t.Fatalf("段数为零应报 ErrBadSegments, got %v", err)
	}
}

// 点到曲线：驻点求解应与暴力扫描一致。
func TestPointSegDistAgainstBruteForce(t *testing.T) {
	c := Cubic{Point{0, 0}, Point{1, 2}, Point{2, -2}, Point{3, 0}}
	p := Point{1.5, 0.7}
	d2, s := pointSegDist(p, c)
	// 暴力扫描验证
	best, bs := math.Inf(1), 0.0
	for i := 0; i <= 20000; i++ {
		v := float64(i) / 20000
		q := c.At(v)
		d := (q.X-p.X)*(q.X-p.X) + (q.Y-p.Y)*(q.Y-p.Y)
		if d < best {
			best, bs = d, v
		}
	}
	if !almost(d2, best, 1e-6) || !almost(s, bs, 1e-3) {
		t.Fatalf("与暴力扫描不符: got (%v,%v) want (%v,%v)", d2, s, best, bs)
	}
}
