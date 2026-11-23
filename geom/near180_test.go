package geom

import (
	"math"
	"testing"
)

func TestNear180Render(t *testing.T) {
	h := 5.0
	v := Point{X: 50, Y: 50}
	t0 := Point{X: 1, Y: 0}
	for _, deg := range []float64{179, 179.9} {
		rad := deg * math.Pi / 180
		t1 := Point{X: math.Cos(rad), Y: math.Sin(rad)}
		for _, st := range []JoinStyle{JoinMiter, JoinRound, JoinBevel} {
			var poly []Point
			switch st {
			case JoinMiter:
				poly = MiterPoly(v, t0, t1, h)
			case JoinRound:
				poly = RoundPoly(v, t0, t1, h, h, 256)
			case JoinBevel:
				poly = BevelPoly(v, t0, t1, h, h/2)
			}
			if !isConvex(poly) {
				t.Fatalf("deg %v style %v 非凸", deg, st)
			}
			r := RasterizePolys([][]Point{poly}, 0.2)
			if r.Count() == 0 {
				t.Fatalf("deg %v style %v 空覆盖", deg, st)
			}
			if r.Components() != 1 {
				t.Fatalf("deg %v style %v 不连通", deg, st)
			}
		}
	}
}
