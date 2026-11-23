package geom

import (
	"math"
	"strconv"
	"testing"
)

// BenchmarkStrokeRaster 两档线宽、三种接头的孤立单接头栅格化耗时基准。
// 与 cmd/strokereport 的口径一致：40° 拐角，圆角 r=h、斜切 r=h/2，p=0.2。
func BenchmarkStrokeRaster(b *testing.B) {
	const pixel = 0.2
	const arcSegs = 128
	v := Point{X: 40, Y: 40}
	t0 := Point{X: 1, Y: 0}
	rad := 40.0 * math.Pi / 180
	t1 := Point{X: math.Cos(rad), Y: math.Sin(rad)}
	for _, w := range []float64{10, 30} {
		h := w / 2
		for _, st := range []struct {
			name string
			poly []Point
		}{
			{"miter", MiterPoly(v, t0, t1, h)},
			{"round", RoundPoly(v, t0, t1, h, h, arcSegs)},
			{"bevel", BevelPoly(v, t0, t1, h, h/2)},
		} {
			b.Run("w"+strconv.Itoa(int(w))+"_"+st.name, func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					RasterizePolys([][]Point{st.poly}, pixel)
				}
			})
		}
	}
}
