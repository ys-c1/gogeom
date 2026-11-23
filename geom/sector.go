package geom

import "math"

// 扇形: A = r²θ/2;弓形 = 扇形 − 三角形: A = r²(θ − sinθ)/2。
// theta 为带符号圆心角(弧度),符号即旋向,|theta| ≤ 2π。

func checkRadiusTheta(r, theta float64) error {
	if !(r > 0) {
		return ErrBadRadius
	}
	if math.Abs(theta) > 2*math.Pi {
		return ErrBadAngle
	}
	return nil
}

// SectorArea 扇形解析面积,带符号(符号随 theta)。
func SectorArea(r, theta float64) (float64, error) {
	if err := checkRadiusTheta(r, theta); err != nil {
		return 0, err
	}
	return r * r * theta / 2, nil
}

// SegmentArea 弓形解析面积,带符号(符号随 theta)。
func SegmentArea(r, theta float64) (float64, error) {
	if err := checkRadiusTheta(r, theta); err != nil {
		return 0, err
	}
	return r * r * (theta - math.Sin(theta)) / 2, nil
}

// SignedArcArea 由圆心 c 与折线顶点 pts 累加带符号面积:
// 每条边贡献 = 弦三角形有向面积 + 弓形有向面积,
// 弓形符号取该边两端顶点相对圆心的夹角符号(atan2 于 (−π, π])。
// 顶点在半径为 r 的圆上时结果精确等于 ½r²·ΣΔθ。
func SignedArcArea(c Point, pts []Point) (float64, error) {
	if len(pts) < 2 {
		return 0, ErrTooFewPts
	}
	var s float64
	for i := 0; i+1 < len(pts); i++ {
		va := Point{pts[i].X - c.X, pts[i].Y - c.Y}
		vb := Point{pts[i+1].X - c.X, pts[i+1].Y - c.Y}
		ra := math.Hypot(va.X, va.Y)
		rb := math.Hypot(vb.X, vb.Y)
		if ra == 0 || rb == 0 {
			return 0, ErrBadRadius
		}
		cross := va.X*vb.Y - va.Y*vb.X
		dot := va.X*vb.X + va.Y*vb.Y
		if cross == 0 && dot < 0 {
			return 0, ErrCollinear
		}
		d := math.Atan2(cross, dot)
		rm := (ra + rb) / 2
		s += cross/2 + rm*rm*(d-math.Sin(d))/2
	}
	return s, nil
}
