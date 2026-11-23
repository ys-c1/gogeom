package geom

func Inside(poly []Point, px, py float64) bool {
	in := false
	n := len(poly)
	for i := 0; i < n; i++ {
		a, b := poly[i], poly[(i+1)%n]
		if (a.Y > py) == (b.Y > py) {
			continue
		}
		t := (py - a.Y) / (b.Y - a.Y)
		if a.X+(b.X-a.X)*t > px {
			in = !in
		}
	}
	return in
}

func Area(poly []Point) (float64, error) {
	if len(poly) < 3 {
		return 0, ErrTooFewPts
	}
	var s float64
	for i := 0; i < len(poly); i++ {
		a, b := poly[i], poly[(i+1)%len(poly)]
		s += a.X*b.Y - b.X*a.Y
	}
	return s / 2, nil
}
