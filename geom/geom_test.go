package geom

import "testing"

func TestBezierEndpoints(t *testing.T) {
	p0, p1 := Point{0, 0}, Point{10, 10}
	got := Bezier(p0, Point{0, 10}, Point{10, 10}, p1, 0)
	if got != p0 {
		t.Fatalf("起点错误 %v", got)
	}
}

func TestFlattenAndLength(t *testing.T) {
	pts, err := Flatten(Point{}, Point{Y: 10}, Point{X: 10, Y: 10}, Point{X: 10}, 8)
	if err != nil || len(pts) != 9 {
		t.Fatalf("细分异常 %v %d", err, len(pts))
	}
	if Length(pts) <= 10 {
		t.Fatal("长度应大于弦长")
	}
	if _, err := Flatten(Point{}, Point{}, Point{}, Point{}, 0); err == nil {
		t.Fatal("段数为零应报错")
	}
}

func TestAreaAndInside(t *testing.T) {
	sq := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	a, err := Area(sq)
	if err != nil || a != 100 {
		t.Fatalf("面积错误 %v %v", a, err)
	}
	if !Inside(sq, 5, 5) || Inside(sq, 15, 5) {
		t.Fatal("内外判定错误")
	}
	if _, err := Area(sq[:2]); err == nil {
		t.Fatal("顶点不足应报错")
	}
}
