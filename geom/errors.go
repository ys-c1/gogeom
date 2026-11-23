package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")

	// 三角化专用错误，均通过 VertexError 携带顶点号
	ErrTooFewVertices = errors.New("顶点不足三")
	ErrAllCollinear   = errors.New("全部顶点共线")
	ErrSelfIntersect  = errors.New("自交多边形无法拆分")
)

// VertexError 携带出错顶点在输入切片中的下标（顶点号）。
// Code 取值：TooFew / AllCollinear / SelfIntersect。
type VertexError struct {
	Code     string
	Msg      string
	Vertices []int
}

func (e *VertexError) Error() string { return e.Msg }

func (e *VertexError) Is(target error) bool {
	if t, ok := target.(*VertexError); ok {
		return e.Code == t.Code
	}
	switch target {
	case ErrTooFewVertices:
		return e.Code == "TooFew"
	case ErrAllCollinear:
		return e.Code == "AllCollinear"
	case ErrSelfIntersect:
		return e.Code == "SelfIntersect"
	}
	return false
}

func errTooFew(verts []int) *VertexError {
	return &VertexError{Code: "TooFew", Msg: "顶点不足三", Vertices: verts}
}

func errCollinear(verts []int) *VertexError {
	return &VertexError{Code: "AllCollinear", Msg: "全部顶点共线（退化输入）", Vertices: verts}
}

func errSelfIntersect(verts []int, msg string) *VertexError {
	if msg == "" {
		msg = "自交多边形无法拆分"
	}
	return &VertexError{Code: "SelfIntersect", Msg: msg, Vertices: verts}
}
