package geom

import (
	"errors"
	"fmt"
)

var (
	ErrBadSegments   = errors.New("segments 必须为正")
	ErrTooFewPts     = errors.New("顶点数不足")
	ErrDegenerate    = errors.New("退化输入")
	ErrSelfIntersect = errors.New("多边形自交未拆分")
	ErrDupVertex     = errors.New("顶点重合")
)

// VertexError 携带出错元素的序号（顶点下标或边下标）。
type VertexError struct {
	Err   error
	Index int
}

func (e *VertexError) Error() string {
	return fmt.Sprintf("%v，序号 %d", e.Err, e.Index)
}

func (e *VertexError) Unwrap() error { return e.Err }
