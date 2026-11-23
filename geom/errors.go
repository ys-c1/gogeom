package geom

import (
	"errors"
	"fmt"
)

var (
	ErrBadSegments  = errors.New("segments 必须为正")
	ErrTooFewPts    = errors.New("顶点数不足")
	ErrDegenerate   = errors.New("退化输入")
	ErrNegativeArea = errors.New("有向面积为负（外环需逆时针）")
	ErrBadWinding   = errors.New("洞环需顺时针（有向面积应为负）")
)

// VertexError 带环号与顶点号的输入错误，便于定位坏点。
type VertexError struct {
	Ring  int // 环号，0 为外环
	Index int // 出问题的顶点号（顶点不足时为缺失点的下标）
	Err   error
}

func (e *VertexError) Error() string {
	return fmt.Sprintf("环%d 顶点%d: %v", e.Ring, e.Index, e.Err)
}

func (e *VertexError) Unwrap() error { return e.Err }
