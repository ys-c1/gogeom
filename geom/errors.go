package geom

import (
	"errors"
	"fmt"
)

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")

	ErrAllCollinear    = errors.New("点集全部共线，无法成层")
	ErrTooFewRemaining = errors.New("剩余点不足，无法剥到该层")
	ErrLayerGap        = errors.New("层号跳号")
	ErrAreaNotMonotone = errors.New("层面积随层号增大")
	ErrBadCoverage     = errors.New("面积覆盖目标必须在 (0,1] 内")
)

// LayerError 把错误归到具体层号，调用方可用 errors.As 取回 Layer。
type LayerError struct {
	Layer int
	Err   error
}

func (e *LayerError) Error() string {
	return fmt.Sprintf("第 %d 层: %v", e.Layer, e.Err)
}

func (e *LayerError) Unwrap() error { return e.Err }
