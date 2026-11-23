package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")

	ErrTooFewSegments     = errors.New("段数不足")
	ErrZeroDerivative     = errors.New("端点导数为零")
	ErrCoincidentControls = errors.New("控制点重合")
)
