package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")
	ErrZeroLength  = errors.New("曲线段长度为零")
)
