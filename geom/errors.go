package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")
	ErrCoincident  = errors.New("曲线重合")
)
