package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")
	ErrBadRadius   = errors.New("radius 必须为正")
	ErrBadAngle    = errors.New("theta 超出 [-2π, 2π]")
	ErrCollinear   = errors.New("pts 顶点与圆心共线,夹角未定")
)
