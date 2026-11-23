package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")
	ErrZeroParts   = errors.New("块数必须为正整数")
	ErrZeroArea    = errors.New("多边形面积为零")
)
