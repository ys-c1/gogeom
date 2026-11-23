package geom

import "errors"

var (
	ErrBadSegments   = errors.New("segments 必须为正")
	ErrTooFewPts     = errors.New("顶点数不足")
	ErrDegenerate    = errors.New("退化输入")
	ErrZeroRadius    = errors.New("接头半径为零")
	ErrBadJoinStyle  = errors.New("非法接头样式")
	ErrCollinearJoin = errors.New("两边共线无法形成接头")
)
