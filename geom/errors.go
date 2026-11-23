package geom

import "errors"

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")

	ErrCollinear      = errors.New("三点共线")
	ErrCoincident     = errors.New("坐标重合，无法定圆")
	ErrRadiusTooLarge = errors.New("半径超出画幅")

	ErrBadRadius       = errors.New("半径必须为正的有限值")
	ErrBadSagitta      = errors.New("弦高误差必须为正的有限值")
	ErrBadSweep        = errors.New("圆弧扫角必须在 (0, 2π] 内")
	ErrBadSamples      = errors.New("采样数必须不小于 2")
	ErrTooManySegments = errors.New("分段数超出上限")
)
