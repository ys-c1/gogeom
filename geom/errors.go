package geom

import (
	"errors"
	"strconv"
)

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")

	// ErrNegativeTolerance 表示抽稀容差为负。
	ErrNegativeTolerance = errors.New("容差不能为负")
	// ErrCollinear 表示所有点共线（或全部重合），没有可抽稀的形状。
	ErrCollinear = errors.New("所有点共线")
	// ErrToleranceExceeded 表示逐点核对时最大偏差超过给定上界。
	ErrToleranceExceeded = errors.New("最大偏差超过容差")
)

// PointError 给几何错误附上点号（0 起），便于定位。
type PointError struct {
	Err   error
	Index int
}

func (e *PointError) Error() string {
	return e.Err.Error() + "（点号 " + strconv.Itoa(e.Index) + "）"
}

func (e *PointError) Unwrap() error { return e.Err }
