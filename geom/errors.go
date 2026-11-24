package geom

import (
	"errors"
	"strconv"
)

var (
	ErrBadSegments = errors.New("segments 必须为正")
	ErrTooFewPts   = errors.New("顶点数不足")
	ErrDegenerate  = errors.New("退化输入")
	ErrCoincident  = errors.New("控制点重合")
	ErrZeroDeriv   = errors.New("导数恒为零")
)

// CurveError 携带出错的段号（控制点段序号或采样段数）。
type CurveError struct {
	Segment int
	Err     error
}

func (e *CurveError) Error() string {
	return "第 " + strconv.Itoa(e.Segment) + " 段: " + e.Err.Error()
}

func (e *CurveError) Unwrap() error { return e.Err }
