package geom

// SamplePair 是一组曲线间距离的样本。
type SamplePair struct {
	Name string
	A, B []Cubic
}

// DistanceSamples 返回平行、相交、交叠、分叉、远离五组样本，
// 演示程序与测试共用，保证输出一致。
func DistanceSamples() []SamplePair {
	return []SamplePair{
		{ // 平行：两条水平线状曲线，间距恒为 1
			Name: "平行",
			A:    []Cubic{{Point{0, 0}, Point{1, 0}, Point{2, 0}, Point{3, 0}}},
			B:    []Cubic{{Point{0, 1}, Point{1, 1}, Point{2, 1}, Point{3, 1}}},
		},
		{ // 相交：水平线与竖直曲线在 (2,0) 处正交穿越
			Name: "相交",
			A:    []Cubic{{Point{0, 0}, Point{1.3333333333333333, 0}, Point{2.6666666666666665, 0}, Point{4, 0}}},
			B:    []Cubic{{Point{2, -2}, Point{2, -0.6666666666666666}, Point{2, 0.6666666666666666}, Point{2, 2}}},
		},
		{ // 交叠：S 形曲线与水平线多次穿越
			Name: "交叠",
			A:    []Cubic{{Point{0, 0}, Point{2, 0}, Point{4, 0}, Point{6, 0}}},
			B:    []Cubic{{Point{0, -1}, Point{2, 3}, Point{4, -3}, Point{6, 1}}},
		},
		{ // 分叉：共享起点 (0,0)，随后各自分开
			Name: "分叉",
			A:    []Cubic{{Point{0, 0}, Point{1, 0}, Point{2, 0}, Point{3, 0}}},
			B:    []Cubic{{Point{0, 0}, Point{1, 0.5}, Point{2, 1.5}, Point{3, 3}}},
		},
		{ // 远离：两条曲线分处远处，最近点是 A 的终点与 B 的起点
			Name: "远离",
			A:    []Cubic{{Point{0, 0}, Point{1, 0}, Point{2, 0}, Point{3, 0}}},
			B:    []Cubic{{Point{10, 4}, Point{11, 4}, Point{12, 4}, Point{13, 4}}},
		},
	}
}
