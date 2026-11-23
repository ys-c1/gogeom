package geom

import "fmt"

// Triangulate 对单条多边形环做耳切三角化，填充接口签名保持稳定：
// 输入顶点序列（顺时针/逆时针均可），输出互不重叠的三角形网格。
// 自交环先按奇偶填充规则拆分为简单面再分别三角化。
func Triangulate(ring []Point) ([]Triangle, error) {
	rep, err := TriangulateReport(ring)
	if err != nil {
		return nil, err
	}
	return rep.Triangles, nil
}

// TriangulateRings 对带洞多边形三角化：outer 为外环，holes 为洞环。
func TriangulateRings(outer []Point, holes [][]Point) ([]Triangle, error) {
	rep, err := TriangulateRingsReport(outer, holes)
	if err != nil {
		return nil, err
	}
	return rep.Triangles, nil
}

// TriangulateReport 与 Triangulate 同流程，额外返回逐块核对报告。
func TriangulateReport(ring []Point) (*TriReport, error) {
	if len(ring) < 3 {
		return nil, errTooFew(allIndices(len(ring)))
	}

	// 全共线（退化）优先判定：共线环的非相邻边在几何上必然“相交”，
	// 必须在自交检测之前报出，否则会被误判为自交。
	if _, _, ok := removeCollinear(ring); !ok {
		return nil, errCollinear(allIndices(len(ring)))
	}

	// 自交：先检测，再拆分；拆不出来的（共线重叠等）带顶点号报错
	if inter, badVerts := hasProperIntersection(ring); inter {
		faces, err := splitSelfIntersecting(ring)
		if err != nil {
			if ve, ok := err.(*VertexError); ok && len(ve.Vertices) == 0 {
				ve.Vertices = badVerts
			}
			return nil, err
		}
		return triangulateFaces(ring, faces, badVerts)
	}

	// 简单环：共线清理 + 耳切
	clean, removed, _ := removeCollinear(ring)
	clean = ensureCCW(clean)
	eps := epsFor(clean)
	origIdx := make([]int, len(clean))
	for i, p := range clean {
		origIdx[i] = indexOfPoint(ring, p, eps)
	}
	tris, rounds, rescueNote := earClip(clean, origIdx, nil)
	est := len(ring) - 2
	note := buildNote(rounds, est, rescueNote, removed, "")

	origArea := ringAbsArea(ring)
	total, overlaps, diff, rel := verifyTriangulation(tris, origArea)
	return &TriReport{
		Triangles:        tris,
		Rounds:           rounds,
		EstRounds:        est,
		RoundsNote:       note,
		CollinearRemoved: uniqSorted(removed),
		SplitFaces:       1,
		TotalTriArea:     total,
		OriginalArea:     origArea,
		AreaDiff:         diff,
		AreaRelDiff:      rel,
		OverlapPairs:     overlaps,
	}, nil
}

// TriangulateRingsReport 带洞多边形的完整报告。
func TriangulateRingsReport(outer []Point, holes [][]Point) (*TriReport, error) {
	if len(outer) < 3 {
		return nil, errTooFew(allIndices(len(outer)))
	}
	for _, h := range holes {
		if len(h) < 3 {
			return nil, errTooFew(allIndices(len(h)))
		}
	}
	// 外环/洞自身均不允许自交
	if inter, badVerts := hasProperIntersection(outer); inter {
		if _, err := splitSelfIntersecting(outer); err != nil {
			if ve, ok := err.(*VertexError); ok && len(ve.Vertices) == 0 {
				ve.Vertices = badVerts
			}
			return nil, err
		}
		return nil, errSelfIntersect(badVerts, "带洞接口的外环不允许自交")
	}
	for _, h := range holes {
		if inter, badVerts := hasProperIntersection(h); inter {
			return nil, errSelfIntersect(badVerts, "带洞接口的洞环不允许自交")
		}
	}

	// 全共线核对
	if _, _, ok := removeCollinear(outer); !ok {
		return nil, errCollinear(allIndices(len(outer)))
	}
	totalInputVerts := len(outer)
	for _, h := range holes {
		totalInputVerts += len(h)
	}

	// 清理各环的近似共线顶点，再桥接合并
	cleanOuter, removedOuter, ok := removeCollinear(outer)
	if !ok {
		return nil, errCollinear(allIndices(len(outer)))
	}
	var cleanHoles [][]Point
	var allRemoved []int
	allRemoved = append(allRemoved, removedOuter...)
	for _, h := range holes {
		ch, rh, okh := removeCollinear(h)
		if !okh {
			return nil, errCollinear(allIndices(len(h)))
		}
		cleanHoles = append(cleanHoles, ch)
		allRemoved = append(allRemoved, rh...)
	}

	merged, bridges, err := mergeHoles(cleanOuter, cleanHoles)
	if err != nil {
		return nil, err
	}
	merged = ensureCCW(merged)
	eps := epsFor(merged)
	origIdx := make([]int, len(merged))
	for i, p := range merged {
		origIdx[i] = indexOfPoint(outer, p, eps)
	}
	// 填充区谓词：在外环内部且不在任何洞内部。用于耳切时排除跨洞假耳。
	fillTest := func(s Point) bool {
		if !Inside(cleanOuter, s.X, s.Y) {
			return false
		}
		for _, h := range cleanHoles {
			if Inside(h, s.X, s.Y) {
				return false
			}
		}
		return true
	}
	tris, rounds, rescueNote := earClip(merged, origIdx, fillTest)

	// 带洞多边形的轮数估算：N - 2*(H+1)，N 为输入顶点总数，H 为洞数
	est := totalInputVerts - 2*(len(holes)+1)
	note := buildNote(rounds, est, rescueNote, allRemoved,
		"带洞桥接会把每个桥接点各重复一次（每洞 +2 顶点），合并环顶点多于原始输入")
	if note == "" {
		note = "实际轮次与估算一致"
	}

	origArea := ringAbsArea(outer)
	for _, h := range holes {
		origArea -= ringAbsArea(h)
	}
	total, overlaps, diff, rel := verifyTriangulation(tris, origArea)
	return &TriReport{
		Triangles:        tris,
		Rounds:           rounds,
		EstRounds:        est,
		RoundsNote:       note,
		CollinearRemoved: uniqSorted(allRemoved),
		Bridges:          bridges,
		SplitFaces:       1,
		TotalTriArea:     total,
		OriginalArea:     origArea,
		AreaDiff:         diff,
		AreaRelDiff:      rel,
		OverlapPairs:     overlaps,
	}, nil
}

// triangulateFaces 处理自交拆分后的多个简单面
func triangulateFaces(ring []Point, faces [][]Point, badVerts []int) (*TriReport, error) {
	var tris []Triangle
	totalRounds := 0
	est := len(ring) - 2
	var allRemoved []int
	var notes []string
	for fi, face := range faces {
		face = ensureCCW(face)
		clean, removed, ok := removeCollinear(face)
		if !ok {
			notes = append(notes, fmt.Sprintf("第 %d 个拆分面退化已跳过", fi))
			continue
		}
		clean = ensureCCW(clean)
		eps := epsFor(clean)
		origIdx := make([]int, len(clean))
		for i, p := range clean {
			origIdx[i] = indexOfPoint(ring, p, eps)
		}
		t, r, rescue := earClip(clean, origIdx, nil)
		tris = append(tris, t...)
		totalRounds += r
		allRemoved = append(allRemoved, removed...)
		if rescue != "" {
			notes = append(notes, rescue)
		}
	}
	if len(tris) == 0 {
		return nil, errSelfIntersect(badVerts, "自交拆分后的面全部退化")
	}

	// 原多边形面积 = 奇偶规则下各填充面面积之和
	var origArea float64
	for _, f := range faces {
		origArea += ringAbsArea(f)
	}
	total, overlaps, diff, rel := verifyTriangulation(tris, origArea)
	note := buildNote(totalRounds, est, joinNotes(notes), allRemoved,
		"自交拆分在交点处引入了新的分割顶点")
	if note == "" {
		note = "自交拆分后逐面耳切，实际轮次与按原顶点数的估算一致"
	}
	return &TriReport{
		Triangles:        tris,
		Rounds:           totalRounds,
		EstRounds:        est,
		RoundsNote:       note,
		CollinearRemoved: uniqSorted(allRemoved),
		SplitFaces:       len(faces),
		TotalTriArea:     total,
		OriginalArea:     origArea,
		AreaDiff:         diff,
		AreaRelDiff:      rel,
		OverlapPairs:     overlaps,
	}, nil
}

// buildNote 汇总耳切轮数与估算不一致的原因。exceedReason 给出实际轮次
// 超出估算时的具体来源（带洞桥接 / 自交拆分引入额外交点等）。
func buildNote(rounds, est int, rescue string, removed []int, exceedReason string) string {
	var note string
	if rescue != "" {
		note = rescue
	}
	if len(removed) > 0 {
		if note != "" {
			note += "；"
		}
		note += "近似共线顶点被清理，实际轮次低于按原始顶点数的估算"
	}
	if rounds > est {
		if note != "" {
			note += "；"
		}
		if exceedReason == "" {
			exceedReason = "拆分引入了额外交点顶点"
		}
		note += "实际轮次超出估算：" + exceedReason
	}
	if rounds == est && note == "" {
		note = "实际轮次与估算一致"
	}
	return note
}

func joinNotes(notes []string) string {
	out := ""
	for i, n := range notes {
		if i > 0 {
			out += "；"
		}
		out += n
	}
	return out
}
