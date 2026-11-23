#!/usr/bin/env python3
"""把 mine/tri_XX_{before,after}.pgm 拼成左右对照 PNG。"""
import os
import sys
from PIL import Image

src = sys.argv[1]
dst = sys.argv[2]
os.makedirs(dst, exist_ok=True)
tags = sorted(
    f.split("_before.pgm")[0]
    for f in os.listdir(src)
    if f.endswith("_before.pgm")
)
gap = 24
for tag in tags:
    b = Image.open(os.path.join(src, f"{tag}_before.pgm")).convert("RGB")
    a = Image.open(os.path.join(src, f"{tag}_after.pgm")).convert("RGB")
    w, h = b.size
    canvas = Image.new("RGB", (w * 2 + gap, h), (255, 255, 255))
    canvas.paste(b, (0, 0))
    canvas.paste(a, (w + gap, 0))
    canvas.save(os.path.join(dst, f"{tag}_pair.png"))
    print(tag, "->", canvas.size)
