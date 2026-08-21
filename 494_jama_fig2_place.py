#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 494 — JAMA Figure 2 (Place): 무제목 지도(412b) 2장 합성 + JAMA 패널 레터.
import os
from PIL import Image, ImageDraw, ImageFont
EH=os.path.join(os.environ.get("EHEC_ROOT","."),"EHEC")
A=Image.open(os.path.join(EH,"JAMA_mapA_RR.png"))
B=Image.open(os.path.join(EH,"JAMA_mapB_EHSA.png"))
H=max(A.height,B.height)
def padH(im):
    if im.height==H: return im
    c=Image.new("RGB",(im.width,H),"white"); c.paste(im,(0,(H-im.height)//2)); return c
A=padH(A.convert("RGB")); B=padH(B.convert("RGB"))
label_h=130; gap=40
canvas=Image.new("RGB",(A.width+gap+B.width,H+label_h),"white")
canvas.paste(A,(0,label_h)); canvas.paste(B,(A.width+gap,label_h))
d=ImageDraw.Draw(canvas)
try: f=ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc",84)
except Exception: f=ImageFont.load_default()
d.text((20,15),"A",font=f,fill="black")
d.text((A.width+gap+20,15),"B",font=f,fill="black")
out=os.path.join(EH,"JAMA_Fig3_place.png")
canvas.save(out,dpi=(300,300))
print("저장:",out,canvas.size)
