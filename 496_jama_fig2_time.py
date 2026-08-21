#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 496 — JAMA Figure 2: Time (A 여름 점유율 안정 · B 연주기 점유율 관측 vs count-matched 영모형 밴드)
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=os.environ.get("EHEC_ROOT",".")
RS=os.path.join(ROOT,"results_ST"); EH=os.path.join(ROOT,"EHEC")
S=pd.read_csv(os.path.join(RS,"MASTER_ehec_summer_share.csv"))
W=pd.read_csv(os.path.join(RS,"MASTER_ehec_wavelet_annshare.csv"))
E=pd.read_csv(os.path.join(RS,"MASTER_ehec_wavelet_nullenvelope.csv"))
plt.rcParams.update({"font.size":10.5,"font.family":["Helvetica","Arial","DejaVu Sans"],"axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
for _a in ax: _a.grid(axis="y",color="0.88",lw=0.6); _a.set_axisbelow(True)
A=ax[0]
A.bar(S.year,S.share_ja,color="#DF8F44",alpha=0.9,edgecolor="white",linewidth=0.4)
A.axhline(S[S.year>=2005].share_ja.mean(),ls="--",color="0.4",lw=1.2)
A.text(2001.2,88,f"2005-2024 mean {S[S.year>=2005].share_ja.mean():.0f}% (P = .86 to .89)",fontsize=8.5,color="0.3")
A.axhline(2/12*100,ls=":",color="0.6",lw=1); A.text(2001.2,18,"uniform 16.7%",fontsize=8,color="0.45",bbox=dict(fc="white",ec="none",alpha=0.85,pad=1.2))
A.set_title("A. July-August share of annual cases",loc="left",fontsize=11,fontweight="bold")
A.set_xlabel("Year"); A.set_ylabel("Percentage of annual cases"); A.set_ylim(0,95)
B=ax[1]
B.fill_between(E.year,E.lo,E.hi,color="#00A1D5",alpha=0.16,label="Constant-seasonality null\n(counts matched to observed; 95% envelope)")
B.plot(E.year,E.mid,"--",color="#00A1D5",lw=1.5)
B.plot(W.year,W.share,"o-",color="#B24745",lw=2,ms=4.5,label="Observed")
B.axvline(2013.5,ls=":",color="0.5",lw=1)
B.text(2013.2,0.205,"change point\n2014",fontsize=8,color="0.45",ha="right")
B.set_title("B. Annual-band share of wavelet spectral power\nvs count-matched null",loc="left",fontsize=11,fontweight="bold")
B.set_xlabel("Year"); B.set_ylabel("Annual-band / total wavelet power")
B.legend(frameon=False,fontsize=8.5,loc="upper left")
plt.tight_layout()
out=os.path.join(EH,"JAMA_Fig2_time.png")
plt.savefig(out,dpi=200,bbox_inches="tight"); print("저장:",out)
