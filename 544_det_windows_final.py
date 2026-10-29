#!/usr/bin/env python3
# 544 — 최종모형 기준 결정요인 창별 그림(eFigure). 코로나 정렬 3창 2011-15/16-19/20-24(542w 산출) + 풀링(540). credible만, 빨강(+)/파랑(−).
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, os
R="results_ST"; S=pd.read_csv(f"{R}/MASTER_final_4yr_split.csv"); T=pd.read_csv(f"{R}/MASTER_final_table2.csv")
F=pd.read_csv(f"{R}/MASTER_final_4yr_fitstat.csv")
groups=[]
for pe in F.period:
    d=S[(S.period==pe)&(S.cred)]; groups.append((f"{pe} (N = {int(F[F.period==pe].N.iloc[0]):,})",d))
d=T[T.cred]; groups.append((f"Pooled 2011-2024 (N = 2,016)",d))
rows=[]; 
for g,d in groups:
    rows.append(("__hdr__",g,None,None,None))
    if len(d)==0: rows.append(("__none__","",None,None,None))
    for _,r in d.sort_values("IRR",ascending=False).iterrows(): rows.append(("pt",r.label,r.IRR,r.lo,r.hi))
fig,ax=plt.subplots(figsize=(7.2,0.30*len(rows)+1.3)); fig.subplots_adjust(left=0.36,right=0.97,top=0.97,bottom=0.10); y=len(rows)
for kind,lab,irr,lo,hi in rows:
    y-=1
    if kind=="__hdr__": ax.text(-0.02,y,lab,fontsize=8.5,fontweight="bold",va="center",ha="right",transform=ax.get_yaxis_transform(),clip_on=False); ax.axhline(y+0.5,color="0.85",lw=0.6)
    elif kind=="__none__": ax.text(1.0,y,"No credible determinant",fontsize=8,style="italic",color="0.5",va="center",ha="center")
    else:
        c="#c0392b" if irr>1 else "#2471a3"; ax.plot([lo,hi],[y,y],color=c,lw=1.3); ax.plot(irr,y,"o",color=c,ms=5); ax.text(-0.02,y,lab,fontsize=8,va="center",ha="right",transform=ax.get_yaxis_transform(),clip_on=False)
ax.set_xscale("log"); ax.axvline(1,color="0.4",lw=0.8,ls="--"); ax.set_ylim(-0.5,len(rows)-0.5); ax.set_yticks([]); ax.set_xlim(0.55,2.3); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks([0.67,1,1.5,2]); ax.set_xticklabels(["0.67","1","1.5","2"],fontsize=8); ax.set_xlabel("Incidence rate ratio per 1 SD (95% credible interval, log scale)",fontsize=8)
for s in ("top","right","left"): ax.spines[s].set_visible(False)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker="o",color="#c0392b",ls="",label="Positive association"),Line2D([],[],marker="o",color="#2471a3",ls="",label="Inverse association")],fontsize=8,loc="lower right",frameon=False)
pass
plt.savefig("EHEC/EHEC_Fig_det_windows.png",dpi=300); print("ok",len(rows))
