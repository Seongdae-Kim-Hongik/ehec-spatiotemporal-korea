#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# fig4_grid_v2d.py — (2026-10-29 /333 15번) 592_ehec_window_grid.py 사본: 라벨 통일("Survey composite (PC1)", "per 1,000"). 출력 fig_src/EHEC_F42b_windows_grid_v2d.png (공유 F42b는 그대로 둠).
import os, numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
_f=lambda v: f"{v:.3f}" if f"{v:.2f}"=="1.00" else f"{v:.2f}"
H=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); ROOT=os.path.dirname(os.path.dirname(H)); RS=os.path.join(ROOT,"results_ST"); OUT=os.path.join(H,"fig_src"); os.makedirs(OUT,exist_ok=True)
mpl.rcParams.update({"font.family":["Helvetica Neue","Arial"],"font.size":10,"pdf.fonttype":42})
s=pd.read_csv(os.path.join(RS,"MASTER_final_4yr_split.csv")); t2=pd.read_csv(os.path.join(RS,"MASTER_final_table2.csv")); f=pd.read_csv(os.path.join(RS,"MASTER_final_4yr_fitstat.csv"))
ixp=os.path.join(RS,"MASTER_final_3win_interaction.csv"); ix=pd.read_csv(ixp)
ixc={(r.term.rsplit("_x_",1)[0],r.term.rsplit("_x_",1)[1]):bool(r.cred) for r in ix.itertuples()}
LAB={"KCHS survey composite (PC1)":"Survey composite (PC1)","Rural-urban axis (survey PC1)":"Survey composite (PC1)"}
norm=lambda l: LAB.get(l,l).replace("per 1000","per 1,000")
cols=list(f.period); ncase={r.period:int(r.N) for r in f.itertuples()}; order=list(t2["cov"]); labs={k:norm(v) for k,v in zip(t2["cov"],t2.label)}
RED="#c0392b"; BLUE="#1f4e79"; RL="#f4d0cb"; BL="#d3dff0"
fig,ax=plt.subplots(figsize=(9.5,9.6)); nr,nc=len(order),len(cols)
for i,v in enumerate(order):
    for j,c in enumerate(cols):
        r=s[(s.period==c)&(s["cov"]==v)]; assert len(r)==1; r=r.iloc[0]; cred=bool(r.cred); pos=r.IRR>1
        fc=(RL if pos else BL) if cred else "white"; ec=(RED if pos else BLUE) if cred else "#cccccc"
        ax.add_patch(Rectangle((j,nr-1-i),1,1,facecolor=fc,edgecolor=ec,linewidth=1.6 if cred else 0.6))
        col=(RED if pos else BLUE) if cred else "#555555"
        ax.text(j+0.5,nr-1-i+0.5,f"{r.IRR:.2f}\n({_f(r.lo)}–{_f(r.hi)})",ha="center",va="center",fontsize=8.6 if cred else 8.0,color=col,fontweight="bold" if cred else "normal")
        k={cols[1]:"W2",cols[2]:"W3"}.get(c)
        if k and ixc.get((v,k)): ax.text(j+0.95,nr-1-i+0.92,"†",ha="right",va="top",fontsize=9,color="#333333")
ax.set_xlim(0,nc); ax.set_ylim(0,nr); ax.set_xticks(np.arange(nc)+0.5); ax.set_xticklabels([f"{c}\n(n = {ncase[c]:,})" for c in cols],fontsize=9.5)
ax.set_yticks(np.arange(nr)+0.5); ax.set_yticklabels([labs[v] for v in order][::-1],fontsize=9.5); ax.xaxis.tick_top()
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(length=0); plt.tight_layout(); fig.savefig(os.path.join(OUT,"EHEC_F42b_windows_grid_v2d.png"),dpi=300,metadata={"Software":None}); plt.close(fig)
print("saved fig_src/EHEC_F42b_windows_grid_v2d.png | labels:",sorted(set(labs.values()))[-3:],"| interaction credible:",int(ix.cred.sum()),"/",len(ix))
