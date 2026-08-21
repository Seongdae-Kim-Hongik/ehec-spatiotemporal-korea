#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 493 — JAMA Figure 1: The 25-year rise and its age structure (3패널 통합).
#  A 전국 crude vs 연령표준화 추세(2001/2008-2024) · B 연령별 발생률(4군, 2008-2024) · C 소아/성인 발생률 비.
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=os.environ.get("EHEC_ROOT",".")
C=os.path.join(ROOT,"data_cache"); RS=os.path.join(ROOT,"results_ST"); EH=os.path.join(ROOT,"EHEC")

tr=pd.read_csv(os.path.join(C,"ST_nat_agestd_trend.csv"))
rr=pd.read_csv(os.path.join(RS,"MASTER_ehec_age_rateratio.csv"))
ac=pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv"))
pop=pd.read_csv(os.path.join(C,"ST_agepop.csv")).rename(columns={"age10":"age_group"})
def band(g):
    if g=="0~9세": return "0-9"
    if g=="10~19세": return "10-19"
    if g in ("60~69세","70세 이상","70세이상"): return "60+"
    return "20-59"
cs=ac.copy(); cs["b"]=cs.age_group.map(band)
pp=pop.copy(); pp["b"]=pp.age_group.map(band)
B=cs.groupby(["year","b"])["cases"].sum().reset_index().merge(
  pp.groupby(["year","b"])["pop"].sum().reset_index(),on=["year","b"])
B=B[B.year.between(2008,2024)]; B["rate"]=B.cases/B["pop"]*1e5

plt.rcParams.update({"font.size":10.5,"font.family":["Helvetica","Arial","DejaVu Sans"],"axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,3,figsize=(15.5,4.6))
for _a in ax: _a.grid(axis="y",color="0.88",lw=0.6); _a.set_axisbelow(True)
A=ax[0]
A.plot(tr.year,tr.crude,"o-",color="#80796B",lw=1.8,ms=4,label="Crude")
t2=tr.dropna(subset=["agestd"])
A.plot(t2.year,t2.agestd,"s-",color="#B24745",lw=2,ms=4,label="Age-standardized (2015 standard population)")
A.fill_between(t2.year,t2.agestd_lo,t2.agestd_hi,color="#B24745",alpha=0.14)
A.set_title("A. National incidence, 2001-2024",loc="left",fontsize=11,fontweight="bold")
A.set_xlabel("Year"); A.set_ylabel("Cases per 100 000"); A.legend(frameon=False,fontsize=9)
Bx=ax[1]
cols={"0-9":"#B24745","10-19":"#DF8F44","20-59":"#79AF97","60+":"#374E55"}
for b,g in B.groupby("b"):
    Bx.plot(g.year,g.rate,"o-",color=cols[b],lw=1.8,ms=3.5,label=b)
Bx.set_yscale("log")
Bx.set_title("B. Age-specific incidence",loc="left",fontsize=11,fontweight="bold")
Bx.set_xlabel("Year"); Bx.set_ylabel("Cases per 100 000 (log scale)"); Bx.legend(frameon=False,fontsize=8.5,ncol=2)
Cx=ax[2]
Cx.fill_between(rr.year,rr.lo1,rr.hi1,color="#6A6599",alpha=0.16)
Cx.plot(rr.year,rr.RR_young_vs_old,"o-",color="#6A6599",lw=2,ms=4)
Cx.axhline(1,ls="--",color="0.6",lw=1)
Cx.set_yscale("log")
Cx.set_title("C. Child (0-9 y) to older (\u226510 y) rate ratio",loc="left",fontsize=11,fontweight="bold")
Cx.set_xlabel("Year"); Cx.set_ylabel("Rate ratio (log scale)")
plt.tight_layout()
out=os.path.join(EH,"JAMA_Fig1_rise_person.png")
plt.savefig(out,dpi=200,bbox_inches="tight"); print("저장:",out)
