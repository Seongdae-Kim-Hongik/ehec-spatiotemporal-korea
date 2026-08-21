#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 492 — JAMA Figure 2-B용: count-matched 영모형의 '연도별' 점유율 궤적 300회 → 밴드(2.5/50/97.5%).
#  (487은 후반-전반 차이만 저장 — 본문 그림엔 관측선 vs 영모형 밴드 겹침이 필요)
import os, numpy as np, pandas as pd, pywt
from statsmodels.tsa.seasonal import STL
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.environ.get("EHEC_ROOT",".")
RS=os.path.join(ROOT,"results_ST")
d0=pd.read_csv(BASE+"/식중독최종.csv",usecols=["disease","year","week","cases"])
d0["cases"]=pd.to_numeric(d0["cases"],errors="coerce").fillna(0); d0=d0[(d0.week>=1)&(d0.week<=52)]
d0=d0[d0.disease=="장출혈성대장균감염증"]
Y0,Y1=2001,2024
idx=pd.MultiIndex.from_product([range(Y0,Y1+1),range(1,53)],names=["year","week"])
w=d0.groupby(["year","week"])["cases"].sum().reindex(idx,fill_value=0).reset_index()
years=w.year.values; nW=len(w)
scales=np.geomspace(4,260,140); periods=pywt.scale2frequency("cmor1.5-1.0",scales)**-1
ann=(periods>=48)&(periods<=56); tt=np.arange(nW)
coi=np.ones((len(scales),nW),bool)
for si,s in enumerate(scales):
    ef=np.sqrt(2)*s; coi[si]=(tt>=ef)&(tt<=nW-1-ef)
def yearly_share(counts):
    x=np.log1p(counts.astype(float))
    x=x-STL(x,period=52,robust=True).fit().trend
    P=np.abs(pywt.cwt(x,scales,"cmor1.5-1.0",sampling_period=1.0)[0])**2
    Pm=np.where(coi,P,np.nan)
    vt=np.isfinite(Pm[ann]).all(axis=0)
    tot=np.nansum(np.where(np.isfinite(Pm),Pm,0),axis=0)
    a=np.nansum(np.where(np.isfinite(Pm[ann]),Pm[ann],0),axis=0)
    st=np.full(nW,np.nan); ok=vt&(tot>0); st[ok]=a[ok]/tot[ok]
    return {y:(np.nanmean(st[years==y]) if vt[years==y].mean()>=0.5 else np.nan) for y in range(Y0,Y1+1)}
prof=w.groupby("week")["cases"].sum(); prof=prof/prof.sum()
Nyr=w.groupby("year")["cases"].sum()
rng=np.random.default_rng(21); NS=300
tra=[]
for k in range(NS):
    xs=np.concatenate([rng.poisson(prof.values*Nyr[y]) for y in range(Y0,Y1+1)])
    tra.append(yearly_share(xs))
T=pd.DataFrame(tra)
env=pd.DataFrame({"year":T.columns,
                  "lo":T.quantile(0.025),"mid":T.quantile(0.5),"hi":T.quantile(0.975)}).dropna()
env.to_csv(os.path.join(RS,"MASTER_ehec_wavelet_nullenvelope.csv"),index=False)
print(env.round(3).to_string(index=False)); print("저장 완료")
