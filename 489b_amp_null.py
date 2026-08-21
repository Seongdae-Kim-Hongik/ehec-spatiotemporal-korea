#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 489b — 방법 C(STL 계절진폭 추세)의 count-matched 영모형: 계절성 불변+N만 관측대로일 때
#  log1p 척도 진폭도 count-floor로 저절로 상승하는지 → 관측 tau(+0.942)의 공정한 p.
import os, numpy as np, pandas as pd
from statsmodels.tsa.seasonal import STL
from scipy.stats import kendalltau
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.environ.get("EHEC_ROOT",".")
RS=os.path.join(ROOT,"results_ST")
d0=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["disease","year","week","cases"])
d0["cases"]=pd.to_numeric(d0["cases"],errors="coerce").fillna(0); d0=d0[(d0.week>=1)&(d0.week<=52)]
d0=d0[d0.disease=="장출혈성대장균감염증"]
Y0,Y1=2001,2024
idx=pd.MultiIndex.from_product([range(Y0,Y1+1),range(1,53)],names=["year","week"])
w=d0.groupby(["year","week"])["cases"].sum().reindex(idx,fill_value=0).reset_index()
years=w.year.values
prof=w.groupby("week")["cases"].sum(); prof=prof/prof.sum()
Nyr=w.groupby("year")["cases"].sum()
def amp_tau(counts):
    x=np.log1p(counts.astype(float))
    seas=STL(x,period=52,robust=True).fit().seasonal
    amps=[np.std(seas[years==y]) for y in range(Y0,Y1+1)]
    return kendalltau(range(Y0,Y1+1),amps)[0]
obs_tau=amp_tau(w.cases.values)
rng=np.random.default_rng(11); NS=300
null=[amp_tau(np.concatenate([rng.poisson(prof.values*Nyr[y]) for y in range(Y0,Y1+1)])) for _ in range(NS)]
null=np.array(null)
p1=(np.sum(null>=obs_tau)+1)/(len(null)+1)   # MC +1 보정
print(f"관측 tau={obs_tau:+.3f} · 영모형(계절불변+N관측) tau 평균 {null.mean():+.3f} (95% [{np.percentile(null,2.5):+.3f},{np.percentile(null,97.5):+.3f}]) · 단측 p={p1:.3f}")
pd.DataFrame({"null_tau":null}).to_csv(os.path.join(RS,"MASTER_ehec_amp_nulltau.csv"),index=False)
