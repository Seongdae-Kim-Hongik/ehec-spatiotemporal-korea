# -*- coding: utf-8 -*-
# patch_na_counts.py — (2026-10-22 v2c) 542·542w·542c·553·553z·546b가 청주 2011-2013 마스킹 칸(NA) 때문에 sum(DD$cases)·mean(cases==0)을 NA로 쓴 서술 열만 채움.
#  모형 추정치는 그대로이고, 값은 우도에 들어간 칸(마스킹 3칸 제외) 기준. 스크립트는 HM 모드와 공유라 고치지 않음(blockers_JKMS.md 9번).
import os, pandas as pd, numpy as np
ROOT=os.environ.get("EHEC_ROOT","."); RS=os.path.join(ROOT,"results_ST"); C=os.path.join(ROOT,"data_cache")
cs=pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv"),dtype={"sigungu_cd":str}).groupby(["sigungu_cd","year"]).cases.sum()
g=pd.MultiIndex.from_product([sorted(pd.read_csv(os.path.join(C,"ST_Estd_panel.csv"),dtype={"sigungu_cd":str}).sigungu_cd.unique()),range(2011,2025)],names=["sigungu_cd","year"])
D=cs.reindex(g).fillna(0).reset_index(); assert len(D)==3206, len(D)
D.loc[(D.sigungu_cd=="33010")&D.year.between(2011,2013),"cases"]=np.nan
obs=D.dropna(); log=[]
def N(y0,y1): return int(obs[obs.year.between(y0,y1)].cases.sum())
def fix(f,col,mapper):
    p=os.path.join(RS,f); d=pd.read_csv(p); b=d[col].copy()
    d[col]=[mapper(r) if pd.isna(r[col]) else r[col] for _,r in d.iterrows()]
    if col!="zero_pct": d[col]=d[col].astype(int)
    d.to_csv(p,index=False); log.append((f,col,b.isna().sum(),list(d[col].unique())[:4]))
per=lambda s: tuple(int(x) for x in str(s).split("-"))
fix("MASTER_final_4yr_fitstat.csv","N",lambda r: N(*per(r.period)))
fix("MASTER_final_4yr_split.csv","N",lambda r: N(*per(r.period)))
fix("MASTER_final_covid_fitstat.csv","N",lambda r: N(*per(r.period)))
fix("MASTER_final_covid_phi_ci.csv","N",lambda r: N(*per(r.period)))
fix("MASTER_monitor_beef_sens.csv","N",lambda r: N(2011,2018) if r.model.endswith("2011_2018") else (N(2019,2024) if r.model.endswith("2019_2024") else N(2011,2024)))   # 546b
zp=round(100*(obs.cases==0).mean(),1); nz=int((obs.cases==0).sum())
fix("MASTER_agestd_supp_zinb.csv","zero_pct",lambda r: zp)
fix("MASTER_supp_nb_vs_zinb_fit.csv","obs_zero",lambda r: nz)
for l in log: print(l)
print("관측 칸",len(obs),"0칸",nz,f"({zp}%)")
