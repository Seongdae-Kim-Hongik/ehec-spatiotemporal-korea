#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 505 — 잔여 오염 공변량 복구 (Z2: 원자료 컬럼 전수확인 → 더 나은 컬럼/병합)
#  g_foreign은 100% 0(완전 사망)이었으나 '인구천명당외국인수'가 15년 실측으로 존재했다.
#  분모 정규화 필요 변수는 인구/면적으로 나눠 밀도화(절대수 = 규모 대리변수 회피, Z5).
import pandas as pd, numpy as np, os
RAW=os.path.join(os.environ.get("EHEC_DATA","data"),"독립변수/국민건강결과.csv")
C="data_cache"
raw=pd.read_csv(RAW,encoding="utf-8-sig",low_memory=False)
raw["region"]=raw.region.astype(str).str.replace(r"\s+","",regex=True)
raw.loc[raw.region=="인천시미추홀구","region"]="인천시남구"
panel=pd.read_csv(os.path.join(C,"ST_curated46_panel.csv"))
panel["region"]=panel.region.astype(str).str.replace(r"\s+","",regex=True)
regs=set(panel.region); raw=raw[raw.region.isin(regs)]
pop=panel[["region","year","pop"]].copy() if "pop" in panel.columns else None

SPEC={  # var: (컬럼 우선순위, 인구정규화 여부, 단위)
 "g_foreign":       (["인구천명당외국인수"],False,"foreign residents per 1000 population"),
 "g_basiclivelihood":(["국민기초생활보장수급자수율","기초생활수급자수율"],False,"% of population receiving basic-livelihood benefits"),
 "e_rural_pop":     (["농촌인구수","농촌인구"],True,"rural residents as a proportion of district population"),
 "a_wastewater":    (["인구천명당폐수방류량"],False,"wastewater discharge per 1000 population"),
 "a_orgload":       (["인구천명당유기물질부하량발생량"],False,"organic pollutant load generated per 1000 population"),
}
out=[]
for var,(cols,need_pop,unit) in SPEC.items():
    cols=[c for c in cols if c in raw.columns]
    if not cols: print(f"  ✗ {var}: 후보 컬럼 없음"); continue
    m=raw[["region","year"]].copy(); m["val"]=np.nan
    for c in cols:
        x=pd.to_numeric(raw[c],errors="coerce"); x=x.where(x>0)
        m["val"]=m.val.fillna(x)
    if need_pop and pop is not None:
        m=m.merge(pop,on=["region","year"],how="left")
        m["val"]=np.where((m["pop"]>0)&m.val.notna(), m.val/m["pop"], np.nan)
    obs=m.dropna(subset=["val"]); yrs=sorted(obs.year.unique())
    cov=obs.groupby("region")["val"].mean()
    cov.rename(var).to_frame().to_csv(os.path.join(C,f"RECOV_{var}.csv"))
    old_zero=(pd.to_numeric(panel[var],errors="coerce").fillna(0)==0).mean()*100
    print(f"✓ {var:18s} {cols[0][:26]:28s} 실측 {len(yrs):2d}년({min(yrs)}-{max(yrs)}) · {cov.notna().sum():3d}/229구 · median {cov.median():>10.3f} · 기존0비율 {old_zero:.1f}%")
    out.append(dict(var=var,unit=unit,n_years=len(yrs),n_regions=int(cov.notna().sum()),median=round(float(cov.median()),4)))
pd.DataFrame(out).to_csv(os.path.join(C,"RECOV_remaining_summary.csv"),index=False)

# b_gw_unfit 값범위(Z4) 재처리: 백분율인데 >100 존재 → 검사합계 0 행의 잔재
g=panel[["region","year","b_gw_unfit"]].copy()
g["b_gw_unfit"]=pd.to_numeric(g.b_gw_unfit,errors="coerce")
bad=(g.b_gw_unfit>100).sum()
g.loc[(g.b_gw_unfit<=0)|(g.b_gw_unfit>100),"b_gw_unfit"]=np.nan
cov=g.dropna().groupby("region")["b_gw_unfit"].mean()
cov.to_frame().to_csv(os.path.join(C,"RECOV_b_gw_unfit.csv"))
print(f"✓ {'b_gw_unfit':18s} {'값범위 [0,100] 강제':28s} >100 제거 {bad}행 · {cov.notna().sum():3d}/229구 · median {cov.median():.2f}%")
