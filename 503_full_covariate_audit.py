#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 503 — 전 공변량 자료품질 전수 감사 (런북 §5-Z, Z1~Z4)
#  목적: 어떤 변수를 '쓸 수 있는지'를 통계와 무관하게 자료 품질만으로 확정한다.
#  ★ 통계적 유의성은 이 단계에서 절대 보지 않는다(사후 선택 방지).
import pandas as pd, numpy as np, os
C="data_cache"
import sys
SRC=sys.argv[1] if len(sys.argv)>1 else "ST_curated46_panel.csv"      # 기본=현행 패널, 인자로 백업 지정 가능
orig=pd.read_csv(os.path.join(C,SRC))
orig=orig[(orig.year>=2011)&(orig.year<=2024)]
skip={"region","year","area","pop"}
rows=[]
for v in [c for c in orig.columns if c not in skip]:
    x=pd.to_numeric(orig[v],errors="coerce")
    zr_by_year=(x.fillna(0)==0).groupby(orig.year).mean()
    bad_years=[int(y) for y in zr_by_year.index if zr_by_year[y]>0.9]
    obs_years=[int(y) for y in zr_by_year.index if zr_by_year[y]<=0.5]
    nz=x[x>0]
    rows.append(dict(var=v,
        zero_pct=round((x.fillna(0)==0).mean()*100,1),
        na_pct=round(x.isna().mean()*100,1),
        n_bad_years=len(bad_years), bad_years=",".join(map(str,bad_years)) if bad_years else "",
        n_obs_years=len(obs_years),
        obs_range=f"{min(obs_years)}-{max(obs_years)}" if obs_years else "없음",
        median=round(float(nz.median()),3) if len(nz) else np.nan,
        max=round(float(x.max()),1) if x.notna().any() else np.nan,
        pct_over100=round(float((x>100).mean()*100),1)))
a=pd.DataFrame(rows)
# 등급: Z1(0비율·결측연도) + Z4(값범위)
def grade(r):
    if r.n_bad_years>=8: return "D_심각오염"
    if r.n_bad_years>=3: return "C_부분오염"
    if r.n_bad_years>=1: return "B_경미"
    return "A_사용가능"
a["grade"]=a.apply(grade,axis=1)
a=a.sort_values(["grade","zero_pct"],ascending=[True,False])
a.to_csv(os.path.join(C,"COVARIATE_QUALITY_AUDIT.csv"),index=False)
print(f"[감사 대상] {SRC}")
print("═══ 등급 분포 ═══")
print(a.grade.value_counts().to_string())
for g in ["D_심각오염","C_부분오염","B_경미"]:
    sub=a[a.grade==g]
    if len(sub):
        print(f"\n─── {g} ({len(sub)}개) ───")
        print(sub[["var","zero_pct","n_bad_years","obs_range","median","max"]].to_string(index=False))
print(f"\n─── A_사용가능 ({len(a[a.grade=='A_사용가능'])}개) ───")
print(a[a.grade=="A_사용가능"][["var","zero_pct","obs_range","median"]].to_string(index=False))
print("\n저장: data_cache/COVARIATE_QUALITY_AUDIT.csv")
