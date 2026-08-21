#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 500 — 결측→0 오염 공변량 복구 (2026-08-21)
#  발견: 국민건강결과.csv 원자료에서 여러 지표가 미수집 연도를 0으로 기록.
#        230 조립 스크립트가 이를 그대로 실값으로 받아 h_pub 64%, g_basiclivelihood 71%가 0.
#  ★ h_pub은 컬럼 선택 자체가 오류 — '주점업수'(절대수, 2016-2019만) 대신
#    '인구천명당주점업수'(2011-2020 전 연도 실측, 진짜 density)를 써야 함.
#  처리: ①정상 컬럼으로 재추출 ②0→결측(비율·밀도는 진짜 0 불가) ③구별 실측평균=시간불변 공변량
#        (축산 2019-2023 스냅샷과 동일한 기존 선례)
import pandas as pd, numpy as np, re, os
RAW=os.path.join(os.environ.get("EHEC_DATA","data"),"독립변수/국민건강결과.csv")
ROOT=os.environ.get("EHEC_ROOT",".")
C=os.path.join(ROOT,"data_cache")

# 복구 대상: (패널 변수명, 원자료 컬럼, 단위)
FIX={"h_pub":("인구천명당주점업수","bars and pubs per 1000 population"),
     "h_fastfood":("인구천명당패스트푸드점수","fast-food outlets per 1000 population"),
     "g_basiclivelihood":("기초생활수급자수율","% of population receiving basic-livelihood benefits"),
     "g_solo_elder":("독거노인가구비율","% of households: older adults living alone"),
     "g_fiscal":("재정자립도","% fiscal independence index"),
     "g_sexratio":("성비","males per 100 females"),
     "h_selfhealth":("주관적건강수준인지율_표준화율","% rating own health as good (age-standardized)")}

raw=pd.read_csv(RAW,encoding="utf-8-sig",low_memory=False)
raw["region"]=raw["region"].astype(str).str.replace(r"\s+","",regex=True)
raw.loc[raw.region=="인천시미추홀구","region"]="인천시남구"
panel=pd.read_csv(os.path.join(C,"ST_curated46_panel.csv"))
panel["region"]=panel["region"].astype(str).str.replace(r"\s+","",regex=True)
regs=sorted(panel.region.unique())
print(f"분석 대상 시군구 {len(regs)}개")

rows=[]
for var,(col,unit) in FIX.items():
    if col not in raw.columns:
        print(f"  ✗ {var}: 원자료에 '{col}' 없음"); continue
    s=raw[["region","year",col]].copy()
    s=s[s.region.isin(regs)]
    s[col]=pd.to_numeric(s[col],errors="coerce")
    s.loc[s[col]<=0,col]=np.nan          # 0 = 미수집(비율·밀도는 진짜 0 불가)
    obs=s.dropna(subset=[col])
    yrs=sorted(obs.year.unique())
    cov=obs.groupby("region")[col].mean()   # 실측 구간 구별 평균 → 시간불변
    rows.append({"var":var,"source_col":col,"unit":unit,
                 "obs_years":f"{min(yrs)}-{max(yrs)}" if yrs else "없음",
                 "n_years":len(yrs),"n_regions":cov.notna().sum(),
                 "median":round(float(cov.median()),3),"min":round(float(cov.min()),3),"max":round(float(cov.max()),3)})
    print(f"  ✓ {var:20s} ← {col:22s} 실측 {len(yrs)}년({min(yrs)}-{max(yrs)}) · {cov.notna().sum()}/{len(regs)}구 · median {cov.median():.3f}")
    cov.rename(var).to_frame().to_csv(os.path.join(C,f"RECOV_{var}.csv"))

pd.DataFrame(rows).to_csv(os.path.join(C,"RECOV_summary.csv"),index=False)
print("\n저장: data_cache/RECOV_*.csv")
