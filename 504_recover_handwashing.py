#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 504 — 손씻기(위생) 공변량 복구 (런북 §5-Z Z2: 원자료 컬럼 전수확인)
#  발견: 동일 지표가 원자료에서 여러 컬럼으로 분산(오타 컬럼 '화장실다년온후' 포함, '_1' 접미사 버전).
#        230 조립은 컬럼 1개만 선택 → 나머지 연도가 통째로 유실되어 c_hw_* 4개가 43~71% 오염.
#  처리: 지표별 컬럼군을 연도 단위로 coalesce(표준화율 우선) → 커버리지 회복 → 구별 평균(시간불변).
import pandas as pd, numpy as np, os
RAW=os.path.join(os.environ.get("EHEC_DATA","data"),"독립변수/국민건강결과.csv")
C="data_cache"
GROUPS={
 "c_hw_meal":   (["식사전손씻기실천율_표준화율","식사전손씻기실천율_표준화율_1","식사전손씻기실천율_조율","식사전손씻기실천율_조율_1"],
                 "% practising handwashing before meals (age-standardized)"),
 "c_hw_toilet": (["화장실다녀온후손씻기실천율_표준화율","화장실다년온후손씻기실천율_표준화율","화장실다녀온후손씻기실천율_표준화율_1",
                  "화장실다녀온후손씻기실천율_조율","화장실다년온후손씻기실천율_조율","화장실다녀온후손씻기실천율_조율_1"],
                 "% practising handwashing after using the toilet (age-standardized)"),
 "c_hw_out":    (["외출후손씻기실천율_표준화율","외출후손씻기실천율_표준화율_1","외출후손씻기실천율_조율","외출후손씻기실천율_조율_1"],
                 "% practising handwashing after going out (age-standardized)"),
 "c_hw_soap":   (["비누손세정제사용률_표준화율","비누,손세정제사용률_표준화율","비누손세정제사용률_조율","비누,손세정제사용률_조율","비누,손세정제사용률_조율_1"],
                 "% using soap or hand sanitiser when washing hands (age-standardized)"),
}
raw=pd.read_csv(RAW,encoding="utf-8-sig",low_memory=False)
raw["region"]=raw.region.astype(str).str.replace(r"\s+","",regex=True)
raw.loc[raw.region=="인천시미추홀구","region"]="인천시남구"
panel=pd.read_csv(os.path.join(C,"ST_curated46_panel.csv"))
regs=set(panel.region.astype(str).str.replace(r"\s+","",regex=True))
raw=raw[raw.region.isin(regs)]
summary=[]
for var,(cols,unit) in GROUPS.items():
    cols=[c for c in cols if c in raw.columns]
    m=raw[["region","year"]].copy(); m["val"]=np.nan
    used=[]
    for c in cols:                                    # 앞선 컬럼(표준화율) 우선, 빈 칸만 채움
        x=pd.to_numeric(raw[c],errors="coerce")
        x=x.where((x>0)&(x<=100))                     # 0=미수집, >100=불가능 → 결측
        before=m.val.notna().sum()
        m["val"]=m.val.fillna(x)
        if m.val.notna().sum()>before: used.append(f"{c}(+{m.val.notna().sum()-before})")
    obs=m.dropna(subset=["val"])
    yrs=sorted(obs.year.unique())
    cov=obs.groupby("region")["val"].mean()
    cov.rename(var).to_frame().to_csv(os.path.join(C,f"RECOV_{var}.csv"))
    summary.append(dict(var=var,unit=unit,n_cols_merged=len(used),
                        obs_years=f"{min(yrs)}-{max(yrs)}",n_years=len(yrs),
                        n_regions=int(cov.notna().sum()),median=round(float(cov.median()),2)))
    print(f"✓ {var:12s} 컬럼 {len(used)}개 병합 → 실측 {len(yrs)}년({min(yrs)}-{max(yrs)}) · {cov.notna().sum()}/229구 · median {cov.median():.1f}%")
    for u in used: print(f"      {u}")
pd.DataFrame(summary).to_csv(os.path.join(C,"RECOV_handwashing_summary.csv"),index=False)
print("\n저장: RECOV_c_hw_*.csv")
