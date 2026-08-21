# -*- coding: utf-8 -*-
# 513 — c_hw_soap 신호 정체 규명. 4종 중 비누만 credible인 이유를 자료 수준에서 검증.
import pandas as pd, numpy as np, os, glob, re
ROOT=os.path.dirname(os.path.abspath(__file__)); C=os.path.join(ROOT,"data_cache")
RAW=os.path.join(os.environ.get("EHEC_DATA","data"),"독립변수/국민건강결과.csv")
d=pd.read_csv(RAW,encoding="utf-8-sig",low_memory=False)
d["region"]=d.region.astype(str).str.replace(r"\s+","",regex=True)
panel=pd.read_csv(os.path.join(C,"ST_curated46_panel.csv"))
panel["_rk"]=panel.region.astype(str).str.replace(r"\s+","",regex=True)
HW=[c for c in panel.columns if c.startswith("c_hw_")]

print("═"*76); print("① 원자료 컬럼 — 손씻기 계열 전수"); print("═"*76)
pat=[c for c in d.columns if ("손씻기" in str(c)) or ("비누" in str(c)) or ("손씻" in str(c))]
for c in pat:
    x=pd.to_numeric(d[c],errors="coerce"); x=x.where((x>0)&(x<=100))
    g=x.notna().groupby(d.year).mean(); obs=[int(y) for y in g.index if g[y]>0.5]
    if obs: print(f"   {str(c)[:44]:46s} {len(obs):2d}년 {min(obs)}-{max(obs)} med={x.median():.1f}")

print("\n"+"═"*76); print("② 패널 손씻기 4종 — 값 분포·상호상관"); print("═"*76)
X=panel.groupby("_rk")[HW].first()
for c in HW:
    v=X[c].dropna(); print(f"   {c:14s} n={len(v):3d} med={v.median():6.2f} IQR[{v.quantile(.25):.2f}-{v.quantile(.75):.2f}] "
        f"range[{v.min():.2f}-{v.max():.2f}] SD={v.std():.2f}")
print("\n   상호상관:")
cc=X.corr()
for i,a in enumerate(HW):
    print("   "+" ".join(f"{a if b==HW[0] else ''}"[:14].ljust(14) if b==HW[0] else "" for b in [HW[0]]) +
          "  ".join(f"{cc.loc[a,b]:+.2f}" for b in HW))

print("\n"+"═"*76); print("③ 사례수와의 원시 상관 (구별 총 발생 vs 각 지표)"); print("═"*76)
case=pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv"))
case["_rk"]=case.region.astype(str).str.replace(r"\s+","",regex=True)
tot=case[case.year.between(2011,2024)].groupby("_rk")["cases"].sum()
est=pd.read_csv(os.path.join(C,"ST_Estd_panel.csv"))
est["_rk"]=est.region.astype(str).str.replace(r"\s+","",regex=True)
E=est[est.year.between(2011,2024)].groupby("_rk")["E_std"].sum()
sir=(tot/E).replace([np.inf,-np.inf],np.nan)
A=X.join(sir.rename("SIR")).join(tot.rename("cases")).dropna()
for c in HW:
    r1=A[c].corr(A.SIR); r2=A[c].corr(np.log1p(A.cases))
    print(f"   {c:14s} r(SIR)={r1:+.3f}   r(log cases)={r2:+.3f}")

print("\n"+"═"*76); print("④ 비누 지표 결측 패턴 — 연도별 실측 여부"); print("═"*76)
soapcols=[c for c in d.columns if "비누" in str(c)]
for c in soapcols[:3]:
    x=pd.to_numeric(d[c],errors="coerce"); x=x.where((x>0)&(x<=100))
    g=x.notna().groupby(d.year).mean()
    yrs="".join("O" if g.get(y,0)>0.5 else "." for y in range(2008,2025))
    print(f"   {str(c)[:40]:42s} 2008→2024: {yrs}")
print("\n   (O=실측, .=미수집) — 실측 연도의 구별 평균을 시간불변 공변량으로 사용")
