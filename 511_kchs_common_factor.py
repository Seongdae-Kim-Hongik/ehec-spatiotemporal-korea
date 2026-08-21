# -*- coding: utf-8 -*-
# 511 — 지역사회건강조사(KCHS) 기반 지표의 공통성분 진단.
#   동기: 흡연율~접종률 r=+0.93은 행태로 설명 불가 → 공통 요인 의심.
import pandas as pd, numpy as np, os
ROOT=os.path.dirname(os.path.abspath(__file__)); C=os.path.join(ROOT,"data_cache")
X=pd.read_csv(os.path.join(C,"D6_health_behaviour_values.csv")).set_index("_rk")
panel=pd.read_csv(os.path.join(C,"ST_curated46_panel.csv"))
panel["_rk"]=panel.region.astype(str).str.replace(r"\s+","",regex=True)
KCHS_HW=[c for c in panel.columns if c.startswith("c_hw_")]
base=panel.groupby("_rk")[[c for c in KCHS_HW+["e_rural_pop","e_urban_pop","e_pop_dens","g_elderly","g_fiscal","h_docs","d_field"] if c in panel.columns]].first()
A=X.join(base,how="inner").dropna()
print(f"분석 대상 {len(A)}개 구 · KCHS 지표 {len(X.columns)+len(KCHS_HW)}개\n")

K=[c for c in list(X.columns)+KCHS_HW if c in A.columns]
Z=(A[K]-A[K].mean())/A[K].std()
u,s,vt=np.linalg.svd(Z.values,full_matrices=False)
ev=s**2/ (s**2).sum()
print("── KCHS 지표 주성분 ──")
for i in range(4): print(f"   PC{i+1}: 분산설명 {ev[i]*100:5.1f}%")
pc1=pd.Series(u[:,0]*s[0],index=A.index)
if np.corrcoef(pc1, A["c_flu_vax"])[0,1]<0: pc1=-pc1
print("\n── PC1 적재량 (부호 = 같은 방향으로 움직임) ──")
ld=pd.Series(vt[0]*(1 if np.corrcoef(u[:,0]*s[0],A['c_flu_vax'])[0,1]>0 else -1),index=K).sort_values()
for k,v in ld.items(): print(f"   {k:16s} {v:+.3f}")
print(f"\n── PC1의 정체: 외부 변수와의 상관 (PC1 단독 분산설명 {ev[0]*100:.1f}%) ──")
for b in ["e_rural_pop","e_urban_pop","e_pop_dens","g_elderly","g_fiscal","h_docs","d_field"]:
    if b in A.columns:
        r=np.corrcoef(pc1,A[b])[0,1]
        flag=" ★" if abs(r)>0.6 else ""
        print(f"   {b:16s} r = {r:+.3f}{flag}")
pd.DataFrame({"_rk":A.index,"kchs_pc1":pc1.values}).to_csv(os.path.join(C,"D6_kchs_pc1.csv"),index=False,encoding="utf-8-sig")
print("\n저장: data_cache/D6_kchs_pc1.csv")
