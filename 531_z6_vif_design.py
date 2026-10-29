# -*- coding: utf-8 -*-
# 531 — Z6(조사지표 주성분) + 도메인 내 공선성 정리 → 모형 투입 사양 확정. IRR 없음.
import pandas as pd, numpy as np, os, json
C=os.environ.get("FBD_CACHE","data_cache")   # (2026-10-21) 출력 위치 지정. kchs_pc1은 530b가 v2 방식으로 다시 계산해 덮어씀
D=pd.read_csv(os.path.join(C,"ST_prereg_panel.csv")); M=pd.read_csv(os.path.join(C,"ST_prereg_panel_META.csv"))
M=M[M.gate=="PASS"]; V=M.key.tolist()
X=D.groupby("region")[V].mean()   # 구별 평균(공선성·PC 진단용)
# ── Z6: 조사(KCHS) 계열 — 국민건강결과 출처 중 % 지표 + 손씻기 ──
survey=[k for k in V if ("국민건강결과" in str(M.set_index("key").loc[k,"source"]) and str(M.set_index("key").loc[k,"source"]).split(":")[1].endswith(("율","률","표준화율"))) or "handwash" in k or "soap" in k]
Xs=X[survey].dropna(axis=1,thresh=int(0.8*len(X))).dropna()
Z=(Xs-Xs.mean())/Xs.std(); u,s_,vt=np.linalg.svd(Z.values,full_matrices=False); ev=s_**2/(s_**2).sum()
pc1=pd.Series(u[:,0]*s_[0],index=Xs.index)
ext=X[["urban_population_share","older_adult_share","fiscal_independence","population_density"]].reindex(Xs.index)
if pc1.corr(ext.urban_population_share)>0: pc1=-pc1   # 부호: 농촌 방향이 +
print(f"═══ Z6: 조사지표 {len(Xs.columns)}개 주성분 ═══")
print(f"  PC1 {ev[0]*100:.1f}% · PC2 {ev[1]*100:.1f}% · PC3 {ev[2]*100:.1f}%")
print("  PC1 ↔ 도시인구 {:+.2f} · 고령 {:+.2f} · 재정 {:+.2f} · 인구밀도 {:+.2f}".format(*[pc1.corr(ext[c]) for c in ext.columns]))
ld=pd.Series(vt[0]*(-1 if pc1.equals(-pd.Series(u[:,0]*s_[0],index=Xs.index)) else 1),index=Xs.columns).sort_values()
print("  최대 적재:",", ".join(f"{k}({v:+.2f})" for k,v in pd.concat([ld.head(3),ld.tail(3)]).items()))
verdict="조사지표 종합축 확인 → 모형에 kchs_pc1 명시 포함, 조사지표는 개별 행태로 해석 금지" if ev[0]>0.35 else "공통성분 약함"
print("  판정:",verdict)
# PC1을 패널에 추가
D["kchs_pc1"]=D.region.map(pc1)
# ── 도메인 내 공선성: |r|>0.85 쌍은 커버리지·실측연수 높은 쪽만 (사전 규칙) ──
print("\n═══ 도메인 내 |r|>0.85 정리 ═══")
drop=[]; meta=M.set_index("key")
for dom,g in M.groupby("domain"):
    ks=[k for k in g.key if k not in drop]; cc=X[ks].corr().abs()
    for i,a in enumerate(ks):
        for b in ks[i+1:]:
            if a in drop or b in drop: continue
            if cc.loc[a,b]>0.85:
                sa=(meta.loc[a,"region_coverage"],meta.loc[a,"obs_years"]); sb=(meta.loc[b,"region_coverage"],meta.loc[b,"obs_years"])
                loser=b if sa>=sb else a; drop.append(loser); print(f"  {dom[:14]:14s} {a} ~ {b} r={cc.loc[a,b]:.2f} → 제외 {loser}")
keep=[k for k in V if k not in drop]
print(f"\n투입 변수: {len(keep)} (제외 {len(drop)})")
# ── 도메인별 투입 세트 저장 ──
design={"pc1":"kchs_pc1","pc1_variance":round(float(ev[0]),3),"domains":{}}
for dom,g in M[M.key.isin(keep)].groupby("domain"):
    design["domains"][dom]={"vars":g.key.tolist(),"direction":dict(zip(g.key,g.direction)),"time_invariant":dict(zip(g.key,g.time_invariant.astype(bool)))}
design["negative_controls"]=M[(M.key.isin(keep))&(M.direction=="0")].key.tolist()
design["detection_controls"]=design["domains"].get("D7 detection",{}).get("vars",[])
# ★2026-10-21: PREREG_DESIGN.json은 사전 지정 설계라 WRITE_DESIGN=1일 때만 덮어씀. 아니면 기존 설계(data_cache)와 비교만 하고, 다르면 _v2_logs에 새 계산본을 남김.
_cur=os.path.join("data_cache","PREREG_DESIGN.json"); _old=json.load(open(_cur)) if os.path.exists(_cur) else None
if os.environ.get("WRITE_DESIGN")=="1":
    json.dump(design,open(os.path.join(C,"PREREG_DESIGN.json"),"w"),ensure_ascii=False,indent=1); print("PREREG_DESIGN.json 덮어씀(WRITE_DESIGN=1)")
else:
    _new=json.loads(json.dumps(design,ensure_ascii=False)); _diff=[] if _old is None else sorted(k for k in set(_new)|set(_old) if _new.get(k)!=_old.get(k))
    print("PREREG_DESIGN.json 덮어쓰지 않음(사전 지정 설계 보존) · 새로 계산한 설계와 기존 설계 차이 항목:",_diff or "없음")
    if _diff:
        _ld=os.path.join(C,"_v2_logs"); os.makedirs(_ld,exist_ok=True); json.dump(design,open(os.path.join(_ld,"PREREG_DESIGN_recomputed.json"),"w"),ensure_ascii=False,indent=1)
D.to_csv(os.path.join(C,"ST_prereg_panel.csv"),index=False,encoding="utf-8-sig")
print("\n도메인별 투입:"); [print(f"  {d:16s} {len(v['vars']):3d}") for d,v in design["domains"].items()]
print(f"  음성대조 {len(design['negative_controls'])} · 검출통제 {len(design['detection_controls'])}")
print("저장: PREREG_DESIGN.json · ST_prereg_panel.csv(+kchs_pc1)")
