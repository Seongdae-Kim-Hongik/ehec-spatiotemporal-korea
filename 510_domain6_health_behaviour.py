# -*- coding: utf-8 -*-
# 510 — D6 건강행태·예방행동 도메인 신설. §5-Z 품질게이트 전면 적용.
#   근거: c_flu_vax·h_selfhealth를 '검출 대리'로만 분류한 사전지정의 배제 근거가 부실했음(개정 R7).
#   IRR 보기 전에 변수·방향가설을 확정한다.
import pandas as pd, numpy as np, os, sys
ROOT=os.path.dirname(os.path.abspath(__file__)); C=os.path.join(ROOT,"data_cache")
RAW=os.path.join(os.environ.get("EHEC_DATA","data"),"독립변수/국민건강결과.csv")

# ── 사전지정: 변수 / 이론적 근거 / 방향가설 (IRR 확인 전 확정) ──
SPEC=[
 # key, 원자료 컬럼, 라벨, 하위축, 방향가설
 ("c_flu_vax",     "연간인플루엔자예방접종률_표준화율","Influenza vaccination coverage",   "예방행동","protective(-)"),
 ("c_infant_chk",  "영유아건강검진율_전체",            "Infant health check-up rate",      "예방행동","양방향: 예방행동(-) vs 소아 진단기회(+)"),
 ("c_nutrilabel",  "영양표시독해율_표준화율",          "Nutrition-label literacy",         "식품안전관심","protective(-)"),
 ("c_walking",     "걷기실천율_표준화율",              "Walking practice rate",            "일반건강행태","음성대조(≈1 기대)"),
 ("c_smoking",     "현재흡연율_표준화율",              "Current smoking rate",             "일반건강행태","risk(+) 건강행태 역지표"),
 ("c_drinking",    "월간음주율_표준화율",              "Monthly drinking rate",            "일반건강행태","risk(+) 외식·회식 노출"),
]
d=pd.read_csv(RAW,encoding="utf-8-sig",low_memory=False)
d["region"]=d.region.astype(str).str.replace(r"\s+","",regex=True)
panel=pd.read_csv(os.path.join(C,"ST_curated46_panel.csv"))
panel["_rk"]=panel.region.astype(str).str.replace(r"\s+","",regex=True)
YEARS=sorted(panel.year.unique()); regs=sorted(panel._rk.unique())

print("═"*78); print("D6 건강행태·예방행동 — §5-Z 품질게이트"); print("═"*78)
out={}; rows=[]
for key,col,lab,axis,hyp in SPEC:
    if col not in d.columns: print(f"✗ {key}: 원자료 컬럼 없음 ({col})"); continue
    sub=d[d.region.isin(regs)][["region","year",col]].copy()
    sub[col]=pd.to_numeric(sub[col],errors="coerce")
    # Z4 값범위 sanity: 백분율 [0,100], 0은 미수집으로 간주(비율지표는 실제 0 불가)
    bad_hi=int((sub[col]>100).sum()); bad_lo=int((sub[col]<0).sum())
    sub[col]=sub[col].where((sub[col]>0)&(sub[col]<=100))
    # Z1 연도별 0/결측 비율 감사
    g=sub.groupby("year")[col].apply(lambda x: x.notna().mean())
    obs_y=[int(y) for y in g.index if g[y]>0.5]
    if not obs_y: print(f"✗ {key}: 실측 연도 없음"); continue
    # 실측 연도만의 구별 평균 → 시간불변 공변량 (기존 D1-D5와 동일 처리)
    m=sub[sub.year.isin(obs_y)].groupby("region")[col].mean()
    cov=float(m.reindex(regs).notna().mean())
    if cov<0.90: print(f"✗ {key}: 지역 커버리지 {cov:.1%} <90% → 제외"); continue
    v=m.reindex(regs)
    out[key]=v
    rows.append(dict(key=key,label=lab,axis=axis,hypothesis=hyp,src_col=col,
        n_years=len(obs_y),year_min=min(obs_y),year_max=max(obs_y),coverage=round(cov,4),
        median=round(float(v.median()),2),p25=round(float(v.quantile(.25)),2),p75=round(float(v.quantile(.75)),2),
        unit="%",range_viol=bad_hi+bad_lo))
    print(f"✓ {key:14s} {len(obs_y):2d}년({min(obs_y)}-{max(obs_y)}) 커버리지{cov:6.1%} "
          f"중앙값{float(v.median()):5.1f}% [{float(v.quantile(.25)):.1f}-{float(v.quantile(.75)):.1f}]")

# Z5 교란상관 검사: 상호 + 기존 검출대리(h_docs)·규모(pop_dens)와의 상관
X=pd.DataFrame(out)
print("\n── Z5 상호상관 (|r|>0.7 = 중복 위험) ──")
cc=X.corr(); hi=[]
for i,a in enumerate(cc.columns):
    for b in cc.columns[i+1:]:
        r=cc.loc[a,b]
        if abs(r)>0.7: hi.append((a,b,r)); print(f"   ⚠️ {a} ~ {b}: r={r:+.3f}")
if not hi: print("   통과 — 모든 쌍 |r|≤0.7")
base=panel.groupby("_rk")[[c for c in ("h_docs","e_pop_dens","g_fiscal") if c in panel.columns]].first().reindex(regs)
print("\n── Z5 기존 변수와의 상관 (검출·규모 대리 여부) ──")
for k in X.columns:
    s=" ".join(f"{b}={X[k].corr(base[b]):+.2f}" for b in base.columns)
    print(f"   {k:14s} {s}")

pd.DataFrame(rows).to_csv(os.path.join(C,"D6_health_behaviour_spec.csv"),index=False,encoding="utf-8-sig")
X.index.name="_rk"; X.reset_index().to_csv(os.path.join(C,"D6_health_behaviour_values.csv"),index=False,encoding="utf-8-sig")
print(f"\n저장: D6 변수 {len(X.columns)}개 · data_cache/D6_health_behaviour_*.csv")
