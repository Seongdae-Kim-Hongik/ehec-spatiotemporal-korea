# -*- coding: utf-8 -*-
# 530 — 사전지정 풀(108) → 분석 패널 조립 + Z1–Z5 게이트. 모형 진입 전 단계. 어떤 IRR도 계산하지 않는다.
#  규칙: (1) 0은 미측정으로 취급(NA) (2) 실측 연도의 구별 평균 = 시간불변 값 (3) 분모: pop=ST_agepop 합, area=지역 중앙값(2026-09-03 정정; 구 하수도 파일 분모는 이중계상)
#        (4) 같은 지표가 여러 컬럼에 흩어진 경우 연도별 coalesce (5) 게이트 실패 변수는 패널에서 제외하고 사유 기록
# ★2026-10-21 v2 정본화(이전 판 530_build_prereg_panel.py.bak_pre_v2): (1) area = final.shp 기하 면적(코드 결합) — 하수도 파일 총면적 중앙값은 2017+ 소계 이중 합산으로 부풀었음
#   (2) pop = FBD_CACHE/ST_agepop.csv(400 v2, 코드 결합, 빈칸 없음 단언 — 옛 보간·중앙값 채움 삭제; 옛 판은 군위 2011-22를 전국 중앙값 인구로 채웠음)
#   (3) 끝에서 545(총조사 농가 밀도, shp 면적) → 548 → 531(kchs_pc1, 설계 파일은 WRITE_DESIGN=1일 때만) → 530b(v2 교정: 세종·제주·비누·정화조 구조적 0·kchs_pc1 v2)를 같은 FBD_CACHE로 이어 실행.
#   원자료 시군구 결합(이름)은 옛 판 그대로 두고(세종 등 빠진 칸은 530b가 코드 결합으로 명시 채움·기록), 출력 위치 = 환경변수 FBD_CACHE(기본 data_cache).
import pandas as pd, numpy as np, os, re, sys, unicodedata as ud, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import v2_fixes as v2
import jkms_yeargate_v2d as yg   # ★2026-10-22 v2d: 563 year_anomaly_gate와 같은 연도 이상 게이트(+원자료 재확인), _v2_promote_261021/19_JKMS패널결함.md
import jkms_rawfix_v2c as jr   # ★2026-10-22 v2c: 주택소유율 2022 척도·격리병상 미발표 연도 교정(_v2_promote_261021/19_JKMS패널결함.md)
RAW=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
CIN="data_cache"; C=v2.cache_dir(); P=pd.read_csv(os.path.join(CIN,"COVARIATE_PREREG_POOL.csv")); I=P[P.status=="INCLUDE"].copy()
nfc=lambda s: ud.normalize("NFC",str(s))
def rd(f):
    p=os.path.join(RAW,f)
    for enc in ("utf-8-sig","cp949"):
        try: d=pd.read_csv(p,encoding=enc,low_memory=False); break
        except Exception: d=None
    if d is None: return None
    d.columns=[nfc(c) for c in d.columns]
    yc=next((c for c in d.columns if c.lower() in ("year","연도","년도")),None); rc="region" if "region" in d.columns else next((c for c in d.columns if c.lower() in ("시군구","지역")),None)
    d=d.rename(columns={yc:"year",rc:"region"})
    d=d.loc[:,~d.columns.duplicated()]                       # 같은 이름 컬럼(예: region 2개) 첫 것만
    if isinstance(d["region"],pd.DataFrame): d["region"]=d["region"].iloc[:,0]
    d["year"]=pd.to_numeric(d["year"],errors="coerce"); d["region"]=d["region"].astype(str).str.replace(r"\s+","",regex=True).map(nfc)
    d["region"]=d.region.replace("인천시미추홀구","인천시남구")
    if f.startswith("국민건강결과"): d=jr.fix_frame(d,tag="530 "+f)   # ★v2c
    return d[d.year.between(2011,2024)]
REG=pd.read_csv(os.path.join(CIN,"ST_curated46_panel.csv"))[["region","year"]]
REG["region"]=REG.region.astype(str).str.replace(r"\s+","",regex=True).map(nfc).replace("인천시미추홀구","인천시남구")
regs=sorted(REG.region.unique()); yrs=list(range(2011,2025)); G=pd.MultiIndex.from_product([regs,yrs],names=["region","year"]).to_frame(index=False)
# 분모
# ★ 2026-09-03 분모 정정(Z10): merged_하수도보급률 '총인구(명)'·'총면적(㎢)'은 2017+ '소계' 이중계상(전국 82M, 지역별 ×2.0; 2010/11/14/16 일부 지역도 중복) → 사용 금지.
#   pop = ST_agepop(연령별 인구 합; 연령표준화 E_std와 동일 분모) · area = 하수도 파일 총면적의 지역별 중앙값(연도불변)
ap=pd.read_csv(os.path.join(C,"ST_agepop.csv"),dtype={"sigungu_cd":str}); assert "sigungu_cd" in ap.columns, "ST_agepop.csv가 v2가 아님 — 400_agestd_assemble.py 먼저"; ap["region"]=ap.region.astype(str).str.replace(r"\s+","",regex=True).map(nfc).replace("인천시미추홀구","인천시남구")
den=ap.groupby(["region","year"])["pop"].sum().reset_index()
_ar=v2.area_by_region(regs); ar=pd.DataFrame({"region":_ar.index,"area":_ar.values})   # ★2026-10-21 shp 면적
G=G.merge(den,on=["region","year"],how="left").merge(ar,on="region",how="left")
_nat=G.groupby("year")["pop"].sum()/1e6; assert _nat.max()<60 and _nat.min()>45, f"Z10 pop 국가합 이상: {_nat.round(1).to_dict()}"
for c in ("pop","area"):
    assert G[c].notna().all() and (G[c]>0).all(), f"분모 {c} 빈칸·0 {int((G[c].isna()|(G[c]<=0)).sum())}칸 — v2는 채우지 않고 멈춤"   # ★2026-10-21
print(f"분모 복구 후 0/NA: pop {int(((G['pop']<=0)|G['pop'].isna()).sum())} · area {int(((G.area<=0)|G.area.isna()).sum())}")
print(f"격자 {len(G)}행 · 분모 결측 pop {G['pop'].isna().sum()} area {G.area.isna().sum()}")

cache={}; out=G.copy(); meta=[]
def getfile(f):
    if f not in cache: cache[f]=rd(f)
    return cache[f]
ZERO_VALID_FILES=("가축두수","국토이용","하수도보급률","하수관로","분뇨","지하수","1인가구","고령인구","종합소득")
ZERO_VALID_COLS=("전문의","요양기관수","의사수","간호사수","약사수","병상","대학","보육시설","유치원","초등학교","사업체","종사자")
def series(d,col,f=""):
    x=pd.to_numeric(d[col].astype(str).str.replace(",",""),errors="coerce")
    # 0 처리 규칙(Z1 정밀화): 조사·비율·행정지표는 0=미측정(NA). 두수·면적·개소·인구 같은 실측 카운트는 0이 유효값.
    if not (any(k in f for k in ZERO_VALID_FILES) or any(k in col for k in ZERO_VALID_COLS)): x=x.where(x!=0)
    if col=="도시지역인구비율": x=x.where(x<=100)   # 원자료 오류: 2020년 일반구 보유 대도시 10행이 분모 중복으로 100 초과 → NA (Z4가 적발)
    return d.assign(_v=x).groupby(["region","year"])["_v"].mean().reset_index()
# 기초생활수급: 6컬럼 coalesce
def basic_livelihood():
    parts=[]
    for f,cols in [("경제수준.csv",["국민기초생활보장수급자수율"]),("국민건강결과_최종.csv",["기초생활수급자수율","국민기초생활수급자수율","국민기초생활보장수급자수율"])]:
        d=getfile(f)
        for c in cols:
            if c in d.columns: parts.append(series(d,c,f))
    m=G[["region","year"]].copy()
    for p in parts: m=m.merge(p.rename(columns={"_v":f"v{len(m.columns)}"}),on=["region","year"],how="left")
    vc=[c for c in m.columns if c.startswith("v")]; m["_v"]=m[vc].bfill(axis=1).iloc[:,0]
    # 절대수(2011-17 국민기초생활수급자수)는 %로 변환해 보충
    d=getfile("국민건강결과_최종.csv")
    if "국민기초생활수급자수" in d.columns:
        a=series(d,"국민기초생활수급자수","국민건강결과_최종.csv").merge(G[["region","year","pop"]],on=["region","year"]); a["_v2"]=a._v/a["pop"]*100
        m=m.merge(a[["region","year","_v2"]],on=["region","year"],how="left"); m["_v"]=m._v.fillna(m._v2)
    return m[["region","year","_v"]]

# Z9 게이트: 이중 정규화 방지 — 원자료 컬럼명이 이미 비율/밀도인데 분모 나눗셈이 지정돼 있으면 중단
_bad=[r.label for _,r in I.iterrows() if any(k in str(r.col) for k in ("천명당","만명당","비율","율","%","률")) and any(k in str(r.normalization) for k in ("per 1000","per 100k","% of pop","% of district","per capita")) and not str(r.normalization).startswith("as-is")]
assert not _bad, f"Z9 이중 정규화 위험: {_bad}"
for _,r in I.iterrows():
    key=re.sub(r"[^a-z0-9]+","_",r.label.lower()).strip("_")[:40]
    f=nfc(r.file); col=nfc(r.col)
    if r.label=="Basic-livelihood recipient rate": s=basic_livelihood()
    else:
        d=getfile(f)
        if d is None or col not in d.columns: meta.append(dict(key=key,label=r.label,status="MISSING",note=f"{f}:{col}")); continue
        if f.startswith("merged_하수도보급률") and "% of district pop" in str(r.normalization):
            # ★ 하수도 파일 인구 계열도 2017+ 이중계상 → 같은 행 '총인구(명)'으로 나눈 내부비율(%)로 상쇄, 0=미측정
            dd=d.assign(_v=pd.to_numeric(d[col].astype(str).str.replace(",",""),errors="coerce")/pd.to_numeric(d["총인구(명)"].astype(str).str.replace(",",""),errors="coerce")*100)
            dd["_v"]=dd._v.where(dd._v.between(0,100)).where(lambda x:x>0); s=dd.groupby(["region","year"])["_v"].mean().reset_index(); r=r.copy(); r["normalization"]="as-is (% of district pop, within-source ratio)"; r["tier"]="B"   # 2018+만 실질 관측 → 구별 평균(시불변) 처리
        else: s=series(d,col,f)
    m=G.merge(s,on=["region","year"],how="left")
    n=str(r.normalization); v=m._v.copy()
    if n.startswith("as-is"): pass                                   # 원자료가 이미 비율/밀도 — 재정규화 금지 (Z9)
    elif "per km2" in n or "share of district" in n:
        if f.startswith("국토"): v=v/1e6                      # 국토이용현황은 ㎡ → ㎢
        v=v/m.area*(100 if "share" in n else 1)                # share = 구 면적 대비 %
    elif "per 1000" in n: v=v/m["pop"]*1000
    elif "per 100k" in n: v=v/m["pop"]*1e5
    elif "% of pop" in n or "% of district pop" in n: v=v/m["pop"].replace(0,np.nan)*100
    elif "% of households" in n or "% of dwellings" in n or "% of older pop" in n: v=v/m["pop"].replace(0,np.nan)*100   # 가구·주택 분모 부재 → 인구 대비로 근사, 라벨에 명시
    elif "per capita" in n and "KRW" in n: v=v/m["pop"]
    elif "persons per household" in n: v=m["pop"]/v
    # Tier B 또는 부분 실측 → 실측 연도 구별 평균을 시간불변으로
    v=v.replace([np.inf,-np.inf],np.nan); m["_n"]=v
    m["_n"]=yg.gate(m,key,r.label,col)   # ★v2d 연도 이상 게이트(보간·시간불변 평균 전)
    obs_years=m.loc[m._n.notna(),"year"].nunique()
    reg_mean=m.groupby("region")._n.mean()
    tinv=(r.tier=="B") or (obs_years<10) or (key in yg.ERRATIC)
    m["_final"]=m.region.map(reg_mean) if tinv else m.groupby("region")._n.transform(lambda s: s.interpolate(limit_direction="both"))
    out[key]=m._final.values
    meta.append(dict(key=key,label=r.label,domain=r.domain,direction=r.direction,tier=r.tier,source=f"{f}:{col}",obs_years=obs_years,time_invariant=tinv,
                     region_coverage=round(reg_mean.notna().mean(),3),median=round(float(np.nanmedian(m._final)),4),status="OK"))
M=pd.DataFrame(meta); V=[k for k in out.columns if k not in ("region","year","pop","area")]
print(f"\n조립: {len(V)}변수 · MISSING {int((M.status=='MISSING').sum())}")
for _,r in M[M.status=="MISSING"].iterrows(): print("   ✗",r.label,"←",r.note)

# ── 게이트 ──
print("\n═══ Z1–Z5 게이트 ═══"); fail=[]
for k in V:
    x=out[k]; mm=M[M.key==k].iloc[0]
    # Z1: 연도별 전무(값 없음) 비율 — 시간불변은 면제
    if not mm.time_invariant:
        zr=out.assign(_z=x.isna()).groupby("year")._z.mean()
        if (zr>0.9).sum()>0: fail.append((k,"Z1",f"{int((zr>0.9).sum())} years >90% missing")); continue
    # Z3 커버리지
    thr=0.60 if mm.time_invariant else 0.80
    if mm.region_coverage<thr: fail.append((k,"Z3",f"region coverage {mm.region_coverage:.2f} < {thr}")); continue
    # Z4 값 범위
    if not np.isfinite(x.dropna()).all(): fail.append((k,"Z4","non-finite values (inf)")); continue
    nz=str(I[I.label==mm.label].normalization.iloc[0])
    if nz.startswith("%") and "change" not in k and (x.max()>100.5 or x.min()<-0.5): fail.append((k,"Z4",f"range [{x.min():.1f},{x.max():.1f}] outside [0,100]")); continue
    if (x<0).any() and "change" not in k and "migration" not in k: fail.append((k,"Z4","negative values")); continue
# Z5 규모대리 상관
scale=out.groupby("region")[["pop","area"]].first(); z5flag={}
for k in V:
    if any(f[0]==k for f in fail): continue
    rm=out.groupby("region")[k].mean()
    for dn in ("pop","area"):
        rr=rm.corr(np.log1p(scale[dn]))
        if abs(rr)>0.7: z5flag[k]=f"|r| with log {dn} = {rr:.2f}"; break
F=pd.DataFrame(fail,columns=["key","gate","reason"]) if fail else pd.DataFrame(columns=["key","gate","reason"])
print(f"통과 {len(V)-len(F)} · 실패 {len(F)} · Z5 플래그(통과하되 도농/규모 축 주의) {len(z5flag)}")
for k,v in z5flag.items(): print(f"   ⚑ [Z5] {k:40s} {v}")
for _,r in F.iterrows(): print(f"   ✗ [{r.gate}] {r.key:40s} {r.reason}")
keep=[k for k in V if k not in set(F.key)]
out[["region","year","pop","area"]+keep].to_csv(os.path.join(C,"ST_prereg_panel.csv"),index=False,encoding="utf-8-sig")
M["gate"]=M.key.map(dict(zip(F.key,F.gate))).fillna("PASS"); M["z5_scale_flag"]=M.key.map(z5flag).fillna(""); M["gate_reason"]=M.key.map(dict(zip(F.key,F.reason))).fillna("")
M.to_csv(os.path.join(C,"ST_prereg_panel_META.csv"),index=False,encoding="utf-8-sig")
os.makedirs(v2.log_dir(),exist_ok=True); jr.write_log(os.path.join(v2.log_dir(),"JKMS_v2c_rawfix_530.csv")); yg.write(v2.log_dir())   # ★v2c 교정 기록
print(f"\n저장: ST_prereg_panel.csv ({len(keep)}변수) · ST_prereg_panel_META.csv")

# ── (2026-08-26 추가) 최종 단계: 한우·젖소를 농림어업총조사 시변 밀도로 교체 (548과 동일 규칙; 재조립 시 자동 적용) ──
# ★2026-10-21 v2: 545(총조사 밀도를 shp 면적으로 다시 계산) → 548 → 531(kchs_pc1) → 530b(v2 교정)를 같은 FBD_CACHE로 이어 실행
import subprocess as _sp
_here=os.path.dirname(os.path.abspath(__file__))
for _s in ("545_census_livestock.py","548_patch_panel_census.py","531_z6_vif_design.py","530b_v2_prereg_fix.py"): _sp.run(["python3",os.path.join(_here,_s)],check=True,cwd=_here)
