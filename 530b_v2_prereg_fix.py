# -*- coding: utf-8 -*-
# 530b — (2026-10-21) JKMS EHEC 원고 패널(ST_prereg_panel.csv)의 v2 교정. 530 → 545 → 548 → 531 다음에 530이 이어 부름(같은 FBD_CACHE).
#  계산은 _v2_refit_261021/11_build_inputs_jkms.py(F·정화조 구조적 0)와 10_build_inputs.py(EHEC kchs_pc1, 531 제조법 + 원자료 보충)를 그대로 옮김.
#  면적·인구 분모는 530이 이미 shp 면적·v2 인구로 계산하므로(옛 11의 P1 비율 환산 불필요) 여기서는 채움·교정만 함.
#   1) kchs_pc1(조사지표 PC1, 부호 = 농촌 +) = 531 제조법 21종 지역 평균, 빈 지역(부산진구·세종·제주시)은 원자료 parquet 코드 결합 평균(0 제외)으로 보충,
#      손씻기·비누는 원자료 코드 재구성 값(전 지역) → 229곳 모두 값(옛 판은 3곳 결측을 셋업이 평균 대입)
#   2) 원고 20변수 중 세종이 14년 모두 빈 10개 = 원자료 parquet 코드 결합(0 제외, 같은 지역 안 보간), 순이동은 /인구×1000; 독거노인(Tier B) = 1인가구60+/계60+ × 인구v2 60+ / 인구 × 100 실측 평균
#   3) 재정자립도·재정자주도 제주시·서귀포시 = KOSIS 제주특별자치도 값(외부 출처 대체)
#   4) 비누 손씻기(soap_sanitizer_use)와 손씻기 4종 = 원자료 코드 재구성(전 지역; 천안·창원·전주의 0 섞인 평균·부산진구·제주시 빈칸 교정)
#   5) 정화조 인구: 원자료 두 판 모두 14년 0이고 하수도 보급률 ≥99%인 48곳 = 구조적 0(explicit_fill zero; 530의 '0=미측정'이 결측으로 바꾼 것),
#      제주시 = 독립변수 3 판 merged_하수도보급률.csv 코드 결합(530 규칙)
#   6) 그 밖에 세종이 14년 모두 빈 변수 = 530과 같은 원자료 파일에서 세종 행을 코드로 찾아 530 규칙(0=미측정·정규화·시불변 평균/보간)으로 채움.
#      같은 계산을 다른 지역에 해서 패널 값을 재현(상대차 1e-6 안 지역 90% 이상, 또는 지역별 상대차 중앙값 1% 이하 = 원자료 판 차이 수준)하고
#      세종 값이 다른 지역 범위(최솟값/2-최댓값×2) 안인 변수만 채우고, 나머지는 결측으로 둠(기록; 세종 원자료 '숫자 두 번 이어붙음' 결함 차단).
#  채움 기록 → FBD_CACHE/_v2_logs/explicit_fill_JKMS_prereg.csv, 남은 결측 → _v2_logs/ST_prereg_panel_remaining_NA.csv
import os, sys, json, warnings, unicodedata as ud
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import v2_fixes as v2
import jkms_rawfix_v2c as jr   # ★2026-10-22 v2c: 530과 같은 원자료 교정(주택소유율 2022 척도·격리병상 미발표 연도)을 parquet에도 적용해 재현 검증을 맞춤
_PQC = jr.fix_frame(v2.parquet_coded(), key_col="sigungu_cd", tag="530b parquet")
import jkms_yeargate_v2d as yg   # ★2026-10-22 v2d: 530이 정한 연도 이상 게이트 결정(_v2_logs/JKMS_v2d_year_gate.json)을 세종 보충·재현 계산에도 같게 적용
YG = yg.load(v2.log_dir())
sk = v2.sk; C = v2.cache_dir(); LG = v2.log_dir(); CIN = os.path.join(v2.ROOT, "data_cache")
CODES = v2.CODES; NM = v2.NM; YY = v2.YY
PF = os.path.join(C, "ST_prereg_panel.csv")
PR0 = pd.read_csv(PF); COLS0 = PR0.columns.tolist()
PR = v2.code(PR0); assert len(PR) == 229 * 14 and not PR.duplicated(["sigungu_cd", "year"]).any()
CODE_OF = dict(zip(PR.region.astype(str), PR.sigungu_cd.astype(str)))   # ★v2d: 게이트 결정의 지역 이름 → 코드
META = pd.read_csv(os.path.join(C, "ST_prereg_panel_META.csv"), keep_default_na=False)
POOL = pd.read_csv(os.path.join(CIN, "COVARIATE_PREREG_POOL.csv")); POOL = POOL[POOL.status == "INCLUDE"]
FILL = []; LOG = []
def log(s): print(s); LOG.append(s)
pq = _PQC.copy().set_index(["sigungu_cd", "year"])   # ★v2c
HW = v2.hw_raw()
before_na = {c: int(PR[c].isna().sum()) for c in COLS0 if c not in ("region", "year") and PR[c].isna().any()}

# ── 1) kchs_pc1 (EHEC·JKMS 조사지표 PC1) ──
PM = META[META.gate == "PASS"]; Vk = [k for k in PM.key if k in PR.columns]; src = PM.set_index("key").source
X = PR.groupby("sigungu_cd")[Vk].mean()
survey = [k for k in Vk if ("국민건강결과" in str(src[k]) and str(src[k]).split(":")[1].endswith(("율", "률", "표준화율"))) or "handwash" in k or "soap" in k]
Xs0 = X[survey].dropna(axis=1, thresh=int(0.8 * len(X))); S21 = Xs0.columns.tolist()
def ehec_pc(Xs):
    p, v, ld = v2.pca1(Xs)
    if p.corr(X.urban_population_share.reindex(Xs.index)) > 0: p = -p
    return p, v
def raw_mean(colspec, key=None):
    c = colspec.split(":")[1]; pqr = pq.reset_index(); x = pd.Series(np.nan, index=pqr.index)
    _d = YG.get(key) if key else None; _ny = set(_d["zero_years"]) | set(_d["ratio_years"]) if _d else set()   # ★v2d: 게이트 NA 연도 제외
    for cc in [c, c.replace("_", "", 1), c.replace("_표준화율", "_조율")]:   # 같은 지표의 표기 변형 열을 순서대로 보충
        if cc in pqr.columns: y = pd.to_numeric(pqr[cc], errors="coerce"); x = x.fillna(y.where(y != 0))
    if _ny: x[pqr.year.isin(list(_ny))] = np.nan
    return x.groupby(pqr.sigungu_cd).mean()
X3 = Xs0.reindex(CODES).copy()
for k in S21:
    if k in v2.HWG: X3[k] = X3[k].fillna(HW[k]); continue
    na = X3[k].isna()
    if na.any():
        rm = raw_mean(src[k], k); X3.loc[na, k] = rm.reindex(X3.index[na]).values
        for cd in X3.index[na]: FILL.append(dict(model="JKMS_prereg", variable=f"kchs_pc1 입력 {k}", sigungu_cd=cd, region=NM[cd], year="2011-2024", value=round(float(X3.at[cd, k]), 3), how="원자료 지역 평균(0 제외), 코드 결합", source=src[k]))
assert X3.notna().all().all(), X3.isna().sum()[lambda s: s > 0]
XF = X3.copy()
for k in S21:
    if k in v2.HWG: XF[k] = HW[k].reindex(XF.index)
pc, pv = ehec_pc(XF)
old = PR.groupby("sigungu_cd").kchs_pc1.first().reindex(pc.index)
log(f"[kchs_pc1] 조사지표 {len(S21)}종 PC1 설명 분산 {pv:.3f} · 229곳 · 531 계산값(226곳)과 상관 {pc.corr(old):.3f}")
assert pc.corr(old) > 0
PR["kchs_pc1"] = PR.sigungu_cd.map(pc).values

# ── 2)-5) 원고 20변수 교정(11_build_inputs_jkms.py F·정화조와 같음) ──
RAW = v2.RAW
jj_indep = v2.kosis_jeju(v2.KOS_INDEP); jj_auto = v2.kosis_jeju(v2.KOS_AUTONOMY)
ap = pd.read_csv(os.path.join(C, "ST_agepop.csv"), dtype={"sigungu_cd": str}); ap = ap[ap.year.between(2011, 2024)]
p60 = ap[ap.age10.isin(["60~69세", "70세 이상"])].groupby(["sigungu_cd", "year"]).pop.sum()
SEJ = {"piped_water_supply_coverage": "상수도보급률", "sewerage_coverage": "하수도보급률", "fast_food_outlets_per_1000_pop": "인구천명당패스트푸드점수",
       "childcare_facilities_per_1000_infants": "유아천명당보육시설수", "sex_ratio": "성비", "current_smoking": "현재흡연율_표준화율",
       "influenza_vaccination_coverage": "연간인플루엔자예방접종률_표준화율", "physicians_per_1000_pop": "인구천명당의료기관종사의사수", "fiscal_independence": "재정자립도", "net_migration": "순이동인구"}
go = v2.code(pd.read_csv(os.path.join(RAW, "고령인구_전처리.csv")))
def interp(s):
    s = s.astype(float); ok = s.notna()
    if ok.sum() == 0: return s
    if ok.sum() == 1: return pd.Series(s[ok].iloc[0], index=s.index)
    return pd.Series(np.interp(np.arange(len(s)), np.where(ok)[0], s[ok].values), index=s.index)
d = PR.copy(); m = d.sigungu_cd == v2.SEJONG; yrs = d.loc[m, "year"].values; assert list(yrs) == YY
for v, col in SEJ.items():
    if not d.loc[m, v].isna().all(): continue
    s = pd.to_numeric(pq.loc[v2.SEJONG][col], errors="coerce").reindex(yrs); s = s.where(s != 0)
    if v == "net_migration": s = s / d.loc[m, "pop"].values * 1000
    s = interp(s.reset_index(drop=True)); d.loc[m, v] = s.values
    FILL.append(dict(model="JKMS_prereg", variable=v, sigungu_cd=v2.SEJONG, region=NM[v2.SEJONG], year="2011-2024", value=round(float(s.mean()), 3), how="원자료 코드 결합, 0 제외, 같은 지역 안 연도 보간(explicit)", source=f"독립변수 3/국민건강결과.parquet:{col}"))
if d.loc[m, "older_adults_living_alone"].isna().all():
    g = go[go.sigungu_cd == v2.SEJONG].set_index("year"); vals = []
    for y in g.index:
        if y in yrs: vals.append(g.at[y, "1인가구_60세 이상 - 계"] / g.at[y, "계_60세 이상 - 계"] * p60[(v2.SEJONG, y)] / float(d.loc[m & (d.year == y), "pop"].iloc[0]) * 100)
    d.loc[m, "older_adults_living_alone"] = float(np.mean(vals))
    FILL.append(dict(model="JKMS_prereg", variable="older_adults_living_alone", sigungu_cd=v2.SEJONG, region=NM[v2.SEJONG], year="2011-2024", value=round(float(np.mean(vals)), 3), how="원자료 비율 코드 결합, Tier B 실측 연도 평균", source="고령인구_전처리.csv + ST_agepop 60세 이상"))
for cd in v2.JEJU:
    mm = d.sigungu_cd == cd
    for v, jj, tb in [("fiscal_independence", jj_indep, "DT_1YL20921 재정자립도"), ("fiscal_autonomy", jj_auto, "DT_1YL20891 재정자주도")]:
        if v not in d.columns: continue
        assert d.loc[mm, v].isna().all(), f"{cd} {v} 이미 값 있음"
        d.loc[mm, v] = d.loc[mm, "year"].map(jj).values
        FILL.append(dict(model="JKMS_prereg", variable=v, sigungu_cd=cd, region=NM[cd], year="2011-2024", value=round(float(jj.mean()), 3), how="외부 출처 대체(explicit): 제주특별자치도 값을 두 행정시에 같게(행정시 값은 원자료에 없음)", source=f"KOSIS {tb} T10(세입과목개편전)"))
# 세종 1인당 GRDP: HM 패널과 같은 외부 출처 대체(KOSIS 세종 GRDP / 인구, 2015-2019 평균)
if "grdp_per_capita" in d.columns and d.loc[d.sigungu_cd == v2.SEJONG, "grdp_per_capita"].isna().all():
    d, gval = v2.fill_grdp_sejong(d)
    FILL.append(dict(model="JKMS_prereg", variable="grdp_per_capita", sigungu_cd=v2.SEJONG, region=NM[v2.SEJONG], year="2011-2024", value=round(gval, 3), how="외부 출처 대체(explicit): KOSIS 시도 GRDP/패널 인구 2015-2019 평균(시간불변, HM 패널과 같음)", source="KOSIS DT_1C81 세종특별자치시 명목 GRDP"))
# 손씻기·비누: 원자료 코드 재구성(전 지역)
d["soap_sanitizer_use"] = d.sigungu_cd.map(HW["handwashing_with_soap"]).values
for v in v2.HWG:
    if v in d.columns: d[v] = d.sigungu_cd.map(HW[v]).values
FILL.append(dict(model="JKMS_prereg", variable="soap_sanitizer_use·handwashing_* 4종", sigungu_cd="(229 전 지역)", region="", year="2011-2024", value="", how="원자료 코드 재구성(전 지역)", source="독립변수 3/국민건강결과.parquet 표준화율 우선·조율 보충, 0·>100 결측, 실측 연도 평균"))
# 정화조 인구: 구조적 0 + 제주시
def sewer_raw(path):
    s_ = pd.read_csv(path, low_memory=False); s_["region"] = s_.region.astype(str).str.replace(r"\s+", "", regex=True).map(v2.nfc).replace("인천시미추홀구", "인천시남구")
    s_ = s_[s_.year.between(2011, 2024)]; s_ = sk.attach_code(s_, "region", "year"); s_["sigungu_cd"] = s_.sigungu_cd.astype(str)
    num = lambda c: pd.to_numeric(s_[c].astype(str).str.replace(",", ""), errors="coerce")
    s_["sep"] = num("하수처리구역외_정화조인구"); r = s_.sep / num("총인구(명)") * 100; s_["r530"] = r.where(r.between(0, 100)).where(lambda x: x > 0)
    return s_.groupby("sigungu_cd").agg(sep_max=("sep", "max"), n_obs=("r530", lambda x: int(x.notna().sum())), r530=("r530", "mean"))
SW_ROOT = sewer_raw(os.path.join(v2.RAW0, "merged_하수도보급률.csv")); SW_IV3 = sewer_raw(os.path.join(RAW, "merged_하수도보급률.csv"))
_pv = d.groupby("sigungu_cd").septic_tank_population.first(); SEPNA = sorted(_pv[_pv.isna()].index); SW_SEW = d.groupby("sigungu_cd").sewerage_coverage.mean()
SEP0 = [c for c in SEPNA if SW_ROOT.sep_max.get(c, np.nan) == 0 and SW_IV3.sep_max.get(c, np.nan) == 0 and SW_SEW.get(c, 0) >= 99]
SEPJ = [c for c in SEPNA if c not in SEP0]; assert SEPJ == ["39010"], SEPJ
JEJU_SEP = float(SW_IV3.r530["39010"]); d.loc[d.sigungu_cd == "39010", "septic_tank_population"] = JEJU_SEP
FILL.append(dict(model="JKMS_prereg", variable="septic_tank_population", sigungu_cd="39010", region=NM["39010"], year="2011-2024", value=round(JEJU_SEP, 3), how="원자료 코드 결합(독립변수 3 판, 530 규칙, 실측 연도 평균) — FBD_DATA_ZIP 루트 판은 제주시 행이 전부 0(병합 누락)", source="독립변수 3/merged_하수도보급률.csv"))
zl = os.path.join(LG, "explicit_fill_JKMS_prereg_septic_zero.csv")
if os.path.exists(zl): os.remove(zl)
_ord = d.index.copy(); d["_ord"] = np.arange(len(d))
d = sk.explicit_fill(d, "septic_tank_population", "zero", "하수처리구역 밖 정화조 인구 0(하수도 완전 보급) — 530 규칙이 0을 결측으로 바꾼 구조적 0", log_path=zl).sort_values("_ord").drop(columns="_ord").reset_index(drop=True)
for c in SEP0: FILL.append(dict(model="JKMS_prereg", variable="septic_tank_population", sigungu_cd=c, region=NM[c], year="2011-2024", value=0.0, how="explicit_fill(how='zero'): 구조적 0", source="merged_하수도보급률.csv 정화조 인구 14년 0(두 판)·하수도 보급률 ≥99%"))
log(f"[septic] 빈칸 {len(SEPNA)}곳 = 구조적 0 {len(SEP0)}곳 + 제주시 {JEJU_SEP:.3f}%")

# ── 6) 그 밖에 세종이 14년 모두 빈 변수: 530 규칙으로 원자료 세종 행을 코드 결합(검증 통과 변수만) ──
CACHE_RAW = {}
def rd(f):
    if f in CACHE_RAW: return CACHE_RAW[f]
    if f == "국민건강결과_최종.csv":   # 루트 csv에는 세종 행이 없음 → 같은 지표의 독립변수 3 parquet(코드 결합; 11의 세종 채움과 같은 출처)
        dd = _PQC.copy(); dd.columns = [v2.nfc(c) for c in dd.columns]; CACHE_RAW[f] = dd; return dd
    p = os.path.join(v2.RAW0, f); dd = None
    for enc in ("utf-8-sig", "cp949"):
        try: dd = pd.read_csv(p, encoding=enc, low_memory=False); break
        except Exception: dd = None
    if dd is not None:
        dd.columns = [v2.nfc(c) for c in dd.columns]
        yc = next((c for c in dd.columns if c.lower() in ("year", "연도", "년도")), None); rc = "region" if "region" in dd.columns else next((c for c in dd.columns if c.lower() in ("시군구", "지역")), None)
        dd = dd.rename(columns={yc: "year", rc: "region"}); dd = dd.loc[:, ~dd.columns.duplicated()]
        dd["year"] = pd.to_numeric(dd["year"], errors="coerce"); dd["region"] = dd["region"].astype(str).str.replace(r"\s+", "", regex=True).map(v2.nfc).replace("인천시미추홀구", "인천시남구")
        dd = dd[dd.year.between(2011, 2024)]
        u = dd.region.unique(); cm = {}
        for nme in u:
            try:
                k = sk.name_kind([nme])[0]
                cm[nme] = sk.to_code([nme])[0] if k == "시군구" else None
            except Exception: cm[nme] = None
        dd["sigungu_cd"] = dd.region.map(cm); dd = dd[dd.sigungu_cd.notna()]
    CACHE_RAW[f] = dd; return dd
ZERO_VALID_FILES = ("가축두수", "국토이용", "하수도보급률", "하수관로", "분뇨", "지하수", "1인가구", "고령인구", "종합소득")
ZERO_VALID_COLS = ("전문의", "요양기관수", "의사수", "간호사수", "약사수", "병상", "대학", "보육시설", "유치원", "초등학교", "사업체", "종사자")
POPD = d.set_index(["sigungu_cd", "year"])["pop"]; AREAD = d.groupby("sigungu_cd").area.first()
def generic(k):
    mrow = META[META.key == k].iloc[0]; f, col = str(mrow.source).split(":", 1)
    if f.startswith("merged_하수도보급률") or mrow.label == "Basic-livelihood recipient rate" or "농림어업총조사" in f: return None, "특수 규칙 변수(하수도 내부 비율·기초생활 coalesce·총조사)"
    pr = POOL[POOL.label == mrow.label]
    if len(pr) != 1: return None, "풀 정의 없음"
    n = str(pr.normalization.iloc[0]); dd = rd(f)
    if dd is None or col not in dd.columns: return None, f"원자료 열 없음 {f}:{col}"
    x = pd.to_numeric(dd[col].astype(str).str.replace(",", ""), errors="coerce")
    if not (any(z in f for z in ZERO_VALID_FILES) or any(z in col for z in ZERO_VALID_COLS)): x = x.where(x != 0)
    if col == "도시지역인구비율": x = x.where(x <= 100)
    s = dd.assign(_v=x).groupby(["sigungu_cd", "year"])["_v"].mean()
    g = pd.MultiIndex.from_product([CODES, YY], names=["sigungu_cd", "year"]); v = s.reindex(g)
    pop = POPD.reindex(g).values; area = AREAD.reindex(g.get_level_values(0)).values
    if n.startswith("as-is"): pass
    elif "per km2" in n or "share of district" in n:
        if f.startswith("국토"): v = v / 1e6
        v = v / area * (100 if "share" in n else 1)
    elif "per 1000" in n: v = v / pop * 1000
    elif "per 100k" in n: v = v / pop * 1e5
    elif "% of pop" in n or "% of district pop" in n: v = v / np.where(pop == 0, np.nan, pop) * 100
    elif "% of households" in n or "% of dwellings" in n or "% of older pop" in n: v = v / np.where(pop == 0, np.nan, pop) * 100
    elif "per capita" in n and "KRW" in n: v = v / pop
    elif "persons per household" in n: v = pop / v
    v = pd.Series(np.asarray(v, dtype=float), index=g).replace([np.inf, -np.inf], np.nan)
    v, _er = yg.apply_coded(v, k, YG, CODE_OF)   # ★v2d
    tinv = str(mrow.time_invariant) in ("True", "true", "1") or _er
    if tinv: fin = v.groupby(level=0).transform("mean")
    else: fin = v.groupby(level=0).transform(lambda s_: s_.interpolate(limit_direction="both"))
    return fin, n
GEN = []
for k in [c for c in COLS0 if c not in ("region", "year", "pop", "area", "kchs_pc1", "soap_sanitizer_use")]:
    if k not in META.key.values: continue
    ms = d.sigungu_cd == v2.SEJONG
    if not d.loc[ms, k].isna().all(): continue
    fin, why = generic(k)
    if fin is None: GEN.append(dict(variable=k, filled=False, reason=why)); continue
    ref = d.set_index(["sigungu_cd", "year"])[k]; cmpi = ref.dropna().index; cmpi = cmpi[cmpi.get_level_values(0) != v2.SEJONG]
    a_ = fin.reindex(cmpi); b_ = ref.reindex(cmpi); rel = ((a_ - b_).abs() / b_.abs().clip(lower=1e-9))
    okreg = (rel.groupby(level=0).max() < 1e-6).mean(); rmed = float(rel.groupby(level=0).mean().median()) if len(rel) else float("nan")
    sj = fin.loc[v2.SEJONG]
    oth = ref.drop(index=v2.SEJONG, level=0, errors="ignore").groupby(level=0).mean(); lo_, hi_ = float(oth.min()), float(oth.max())
    plaus = sj.notna().all() and (sj.mean() >= (lo_ / 2 if lo_ > 0 else lo_ * 2 - 1e-9)) and (sj.mean() <= (hi_ * 2 if hi_ > 0 else hi_ / 2 + 1e-9))   # 세종 원자료 '숫자 두 번 이어붙음' 결함 등 차단
    accept = (okreg >= 0.90 or rmed <= 0.01)   # 다른 지역 재현: 상대차 1e-6 안 지역 90% 이상, 또는 지역별 평균 상대차의 중앙값 1% 이하(원자료 판 차이 수준)
    if accept and plaus:
        d.loc[ms, k] = sj.reindex(d.loc[ms, "year"]).values
        GEN.append(dict(variable=k, filled=True, reason=f"다른 지역 재현 {okreg:.3f}(상대차 중앙 {rmed:.2e})", value_mean=round(float(sj.mean()), 4), normalization=why))
        FILL.append(dict(model="JKMS_prereg", variable=k, sigungu_cd=v2.SEJONG, region=NM[v2.SEJONG], year="2011-2024", value=round(float(sj.mean()), 4), how=f"530 규칙 재계산, 원자료 세종 행 코드 결합(다른 지역 재현 {okreg:.3f}, 상대차 중앙 {rmed:.1e}; 국민건강결과는 독립변수 3 parquet 판)", source=str(META[META.key == k].source.iloc[0])))
    else:
        GEN.append(dict(variable=k, filled=False, reason=f"다른 지역 재현 {okreg:.3f}(상대차 중앙 {rmed:.2e}), 세종 원자료 빈칸 {int(sj.isna().sum())}, 세종 평균 {float(sj.mean()) if sj.notna().any() else float('nan'):.4g} (다른 지역 범위 {lo_:.4g}-{hi_:.4g})"))
pd.DataFrame(GEN).to_csv(os.path.join(LG, "ST_prereg_panel_sejong_generic.csv"), index=False, encoding="utf-8-sig")
log(f"[세종 일반 규칙] 대상 {len(GEN)}개 · 채움 {sum(g['filled'] for g in GEN)}개 · 남김 {sum(not g['filled'] for g in GEN)}개 → _v2_logs/ST_prereg_panel_sejong_generic.csv")

# ── 7) v2b(2026-10-21 정지 지점 해결): 원고 다중우주(533)·1단계(532) 풀의 남은 결측 14개 ──
#  결정 순서 (1) 원자료(또는 그 원출처)에서 sigungu_cd 결합으로 되살림 → (2) 구조적 결측은 외부 출처 대체 → (3) 그래도 남으면 채우지 않고 풀에서 뺌(_v2_promote_261021/15_세트풀_제외변수.csv). 평균 채움 금지.
#  되살린 값은 530 규칙(0=미측정·정규화·시불변 평균/같은 지역 보간)을 그대로 적용하고, 같은 계산이 다른 지역의 패널 값을 재현하는지(상대차 1e-6 안 지역 비율) 확인한 뒤에만 넣음.
#  국민건강결과는 루트 csv에 세종 행이 없어 독립변수 3 parquet 세종 행을 씀. 루트 csv와 parquet은 2011-2023 값이 같고(10/21 대조) 루트 csv에 2024가 없으므로 parquet도 2023까지만 씀.
V2B = []
YY23 = [y for y in YY if y <= 2023]
def _rule530(k, raw, keep_regions=False):
    """raw: (sigungu_cd, year) 인덱스의 원자료 값 → 530 규칙 최종값(229 × 14)."""
    mrow = META[META.key == k].iloc[0]; f, col = str(mrow.source).split(":", 1); n = str(POOL[POOL.label == mrow.label].normalization.iloc[0])
    x = pd.to_numeric(raw, errors="coerce")
    if not (any(z in f for z in ZERO_VALID_FILES) or any(z in col for z in ZERO_VALID_COLS)): x = x.where(x != 0)
    g = pd.MultiIndex.from_product([CODES, YY], names=["sigungu_cd", "year"]); v = x.groupby(level=[0, 1]).mean().reindex(g)
    pop = POPD.reindex(g).values; area = AREAD.reindex(g.get_level_values(0)).values
    if n.startswith("as-is"): pass
    elif "per km2" in n or "share of district" in n:
        if f.startswith("국토"): v = v / 1e6
        v = v / area * (100 if "share" in n else 1)
    elif "per 1000" in n: v = v / pop * 1000
    elif "per 100k" in n: v = v / pop * 1e5
    elif "% of pop" in n or "% of district pop" in n or "% of households" in n or "% of dwellings" in n or "% of older pop" in n: v = v / np.where(pop == 0, np.nan, pop) * 100
    elif "per capita" in n and "KRW" in n: v = v / pop
    elif "persons per household" in n: v = pop / v
    v = pd.Series(np.asarray(v, dtype=float), index=g).replace([np.inf, -np.inf], np.nan)
    v, _er = yg.apply_coded(v, k, ({**YG, k: {**YG[k], "na_regions": []}} if (keep_regions and k in YG) else YG), CODE_OF)   # ★v2d(keep_regions: 게이트 (c)로 비운 지역을 다른 판 원자료로 되살릴 때)
    tinv = str(mrow.time_invariant) in ("True", "true", "1") or _er
    return v.groupby(level=0).transform("mean") if tinv else v.groupby(level=0).transform(lambda s_: s_.interpolate(limit_direction="both"))
def _repro(k, fin, who):
    ref = d.set_index(["sigungu_cd", "year"])[k]; idx = ref.dropna().index; idx = idx[~idx.get_level_values(0).isin(who)]
    rel = ((fin.reindex(idx) - ref.reindex(idx)).abs() / ref.reindex(idx).abs().clip(lower=1e-9)).groupby(level=0).max()
    oth = ref.drop(index=list(who), level=0, errors="ignore").groupby(level=0).mean()
    return float((rel < 1e-6).mean()), float(oth.min()), float(oth.max())
def _put(k, fin, who, rule, how, src, min_repro=0.90):   # 530b 6단계와 같은 재현 기준(지역 90% 이상; 나머지는 루트 csv·parquet 판 차이)
    rp, lo_, hi_ = _repro(k, fin, who)
    vals = {cd: fin.loc[cd].reindex(YY) for cd in who}
    ok = rp >= min_repro and all(s_.notna().all() for s_ in vals.values()) and all((lo_ / 2 if lo_ > 0 else lo_ * 2 - 1e-9) <= s_.mean() <= (hi_ * 2 if hi_ > 0 else hi_ / 2 + 1e-9) for s_ in vals.values())
    for cd, s_ in vals.items():
        mm = d.sigungu_cd == cd; assert d.loc[mm, k].isna().all(), (k, cd)
        if ok:
            d.loc[mm, k] = s_.reindex(d.loc[mm, "year"]).values
            FILL.append(dict(model="JKMS_prereg", variable=k, sigungu_cd=cd, region=NM[cd], year="2011-2024", value=round(float(s_.mean()), 4), how=f"[v2b 규칙 ({rule})] {how} — 같은 계산의 다른 지역 재현 {rp:.3f}", source=src))
    V2B.append(dict(variable=k, rule=rule if ok else "3", method=how, regions=", ".join(NM[c] for c in who), cells=14 * len(who), filled=ok, other_region_reproduction=round(rp, 4),
                    value_mean="; ".join(f"{NM[c]} {vals[c].mean():.4g}" for c in who), other_range=f"{lo_:.4g}-{hi_:.4g}", source=src))
    log(f"[v2b] {k}: {'채움' if ok else '채우지 않음(검증 실패)'} ({rule}) 재현 {rp:.3f} · " + "; ".join(f"{NM[c]} {vals[c].mean():.4g}" for c in who) + f" (다른 지역 {lo_:.4g}-{hi_:.4g})")
    return ok
PQ23 = _PQC.copy(); PQ23 = PQ23[PQ23.year <= 2023].set_index(["sigungu_cd", "year"])
def _pq(col, alias=None, who=(v2.SEJONG,)):
    """parquet ≤2023 원자료 열. alias가 있으면 세종(who) 행에 한해, 다른 지역이 그 해 본 열에 실측을 가진 해에만 표기 변형 열로 보충.
       본 열이 그 해 다른 지역 모두 0(연도 전체 0, 미발표 해)이면 세종도 같은 해 0(다른 지역과 같은 처리)."""
    x = pd.to_numeric(PQ23[col], errors="coerce").copy()
    if alias is None: return x
    y = pd.to_numeric(PQ23[alias], errors="coerce"); oth = x[~x.index.get_level_values(0).isin(who)]
    for yr in YY23:
        o = oth.xs(yr, level=1); nz = o.notna() & (o != 0)
        for cd in who:
            if not pd.isna(x.get((cd, yr), np.nan)): continue
            if o.notna().mean() >= 0.98 and (o.fillna(0) == 0).mean() >= 0.98: x.loc[(cd, yr)] = 0.0          # 연도 전체 0(미발표)인 해
            elif nz.mean() >= 0.5 and (cd, yr) in y.index: x.loc[(cd, yr)] = y.loc[(cd, yr)]                  # 본 열이 자료를 담은 해
    return x
SEJ1 = (v2.SEJONG,)
# (1) 국민건강결과(parquet ≤2023) 세종 행, 530 규칙 그대로
for k, col in [("hospital_beds_per_1000_pop", "인구천명당의료기관병상수"), ("businesses_per_1000_pop", "인구천명당사업체수"), ("public_health_center_use", "연간보건기관이용률_표준화율"), ("home_ownership_rate", "주택소유율")]:
    if k in d.columns and d.loc[d.sigungu_cd == v2.SEJONG, k].isna().all():
        _put(k, _rule530(k, _pq(col)), SEJ1, "1", f"원자료 세종 행 코드 결합(독립변수 3 parquet '{col}', 2011-2023), 530 규칙" + (" — 2022년 원자료는 전 지역이 0-1 비율 척도(다른 지역과 같은 처리)" if k == "home_ownership_rate" else ""), f"독립변수 3/국민건강결과.parquet:{col}")
# (1) 같은 지표의 표기 변형 열(세종 행만 변형 열에 값이 있음)
for k, col, al in [("obesity_self_reported", "비만율_자가보고_표준화율", "비만율자가보고_표준화율"), ("isolation_beds", "격리병실_병상수", "격리병실병상수")]:
    if k in d.columns and d.loc[d.sigungu_cd == v2.SEJONG, k].isna().all():
        if k == "isolation_beds":   # ★v2c: 행 coalesce·미발표 연도 처리는 jr.fix_frame이 이미 함(2019 변형 열, 2018·2023 결측 → 530 보간)
            _put(k, _rule530(k, _pq(col)), SEJ1, "1", f"원자료 세종 행 코드 결합: 본 열 '{col}'이 빈 행은 표기 변형 열 '{al}'로 보충, 본 열 전 지역 0인 해는 v2c 규칙(2019 변형 열, 2018·2023 결측 → 같은 지역 보간), 530 규칙", f"독립변수 3/국민건강결과.parquet:{col}|{al}"); continue
        _put(k, _rule530(k, _pq(col, al)), SEJ1, "1", f"원자료 세종 행 코드 결합: 본 열 '{col}'이 비어 같은 지표의 표기 변형 열 '{al}'로 보충(다른 지역이 본 열에 실측을 가진 해만, 연도 전체 0인 해는 다른 지역과 같이 0), 530 규칙", f"독립변수 3/국민건강결과.parquet:{col}|{al}")
# (1) 1인가구: 루트 원자료 세종 행의 '숫자 두 번 이어붙음' 결함을 563과 같은 규칙(짝수 자릿수·앞뒤 절반 같음 → 절반)으로 복원
k = "single_person_households"
if k in d.columns and d.loc[d.sigungu_cd == v2.SEJONG, k].isna().all():
    dd1 = rd("1인가구_전처리.csv").copy(); x = pd.to_numeric(dd1["1인가구"].astype(str).str.replace(",", ""), errors="coerce")
    def _half(val):
        if pd.isna(val) or abs(val) < 1e7 or val != round(val): return val
        s_ = str(int(abs(val))); return float(s_[:len(s_) // 2]) if len(s_) % 2 == 0 and s_[:len(s_) // 2] == s_[len(s_) // 2:] else val
    nh = int(((x != x.map(_half)) & x.notna()).sum()); x = x.map(_half)
    raw = pd.Series(x.values, index=pd.MultiIndex.from_arrays([dd1.sigungu_cd.astype(str).values, dd1.year.astype(int).values], names=["sigungu_cd", "year"]))
    _put(k, _rule530(k, raw), SEJ1, "1", f"루트 원자료 1인가구_전처리.csv 세종 행 코드 결합, 두 번 이어붙은 값 {nh}칸을 절반으로 복원(563과 같은 규칙), 530 규칙", "FBD_DATA_ZIP/1인가구_전처리.csv:1인가구")
# (1) 기초생활수급자 비율: 530 coalesce 규칙을 코드 결합으로(경제수준.csv 세종 행은 이름 표기 '세종특별자치시세종특별자치시'로 빠졌음; 국민건강결과는 세종만 parquet ≤2023)
k = "basic_livelihood_recipient_rate"
if k in d.columns and d.loc[d.sigungu_cd == v2.SEJONG, k].isna().all():
    ec = rd("경제수준.csv"); hroot = None
    for enc in ("utf-8-sig", "cp949"):
        try: hroot = pd.read_csv(os.path.join(v2.RAW0, "국민건강결과_최종.csv"), encoding=enc, low_memory=False); break
        except Exception: hroot = None
    hroot.columns = [v2.nfc(c) for c in hroot.columns]; hroot = hroot[hroot.year.between(2011, 2024)]; hroot = v2.code(hroot)
    hs = PQ23.reset_index(); hs = hs[hs.sigungu_cd == v2.SEJONG]; hh = pd.concat([hroot, hs], ignore_index=True)
    def _ser(dd_, col):
        x_ = pd.to_numeric(dd_[col].astype(str).str.replace(",", ""), errors="coerce"); x_ = x_.where(x_ != 0)
        return pd.Series(x_.values, index=pd.MultiIndex.from_arrays([dd_.sigungu_cd.astype(str).values, dd_.year.astype(int).values])).groupby(level=[0, 1]).mean()
    g_ = pd.MultiIndex.from_product([CODES, YY], names=["sigungu_cd", "year"]); parts = []
    if "국민기초생활보장수급자수율" in ec.columns: parts.append(_ser(ec, "국민기초생활보장수급자수율").reindex(g_))
    for c_ in ["기초생활수급자수율", "국민기초생활수급자수율", "국민기초생활보장수급자수율"]:
        if c_ in hh.columns: parts.append(_ser(hh, c_).reindex(g_))
    v_ = pd.concat(parts, axis=1).bfill(axis=1).iloc[:, 0]
    if "국민기초생활수급자수" in hh.columns: v_ = v_.fillna(_ser(hh, "국민기초생활수급자수").reindex(g_) / POPD.reindex(g_) * 100)
    fin = v_.groupby(level=0).transform(lambda s_: s_.interpolate(limit_direction="both"))   # 530: as-is, 시변(보간)
    _put(k, fin, SEJ1, "1", "530 기초생활 coalesce 규칙을 코드 결합으로 재계산(경제수준.csv 세종 행 '세종특별자치시세종특별자치시' 이름 결합 누락, 국민건강결과 세종 행은 parquet ≤2023)", "FBD_DATA_ZIP/경제수준.csv:국민기초생활보장수급자수율 + 국민건강결과(기초생활 열 3종·수급자수/인구)")
# (1) 종합소득세: 원자료 종합소득세.csv의 원출처(국세청 KOSIS DT_133N_A3212, 다른 지역 값 일치)에서 세종 행을 코드 결합
kt = pd.read_csv(v2.KOS_INCOMETAX, dtype={"C1": str}); sdn = kt[kt.C1.str.len() == 3].drop_duplicates("C1").set_index("C1").C1_NM
kt["region"] = np.where(kt.C1.str.len() == 5, kt.C1.str[:3].map(sdn) + " " + kt.C1_NM, np.where(kt.C1 == "A09", "세종특별자치시", None)); kt = kt[kt.region.notna()]
import contextlib as _cl, io as _io
_u = kt[["region"]].drop_duplicates(); _bad = []
for _r in _u.region:
    try:
        with _cl.redirect_stdout(_io.StringIO()): sk.to_code([_r])
    except Exception: _bad.append(_r)
kt = kt[~kt.region.isin(_bad)]
with _cl.redirect_stdout(_io.StringIO()): _cc = sk.attach_code(kt[["region"]].drop_duplicates(), "region", None)
kt["sigungu_cd"] = kt.region.map(dict(zip(_cc.region, _cc.sigungu_cd.astype(str))))
for k, item in [("income_tax_filers", "신고인원"), ("aggregate_taxable_income", "종합소득금액")]:
    if k in d.columns and d.loc[d.sigungu_cd == v2.SEJONG, k].isna().all():
        s_ = kt[kt.C2_NM == item]; raw = pd.to_numeric(s_.DT, errors="coerce"); raw = pd.Series(raw.values, index=pd.MultiIndex.from_arrays([s_.sigungu_cd.values, s_.PRD_DE.astype(int).values])).groupby(level=[0, 1]).sum(min_count=1)
        _put(k, _rule530(k, raw), SEJ1, "1", f"원자료 종합소득세.csv의 원출처 KOSIS DT_133N_A3212 '{item}' 세종 행 코드 결합(원자료에는 세종 행 없음; 다른 지역 값은 원자료와 같음), 530 규칙(Tier B 2018-2024 평균)" + (f"; 코드화 못 한 KOSIS 표기 {len(_bad)}개는 대조에서만 제외" if _bad else ""), f"KOSIS DT_133N_A3212(국세청) {item}")
# (1) 내과 의원 수: 옹진군 원자료 '의원_내과'는 14년 모두 0(의원 합계 2-3곳 모두 내과 아님) → 530의 '0=미측정' 규칙이 결측으로 바꾼 실측 0. 구조적 0으로 명시(explicit_fill zero)
k = "internal_medicine_clinics"; ONGJIN = sk.to_code(["인천광역시 옹진군"])[0]
if k in d.columns and d.loc[d.sigungu_cd == ONGJIN, k].isna().all():
    o_ = PQ23.xs(ONGJIN, level=0); im = pd.to_numeric(o_["의원_내과"], errors="coerce"); tot = pd.to_numeric(o_["의원_합계"], errors="coerce")
    ok0 = bool((im.fillna(-1) == 0).all() and (tot > 0).all())
    if ok0:
        mm = d.sigungu_cd == ONGJIN; d.loc[mm, k] = 0.0
        FILL.append(dict(model="JKMS_prereg", variable=k, sigungu_cd=ONGJIN, region=NM[ONGJIN], year="2011-2024", value=0.0, how="[v2b 규칙 (1)] explicit_fill(how='zero'): 원자료 '의원_내과' 14년 0이고 '의원_합계'는 " + ",".join(str(int(t)) for t in tot.values) + "곳(내과 아닌 의원만) → 530 '0=미측정' 규칙이 결측으로 바꾼 실측 0", source="독립변수 3/국민건강결과.parquet:의원_내과·의원_합계"))
    V2B.append(dict(variable=k, rule="1" if ok0 else "3", method="원자료 실측 0 복원(explicit_fill zero; 의원 합계 2-3곳 중 내과 0)", regions=NM[ONGJIN], cells=14, filled=ok0, other_region_reproduction="", value_mean=f"{NM[ONGJIN]} 0", other_range="", source="독립변수 3/국민건강결과.parquet:의원_내과·의원_합계"))
    log(f"[v2b] {k}: 옹진군 {'구조적 0 명시' if ok0 else '전제 불일치 — 채우지 않음'}")
# ── 8) v2d(2026-10-22 연도 이상 게이트 뒤처리): 게이트가 비운 칸 중 원자료로 되살릴 수 있는 것 ──
#  (가) 게이트 (c)(비율형 전 연도 0 지역 → NA) 지역은 같은 이름의 독립변수 3 판 원자료를 sigungu_cd로 결합해 530 규칙으로 다시 계산(다른 지역 재현 ≥ 0.90일 때만)
#  (나) 내과 의원 수: 게이트가 비운 연도(2018·2019 지역 어긋남)를 뺀 모든 해에 '의원_내과' 0이고 '의원_합계' > 0인 지역 = 실측 0(옹진과 같은 규칙)
def rd3(f):
    p = os.path.join(v2.RAW, f); dd = None
    for enc in ("utf-8-sig", "cp949"):
        try: dd = pd.read_csv(p, encoding=enc, low_memory=False); break
        except Exception: dd = None
    if dd is None: return None
    dd.columns = [v2.nfc(c) for c in dd.columns]
    yc = next((c for c in dd.columns if c.lower() in ("year", "연도", "년도")), None); dd = dd.rename(columns={yc: "year"}); dd = dd.loc[:, ~dd.columns.duplicated()]
    dd["year"] = pd.to_numeric(dd["year"], errors="coerce"); dd["region"] = dd["region"].astype(str).str.replace(r"\s+", "", regex=True).map(v2.nfc).replace("인천시미추홀구", "인천시남구")
    dd = dd[dd.year.between(2011, 2024)]; cm = {}
    for nme in dd.region.unique():
        try: cm[nme] = sk.to_code([nme])[0] if sk.name_kind([nme])[0] == "시군구" else None
        except Exception: cm[nme] = None
    dd["sigungu_cd"] = dd.region.map(cm); return dd[dd.sigungu_cd.notna()]
V2D = []
for k, dk in YG.items():
    if not dk.get("na_regions") or k not in d.columns: continue
    who = tuple(c for c in (CODE_OF.get(r) for r in dk["na_regions"]) if c and d.loc[d.sigungu_cd == c, k].isna().all())
    if not who: continue
    f, col = str(META[META.key == k].source.iloc[0]).split(":", 1); dd = rd3(f)
    if dd is None or col not in dd.columns:
        V2D.append(dict(variable=k, regions=", ".join(NM[c] for c in who), filled=False, how=f"독립변수 3 판에 {f}:{col} 없음")); continue
    raw = pd.Series(pd.to_numeric(dd[col].astype(str).str.replace(",", ""), errors="coerce").values, index=pd.MultiIndex.from_arrays([dd.sigungu_cd.astype(str).values, dd.year.astype(int).values]))
    raw = raw.where(raw != 0)   # 이 지역들은 0이 미측정(루트 판 병합 누락)이라 0 제외
    fin = _rule530(k, raw, keep_regions=True)
    ref = d.set_index(["sigungu_cd", "year"])[k]; ci_ = ref.dropna().index; ci_ = ci_[~ci_.get_level_values(0).isin(who)]
    rel = ((fin.reindex(ci_) - ref.reindex(ci_)).abs() / ref.reindex(ci_).abs().clip(lower=1e-9)); rmed = float(rel.groupby(level=0).mean().median())
    # 재현 기준은 530b 6단계와 같음: 상대차 1e-6 안 지역 90% 이상, 또는 지역별 평균 상대차의 중앙값 1% 이하(원자료 판 차이 수준)
    ok = _put(k, fin, who, "1", f"[v2d] 게이트 (c)로 비운 지역(루트 판 전 연도 0 = 병합 누락)을 독립변수 3 판 {f} '{col}' 코드 결합(0 제외)으로 530 규칙 재계산(지역별 상대차 중앙 {rmed:.2e})", f"독립변수 3/{f}:{col}", min_repro=(0.0 if rmed <= 0.01 else 0.90))
    V2D.append(dict(variable=k, regions=", ".join(NM[c] for c in who), filled=ok, how=f"독립변수 3 판 코드 결합, 다른 지역 상대차 중앙 {rmed:.3g}"))
k = "internal_medicine_clinics"
if k in d.columns:
    ny = set(YG.get(k, {}).get("zero_years", [])) | set(YG.get(k, {}).get("ratio_years", []))
    for cd in sorted(d.loc[d[k].isna(), "sigungu_cd"].unique()):
        if not d.loc[d.sigungu_cd == cd, k].isna().all(): continue
        o_ = PQ23.xs(cd, level=0); o_ = o_[~o_.index.isin(list(ny))]; im = pd.to_numeric(o_["의원_내과"], errors="coerce"); tot = pd.to_numeric(o_["의원_합계"], errors="coerce")
        ok0 = bool(len(o_) and (im.fillna(-1) == 0).all() and (tot > 0).all())
        ok00 = bool(len(o_) and (im.fillna(-1) == 0).all() and (tot.fillna(-1) == 0).all())   # 의원 자체가 매년 0으로 기록된 지역(울릉군 등): 내과 의원 0 of 0
        if ok0 or ok00:
            d.loc[d.sigungu_cd == cd, k] = 0.0
            FILL.append(dict(model="JKMS_prereg", variable=k, sigungu_cd=cd, region=NM[cd], year="2011-2024", value=0.0, how=f"[v2d] explicit_fill(how='zero'): 게이트 연도 {sorted(ny)}를 뺀 모든 해 '의원_내과' 0, '의원_합계' " + ("> 0 → 실측 0(옹진과 같은 규칙)" if ok0 else "도 매년 0(의원 없음, 보건의료원만) → 내과 의원 0 of 0"), source="독립변수 3/국민건강결과.parquet:의원_내과·의원_합계"))
        V2D.append(dict(variable=k, regions=NM[cd], filled=ok0 or ok00, how="실측 0 확인" if ok0 else ("의원 합계도 매년 0 → 0 of 0" if ok00 else "전제 불일치 — 채우지 않음")))
pd.DataFrame(V2D).to_csv(os.path.join(LG, "ST_prereg_panel_v2d_gate_restore.csv"), index=False, encoding="utf-8-sig")
for r_ in V2D: log(f"[v2d] {r_['variable']} {r_['regions']}: {'채움' if r_['filled'] else '채우지 않음'} ({r_['how']})")
# (2)·(3) 원자료에도 원출처에도 같은 정의의 값이 없는 변수 → 채우지 않음(풀에서 뺌, 15_세트풀_제외변수.csv)
for k, why in [("unemployment_rate", "원자료 경제수준.csv 실업률은 2012-2019 시·군(155곳)만 있고 광역시 자치구·군 74곳은 전 기간 없음(지역별고용조사 시군 단위). KOSIS 시군구 실업률(INH_1ES3A01S)은 자치구를 2021년 이후만 제공해 같은 기간 정의로 되살릴 수 없음. 창원·군위·제주시·서귀포시 4곳은 원자료에 있으나(이름 결합 누락) 74곳이 남아 (3)"),
               ("employment_rate", "실업률과 같은 원자료·같은 74곳 구조적 결측(지역별고용조사 시군 단위) → (3)"),
               ("welfare_share_of_general_budget", "제주시·서귀포시는 행정시라 원자료(일반회계 중 사회복지예산 비중)에 값이 없고, KOSIS에 같은 정의의 제주특별자치도 예산 비중 표가 없음(있는 표 DT_1YL20912는 '결산' 비중·2020년부터라 정의·기간이 다름) → (2) 불가, (3)")]:
    if k in d.columns and d[k].isna().any():
        V2B.append(dict(variable=k, rule="3", method="채우지 않음(풀에서 제외)", regions=", ".join(sorted(set(d.loc[d[k].isna(), "region"])))[:200], cells=int(d[k].isna().sum()), filled=False, other_region_reproduction="", value_mean="", other_range="", source=why))
        log(f"[v2b] {k}: 결측 {int(d[k].isna().sum())}칸 남김 → 풀에서 제외(3)")
pd.DataFrame(V2B).to_csv(os.path.join(LG, "ST_prereg_panel_v2b_pool_decisions.csv"), index=False, encoding="utf-8-sig")

# ── 출력 ──
cols =[c for c in COLS0 if c != "soap_sanitizer_use"] + ["soap_sanitizer_use"]
out = d[cols]
rem = {c: int(out[c].isna().sum()) for c in cols if out[c].isna().any()}
pd.DataFrame([dict(variable=k, na_cells=v, regions=", ".join(sorted(set(out.loc[out[k].isna(), "region"])))[:300], na_before=before_na.get(k, 0)) for k, v in rem.items()]).to_csv(os.path.join(LG, "ST_prereg_panel_remaining_NA.csv"), index=False, encoding="utf-8-sig")
log(f"[남은 결측] {len(rem)}개 변수: " + ", ".join(f"{k} {v}" for k, v in rem.items()))
g = v2.panel_gate(v2.code(out[["region", "year", "pop", "area"]]), cols=["pop", "area"], positive=["pop", "area"], years=range(2011, 2025), code_col="sigungu_cd", label="ST_prereg_panel v2 분모", verbose=False, area_col="area")
out.to_csv(PF, index=False, encoding="utf-8-sig")
v2.write_fill_log(FILL, "explicit_fill_JKMS_prereg.csv")
open(os.path.join(LG, "530b_v2_prereg_fix.log"), "w").write("\n".join(LOG))
jr.write_log(os.path.join(LG, "JKMS_v2c_rawfix_530b.csv"))   # ★v2c
print(f"저장: {PF} ({out.shape})")
