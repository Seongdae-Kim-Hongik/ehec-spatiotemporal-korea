# -*- coding: utf-8 -*-
# v2_fixes.py — (2026-10-21) FBD 시공간 캐시 빌더 공용 v2 교정 모듈. 400·450·450b·563·565·530b가 import.
#  출처: _v2_refit_261021/10_build_inputs.py(손씻기 원자료 코드 재구성·제주 재정자주도·세종 GRDP·kchs_pc1)와
#        11_build_inputs_jkms.py(JKMS 세종 J20·독거노인·제주 재정자립도·비누 손씻기·정화조 구조적 0)의 계산을 그대로 옮김.
#  규칙: 시군구 결합은 sigungu_key.attach_code(코드), 결측 채움은 명시 기록(explicit_fill 또는 FILL 로그)만.
#  출력 위치: 환경변수 FBD_CACHE(기본 data_cache). 검증 때는 FBD_CACHE=<임시 폴더>로 돌려 정본을 건드리지 않음.
import os, sys, json, warnings, unicodedata as ud
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
GD = os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
ROOT = os.path.join(GD, "FBD_SpatioTemporal")
GATE = os.path.join(GD, "논문_마스터자료", "00_분석전_데이터게이트")
if GATE not in sys.path: sys.path.insert(0, GATE)
import sigungu_key as sk                     # noqa: E402
from panel_gate import panel_gate, report_text  # noqa: E402,F401
RAW0 = os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
RAW = os.path.join(RAW0, "독립변수 3")
SRC = os.path.join(ROOT, "data_cache", "_v2_sources")          # 외부 출처(KOSIS) 사본 — 정본 위치
KOS_AUTONOMY = os.path.join(SRC, "DT_1YL20891.csv")             # 재정자주도(세입과목개편전) T10
KOS_INDEP = os.path.join(SRC, "DT_1YL20921.csv")                # 재정자립도(세입과목개편전) T10
KOS_GRDP = os.path.join(SRC, "KOSIS_DT_1C81_GRDP_sido_2015_2019.csv")
YY = list(range(2011, 2025))
nfc = lambda s: ud.normalize("NFC", str(s))
KEY = pd.read_csv(os.path.join(GATE, "sigungu_key.csv"), dtype={"sigungu_cd": str})
CODES = KEY.sigungu_cd.tolist(); assert len(CODES) == 229
NM = dict(zip(KEY.sigungu_cd, KEY.canonical_region))
AREA = dict(zip(KEY.sigungu_cd, KEY.area_km2.astype(float)))   # final.shp 기하 면적(EPSG:5179)
JEJU = ["39010", "39020"]; SEJONG = "29010"; CHEONGJU = "33010"


def cache_dir():
    d = os.environ.get("FBD_CACHE", os.path.join(ROOT, "data_cache")); os.makedirs(d, exist_ok=True); return d


def log_dir():
    d = os.path.join(cache_dir(), "_v2_logs"); os.makedirs(d, exist_ok=True); return d


def code(df, col="region", ycol="year"):
    """region 표기(공백 제거·NFC)를 sigungu_cd로. 못 찾으면 멈춤."""
    d = df.copy(); d[col] = d[col].astype(str).str.replace(r"\s+", "", regex=True).map(nfc)
    d = sk.attach_code(d, col, ycol if (ycol and ycol in d.columns) else None)
    d["sigungu_cd"] = d.sigungu_cd.astype(str); return d


def area_by_region(regs):
    """region 목록 → shp 면적(km²), 코드 결합."""
    a = code(pd.DataFrame({"region": list(regs)}), "region", None)
    s = a.sigungu_cd.map(AREA); assert s.notna().all(); return pd.Series(s.values, index=a.region.values)


# ── 손씻기 4종: 원자료 parquet 코드 결합(표준화율 우선·조율 보충, 0·>100 → 결측, 2011-2024 지역 평균) ──
HWG = {"handwashing_before_meals": ["식사전손씻기실천율_표준화율", "식사전손씻기실천율_조율"],
       "handwashing_after_toilet": ["화장실다녀온후손씻기실천율_표준화율", "화장실다년온후손씻기실천율_표준화율", "화장실다녀온후손씻기실천율_조율", "화장실다년온후손씻기실천율_조율"],
       "handwashing_after_going_out": ["외출후손씻기실천율_표준화율", "외출후손씻기실천율_조율"],
       "handwashing_with_soap": ["비누손세정제사용률_표준화율", "비누,손세정제사용률_표준화율", "비누손세정제사용률_조율", "비누,손세정제사용률_조율"]}
_PQ = None


def parquet_coded():
    global _PQ
    if _PQ is None:
        pq = code(pd.read_parquet(os.path.join(RAW, "국민건강결과.parquet")))
        assert not pq.duplicated(["sigungu_cd", "year"]).any()
        _PQ = pq[pq.year.between(2011, 2024)].copy()
    return _PQ


def hw_raw():
    pq = parquet_coded(); HW = {}
    for v, cols in HWG.items():
        x = pd.Series(np.nan, index=pq.index)
        for c in cols:
            if c in pq.columns:
                y = pd.to_numeric(pq[c], errors="coerce"); x = x.fillna(y.where((y > 0) & (y <= 100)))
        HW[v] = pd.Series(x.values, index=pq.index).groupby(pq.sigungu_cd).mean()
    HW = pd.DataFrame(HW).reindex(CODES); assert HW.notna().all().all(); return HW


def kosis_jeju(path):
    ko = pd.read_csv(path); ko = ko[ko.ITM_ID == "T10"]
    jj = ko[ko.C1_NM == "제주특별자치도"].set_index("PRD_DE").DT.astype(float).reindex(YY); assert jj.notna().all(); return jj


def apply_hm_fixes(d, fill_log=None):
    """HM 패널(region, year, pop, …)에 v2 교정(F1)을 얹음: 손씻기 4종 전 지역 원자료, 제주 재정자주도 KOSIS, 세종 1인당 GRDP KOSIS.
       d는 sigungu_cd 열이 있어야 함. 반환: 교정된 d. fill_log(list)에 채움 기록을 덧붙임."""
    d = d.copy(); HW = hw_raw()
    for v in HWG: d[v] = d.sigungu_cd.map(HW[v]).values
    jj = kosis_jeju(KOS_AUTONOMY)
    for cd in JEJU:
        m = d.sigungu_cd == cd; assert d.loc[m, "fiscal_autonomy"].isna().all(), f"{cd} 재정자주도가 이미 있음"
        d.loc[m, "fiscal_autonomy"] = d.loc[m, "year"].map(jj).values
    d, val = fill_grdp_sejong(d)
    if fill_log is not None:
        for v in HWG:
            fill_log.append(dict(model="HM", variable=v, sigungu_cd="(229 전 지역)", region="", year="2011-2024", value="", how="원자료 코드 재구성(전 지역): 손씻기_전처리.csv의 0 섞인 평균·빈칸 대신", source="독립변수 3/국민건강결과.parquet 표준화율 우선·조율 보충, 0·>100 결측, 실측 연도 평균"))
        for cd in JEJU:
            for y in YY: fill_log.append(dict(model="HM", variable="fiscal_autonomy", sigungu_cd=cd, region=NM[cd], year=y, value=jj[y], how="외부 출처 대체(explicit): 제주특별자치도 값을 두 행정시에 같게(행정시 값은 원자료에 없음)", source="KOSIS DT_1YL20891 T10 재정자주도(세입과목개편전)"))
        for y in YY: fill_log.append(dict(model="HM", variable="grdp_per_capita", sigungu_cd=SEJONG, region=NM[SEJONG], year=y, value=val, how="외부 출처 대체(explicit): KOSIS 시도 GRDP/패널 인구 2015-2019 평균(시간불변, 패널 규칙과 같음)", source="KOSIS DT_1C81 세종특별자치시 명목 GRDP"))
    return d


KOS_INCOMETAX = os.path.join(SRC, "DT_133N_A3212.csv")         # 국세청 종합소득세 주요항목 신고 현황Ⅱ(시·군·구) 2018~ — 원자료 종합소득세.csv의 원출처(값 일치 확인, 세종 행만 원자료에서 빠짐)


def kosis_sigungu(path):
    """KOSIS 시군구 표(C1 5자리 + 시도 2자리 행)를 '시도 시군구' 표기로 코드화(sigungu_cd 열 추가). 통합 전 단위(청원·마산 등)는 같은 코드로 붙으므로 쓰는 쪽에서 거름."""
    import contextlib, io
    k = pd.read_csv(path, dtype={"C1": str})
    sd = k[k.C1.str.len() == 2].drop_duplicates("C1").set_index("C1").C1_NM
    k5 = k[k.C1.str.len() == 5].copy(); k5["region"] = k5.C1.str[:2].map(sd) + " " + k5.C1_NM
    u = k5[["C1", "region"]].drop_duplicates()
    with contextlib.redirect_stdout(io.StringIO()): cc = sk.attach_code(u, "region", None)
    k5["sigungu_cd"] = k5.C1.map(dict(zip(cc.C1, cc.sigungu_cd.astype(str))))
    return k5


# ── v2b(2026-10-21 정지 지점 해결): 평균 채움 금지 뒤 멈춘 비교 세트의 재정 지표 ──
#  결정 순서 (1) 원자료 출처에서 sigungu_cd 결합으로 되살림 → (2) 구조적 결측은 주모형과 같은 외부 출처 대체(제주특별자치도 값) → (3) 그래도 남으면 세트·풀에서 빼고 기록(평균 채움 금지)
HM_FIXED_V2B = ["fiscal_independence", "fiscal_self_reliance"]


def apply_hm_fixes_v2b(d, fill_log=None):
    """HM 패널 v2b 교정(게이트 뒤, apply_hm_fixes 다음):
       fiscal_independence(국민건강결과 '재정자립도' = KOSIS DT_1YL20921 T10 소수 1자리) 제주시·서귀포시 → (2) KOSIS 제주특별자치도 T10(JKMS 530b와 같은 값),
       fiscal_self_reliance(재정자립도_전처리.csv '재정자립도(세입과목개편전)' = 같은 T10) 부산진구 → (1) 원출처 KOSIS T10 코드 결합(원자료 정밀도 소수 1자리; 원자료 행은 0=결측),
       제주시·서귀포시 → (2) KOSIS 제주특별자치도 T10. 이미 값이 있는 칸은 건드리지 않고(단언) 결측 칸만 채움."""
    d = d.copy(); jj = kosis_jeju(KOS_INDEP)
    kk = kosis_sigungu(KOS_INDEP); kk = kk[(kk.ITM_ID == "T10") & kk.sigungu_cd.notna()]
    BUSANJIN = sk.to_code(["부산광역시 부산진구"])[0]
    bj = kk[(kk.sigungu_cd == BUSANJIN) & (kk.C1_NM == "부산진구")].set_index("PRD_DE").DT.astype(float).reindex(YY); assert bj.notna().all()
    plan = [("fiscal_independence", cd, "2", jj, "외부 출처 대체(explicit): 제주특별자치도 값을 두 행정시에 같게(행정시 값은 원자료·KOSIS에 없음)", "KOSIS DT_1YL20921 T10 재정자립도(세입과목개편전) 제주특별자치도") for cd in JEJU] + \
           [("fiscal_self_reliance", BUSANJIN, "1", bj.round(1), "원자료 원출처 코드 결합(explicit): 재정자립도_전처리.csv의 부산진구 행이 전 연도 0(이름 결합 누락), 원출처 KOSIS 값을 원자료 정밀도(소수 1자리)로", "KOSIS DT_1YL20921 T10 부산진구(sigungu_cd 결합)")] + \
           [("fiscal_self_reliance", cd, "2", jj, "외부 출처 대체(explicit): 제주특별자치도 값을 두 행정시에 같게(행정시 값은 원자료·KOSIS에 없음)", "KOSIS DT_1YL20921 T10 재정자립도(세입과목개편전) 제주특별자치도") for cd in JEJU]
    for v, cd, rule, ser, how, src in plan:
        m = d.sigungu_cd == cd; assert m.sum() == len(YY), (v, cd)
        assert d.loc[m, v].isna().all(), f"{cd} {v}에 이미 값이 있음 — v2b 교정 전제 깨짐"
        d.loc[m, v] = d.loc[m, "year"].map(ser).values
        if fill_log is not None:
            for y in YY: fill_log.append(dict(model="HM", variable=v, sigungu_cd=cd, region=NM[cd], year=y, value=float(ser[y]), how=f"[v2b 규칙 ({rule})] " + how, source=src))
    return d


def fill_grdp_sejong(d):
    """세종 1인당 GRDP(원자료 parquet에 세종 GRDP 없음) = KOSIS DT_1C81 세종 명목 GRDP / 패널 인구, 2015-2019 평균(시간불변, 패널 Tier B 규칙과 같음).
       같은 규칙이 종로 패널 값을 재현하는지 단언. d는 sigungu_cd·year·pop·grdp_per_capita 열 필요."""
    d = d.copy()
    gk = pd.read_csv(KOS_GRDP); gs = gk[gk.C1_NM == "세종특별자치시"].set_index("PRD_DE").DT.astype(float); assert len(gs) == 5
    gr = parquet_coded().set_index(["sigungu_cd", "year"])["GRDP_당해년가격"]
    m = d.sigungu_cd == SEJONG; assert d.loc[m, "grdp_per_capita"].isna().all()
    pp = d[m].set_index("year")["pop"]; val = float(np.mean([gs[y] / pp[y] for y in range(2015, 2020)]))
    z = d[d.sigungu_cd == "11010"].set_index("year"); rep = float(np.mean([gr[("11010", y)] / z["pop"][y] for y in range(2015, 2020)]))
    assert abs(rep - z.grdp_per_capita.iloc[0]) < 1e-6 * max(1, abs(rep)), (rep, z.grdp_per_capita.iloc[0])   # 같은 규칙 재현(종로)
    d.loc[m, "grdp_per_capita"] = val
    return d, val


# ── kchs_pc1 ──
def pca1(R):
    Z = (R - R.mean()) / R.std(); u, s, vt = np.linalg.svd(Z.values, full_matrices=False)
    return pd.Series(u[:, 0] * s[0], index=R.index), float(s[0] ** 2 / np.sum(s ** 2)), pd.Series(vt[0], index=R.columns)


def kchs_handwash(R):
    """손씻기 종합지표: 지역 평균 손씻기 지표들의 PC1. 부호 = 손씻기 실천이 많을수록 +(적재 합 > 0)."""
    assert R.notna().all().all(), R.isna().sum()[lambda s: s > 0]
    pc, var, ld = pca1(R)
    if ld.sum() < 0: pc = -pc; ld = -ld
    return pc, var, ld


K5 = ["handwashing_before_meals", "handwashing_after_toilet", "handwashing_after_going_out", "handwashing_with_soap", "handwashing_after_outing"]   # 노로·살모넬라(565 v2 제조법)
K4 = ["handwashing_before_meals", "handwashing_after_going_out", "handwashing_with_soap", "handwashing_after_toilet"]                              # 이질·장티푸스(원자료 손씻기 4종)


def write_fill_log(rows, name):
    if not rows: return None
    p = os.path.join(log_dir(), name); pd.DataFrame(rows).to_csv(p, index=False, encoding="utf-8-sig"); return p
