#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 400 — EHEC 연령표준화 자료 조립 (2026-10-21 v2 정본화: 400v2_agestd_assemble.py의 시군구 코드 결합판을 정본 이름으로 씀).
#  이전 판(이름 결합, 미추홀 2018-24·군위·여주군·당진군·연기군·청원군 탈락, 연령 칸 0 채움)은 400_agestd_assemble.py.bak_pre_v2.
#  산출(FBD_CACHE, 기본 data_cache/):
#   ST_agepop.csv(+ 같은 내용 ST_agepop_mfds.csv — 식약처형 옛 경계 보정본 이름을 쓰는 스크립트용), ST_EHEC_agecase.csv, ST_natrate.csv, ST_Estd_panel.csv
#   (region = sigungu_key canonical_region, sigungu_cd 열 포함) · 로그 _v2_logs/agepop_v2_refmonth.csv·agepop_v2_dropped_names.csv
#  쓰지 않는 것: ST_stdpop2015.csv(2015 표준인구는 고정 표준으로 유지), ST_nat_agestd_trend.csv(400c가 포털 완전집계로 만듦).
#  아래 바뀐 점 1)-7)은 400v2 설명 그대로.
import pandas as pd, numpy as np, os, sys, warnings
GD = os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
GATE = os.path.join(GD, "논문_마스터자료/00_분석전_데이터게이트"); sys.path.insert(0, GATE)
import sigungu_key as sk
from panel_gate import panel_gate
ROOT = os.path.join(GD, "FBD_SpatioTemporal"); OUT = os.environ.get("FBD_CACHE", os.path.join(ROOT, "data_cache")); LOG = os.path.join(OUT, "_v2_logs")
os.makedirs(LOG, exist_ok=True)
KDCA = os.path.join(GD, "건상/성별_종합 (3).csv")
POP = os.path.join(GD, "인구데이터_통합.csv")
DZ = "장출혈성대장균감염증"
T = sk.load_tables(); KEY = T["key"].copy(); KEY["sigungu_cd"] = KEY.sigungu_cd.astype(str).str.zfill(5)
NAME = dict(zip(KEY.sigungu_cd, KEY.canonical_region)); CODES = sorted(NAME)
assert len(CODES) == 229
AGE10 = ["0~9세", "10~19세", "20~29세", "30~39세", "40~49세", "50~59세", "60~69세", "70세 이상"]


def pop5to10(a):
    m = {"0~4세": "0~9세", "5~9세": "0~9세", "10~14세": "10~19세", "15~19세": "10~19세",
         "20~24세": "20~29세", "25~29세": "20~29세", "30~34세": "30~39세", "35~39세": "30~39세",
         "40~44세": "40~49세", "45~49세": "40~49세", "50~54세": "50~59세", "55~59세": "50~59세",
         "60~64세": "60~69세", "65~69세": "60~69세"}
    if a in m: return m[a]
    if a in ("70~74세", "75~79세", "80~84세", "85~89세", "90~94세", "95~99세", "100세 이상"): return "70세 이상"
    return None


# ================= S1b: 인구(주민등록) → (sigungu_cd, year, age10) =================
def build_agepop():
    parts = []
    for ch in pd.read_csv(POP, usecols=["행정구역", "날짜", "구분", "연령대", "인구수"], dtype=str, chunksize=3_000_000):
        ch = ch[(ch["구분"] == "계") & ch["연령대"].map(lambda a: pop5to10(a) is not None)]
        parts.append(ch)
    pp = pd.concat(parts).drop_duplicates(subset=["행정구역", "날짜", "구분", "연령대"])
    n0 = pp["행정구역"].str.contains(r"\(3611", na=False).sum()
    pp = pp[~pp["행정구역"].str.contains(r"\(3611", na=False)].copy()          # 400과 같게 세종 중복코드(3611) 제외
    print(f"[인구] 세종 중복코드(3611) {n0}행 제외(400과 동일)")
    pp["year"] = pp["날짜"].str[:4].astype(int); pp["mo"] = pp["날짜"].str[5:7].astype(int)
    pp["인구수"] = pd.to_numeric(pp["인구수"].str.replace(",", ""), errors="coerce")
    names = pp["행정구역"].drop_duplicates().tolist()
    kinds = pd.Series(sk.name_kind(names), index=names)
    print("[인구] 행정구역 표기", len(names), "개, 종류:", kinds.value_counts().to_dict())
    drop = kinds[~kinds.isin(["시군구", "통합전"])]
    dropped = pp[pp["행정구역"].isin(drop.index)].groupby("행정구역").size()
    pd.DataFrame({"행정구역": dropped.index, "종류": drop.reindex(dropped.index).values, "행수": dropped.values}).to_csv(
        os.path.join(LOG, "agepop_v2_dropped_names.csv"), index=False, encoding="utf-8-sig")
    print(f"[인구] 시도·일반구·검토필요(출장소) {len(drop)}개 표기 {int(dropped.sum())}행 제외: 검토필요={drop[drop=='검토필요'].index.tolist()}")
    keep = pp[pp["행정구역"].map(kinds).isin(["시군구", "통합전"])].copy()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        keep = sk.attach_code(keep, "행정구역", "year")
    for x in w: print("[attach_code 경고]", str(x.message)[:600])
    # 기준 월: (코드, 연도)마다 하나 — 12월, 없으면 그 코드의 마지막 달
    ref = keep.groupby(["sigungu_cd", "year"]).mo.transform(lambda s: 12 if (s == 12).any() else s.max())
    good = keep[keep.mo == ref].copy()
    # 로그: (코드, 연도) 기준 월과 합산한 원자료 이름(통합·개명 해 확인용)
    lg = good.groupby(["sigungu_cd", "year"]).agg(ref_month=("mo", "first"), names=("행정구역", lambda s: " + ".join(sorted(set(s))))).reset_index()
    lg["n_names"] = lg.names.str.count(r"\+") + 1
    lg.insert(1, "region", lg.sigungu_cd.map(NAME))
    lg.to_csv(os.path.join(LOG, "agepop_v2_refmonth.csv"), index=False, encoding="utf-8-sig")
    print("[인구] 기준 월이 12월이 아닌 (코드, 연도):", lg[lg.ref_month != 12][["region", "year", "ref_month"]].to_dict("records")[:10])
    print("[인구] 둘 이상 이름을 합산한 (코드, 연도):\n", lg[lg.n_names > 1][["region", "year", "names"]].to_string(index=False))
    good["age10"] = good["연령대"].map(pop5to10)
    pop = good.groupby(["sigungu_cd", "year", "age10"], as_index=False)["인구수"].sum().rename(columns={"인구수": "pop"})
    pop = pop[pop.year.between(2008, 2024)]
    # 8개 연령 칸이 모두 있는지(빈 칸 0 채움 금지)
    cnt = pop.groupby(["sigungu_cd", "year"]).age10.nunique()
    assert (cnt == 8).all(), f"연령 칸 누락: {cnt[cnt != 8].head()}"
    assert pop["pop"].notna().all() and (pop["pop"] >= 0).all()
    pop["region"] = pop.sigungu_cd.map(NAME)
    return pop[["region", "year", "age10", "pop", "sigungu_cd"]].sort_values(["sigungu_cd", "year", "age10"]).reset_index(drop=True)


# ================= S1a: KDCA 연령별 사례 → (sigungu_cd, year, age10) =================
def build_kdca_case(dz, kd=None):
    if kd is None:
        kd = pd.read_csv(KDCA, encoding="utf-8-sig")
    kd = kd.drop_duplicates(["year", "sido", "sigungu", "disease", "sex", "age_group"])
    kd = kd[(kd.disease == dz) & (kd.sex == "계")].copy()
    kd["nm"] = kd.sido.astype(str).str.strip() + " " + kd.sigungu.astype(str).str.strip()
    u = kd[["nm"]].drop_duplicates(); u["kind"] = sk.name_kind(u.nm.tolist())
    u["sigungu_cd"] = sk.to_code(u.nm.tolist(), subdistrict="city")     # KDCA 교차표는 현행 경계(일반구 단위) 표기
    # 시 합계 행과 일반구 행이 같은 코드에 함께 있으면 이중 합산 위험 → 부천시(31050)만 예외: '부천시' 행은 2024 일반구 재설치 뒤 원미구 몫(부천시 < 소사+오정)이라 합산이 맞음
    mix = u.groupby("sigungu_cd").kind.agg(lambda s: set(s))
    mixed = [c for c, s in mix.items() if s == {"시군구", "일반구"}]
    assert set(mixed) <= {"31050"}, f"시 합계+일반구 혼재 코드: {mixed}"
    kd = kd.merge(u[["nm", "sigungu_cd"]], on="nm", how="left")
    mis = int(kd[kd.age_group == "연령미상"].cases.sum())
    kd = kd[kd.age_group.isin(AGE10)]
    case = kd.groupby(["sigungu_cd", "year", "age_group"], as_index=False)["cases"].sum()
    case = case[case.year.between(2001, 2024)]
    case["region"] = case.sigungu_cd.map(NAME)
    return case[["region", "year", "age_group", "cases", "sigungu_cd"]], mis


def estd_from(pop, case):
    c_sp = case[case.year.between(2011, 2024)].groupby("age_group")["cases"].sum().reindex(AGE10).fillna(0)  # 그 연령대 사례 0건 = 실제 0
    p_sp = pop[pop.year.between(2011, 2024)].groupby("age10")["pop"].sum().reindex(AGE10)
    natrate = (c_sp / p_sp); natrate.name = "rate"
    pw = pop.pivot_table(index=["sigungu_cd", "year"], columns="age10", values="pop", aggfunc="sum").reindex(columns=AGE10)
    assert pw.notna().all().all(), "연령 칸 결측(채우지 않음)"
    E = (pw * natrate.values).sum(axis=1).reset_index(); E.columns = ["sigungu_cd", "year", "E_std"]
    E["region"] = E.sigungu_cd.map(NAME)
    return E[["region", "year", "E_std", "sigungu_cd"]], natrate


if __name__ == "__main__":
    POPF = os.path.join(OUT, "ST_agepop.csv")
    pop = build_agepop(); pop.to_csv(POPF, index=False, encoding="utf-8-sig")
    pop.to_csv(os.path.join(OUT, "ST_agepop_mfds.csv"), index=False, encoding="utf-8-sig")   # 옛 식약처형 경계 보정본 이름(같은 내용)
    tot = pop.groupby(["sigungu_cd", "year"], as_index=False)["pop"].sum()
    print("\n[게이트] ST_agepop (연령 합) 2008-2024")
    panel_gate(tot, cols=["pop"], positive=["pop"], years=range(2008, 2025), code_col="sigungu_cd", label="ST_agepop v2 2008-2024")
    case, mis = build_kdca_case(DZ)
    case.to_csv(os.path.join(OUT, "ST_EHEC_agecase.csv"), index=False, encoding="utf-8-sig")
    E, natrate = estd_from(pop, case)
    natrate.to_csv(os.path.join(OUT, "ST_natrate.csv"))
    E.to_csv(os.path.join(OUT, "ST_Estd_panel.csv"), index=False, encoding="utf-8-sig")
    print("\n[게이트] ST_Estd_panel 2011-2024")
    panel_gate(E, cols=["E_std"], positive=["E_std"], years=range(2011, 2025), code_col="sigungu_cd", label="EHEC E v2")
    n = int(case[case.year.between(2011, 2024)].cases.sum()); es = E[E.year.between(2011, 2024)].E_std.sum()
    print(f"[EHEC] 연령미상 {mis} · 2011-24 사례 {n} · E_std 합 {es:.1f} (SIR≈{n/es:.4f}) · 사례 코드 {case.sigungu_cd.nunique()}/229 · 인구 코드 {pop.sigungu_cd.nunique()}/229")
    print(f"완료 → {OUT}/ST_{{agepop,agepop_mfds,EHEC_agecase,natrate,Estd_panel}}.csv")
