#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""401b_children_pct_v2.py — (2026-10-29 저자 결정 4) 0-4세·0-9세 인구 비율 패널을 v2 규칙(시군구 코드 결합·채움 금지)으로 만든다.
401_children_pct.py(이름 결합, v1 인구, 연도별 225-228곳)를 대체하는 정본 빌더.
입력: G.Downloads/인구데이터_통합.csv(행안부 월별 5세 단위 인구, 2008-2024). (코드, 연도)마다 12월(없으면 그 해 마지막 달) 값.
출력: data_cache/ST_children_pct.csv — sigungu_cd, region(canonical), year 2008-2024, children_0_4, children_0_9 (229 × 17 완비 단언)
     구조적 결측(세종 2012 출범 전 등)만 가장 가까운 관측 연도로 채우고 data_cache/_v2_logs/explicit_fill_children_pct.csv에 기록(explicit fill).
검증: 2011-2024 children_0_4를 JKMS 표 S3 1차 재적합 입력(JKMS_EHEC_v2d_261029/_internal/s3_refit/ST_children_pct_v2d.csv, 통합 전 단위 미합산 판)과 대조 — 청주 2011-2013(청원군 합산, 모형에서는 마스킹 칸)만 달라야 함. 표 S3 최종값은 이 정본으로 s3_child04_refit_v2d.R을 다시 돌려 확인.
사용: python3 401b_children_pct_v2.py   (기존 정본은 먼저 ST_children_pct.csv.bak_pre_v2d_261029로 백업; 이미 있으면 덮어쓰지 않음)
"""
import os, sys, shutil, pandas as pd, numpy as np
GD = os.path.join(os.environ.get("EHEC_ROOT", "."), "G.Downloads")
ROOT = os.path.join(GD, "FBD_SpatioTemporal"); CACHE = os.path.join(ROOT, "data_cache"); LOGS = os.path.join(CACHE, "_v2_logs"); os.makedirs(LOGS, exist_ok=True)
GATE = os.environ.get("EHEC_GATE", os.path.join(os.environ.get("EHEC_ROOT", "."), "gate")); sys.path.insert(0, GATE)
from sigungu_key import to_code, name_kind, KEY_CSV
from panel_gate import panel_gate
Y0, Y1 = 2008, 2024
p = pd.read_csv(os.path.join(GD, "인구데이터_통합.csv"), dtype=str, usecols=["행정구역", "날짜", "구분", "연령대", "인구수"])
p = p[(p["구분"] == "계") & (p["연령대"].isin(["총인구수", "0~4세", "5~9세"]))].drop_duplicates(subset=["행정구역", "날짜", "연령대"])
p = p[~p["행정구역"].str.contains(r"\(3611", na=False)]   # 세종 중복 코드
p["name"] = p["행정구역"].str.replace(r"\s*\(\d+\)\s*$", "", regex=True).str.strip()
p["yr"] = p["날짜"].str[:4].astype(int); p["mo"] = p["날짜"].str[5:7].astype(int)
p = p[p.yr.between(Y0, Y1)]
names = sorted(p["name"].unique()); kinds = dict(zip(names, name_kind(names)))
codes = {}
for n in names:
    if kinds[n] in ("시도", "전국", "일반구", "모호", "미확인", "검토필요"): continue
    try: codes[n] = to_code(n, subdistrict="error")
    except Exception as e: print("skip", n, kinds[n], str(e)[:80])
print("이름 종류:", pd.Series(list(kinds.values())).value_counts().to_dict(), "| 코드화", len(codes))
p = p[p["name"].isin(codes)].copy(); p["sigungu_cd"] = p["name"].map(codes)
p["pop"] = pd.to_numeric(p["인구수"].str.replace(",", ""), errors="coerce")
# 기준 월: (코드, 연도)마다 하나 — 12월, 없으면 그 코드의 마지막 달(400_agestd_assemble.py와 같은 규칙). 통합 전 단위(청주+청원, 창원+마산+진해, 연기)는 같은 코드로 합산.
ref = p.groupby(["sigungu_cd", "yr"]).mo.transform(lambda s: 12 if (s == 12).any() else s.max()); p = p[p.mo == ref]
lg = p.groupby(["sigungu_cd", "yr"])["name"].agg(lambda s: " + ".join(sorted(set(s)))).reset_index(); print("[합산 (코드, 연도)]", lg[lg["name"].str.contains(r"\+")].to_dict("records"))
w = p.pivot_table(index=["sigungu_cd", "yr"], columns="연령대", values="pop", aggfunc="sum").reset_index()
w = w.rename(columns={"yr": "year", "총인구수": "tot", "0~4세": "p04", "5~9세": "p59"}); w = w[w.tot > 0]
w["children_0_4"] = w.p04 / w.tot * 100; w["children_0_9"] = (w.p04 + w.p59) / w.tot * 100
key = pd.read_csv(KEY_CSV, dtype=str)
full = pd.MultiIndex.from_product([key.sigungu_cd, range(Y0, Y1 + 1)], names=["sigungu_cd", "year"]).to_frame(index=False)
m = full.merge(w[["sigungu_cd", "year", "children_0_4", "children_0_9"]], how="left"); miss = m[m.children_0_4.isna()]
print("결측 칸", len(miss), miss.merge(key[["sigungu_cd", "canonical_region"]]).groupby("canonical_region").year.apply(list).to_dict())
fills = []
for cd, g in miss.groupby("sigungu_cd"):
    obs = m[(m.sigungu_cd == cd) & m.children_0_4.notna()]; assert len(obs) >= 10, (cd, len(obs))
    for y in g.year:
        near = obs.iloc[(obs.year - y).abs().argsort().iloc[0]]
        for c in ("children_0_4", "children_0_9"): m.loc[(m.sigungu_cd == cd) & (m.year == y), c] = near[c]
        fills.append((cd, y, int(near.year), round(near.children_0_4, 4), round(near.children_0_9, 4)))
pd.DataFrame(fills, columns=["sigungu_cd", "year", "from_year", "children_0_4", "children_0_9"]).to_csv(os.path.join(LOGS, "explicit_fill_children_pct.csv"), index=False)
assert m[["children_0_4", "children_0_9"]].notna().all().all() and len(m) == 229 * (Y1 - Y0 + 1)
m = m.merge(key[["sigungu_cd", "canonical_region"]]).rename(columns={"canonical_region": "region"})
out = m[["sigungu_cd", "region", "year", "children_0_4", "children_0_9"]].sort_values(["sigungu_cd", "year"])
# 표 S3 재적합 입력(2011-2024, 0-4세)과 동일한지
ref = os.path.join(ROOT, "JKMS_FBD", "JKMS_EHEC_v2d_261029", "_internal", "s3_refit", "ST_children_pct_v2d.csv")
if os.path.exists(ref):
    r = pd.read_csv(ref, dtype={"sigungu_cd": str}, encoding="utf-8-sig"); j = r.merge(out, on=["sigungu_cd", "year"], suffixes=("_ref", ""))
    j["d"] = (j.children_0_4_ref - j.children_0_4).abs(); dif = j[j.d > 1e-9]; assert len(j) == 229 * 14
    print("표 S3 재적합 입력(s3_children_pct_v2d, 통합 전 단위 미합산)과 다른 칸:", len(dif), dif[["sigungu_cd", "year", "children_0_4_ref", "children_0_4"]].round(4).to_dict("records"))
    assert set(dif.sigungu_cd) <= {"33010"}, "청주(2011-2013, 모형에서 마스킹) 외 차이"
# 게이트
rep = panel_gate(out, cols=["children_0_4", "children_0_9"], positive=["children_0_4", "children_0_9"], years=range(Y0, Y1 + 1), code_col="sigungu_cd", label="ST_children_pct_v2", verbose=True, raise_on_fail=True)
dst = os.path.join(CACHE, "ST_children_pct.csv"); bak = dst + ".bak_pre_v2d_261029"
if os.path.exists(dst) and not os.path.exists(bak): shutil.copy2(dst, bak); print("백업", bak)
out.to_csv(dst, index=False, encoding="utf-8-sig")
print("saved", dst, len(out), "행 | 채움", len(fills), "| 전국 평균 0-4세 2008 %.2f 2011 %.2f 2024 %.2f" % tuple(out[out.year == y].children_0_4.mean() for y in (2008, 2011, 2024)))
