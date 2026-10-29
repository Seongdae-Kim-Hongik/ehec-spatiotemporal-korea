# -*- coding: utf-8 -*-
"""s3_children_pct_v2d.py — (2026-10-29 /333 반영 8번) 표 S3 재적합용 0-4세 인구 비율을 v2 규칙(시군구 코드 결합)으로 다시 만든다.
입력: G.Downloads/인구데이터_통합.csv(행안부 월별 5세 단위 인구, 2008-2024). 각 연도 12월(없으면 그 해 마지막 달) 값.
출력: _internal/s3_refit/ST_children_pct_v2d.csv (sigungu_cd, region, year, children_0_4) — 229 × 2011-2024 완비 단언. 공유 data_cache는 건드리지 않음.
"""
import os, re, sys, pandas as pd
GD = os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
GATE = os.path.join(GD, "논문_마스터자료", "00_분석전_데이터게이트"); sys.path.insert(0, GATE)
from sigungu_key import to_code, name_kind, load_tables, KEY_CSV
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(H, "_internal", "s3_refit"); os.makedirs(OUT, exist_ok=True)
p = pd.read_csv(os.path.join(GD, "인구데이터_통합.csv"), dtype=str, usecols=["행정구역", "날짜", "구분", "연령대", "인구수"])
p = p[(p["구분"] == "계") & (p["연령대"].isin(["총인구수", "0~4세"]))].drop_duplicates(subset=["행정구역", "날짜", "연령대"])
p = p[~p["행정구역"].str.contains(r"\(3611", na=False)]   # 세종 중복 코드(401과 같은 처리)
p["name"] = p["행정구역"].str.replace(r"\s*\(\d+\)\s*$", "", regex=True).str.strip()
p["yr"] = p["날짜"].str[:4].astype(int); p["mo"] = p["날짜"].str[5:7].astype(int)
p = p[p.yr.between(2011, 2024)]
names = sorted(p["name"].unique()); kinds = dict(zip(names, name_kind(names)))
codes = {}
for n in names:
    k = kinds[n]
    if k in ("시도", "전국", "일반구", "모호", "미확인", "검토필요"): continue
    try: codes[n] = to_code(n, subdistrict="error")
    except Exception as e: print("skip", n, k, str(e)[:80])
print("이름 종류:", pd.Series(list(kinds.values())).value_counts().to_dict(), "| 코드화", len(codes))
p = p[p["name"].isin(codes)].copy(); p["sigungu_cd"] = p["name"].map(codes)
p["pop"] = pd.to_numeric(p["인구수"].str.replace(",", ""), errors="coerce")
# 연도별 12월 우선
p["pref"] = (p.mo == 12).astype(int)
p = p.sort_values(["sigungu_cd", "yr", "연령대", "pref", "mo"]).groupby(["sigungu_cd", "yr", "연령대"], as_index=False).last()
w = p.pivot_table(index=["sigungu_cd", "yr"], columns="연령대", values="pop", aggfunc="first").reset_index()
w = w.rename(columns={"yr": "year", "총인구수": "tot", "0~4세": "p04"}); w = w[w.tot > 0]
w["children_0_4"] = w.p04 / w.tot * 100
key = pd.read_csv(KEY_CSV, dtype=str); full = pd.MultiIndex.from_product([key.sigungu_cd, range(2011, 2025)], names=["sigungu_cd", "year"]).to_frame(index=False)
m = full.merge(w[["sigungu_cd", "year", "children_0_4", "tot", "p04"]], how="left"); miss = m[m.children_0_4.isna()]
print("결측 칸", len(miss)); print(miss.merge(key[["sigungu_cd", "canonical_region"]]).groupby("canonical_region").year.apply(list).to_dict())
fills = []
if len(miss):   # 구조적 결측(세종 2011 출범 전 등)만 가장 가까운 관측 연도로 채우고 기록(explicit fill)
    for cd, g in miss.groupby("sigungu_cd"):
        obs = m[(m.sigungu_cd == cd) & m.children_0_4.notna()]
        assert len(obs) >= 10, (cd, len(obs))
        for y in g.year:
            near = obs.iloc[(obs.year - y).abs().argsort().iloc[0]]
            m.loc[(m.sigungu_cd == cd) & (m.year == y), "children_0_4"] = near.children_0_4; fills.append((cd, y, int(near.year), round(near.children_0_4, 4)))
    pd.DataFrame(fills, columns=["sigungu_cd", "year", "from_year", "value"]).to_csv(os.path.join(OUT, "children_explicit_fill_log.csv"), index=False)
assert m.children_0_4.notna().all() and len(m) == 229 * 14
m = m.merge(key[["sigungu_cd", "canonical_region"]]).rename(columns={"canonical_region": "region"})
m[["sigungu_cd", "region", "year", "children_0_4"]].to_csv(os.path.join(OUT, "ST_children_pct_v2d.csv"), index=False, encoding="utf-8-sig")
print("saved", len(m), "| 채움", len(fills), "| 전국 평균 2011 %.2f, 2024 %.2f" % (m[m.year == 2011].children_0_4.mean(), m[m.year == 2024].children_0_4.mean()))
