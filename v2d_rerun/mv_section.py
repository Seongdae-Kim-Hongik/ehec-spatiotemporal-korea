# -*- coding: utf-8 -*-
# mv_section.py — (2026-10-22 v2c) 다중우주 끝난 뒤 00_결과요약.md의 다중우주 절을 채움(전 = 재실행 전 PREREG_MV_summary.csv·보충표 S15 정본, 후 = 이번).
import os, glob, re, pandas as pd
from docx import Document
ROOT = os.environ.get("EHEC_ROOT",".")
RS = os.path.join(ROOT, "results_ST"); RR = os.path.join(ROOT, "_v2_rerun_261022", "JKMS"); BK = os.path.join(ROOT, "_pre_v2_rerun_261022", "JKMS")
new = pd.read_csv(os.path.join(RS, "PREREG_MV_summary.csv")).set_index("term")
old = pd.read_csv(os.path.join(BK, "results_ST", "PREREG_MV_summary.csv")).set_index("term")
res = sorted(glob.glob(os.path.join(RS, "PREREG_MV_results_seed1[01][0-9].csv")))
nspec = sum(pd.read_csv(f).spec.nunique() for f in res)
LAB = {"beef_cattle_farms": "한육우 농가 밀도", "reservoir_pond_land_area": "저수지·연못 면적", "current_smoking": "현재흡연", "fiscal_independence": "재정자립도",
       "fiscal_autonomy": "재정자주도", "residential_building_site_area": "대지 면적(택지)", "net_migration": "순이동", "soap_sanitizer_use": "비누 손씻기",
       "swine_density": "돼지 밀도(음성대조)", "poultry_density": "가금 밀도(음성대조)", "fish_farm_land_area": "양어장 면적(음성대조)", "factory_land_area": "공장 용지(음성대조)", "walking_practice": "걷기 실천(음성대조)"}
keys = list(LAB) + [k for k in new.index if new.loc[k, "cred_pct"] >= 50 and k not in LAB and k not in ("sex_ratio", "kchs_pc1")] \
       + [k for k in old.index if old.loc[k, "cred_pct"] >= 50 and k not in LAB and k not in ("sex_ratio", "kchs_pc1") and k in new.index]
seen = set(); keys = [k for k in keys if not (k in seen or seen.add(k))]
def cell(df, k):
    if k not in df.index: return "-"
    r = df.loc[k]; return f"{r.cred_pct:.1f}% (n {int(r.n_spec)}, 중앙 {r.median_IRR:.2f}, {r.verdict if isinstance(r.verdict, str) else 'null(<20%)'})"
L = [f"- 수렴 사양 {nspec:,}개(씨드 101-112 × 250, 결과 파일 {len(res)}개). 판정 규칙은 사전 지정 그대로(≥ 50% credible + 방향 일관 = robust).", ""]
L += ["| 변수 | 전(원고) | 후 | 판정 변화 |", "|---|---|---|---|"]
for k in keys:
    ov = str(old.loc[k, "verdict"]) if k in old.index else "-"; nv = str(new.loc[k, "verdict"]) if k in new.index else "-"
    L.append(f"| {LAB.get(k, k)} | {cell(old, k)} | {cell(new, k)} | {'◆ ' + ov + ' → ' + nv if ov != nv else '같음'} |")
rob_old = set(old[(old.cred_pct >= 50)].index) - {"sex_ratio", "kchs_pc1"}; rob_new = set(new[(new.cred_pct >= 50)].index) - {"sex_ratio", "kchs_pc1"}
L += ["", f"- 50% 이상 credible 변수: 전 {len(rob_old)}개 → 후 {len(rob_new)}개. 새로 넘은 것 {', '.join(LAB.get(k, k) for k in sorted(rob_new - rob_old)) or '없음'}, 빠진 것 {', '.join(LAB.get(k, k) for k in sorted(rob_old - rob_new)) or '없음'}."]
# S15 축별
s15 = pd.read_csv(os.path.join(RR, "S15_axis_robustness_v2c.csv"))
d = Document(os.path.join(BK, "JKMS_FBD/JKMS_EHEC/Supplementary/Supplementary_Tables_EHEC.docx")); t = d.tables[14]
o15 = {r.cells[0].text.strip(): [c.text for c in r.cells[1:]] for r in t.rows[1:]}
L += ["", "표 S15 축별 credible 비율(%)(전 → 후)", "", "| 변수 | 전체 | PC1 축 없음/있음 | 확인 보정 없음/의사/검사/전부 | 구조 BYM2/+상호작용 | 그래프 Q/R/kNN5 | 변수 수 1-4/5-8/9+ |", "|---|---|---|---|---|---|---|"]
for r in s15.itertuples():
    o = o15.get(r.label, ["-"] * 5)
    L.append(f"| {r.label} | {o[0]} → {r.all} | {o[1]} → {r.axis_no_yes} | {o[2]} → {r.ascertainment_none_phys_lab_all} | {o[3]} → {r.structure_bym2_int} | {o[4]} → {r.graph_q_r_knn5} | {r.setsize_1_4_5_8_9p} |")
L += ["", "원고 줄 대조: 초록 41-42행(한육우 100%·재정자립도 70%), 결과 220-230행, 고찰 236-237·257-259·268행, S12 각주(수렴 사양 수), S13·S14 아래·S15, Fig S17. 사양 곡선 그림은 `EHEC/EHEC_Fig_specification_curve.png`(543)."]
F = os.path.join(RR, "00_결과요약.md"); s = open(F, encoding="utf-8").read()
s = re.sub(r"<!-- MV_SECTION_START -->.*?<!-- MV_SECTION_END -->", "<!-- MV_SECTION_START -->\n" + "\n".join(L) + "\n<!-- MV_SECTION_END -->", s, flags=re.S)
open(F, "w", encoding="utf-8").write(s); print("\n".join(L))
