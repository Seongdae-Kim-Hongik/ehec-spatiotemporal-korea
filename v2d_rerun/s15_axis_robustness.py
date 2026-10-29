# -*- coding: utf-8 -*-
# s15_axis_robustness.py — (2026-10-22 v2c) Table S15 축별 강건성. 551_appendix_jmir_layout.py의 같은 계산(씨드 101-112 결과·사양 결합, 사양당 첫 적합)을
#  그대로 옮겨 CSV로만 씀(551은 출력이 jmir_EHEC 부록 docx로 고정돼 있어 돌리지 않음). 추가로 사양 변수 수(n_vars)별 credible 비율.
import os, re, glob, pandas as pd
ROOT = os.environ.get("EHEC_ROOT",".")
R = os.path.join(ROOT, "results_ST"); OUT = os.path.dirname(os.path.abspath(__file__))
def loadmv(pat, idc):
    out = []
    for f in sorted(glob.glob(os.path.join(R, pat))):
        seed = int(re.search(r"seed(\d+)", f).group(1)); d = pd.read_csv(f); d["uid"] = seed * 10000 + d[idc]; out.append(d)
    return pd.concat(out, ignore_index=True)
RR = loadmv("PREREG_MV_results_seed1[01][0-9].csv", "spec"); RR = RR[RR.DIC == RR.groupby("uid").DIC.transform("first")].drop_duplicates(["uid", "term"])
SS = loadmv("PREREG_MV_specs_seed1[01][0-9].csv", "id"); RR = RR.merge(SS[["uid", "pc1", "det", "rum", "st", "graph", "n_vars"]], on="uid")
for c in ("pc1", "cred"): RR[c] = RR[c].astype(str).str.upper().eq("TRUE")
rows = []
VARS = [("beef_cattle_farms", "Beef-cattle farm density"), ("fiscal_independence", "Fiscal independence"), ("fiscal_autonomy", "Fiscal autonomy"), ("net_migration", "Net migration"),
        ("current_smoking", "Current smoking"), ("reservoir_pond_land_area", "Reservoir and pond area"), ("soap_sanitizer_use", "Handwashing with soap")]
for v, lab in VARS:
    r = RR[RR.term == v]
    if not len(r): continue
    f = lambda m: f"{r[m].cred.mean()*100:.0f}" if m.sum() >= 10 else "-"
    nv = r.n_vars.clip(upper=12); bins = pd.cut(nv, [0, 4, 8, 99], labels=["1-4", "5-8", "9+"])
    rows.append(dict(term=v, label=lab, n=len(r), all=f"{r.cred.mean()*100:.0f}", axis_no_yes=f"{f(~r.pc1)} / {f(r.pc1)}",
                     ascertainment_none_phys_lab_all=f"{f(r.det=='none')} / {f(r.det=='docs')} / {f(r.det=='lab')} / {f(r.det=='full')}",
                     structure_bym2_int=f"{f(r.st=='bym2')} / {f(r.st=='bym2_int')}", graph_q_r_knn5=f"{f(r.graph=='queen')} / {f(r.graph=='rook')} / {f(r.graph=='knn5')}",
                     setsize_1_4_5_8_9p=" / ".join(f(bins == b) for b in ["1-4", "5-8", "9+"]),
                     min_level=min((float(x) for x in [f(~r.pc1), f(r.pc1), f(r.det=='none'), f(r.det=='docs'), f(r.det=='lab'), f(r.det=='full'), f(r.st=='bym2'), f(r.st=='bym2_int'), f(r.graph=='queen'), f(r.graph=='rook'), f(r.graph=='knn5')] if x != "-"), default=float("nan"))))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "S15_axis_robustness_v2c.csv"), index=False, encoding="utf-8-sig")
print(pd.DataFrame(rows).to_string())
print("사양 수", RR.uid.nunique())
