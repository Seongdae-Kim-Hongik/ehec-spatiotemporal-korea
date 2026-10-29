#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 450 — 연령표준화 자료조립(일반화). 인자: 질병한글명 영문명   예) python3 450_agestd_assemble_generic.py 세균성이질 Shigellosis
#  (2026-10-21 v2 정본화) 450v2_agestd_assemble_generic.py의 시군구 코드 결합판을 정본 이름으로 씀. 이전 판은 450_agestd_assemble_generic.py.bak_pre_v2.
#   - 인구: FBD_CACHE/ST_agepop.csv(400 v2 산출, sigungu_cd 포함)만 씀. 없거나 sigungu_cd가 없으면 멈춤(옛 판처럼 원자료를 이름 결합으로 다시 만들지 않음).
#   - KDCA 사례: 400.build_kdca_case(drop_duplicates — 원자료 A형간염 113,400행·E형간염 45,360행이 통째로 두 번 들어 있음, to_code, 일반구→시 합산).
#   - 연령 칸 0 채움(pivot fill) 삭제, 8개 칸 단언. panel_gate로 E가 229 × 2011-2024 모두 있는지 확인.
#   - 전국 직접표준화 추세(ST_<eng>_trend.csv)는 옛 판과 같은 정의(표준인구 = 그해 2015 전국 연령 인구, 사례 = 시군구 교차표 합)로 v2 인구·중복 제거 사례에서 다시 계산.
#     이질·장티푸스 추세는 이어서 450c(포털 전국 완전집계)가 덮어씀.
#  산출(FBD_CACHE, 기본 data_cache/): ST_<eng>_agecase.csv · ST_<eng>_Estd.csv · ST_<eng>_natrate.csv · ST_<eng>_trend.csv
import os, sys, importlib.util
import numpy as np, pandas as pd
DZ = sys.argv[1]; ENG = sys.argv[2]
ROOT = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("a400", os.path.join(ROOT, "400_agestd_assemble.py")); a4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a4)
OUT = a4.OUT; AGE10 = a4.AGE10
POPF = os.path.join(OUT, "ST_agepop.csv")
if not os.path.exists(POPF): raise SystemExit(f"{POPF} 없음 — 400_agestd_assemble.py를 먼저 실행")
pop = pd.read_csv(POPF, dtype={"sigungu_cd": str})
if "sigungu_cd" not in pop.columns: raise SystemExit("ST_agepop.csv에 sigungu_cd 없음 — v2 이전 판. 400_agestd_assemble.py를 다시 실행")
case, mis = a4.build_kdca_case(DZ)
case.to_csv(os.path.join(OUT, f"ST_{ENG}_agecase.csv"), index=False, encoding="utf-8-sig")
E, natrate = a4.estd_from(pop, case)
natrate.to_csv(os.path.join(OUT, f"ST_{ENG}_natrate.csv"))
E.to_csv(os.path.join(OUT, f"ST_{ENG}_Estd.csv"), index=False, encoding="utf-8-sig")
print(f"\n[게이트] {ENG} E 2011-2024")
a4.panel_gate(E, cols=["E_std"], positive=["E_std"], years=range(2011, 2025), code_col="sigungu_cd", label=f"{ENG} E v2")
# 전국 직접표준화 추세(옛 판과 같은 정의)
FS = os.path.join(os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP"), "식중독최종.csv")
fs = pd.read_csv(FS, usecols=["region", "year", "population"])
nat_totpop = fs.groupby(["year", "region"], as_index=False)["population"].max().groupby("year")["population"].sum()
std2015 = pop[pop.year == 2015].groupby("age10")["pop"].sum().reindex(AGE10)
tot_case = case.groupby("year")["cases"].sum(); rows = []
for y in range(2001, 2025):
    O = int(tot_case.get(y, 0)); Ncr = nat_totpop.get(y, np.nan)
    crude = O / Ncr * 1e5 if Ncr == Ncr and Ncr > 0 else np.nan; astd = asl = asu = np.nan
    if 2008 <= y <= 2024:
        cy = case[case.year == y].groupby("age_group")["cases"].sum().reindex(AGE10, fill_value=0)   # 그 연령 사례 없음 = 0건(교차표에 행이 없는 연령)
        py = pop[pop.year == y].groupby("age10")["pop"].sum().reindex(AGE10); assert py.notna().all()
        ry = (cy / py).values; w = std2015.values; W = w.sum()
        astd = np.nansum(ry * w) / W * 1e5; var = np.nansum((w ** 2) * (cy.values) / (py.values ** 2)); se = np.sqrt(var) / W * 1e5
        asl = astd - 1.96 * se; asu = astd + 1.96 * se
    rows.append(dict(year=y, cases=O, crude=crude, agestd=astd, agestd_lo=asl, agestd_hi=asu))
pd.DataFrame(rows).to_csv(os.path.join(OUT, f"ST_{ENG}_trend.csv"), index=False)
n = int(case[case.year.between(2011, 2024)].cases.sum()); es = E[E.year.between(2011, 2024)].E_std.sum()
print(f"[{ENG}] 연령미상 {mis} · 2011-24 사례 {n} · E_std 합 {es:.1f} (SIR≈{n/max(es,1e-9):.4f}) · 사례 코드 {case.sigungu_cd.nunique()}/229")
