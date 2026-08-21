#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 507 — Specification curve 작도 + 사전지정 판정 규칙 적용 (506 산출 사용)
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=os.environ.get("EHEC_ROOT",".")
RS=os.path.join(ROOT,"results_ST"); EH=os.path.join(ROOT,"EHEC")
R=pd.read_csv(os.path.join(RS,"MASTER_multiverse_results.csv"))
S=pd.read_csv(os.path.join(RS,"MASTER_multiverse_summary.csv"))
specs=pd.read_csv(os.path.join(RS,"MASTER_multiverse_specs.csv"))
R=R[R.term!="(Intercept)"]

# ── 사전지정 판정 (PREREGISTERED §4) ──
def verdict(r):
    if r.cred_pct>=50 and (r.pct_above1>=90 or r.pct_above1<=10): return "강건(robust)"
    if r.cred_pct>=50: return "강건하나 방향불안정"
    if r.cred_pct>=20: return "시사적(suggestive)"
    return "null"
S["verdict"]=S.apply(verdict,axis=1)
S["dir"]=np.where(S.median_IRR>1,"위험(+)","보호(-)")
S=S.sort_values("cred_pct",ascending=False)
S.to_csv(os.path.join(RS,"MASTER_multiverse_summary.csv"),index=False)

LAB={"d_cattle":"Cattle density","d_dairy":"Dairy density","d_pasture":"Pasture area","rum_TMid":"Ruminant mid tertile",
"rum_THigh":"Ruminant high tertile","d_field":"Upland-field area","d_paddy":"Paddy area","e_forest":"Forest area",
"a_sewer":"Public sewerage","a_septic":"Septic-tank population","a_water_supply":"Piped water supply",
"a_advtreat":"Advanced wastewater treatment","b_gw_use":"Groundwater use","b_gw_unfit":"Unfit groundwater",
"b_gw_home":"Household wells","a_wastewater":"Wastewater discharge","a_orgload":"Organic pollutant load",
"c_hw_toilet":"Handwashing after toilet","c_hw_meal":"Handwashing before meals","c_hw_out":"Handwashing after going out",
"c_hw_soap":"Soap/sanitiser use","h_fastfood":"Fast-food outlets","h_pub":"Bars and pubs","g_fiscal":"Fiscal independence",
"g_basiclivelihood":"Basic-livelihood recipients","g_solo_elder":"Older adults living alone","g_sexratio":"Sex ratio",
"g_foreign":"Foreign residents","e_urban_pop":"Urban population","e_rural_pop":"Rural population","e_pop_dens":"Population density"}

print("═══ 사전지정 판정 결과 ═══")
for v in ["강건(robust)","강건하나 방향불안정","시사적(suggestive)","null"]:
    sub=S[S.verdict==v]
    if not len(sub): continue
    print(f"\n[{v}] {len(sub)}개")
    for _,r in sub.iterrows():
        print(f"   {LAB.get(r.term,r.term):32s} credible {r.cred_pct:5.1f}% ({r.n_cred}/{r.n_spec}) · median IRR {r.median_IRR:.3f} · {r.dir} · IRR>1 {r.pct_above1:.0f}%")

# ── 그림: 상위 12변수 specification curve ──
top=S.head(12).term.tolist()
from matplotlib.ticker import NullFormatter
plt.rcParams.update({"font.family":["Helvetica","Arial","DejaVu Sans"]})
fig,ax=plt.subplots(len(top),1,figsize=(13.5,1.35*len(top)),sharex=True)
NMAX=max(len(R[R.term==t]) for t in top)
for i,t in enumerate(top):
    d=R[R.term==t].sort_values("IRR").reset_index(drop=True)
    a=ax[i]; x=np.arange(len(d))
    a.axhline(1,color="0.6",lw=1.2,ls="--")
    a.vlines(x,d.lo,d.hi,color=np.where(d.cred,"#B24745","0.75"),lw=1.1,alpha=0.85)
    a.plot(x,d.IRR,".",ms=4.0,color="#374E55")
    a.set_yscale("log")
    a.set_ylabel(LAB.get(t,t),fontsize=12,rotation=0,ha="right",va="center")
    a.set_ylim(0.4,2.6); a.set_yticks([0.5,1,2]); a.set_yticklabels(["0.5","1","2"],fontsize=11)
    a.yaxis.set_minor_formatter(NullFormatter())          # 로그 minor 눈금 라벨(6x10^-1 등) 제거
    a.tick_params(axis="y",which="minor",length=0)
    a.set_xlim(-2,NMAX+2)                                  # 패널 간 x 범위 통일
    row=S[S.term==t].iloc[0]
    # credible 라벨은 축 바깥에 두어 데이터와 겹치지 않게 한다
    a.text(1.012,0.5,f"{row.cred_pct:.0f}%",transform=a.transAxes,ha="left",va="center",
           fontsize=12,fontweight="bold",color="#B24745")
    for s in ("top","right"): a.spines[s].set_visible(False)
ax[-1].set_xlabel("Specifications, ordered by incidence rate ratio",fontsize=13)
ax[-1].tick_params(axis="x",labelsize=11)
ax[0].text(1.012,1.55,"Credible",transform=ax[0].transAxes,ha="left",va="center",
           fontsize=11,color="#B24745")
plt.tight_layout()
plt.subplots_adjust(right=0.93)
out=os.path.join(EH,"EHEC_Fig_specification_curve.png")
plt.savefig(out,dpi=200,bbox_inches="tight"); print("\n저장:",out)
