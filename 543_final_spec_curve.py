# -*- coding: utf-8 -*-
# 543 — 사전지정 멀티버스 v2 specification curve (상위 12변수)
import pandas as pd, numpy as np, glob, re, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
RS="results_ST"
def load(pat,idc):
    out=[]
    for f in sorted(glob.glob(os.path.join(RS,pat))):
        seed=int(re.search(r"seed(\d+)",f).group(1)); d=pd.read_csv(f); d["uid"]=seed*10000+d[idc]; out.append(d)
    return pd.concat(out,ignore_index=True)
R=load("PREREG_MV_results_seed1[01][0-9].csv","spec"); R=R[R.DIC==R.groupby("uid").DIC.transform("first")].drop_duplicates(["uid","term"])
S=pd.read_csv(os.path.join(RS,"PREREG_MV_summary.csv")); S=S[~S.term.isin(["sex_ratio","kchs_pc1"])].sort_values("cred_pct",ascending=False)
LAB={"beef_cattle_farms":"Beef-cattle farm density","current_smoking":"Current smoking","physicians_count":"Physicians (ascertainment)","fiscal_autonomy":"Fiscal autonomy","fiscal_independence":"Fiscal independence","net_migration":"Net migration","reservoir_pond_land_area":"Reservoir and pond area","older_adults_living_alone":"Older adults living alone","households":"Households","single_person_households":"Single-person households","unemployment_rate":"Unemployment rate","welfare_share_of_general_budget":"Welfare budget share","dairy_farms":"Dairy farm density","grdp_per_capita":"GRDP per capita","drinking_water_groundwater_use":"Groundwater use","soap_sanitizer_use":"Handwashing with soap","bars_pubs_per_1000_pop":"Bars and pubs"}
top=S.head(12).term.tolist()
plt.rcParams.update({"font.family":["Helvetica","Arial","DejaVu Sans"]})
fig,ax=plt.subplots(len(top),1,figsize=(13.5,1.35*len(top)),sharex=True); NMAX=max(len(R[R.term==t]) for t in top)
for i,t in enumerate(top):
    d=R[R.term==t].sort_values("IRR").reset_index(drop=True); a=ax[i]; x=np.arange(len(d))
    a.axhline(1,color="0.6",lw=1.2,ls="--"); a.vlines(x,d.lo,d.hi,color=np.where(d.cred,"#B24745","0.75"),lw=0.9,alpha=0.85); a.plot(x,d.IRR,".",ms=3.2,color="#374E55")
    a.set_yscale("log"); a.set_ylabel(LAB.get(t,t),fontsize=12,rotation=0,ha="right",va="center"); a.set_ylim(0.4,2.6); a.set_yticks([0.5,1,2]); a.set_yticklabels(["0.5","1","2"],fontsize=11)
    a.yaxis.set_minor_formatter(NullFormatter()); a.tick_params(axis="y",which="minor",length=0); a.set_xlim(-5,NMAX+5)
    row=S[S.term==t].iloc[0]; a.text(1.012,0.5,f"{row.cred_pct:.0f}%",transform=a.transAxes,ha="left",va="center",fontsize=12,fontweight="bold",color="#B24745")
    for sp in ("top","right"): a.spines[sp].set_visible(False)
ax[-1].set_xlabel("Specifications containing the covariate, ordered by incidence rate ratio",fontsize=13); ax[-1].tick_params(axis="x",labelsize=11)
ax[0].text(1.012,1.55,"Credible",transform=ax[0].transAxes,ha="left",va="center",fontsize=11,color="#B24745")
plt.tight_layout(); plt.subplots_adjust(right=0.93)
out="EHEC/EHEC_Fig_specification_curve.png"; plt.savefig(out,dpi=200,bbox_inches="tight"); print("저장:",out)
