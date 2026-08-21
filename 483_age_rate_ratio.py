#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 483 - Age-specific incidence rate ratios over time.
#  구성비(%)는 인구구조 변화에 오염 → 연령별 10만명당 발생률의 '비'로 보면 소아 상대위험 추이가 드러남.
#  0-9세 발생률 ÷ 10세이상 발생률, 그리고 ÷ 60세이상 발생률. 연도별 + 95%CI(포아송) + 추세검정.
#  Note: KDCA reports 10-year age bands (60-69, 70+), so 65+ cannot be derived; 60+ is used and stated in the paper.
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kendalltau
ROOT=os.environ.get("EHEC_ROOT",".")
C=os.path.join(ROOT,"data_cache"); EH=os.path.join(ROOT,"EHEC"); R=os.path.join(ROOT,"results_ST")

case=pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv"))
pop =pd.read_csv(os.path.join(C,"ST_agepop.csv")).rename(columns={"age10":"age_group"})
Y0,Y1=2008,2024   # 연령별 인구 커버리지
case=case[case.year.between(Y0,Y1)]
pop = pop[pop.year.between(Y0,Y1)]

YOUNG="0~9세"
def grp(g): return "young" if g==YOUNG else "old"
cs=case.copy(); cs["grp"]=cs.age_group.map(grp)
pp=pop.copy();  pp["grp"]=pp.age_group.map(grp)
# 60+ 별도
ELDER={"60~69세","70세 이상","70세이상"}
cs["is_eld"]=cs.age_group.isin(ELDER); pp["is_eld"]=pp.age_group.isin(ELDER)

C2=cs.groupby(["year","grp"])["cases"].sum().unstack(fill_value=0)
P2=pp.groupby(["year","grp"])["pop"].sum().unstack(fill_value=0)
Ce=cs[cs.is_eld].groupby("year")["cases"].sum(); Pe=pp[pp.is_eld].groupby("year")["pop"].sum()
Cy=cs[cs.age_group==YOUNG].groupby("year")["cases"].sum(); Py=pp[pp.age_group==YOUNG].groupby("year")["pop"].sum()

def rr_ci(c1,n1,c2,n2):
    """발생률 비 (그룹1/그룹2) + 95% CI. 0 보호 위해 Haldane 0.5 보정."""
    r1=(c1)/n1; r2=(c2)/n2
    rr=r1/r2
    se=np.sqrt(1/np.maximum(c1,0.5)+1/np.maximum(c2,0.5))
    lo=rr*np.exp(-1.96*se); hi=rr*np.exp(1.96*se)
    return rr,lo,hi

yrs=sorted(set(C2.index)&set(P2.index))
rows=[]
for y in yrs:
    cyoung=C2.loc[y,"young"]; nyoung=P2.loc[y,"young"]
    cold=C2.loc[y,"old"];     nold=P2.loc[y,"old"]
    rr1,lo1,hi1=rr_ci(cyoung,nyoung,cold,nold)                 # 0-9 vs 10+
    rr2,lo2,hi2=rr_ci(Cy.get(y,0),Py.get(y,1),Ce.get(y,0),Pe.get(y,1))  # 0-9 vs 60+
    r_y=cyoung/nyoung*1e5; r_o=cold/nold*1e5; r_e=Ce.get(y,0)/Pe.get(y,1)*1e5
    rows.append(dict(year=y,c_young=int(cyoung),c_old=int(cold),
                     rate_young=r_y,rate_old=r_o,rate_eld=r_e,
                     RR_young_vs_old=rr1,lo1=lo1,hi1=hi1,
                     RR_young_vs_eld=rr2,lo2=lo2,hi2=hi2))
df=pd.DataFrame(rows)
df.to_csv(os.path.join(R,"MASTER_ehec_age_rateratio.csv"),index=False)

# ── 추세검정 (RR + log-linear slope) ──
def trend(col):
    t,p=kendalltau(df.year,df[col])
    b=np.polyfit(df.year-df.year.min(),np.log(df[col]),1)[0]  # 연간 log 변화
    return t,p,(np.exp(b)-1)*100
print("=== 연령별 발생률 비(RR) 시계열 2008-2024 ===")
print(df[["year","c_young","c_old","rate_young","rate_old","rate_eld","RR_young_vs_old","RR_young_vs_eld"]].round(2).to_string(index=False))
for col,lab in [("RR_young_vs_old","0-9 vs 10+"),("RR_young_vs_eld","0-9 vs 60+")]:
    t,p,pct=trend(col)
    print(f"\n[추세] {lab}: Kendall tau={t:+.3f} p={p:.3f} · 연평균 {pct:+.1f}%/yr · {df[col].iloc[0]:.2f}→{df[col].iloc[-1]:.2f}")

# ── 그림 ──
plt.rcParams.update({"font.size":11,"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(13,5))
fig.suptitle("EHEC — age-specific incidence RATE RATIO over time (population-adjusted)",
             fontsize=13,fontweight="bold",x=0.02,ha="left")
fig.text(0.02,0.925,"Rate ratio removes the population-structure artefact in case-share plots. RR>1 = children have higher incidence per capita.",
         fontsize=8.5,color="0.35")
for a,(col,lo,hi,lab,c) in zip(ax,[("RR_young_vs_old","lo1","hi1","0–9y vs ≥10y","#c0392b"),
                                   ("RR_young_vs_eld","lo2","hi2","0–9y vs ≥60y","#e67e22")]):
    a.fill_between(df.year,df[lo],df[hi],color=c,alpha=0.15)
    a.plot(df.year,df[col],"o-",color=c,lw=2,ms=5)
    a.axhline(1,ls="--",color="0.5",lw=1)
    a.set_yscale("log"); a.set_xlabel("Year"); a.set_ylabel("Rate ratio (log scale)")
    t,p,pct=trend(col)
    a.set_title(f"{lab}\nKendall τ={t:+.2f}, p={p:.3f}  ·  {pct:+.1f}%/yr",fontsize=11,loc="left")
    for yv in [2020]: a.axvline(yv,ls=":",color="0.7",lw=0.8)
plt.tight_layout(rect=[0,0,1,0.9])
out=os.path.join(EH,"EHEC_Fig_age_rateratio.png")
plt.savefig(out,dpi=150,bbox_inches="tight"); print("\n저장:",out)

# ── [/333 보강] 정식 주검정: NB GLM 연령군×연도 상호작용 (offset=log pop) ──
#  Kendall tau(비가중·순위)는 보조로 강등. 이것이 "두 연령군 증가율이 다른가"의 직접 검정.
import statsmodels.api as sm, statsmodels.formula.api as smf
cs2=case.copy(); cs2["grp"]=np.where(cs2.age_group==YOUNG,"young","old")
pp2=pop.copy();  pp2["grp"]=np.where(pp2.age_group==YOUNG,"young","old")
Cg=cs2.groupby(["year","grp"])["cases"].sum().reset_index()
Pg=pp2.groupby(["year","grp"])["pop"].sum().reset_index()
mg=Cg.merge(Pg,on=["year","grp"]); mg["t"]=mg.year-mg.year.min(); mg["young"]=(mg.grp=="young").astype(int)
nbfit=smf.glm("cases ~ young*t",data=mg,family=sm.families.NegativeBinomial(alpha=1.0),offset=np.log(mg["pop"])).fit()
bI=nbfit.params["young:t"]; pI=nbfit.pvalues["young:t"]; ciI=nbfit.conf_int().loc["young:t"]
print(f"\n[정식검정 NB GLM] young:t = {bI:+.4f} (95%CI {ciI[0]:+.4f}~{ciI[1]:+.4f}) p={pI:.3f}"
      f" → 증가율 차이 {'유의' if pI<0.05 else '유의하지 않음'} (소아 기울기 {(np.exp(nbfit.params['t']+bI)-1)*100:+.1f} vs 10+ {(np.exp(nbfit.params['t'])-1)*100:+.1f}%/yr)")
