#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 462 — EHEC 시간축 트랜지션 (MAIN, /333 만장일치 개정).
#  핵심: share(구성비)는 인구고령화 교란 → age-specific RATE(발생률)로 전환(소아 실제감소 vs 희석 분리).
#  ① 연령별 발생률 4그룹(0-9·10-19·20-59·60+), 2008-2024(연령인구 가용), 각 MK. + 평균연령 추세(MK).
#  ② 계절성: 월별 Gini(집중도, 유지) + 원형통계 평균각(피크 시점 phase) + 벡터길이 r. (raw peak month·조화회귀 아님-저N)
#  DON'T(저N): joinpoint·성×연령·onset/offset. 8밴드 share는 부록(460).
import os, numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import kendalltau
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); C=os.path.join(ROOT,"data_cache"); EH=os.path.join(ROOT,"EHEC")
mpl.rcParams.update({"font.family":"Helvetica Neue","font.size":10.5,"axes.edgecolor":"#444","axes.linewidth":0.8})
COVID="#d00000"
def mk(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);m=np.isfinite(x)&np.isfinite(y)
    if m.sum()<4: return np.nan,np.nan
    return kendalltau(x[m],y[m])
def gini(v):
    v=np.sort(np.asarray(v,float)); n=len(v); s=v.sum()
    return np.nan if s==0 else (2*np.sum(np.arange(1,n+1)*v)/(n*s))-(n+1)/n

# ===== 데이터 =====
k=pd.read_csv(os.path.join(GD,"건상/성별_종합 (3).csv"),encoding="utf-8-sig")
k=k[(k.disease=="장출혈성대장균감염증")&(k.sex=="계")&(~k.age_group.isin(["계","연령미상"]))]
ap=pd.read_csv(os.path.join(C,"ST_agepop.csv"))
G4={"0~9세":"Children (0–9)","10~19세":"School-age (10–19)","20~29세":"Adult (20–59)","30~39세":"Adult (20–59)",
    "40~49세":"Adult (20–59)","50~59세":"Adult (20–59)","60~69세":"Elderly (60+)","70세 이상":"Elderly (60+)"}
ORD=[("Children (0–9)","#c1121f"),("School-age (10–19)","#f4a261"),("Adult (20–59)","#2a9d8f"),("Elderly (60+)","#1a1a2e")]
# ① 연령별 발생률 (2008-2024)
kc=k[k.year.between(2008,2024)].copy(); kc["g"]=kc.age_group.map(G4)
cg=kc.groupby(["year","g"])["cases"].sum()
ap["g"]=ap.age10.map(G4); pg=ap.groupby(["year","g"])["pop"].sum()
rate=(cg/pg*1e5).unstack()  # year×group
# 평균연령 (2001-2024, share 기반이라 전기간 가능)
mid={"0~9세":5,"10~19세":15,"20~29세":25,"30~39세":35,"40~49세":45,"50~59세":55,"60~69세":65,"70세 이상":75}
k["mid"]=k.age_group.map(mid)
ma=k.groupby("year").apply(lambda d:(d.cases*d.mid).sum()/max(d.cases.sum(),1)).rename("meanage")
# ② 계절성 (식약처 월별, 2001-2024)
m=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["disease","year","month","cases"])
m=m[(m.disease=="장출혈성대장균감염증")&(m.year.between(2001,2024))]
piv=m.pivot_table(index="year",columns="month",values="cases",aggfunc="sum").reindex(columns=range(1,13),fill_value=0).fillna(0)
th=2*np.pi*(np.arange(1,13)-0.5)/12
rows=[]
for y in piv.index:
    v=piv.loc[y].values; T=v.sum()
    if T==0: continue
    Cc=(v*np.cos(th)).sum()/T; Ss=(v*np.sin(th)).sum()/T
    ang=np.arctan2(Ss,Cc); pm=(ang/(2*np.pi)*12)%12; pm=pm+12 if pm<0 else pm  # 원형평균 '월'(0=1월초)
    rows.append(dict(year=y,N=int(T),gini=gini(v)*100,peak_month=pm+0.5,r=np.hypot(Cc,Ss)))
S=pd.DataFrame(rows).set_index("year")
Sr=S[S.N>=20]  # 저N 민감도

# ===== 그림 (2패널) =====
fig,ax=plt.subplots(1,2,figsize=(13.6,5.6)); covid=2020
# Panel A: 연령별 발생률 4그룹 + 평균연령(우축)
for name,col in ORD:
    if name in rate.columns:
        v=rate[name]; ta,pa=mk(v.index.values,v.values)
        ax[0].plot(v.index,v.values,color=col,lw=2,marker="o",ms=3.5,label=f"{name}: τ={ta:+.2f}{'*' if pa<0.05 else ''}")
ax[0].axvline(covid,color=COVID,ls=(0,(4,3)),lw=1.0)
ax[0].set_title("A. Age-specific incidence rate (4 groups)  ·  KDCA / KOSIS population (2008–2024)",fontweight="bold",loc="left",fontsize=10)
ax[0].set_ylabel("Incidence per 100,000"); ax[0].grid(axis="y",color="#eee",lw=0.6)
ax2=ax[0].twinx(); tam,pam=mk(ma.index.values,ma.values)
ax2.plot(ma.index,ma.values,color="#7a5195",lw=1.3,ls=(0,(3,2)),label=f"Mean age: τ={tam:+.2f}{'*' if pam<0.05 else ''}")
ax2.set_ylabel("Mean age of cases (yr)",color="#7a5195",fontsize=9); ax2.tick_params(labelcolor="#7a5195",labelsize=8)
l1,la1=ax[0].get_legend_handles_labels(); l2,la2=ax2.get_legend_handles_labels()
ax[0].legend(l1+l2,la1+la2,fontsize=7.8,frameon=False,loc="upper left",title="Age-specific rate: Mann–Kendall τ (* p<0.05)",title_fontsize=7.8)
for s in ("top",): ax[0].spines[s].set_visible(False)
ax[0].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=6))
# Panel B: 계절 Gini + 원형 피크위상/r
tg,pg2=mk(S.index.values,S["gini"].values); tp,pp=mk(S.index.values,S["peak_month"].values); tr,pr=mk(S.index.values,S["r"].values)
z=np.polyfit(S.index.values,S["gini"].values,1); ax[1].plot(S.index,np.polyval(z,S.index.values),color="#c1121f",lw=1.3,ls=(0,(4,3)))
ax[1].plot(S.index,S["gini"].values,color="#e76f51",lw=2,marker="o",ms=4)
ax[1].axvline(covid,color=COVID,ls=(0,(4,3)),lw=1.0)
ax[1].set_title("B. Seasonal concentration & peak timing  ·  source: MFDS monthly",fontweight="bold",loc="left",fontsize=10)
ax[1].set_ylabel("Monthly Gini concentration (×100)"); ax[1].grid(axis="y",color="#eee",lw=0.6)
for s in ("top","right"): ax[1].spines[s].set_visible(False)
ax[1].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=6))
def mo(x):
    mn=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]; i=int(x)%12; return f"{mn[i-1]} (~mo {x:.1f})"
box=(f"Seasonality — Mann–Kendall (all yrs):\n"
     f"  Gini concentration   τ={tg:+.2f} p={pg2:.3f}\n"
     f"  Circular peak month  τ={tp:+.2f} p={pp:.2f}  ({S['peak_month'].iloc[0]:.1f}→{S['peak_month'].iloc[-1]:.1f})\n"
     f"  Vector length r      τ={tr:+.2f} p={pr:.2f}\n"
     f"→ concentration {'DECREASED' if tg<0 and pg2<0.05 else 'n.s.'}; peak timing {'shifted' if pp<0.05 else 'stable'}")
ax[1].text(0.02,0.03,box,transform=ax[1].transAxes,fontsize=7.3,va="bottom",ha="left",family="monospace",
           bbox=dict(boxstyle="round,pad=0.4",fc="#fbfbfb",ec="#ccc"))
fig.suptitle("EHEC — temporal transition: age-specific incidence (rate, 2008–2024) and seasonality (2001–2024)",x=0.012,y=0.995,ha="left",fontsize=13,fontweight="bold")
fig.text(0.012,0.955,"Panel A uses age-specific RATES (cases ÷ age-specific population) to separate real age-group changes from population-ageing dilution; four epidemiologically meaningful bands. "
    "8-band composition shares are in the Supplement. Seasonality (Panel B) is MFDS EHEC-specific monthly counts; peak timing via circular mean.",fontsize=7.4,color="#777")
plt.tight_layout(rect=[0,0,1,0.93])
plt.savefig(os.path.join(EH,"EHEC_Fig_transition_rate.png"),dpi=300,facecolor="white"); plt.close()
print("=== 연령별 발생률 4그룹 MK (2008-2024) ===")
for name,_ in ORD:
    if name in rate.columns:
        v=rate[name]; t,p=mk(v.index.values,v.values); print(f"  {name:<18} {v.iloc[0]:.3f}→{v.iloc[-1]:.3f}/10만  τ={t:+.2f} p={p:.3f}")
print(f"평균연령: {ma.iloc[0]:.1f}→{ma.iloc[-1]:.1f}세  τ={tam:+.2f} p={pam:.3f}")
print(f"계절: Gini τ={tg:+.2f}(p{pg2:.3f}) · 피크월 τ={tp:+.2f}(p{pp:.2f}) {S['peak_month'].iloc[0]:.1f}→{S['peak_month'].iloc[-1]:.1f} · r τ={tr:+.2f}(p{pr:.2f})")
print("저장: EHEC_Fig_transition_rate.png")
