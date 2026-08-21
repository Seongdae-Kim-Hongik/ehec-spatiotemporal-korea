#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 법정감염병 8종 성·연령 descriptive (KDCA 성별_종합, 2001-2024, 17시도 집계=전국)
#  Descriptive sex- and age-specific epidemiology, restricted to the eight outcomes with age data available.
#  산출: 각 질병폴더 <eng>_Fig_demographics.png (3패널: 성별추세·연령분포·연령x연도) + 마스터 CSV.
import os
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
SRC=os.path.join(GD,"건상","성별_종합 (3).csv")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); RS=os.path.join(ROOT,"results_ST")
EN={"A형간염":"HepatitisA","E형간염":"HepatitisE","비브리오패혈증":"V_vulnificus","세균성이질":"Shigellosis","장출혈성대장균감염증":"EHEC","장티푸스":"Typhoid","콜레라":"Cholera","파라티푸스":"Paratyphoid"}
DISP={"HepatitisA":"Hepatitis A","HepatitisE":"Hepatitis E","V_vulnificus":"Vibrio vulnificus sepsis","Shigellosis":"Shigellosis","EHEC":"Enterohaemorrhagic E. coli (EHEC)","Typhoid":"Typhoid fever","Cholera":"Cholera","Paratyphoid":"Paratyphoid fever"}
AGES=["0~9세","10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70세 이상"]
AGE_EN=["0–9","10–19","20–29","30–39","40–49","50–59","60–69","≥70"]

try: df=pd.read_csv(SRC,encoding="utf-8-sig")
except: df=pd.read_csv(SRC,encoding="cp949")
df["cases"]=pd.to_numeric(df["cases"],errors="coerce").fillna(0).astype(int)
df["year"]=df["year"].astype(int)
df=df[(df["year"]>=2001)&(df["year"]<=2024)]           # 2025 부분연도 제외

C_M="#2166AC"; C_F="#C1121F"; COVID="#d00000"
mpl.rcParams.update({"font.family":"Helvetica Neue","font.size":10.5,"axes.edgecolor":"#444","axes.linewidth":0.8})

summ=[]
for kdz,eng in EN.items():
    d=df[df["disease"]==kdz]
    tot=int(d[(d.sex=="계")&(d.age_group=="계")]["cases"].sum())
    if tot==0: continue
    # 성별
    m=int(d[(d.sex=="남")&(d.age_group=="계")]["cases"].sum()); f=int(d[(d.sex=="여")&(d.age_group=="계")]["cases"].sum())
    # 성별 연도추세
    sx=d[(d.sex.isin(["남","여"]))&(d.age_group=="계")].groupby(["year","sex"])["cases"].sum().unstack(fill_value=0).reindex(range(2001,2025),fill_value=0)
    # 연령분포(sex=계), 미상 별도
    ag=d[(d.sex=="계")&(d.age_group.isin(AGES))].groupby("age_group")["cases"].sum().reindex(AGES,fill_value=0)
    unk=int(d[(d.sex=="계")&(d.age_group=="연령미상")]["cases"].sum())
    known=int(ag.sum())
    # 연령x연도 (sex=계, % within year)
    ay=d[(d.sex=="계")&(d.age_group.isin(AGES))].groupby(["year","age_group"])["cases"].sum().unstack(fill_value=0).reindex(index=range(2001,2025),columns=AGES,fill_value=0)
    aypct=ay.div(ay.sum(axis=1).replace(0,np.nan),axis=0)*100
    # 대표 연령대
    domage=AGE_EN[int(np.argmax(ag.values))] if known>0 else "—"
    med_share70=ag.values[-1]/known*100 if known>0 else 0
    summ.append(dict(disease=kdz,eng=eng,total=tot,male=m,female=f,MF_ratio=round(m/max(f,1),2),
        pct_male=round(m/tot*100,1),age_unknown=unk,pct_unknown=round(unk/tot*100,1),
        dominant_age=domage,pct_ge70=round(med_share70,1)))
    # ---- Figure: 3패널 ----
    fig,ax=plt.subplots(1,3,figsize=(15,4.6),gridspec_kw={"width_ratios":[1.15,0.9,1.35],"wspace":0.28})
    # A) 성별 연도추세
    yrs=np.arange(2001,2025)
    ax[0].plot(yrs,sx["남"],color=C_M,lw=1.8,marker="o",ms=3,label="Male")
    ax[0].plot(yrs,sx["여"],color=C_F,lw=1.8,marker="o",ms=3,label="Female")
    ax[0].axvline(2020,color=COVID,ls=(0,(4,3)),lw=1);ax[0].annotate("2020 COVID-19",xy=(2020,ax[0].get_ylim()[1]),xytext=(3,-2),textcoords="offset points",color=COVID,fontsize=8,va="top")
    ax[0].set_title("Annual cases by sex",fontsize=11,fontweight="bold");ax[0].set_xlabel("Year");ax[0].set_ylabel("Cases")
    ax[0].legend(frameon=False,fontsize=9);ax[0].xaxis.set_major_locator(MaxNLocator(6,integer=True))
    for sp in ("top","right"):ax[0].spines[sp].set_visible(False)
    ax[0].grid(axis="y",color="#eee",lw=0.6)
    # B) 연령분포 (수평 막대, % of known)
    pct=ag.values/max(known,1)*100
    ax[1].barh(range(len(AGE_EN)),pct,color="#4C72B0",edgecolor="white")
    ax[1].set_yticks(range(len(AGE_EN)));ax[1].set_yticklabels(AGE_EN,fontsize=9);ax[1].invert_yaxis()
    ax[1].set_title("Age distribution (%)",fontsize=11,fontweight="bold");ax[1].set_xlabel("% of cases (known age)")
    for sp in ("top","right"):ax[1].spines[sp].set_visible(False)
    ax[1].grid(axis="x",color="#eee",lw=0.6)
    # C) 연령x연도 히트맵 (% within year)
    im=ax[2].imshow(aypct.T.values,aspect="auto",cmap="YlOrRd",vmin=0,vmax=np.nanmax(aypct.values) if known>0 else 1,origin="upper",
        extent=[2001,2024,len(AGE_EN)-0.5,-0.5])
    ax[2].set_yticks(range(len(AGE_EN)));ax[2].set_yticklabels(AGE_EN,fontsize=8)
    ax[2].set_title("Age composition over time (% within year)",fontsize=11,fontweight="bold");ax[2].set_xlabel("Year")
    ax[2].axvline(2020,color="#111",ls=(0,(3,2)),lw=0.8)
    cb=fig.colorbar(im,ax=ax[2],fraction=0.046,pad=0.02);cb.ax.tick_params(labelsize=7);cb.set_label("%",fontsize=8)
    ttl=f"{DISP[eng]} — sex & age descriptive, 2001–2024 (n={tot:,}; {round(m/tot*100)}% male; unknown age {round(unk/tot*100,1)}%)"
    fig.suptitle(ttl,x=0.01,y=0.985,ha="left",fontsize=12.5,fontweight="bold",color="#1a1a2e")
    fig.text(0.01,0.915,"KDCA statutory notifiable-disease notifications (national, 17 provinces aggregated). Panel C column-normalised to show shifting age composition.",ha="left",fontsize=8,color="#666")
    fig.subplots_adjust(top=0.82,bottom=0.13,left=0.055,right=0.995)
    out=os.path.join(ROOT,eng,f"{eng}_Fig_demographics.png")
    os.makedirs(os.path.dirname(out),exist_ok=True)
    fig.savefig(out,dpi=250,facecolor="white");plt.close(fig)
    print(f"[{kdz}] {eng}_Fig_demographics.png (n={tot:,}, M/F {m}/{f}, dom {domage})")

S=pd.DataFrame(summ)
S.to_csv(os.path.join(RS,"MASTER_sexage_descriptive.csv"),index=False,encoding="utf-8-sig")
print("\n마스터: results_ST/MASTER_sexage_descriptive.csv")
print(S[["disease","total","male","female","MF_ratio","dominant_age","pct_ge70","pct_unknown"]].to_string(index=False))
