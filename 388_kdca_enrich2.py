#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# KDCA 8종 풍부화 2: (3)연령특이 발생률(rate) + rate기준 연령대 추세 (4)성×연령 부담 히트맵
#  연령별 인구=인구데이터_통합.csv(2008~). 산출: 그림 2종/질병 + MASTER_kdca_agerate.csv
import os, numpy as np, pandas as pd
import pymannkendall as mk
from scipy.stats import theilslopes
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
SRC=os.path.join(GD,"건상","성별_종합 (3).csv"); POPF=os.path.join(GD,"인구데이터_통합.csv")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); RS=os.path.join(ROOT,"results_ST")
EN={"A형간염":"HepatitisA","E형간염":"HepatitisE","비브리오패혈증":"V_vulnificus","세균성이질":"Shigellosis","장출혈성대장균감염증":"EHEC","장티푸스":"Typhoid","콜레라":"Cholera","파라티푸스":"Paratyphoid"}
DISP={"HepatitisA":"Hepatitis A","HepatitisE":"Hepatitis E","V_vulnificus":"Vibrio vulnificus sepsis","Shigellosis":"Shigellosis","EHEC":"Enterohaemorrhagic E. coli (EHEC)","Typhoid":"Typhoid fever","Cholera":"Cholera","Paratyphoid":"Paratyphoid fever"}
AGES=["0~9세","10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70세 이상"]; AGE_EN=["0–9","10–19","20–29","30–39","40–49","50–59","60–69","≥70"]
B5={"0–9":["0~4세","5~9세"],"10–19":["10~14세","15~19세"],"20–29":["20~24세","25~29세"],"30–39":["30~34세","35~39세"],"40–49":["40~44세","45~49세"],"50–59":["50~54세","55~59세"],"60–69":["60~64세","65~69세"],"≥70":["70~74세","75~79세","80~84세","85~89세","90~94세","95~99세","100세 이상"]}
YRS=list(range(2008,2025))
mpl.rcParams.update({"font.family":"Helvetica Neue","font.size":10,"axes.edgecolor":"#444"})

# 연도×성×연령밴드 국가인구 (계/남/여)
pop=pd.read_csv(POPF,usecols=["행정구역","날짜","구분","연령대","인구수"])
pop["code"]=pop["행정구역"].str.extract(r"\((\d+)\)")
pop=pop[pop["code"].notna()].copy(); pop["code"]=pop["code"].astype(np.int64)
pop=pop[pop["code"]%100000000==0]      # 시도만→전국합
pop["year"]=pop["날짜"].str[:4].astype(int); pop=pop[(pop.year>=2008)&(pop.year<=2024)]
def popband(year,sex):
    py=pop[(pop.year==year)&(pop["구분"]==sex)]
    return {ae:float(py[py["연령대"].isin(parts)]["인구수"].sum()) for ae,parts in B5.items()}
POPB={(y,sx):popband(y,sx) for y in YRS for sx in ["계","남","여"]}

try: df=pd.read_csv(SRC,encoding="utf-8-sig")
except: df=pd.read_csv(SRC,encoding="cp949")
df["cases"]=pd.to_numeric(df["cases"],errors="coerce").fillna(0).astype(int); df["year"]=df.year.astype(int)
df=df[(df.year>=2008)&(df.year<=2024)]

rows=[]
for kdz,eng in EN.items():
    d=df[df.disease==kdz]
    if int(d[(d.sex=="계")&(d.age_group=="계")]["cases"].sum())<30: continue
    # (3) 연령특이 발생률 행렬 (age × year), sex=계
    R=np.zeros((8,len(YRS)))
    for j,y in enumerate(YRS):
        pb=POPB[(y,"계")]
        for i,(a,ae) in enumerate(zip(AGES,AGE_EN)):
            c=int(d[(d.sex=="계")&(d.age_group==a)&(d.year==y)]["cases"].sum()); p=pb[ae]
            R[i,j]=c/p*1e5 if p>0 else 0
    # rate기준 연령대 추세 (MK/Sen)
    for i,ae in enumerate(AGE_EN):
        v=R[i];
        if (v>0).sum()<5: rows.append(dict(disease=kdz,eng=eng,age=ae,rate_tau=None,rate_p=None,rate_Sen=None,note="sparse")); continue
        lr=np.log(v+1e-6); r=mk.hamed_rao_modification_test(lr); sl=theilslopes(lr,np.array(YRS,float))[0]
        rows.append(dict(disease=kdz,eng=eng,age=ae,rate_tau=round(r.Tau,3),rate_p=r.p,rate_Sen=round((np.exp(sl)-1)*100,2),note=""))
    # 그림3: 연령특이 발생률 시계열 (8밴드 라인)
    fig,ax=plt.subplots(figsize=(6.8,4.2)); cmap=plt.cm.viridis(np.linspace(0,0.92,8))
    for i,ae in enumerate(AGE_EN): ax.plot(YRS,R[i],color=cmap[i],lw=1.5,marker="o",ms=2.6,label=ae)
    ax.set_yscale("symlog",linthresh=0.01);ax.axvline(2020,color="#d00000",ls=(0,(4,3)),lw=0.8)
    ax.set_xlabel("Year");ax.set_ylabel("Age-specific incidence /100,000 (symlog)")
    ax.set_title(f"{DISP[eng]} — age-specific incidence rates, 2008–2024",fontsize=11,fontweight="bold",loc="left")
    ax.legend(frameon=False,fontsize=7.5,ncol=2,title="Age band",title_fontsize=8,loc="best")
    ax.grid(color="#eee",lw=0.5); [ax.spines[s].set_visible(False) for s in ("top","right")]
    fig.tight_layout();fig.savefig(os.path.join(ROOT,eng,f"{eng}_Fig_agerate.png"),dpi=250,facecolor="white");plt.close(fig)
    # (4) 성×연령 부담 히트맵 (전기간 누적 사례수 %, 남/여 × 8밴드)
    H=np.zeros((2,8))
    for si,sx in enumerate(["남","여"]):
        for i,a in enumerate(AGES): H[si,i]=int(d[(d.sex==sx)&(d.age_group==a)]["cases"].sum())
    tot=H.sum(); Hp=H/tot*100 if tot>0 else H
    fig,ax=plt.subplots(figsize=(8.2,2.6)); im=ax.imshow(Hp,cmap="YlOrRd",aspect="auto")
    ax.set_xticks(range(8));ax.set_xticklabels(AGE_EN);ax.set_yticks([0,1]);ax.set_yticklabels(["Male","Female"])
    for si in range(2):
        for i in range(8): ax.text(i,si,f"{Hp[si,i]:.1f}",ha="center",va="center",fontsize=8,color="#111" if Hp[si,i]<Hp.max()*0.6 else "white")
    cb=fig.colorbar(im,fraction=0.025,pad=0.02);cb.set_label("% of all cases",fontsize=8)
    ax.set_title(f"{DISP[eng]} — sex × age burden (% of cumulative cases, 2001–2024)",fontsize=10.5,fontweight="bold",loc="left")
    ax.set_xlabel("Age band (years)",fontsize=9)
    fig.tight_layout();fig.savefig(os.path.join(ROOT,eng,f"{eng}_Fig_sexage_heat.png"),dpi=250,facecolor="white");plt.close(fig)
    print(f"[{kdz}] age-rate trend + sex×age heat OK")

pd.DataFrame(rows).to_csv(os.path.join(RS,"MASTER_kdca_agerate.csv"),index=False,encoding="utf-8-sig")
print("\n마스터: MASTER_kdca_agerate.csv")
