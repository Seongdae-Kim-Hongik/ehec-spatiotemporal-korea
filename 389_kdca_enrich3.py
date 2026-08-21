#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# KDCA 성·연령 8종 심화3 (장티푸스 논문 깊이): (1)연령대별 남:여 발생률비 (2)성별 연령표준화 추세
#  (3)성×연령 상호작용 NB GLM LRT (4)환자 평균연령 시간추세.  산출: 그림 2종/질병 + 마스터 CSV 4종
import os, numpy as np, pandas as pd
import pymannkendall as mk
from scipy.stats import theilslopes, norm, chi2
import statsmodels.api as sm
import patsy
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
SRC=os.path.join(GD,"건상","성별_종합 (3).csv"); POPF=os.path.join(GD,"인구데이터_통합.csv")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); RS=os.path.join(ROOT,"results_ST")
EN={"A형간염":"HepatitisA","E형간염":"HepatitisE","비브리오패혈증":"V_vulnificus","세균성이질":"Shigellosis","장출혈성대장균감염증":"EHEC","장티푸스":"Typhoid","콜레라":"Cholera","파라티푸스":"Paratyphoid"}
DISP={"HepatitisA":"Hepatitis A","HepatitisE":"Hepatitis E","V_vulnificus":"Vibrio vulnificus sepsis","Shigellosis":"Shigellosis","EHEC":"Enterohaemorrhagic E. coli (EHEC)","Typhoid":"Typhoid fever","Cholera":"Cholera","Paratyphoid":"Paratyphoid fever"}
AGES=["0~9세","10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70세 이상"]
AGE_EN=["0–9","10–19","20–29","30–39","40–49","50–59","60–69","≥70"]
AGE_MID=[5,15,25,35,45,55,65,77]  # 대표연령(평균연령 계산용; 70+ → 77)
B5={"0–9":["0~4세","5~9세"],"10–19":["10~14세","15~19세"],"20–29":["20~24세","25~29세"],"30–39":["30~34세","35~39세"],"40–49":["40~44세","45~49세"],"50–59":["50~54세","55~59세"],"60–69":["60~64세","65~69세"],"≥70":["70~74세","75~79세","80~84세","85~89세","90~94세","95~99세","100세 이상"]}
YRS=list(range(2008,2025))
C_M="#2166AC"; C_F="#C1121F"; COVID="#d00000"
mpl.rcParams.update({"font.family":"Helvetica Neue","font.size":10,"axes.edgecolor":"#444"})

# --- 연도×성×연령밴드 국가인구 ---
pop=pd.read_csv(POPF,usecols=["행정구역","날짜","구분","연령대","인구수"])
pop["code"]=pop["행정구역"].str.extract(r"\((\d+)\)")
pop=pop[pop["code"].notna()].copy(); pop["code"]=pop["code"].astype(np.int64)
pop=pop[pop["code"]%100000000==0]      # 시도만→전국합
pop["year"]=pop["날짜"].str[:4].astype(int); pop=pop[(pop.year>=2008)&(pop.year<=2024)]
def popband(year,sex):
    py=pop[(pop.year==year)&(pop["구분"]==sex)]
    return {ae:float(py[py["연령대"].isin(parts)]["인구수"].sum()) for ae,parts in B5.items()}
POPB={(y,sx):popband(y,sx) for y in YRS for sx in ["계","남","여"]}
# 2015 표준인구(전국 계, 연령밴드) — 직접 표준화 가중
STDPOP={ae:POPB[(2015,"계")][ae] for ae in AGE_EN}; STDTOT=sum(STDPOP.values())
STDW={ae:STDPOP[ae]/STDTOT for ae in AGE_EN}

try: df=pd.read_csv(SRC,encoding="utf-8-sig")
except: df=pd.read_csv(SRC,encoding="cp949")
df["cases"]=pd.to_numeric(df["cases"],errors="coerce").fillna(0).astype(int); df["year"]=df.year.astype(int)
df=df[(df.year>=2008)&(df.year<=2024)]

def cell(d,sx,ae_kr,y): return int(d[(d.sex==sx)&(d.age_group==ae_kr)&(d.year==y)]["cases"].sum())

rr_rows=[]; astd_rows=[]; inter_rows=[]; mage_rows=[]
for kdz,eng in EN.items():
    d=df[df.disease==kdz]
    if int(d[(d.sex=="계")&(d.age_group=="계")]["cases"].sum())<30: continue
    OUT=os.path.join(ROOT,eng); os.makedirs(OUT,exist_ok=True)

    # ===== (1) 연령대별 남:여 발생률비 (전기간 누적, Katz 로그 CI) =====
    rr=[]; rr_lo=[]; rr_hi=[]
    for a_kr,ae in zip(AGES,AGE_EN):
        cm_=sum(cell(d,"남",a_kr,y) for y in YRS); cf_=sum(cell(d,"여",a_kr,y) for y in YRS)
        pm=sum(POPB[(y,"남")][ae] for y in YRS); pf=sum(POPB[(y,"여")][ae] for y in YRS)
        if cm_>0 and cf_>0 and pm>0 and pf>0:
            r=(cm_/pm)/(cf_/pf); se=np.sqrt(1/cm_+1/cf_)
            lo=r*np.exp(-1.96*se); hi=r*np.exp(1.96*se)
        else: r=lo=hi=np.nan
        rr.append(r); rr_lo.append(lo); rr_hi.append(hi)
        rr_rows.append(dict(disease=kdz,eng=eng,age=ae,M=cm_,F=cf_,MF_RR=None if np.isnan(r) else round(r,3),
                            lo=None if np.isnan(lo) else round(lo,3),hi=None if np.isnan(hi) else round(hi,3)))
    fig,ax=plt.subplots(figsize=(6.4,4.0)); yv=np.arange(len(AGE_EN))
    for i in range(len(AGE_EN)):
        if not np.isnan(rr[i]):
            ax.plot([rr_lo[i],rr_hi[i]],[yv[i],yv[i]],color="#888",lw=1.2,zorder=1)
            col=C_M if rr[i]>1 else C_F
            ax.scatter(rr[i],yv[i],s=42,color=col,zorder=2,edgecolor="white",linewidth=0.6)
    ax.axvline(1,color="#333",ls=(0,(4,3)),lw=0.9)
    ax.set_yticks(yv);ax.set_yticklabels(AGE_EN);ax.invert_yaxis();ax.set_xscale("log")
    ax.set_xlabel("Male : female incidence rate ratio (log scale)");ax.set_ylabel("Age band (years)")
    ax.set_title(f"{DISP[eng]} — sex incidence rate ratio by age, 2008–2024",fontsize=11,fontweight="bold",loc="left")
    ax.grid(axis="x",color="#eee",lw=0.5); [ax.spines[s].set_visible(False) for s in ("top","right")]
    fig.tight_layout();fig.savefig(os.path.join(OUT,f"{eng}_Fig_sexage_rr.png"),dpi=250,facecolor="white");plt.close(fig)

    # ===== (2) 성별 연령표준화 연간율 + MK/Sen (2015 표준) =====
    def astd_series(sx):
        out=[]
        for y in YRS:
            rate=0.0
            for a_kr,ae in zip(AGES,AGE_EN):
                p=POPB[(y,sx)][ae]; c=cell(d,sx,a_kr,y)
                if p>0: rate+=STDW[ae]*(c/p*1e5)
            out.append(rate)
        return np.array(out)
    fig,ax=plt.subplots(figsize=(6.6,4.0))
    for sx,lab,col,mk_ in [("남","Male",C_M,"o"),("여","Female",C_F,"s")]:
        v=astd_series(sx)
        if (v>0).sum()>=5:
            lr=np.log(v+1e-6); r=mk.hamed_rao_modification_test(lr); sl=theilslopes(lr,np.array(YRS,float))[0]
            sen=(np.exp(sl)-1)*100
            astd_rows.append(dict(disease=kdz,eng=eng,sex=lab,MK_tau=round(r.Tau,3),MK_p=r.p,astd_Sen_pct=round(sen,2),trend=r.trend))
            ax.plot(YRS,v,color=col,lw=1.7,marker=mk_,ms=3.4,label=f"{lab} (Sen {sen:+.1f}%/yr)")
        else:
            ax.plot(YRS,v,color=col,lw=1.7,marker=mk_,ms=3.4,label=lab)
    ax.axvline(2020,color=COVID,ls=(0,(4,3)),lw=0.8)
    ax.set_xlabel("Year");ax.set_ylabel("Age-standardised incidence /100,000 (2015 ref)")
    ax.set_title(f"{DISP[eng]} — sex-specific age-standardised trend, 2008–2024",fontsize=11,fontweight="bold",loc="left")
    ax.legend(frameon=False,fontsize=8.5);ax.grid(color="#eee",lw=0.6)
    [ax.spines[s].set_visible(False) for s in ("top","right")];ax.xaxis.set_major_locator(MaxNLocator(6,integer=True))
    fig.tight_layout();fig.savefig(os.path.join(OUT,f"{eng}_Fig_agestd_sex.png"),dpi=250,facecolor="white");plt.close(fig)

    # ===== (3) 성×연령 상호작용 NB2 LRT (연도×성×연령 반복측정으로 과대분산 추정) =====
    rows=[]
    for y in YRS:
        for sx,sl_ in [("남","Male"),("여","Female")]:
            for a_kr,ae in zip(AGES,AGE_EN):
                py=POPB[(y,sx)][ae]
                if py>0: rows.append(dict(year=y,sex=sl_,age=ae,cases=cell(d,sx,a_kr,y),logpy=np.log(py)))
    md=pd.DataFrame(rows)
    try:
        yf,Xf=patsy.dmatrices("cases ~ C(sex)*C(age) + year",md,return_type="dataframe")
        _,Xr=patsy.dmatrices("cases ~ C(sex)+C(age) + year",md,return_type="dataframe")
        off=md["logpy"].values
        ff=sm.NegativeBinomial(yf,Xf,offset=off).fit(disp=0,maxiter=300)
        fr=sm.NegativeBinomial(yf,Xr,offset=off).fit(disp=0,maxiter=300)
        dfree=int(Xf.shape[1]-Xr.shape[1])
        lr_stat=2*(ff.llf-fr.llf); pval=chi2.sf(max(lr_stat,0),dfree)
        inter_rows.append(dict(disease=kdz,eng=eng,LRT_chi2=round(lr_stat,2),df=dfree,p=pval,
                               interaction="yes" if pval<0.05 else "no"))
    except Exception as e:
        inter_rows.append(dict(disease=kdz,eng=eng,LRT_chi2=None,df=None,p=None,interaction=f"NA ({type(e).__name__})"))

    # ===== (4) 환자 평균연령 시간추세 (성별) =====
    for sx,lab in [("계","Total"),("남","Male"),("여","Female")]:
        ma=[]
        for y in YRS:
            cs=np.array([cell(d,sx,a_kr,y) for a_kr in AGES],float); tot=cs.sum()
            ma.append(float(np.dot(cs,AGE_MID)/tot) if tot>0 else np.nan)
        ma=np.array(ma); ok=~np.isnan(ma)
        if ok.sum()>=6:
            yy=np.array(YRS)[ok]; vv=ma[ok]; r=mk.hamed_rao_modification_test(vv); sl=theilslopes(vv,yy.astype(float))[0]
            mage_rows.append(dict(disease=kdz,eng=eng,sex=lab,mean_age_first=round(vv[0],1),mean_age_last=round(vv[-1],1),
                                  Sen_yr=round(sl,3),MK_tau=round(r.Tau,3),MK_p=r.p,trend=r.trend))
    print(f"[{kdz}] sexage-RR + agestd-sex + interaction + mean-age OK")

pd.DataFrame(rr_rows).to_csv(os.path.join(RS,"MASTER_kdca_agesex_rr.csv"),index=False,encoding="utf-8-sig")
pd.DataFrame(astd_rows).to_csv(os.path.join(RS,"MASTER_kdca_agestd_sex.csv"),index=False,encoding="utf-8-sig")
pd.DataFrame(inter_rows).to_csv(os.path.join(RS,"MASTER_kdca_sexageinter.csv"),index=False,encoding="utf-8-sig")
pd.DataFrame(mage_rows).to_csv(os.path.join(RS,"MASTER_kdca_meanage.csv"),index=False,encoding="utf-8-sig")
print("\n마스터: MASTER_kdca_agesex_rr / agestd_sex / sexageinter / meanage")
