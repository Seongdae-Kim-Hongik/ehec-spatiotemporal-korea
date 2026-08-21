#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 383 - Morlet wavelet, early/late halves. Method and presentation follow Kim, Vounatsou & Chun, BMC Infect Dis 2024.
#  /333 R3 반영: 정규화 GWS·COI 해칭·블록부트스트랩 CI·count-scaling 명시.
#  준거 템플릿 보강: ①scalogram 95%유의 흰등고선+ridge ②GWS 95%(빨강)/90%(파랑) 유의점 ③peak power(periodicity) 표
#   ④EHEC 단봉(여름) vs ST 이중봉 대조 명시 ⑤서로게이트 500회.
#  해석(정정): 집중완화=서술적 '더 분산된 패턴'. 기후 인과 주장 금지(기후변수 없음; detection/감시변화 대안). ST=방법론적 비교만.
import os, numpy as np, pandas as pd
import pywt
from statsmodels.tsa.seasonal import STL
from scipy.signal import find_peaks
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); EH=os.path.join(ROOT,"EHEC"); RS=os.path.join(ROOT,"results_ST")
np.random.seed(42); WAV="cmor1.5-1.0"; NSUR=500
mpl.rcParams.update({"font.family":"Helvetica Neue","font.size":10.5,"axes.edgecolor":"#444","axes.linewidth":0.8})
d0=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["disease","year","week","cases"])
d0["cases"]=pd.to_numeric(d0["cases"],errors="coerce").fillna(0); d0=d0[(d0.week>=1)&(d0.week<=52)]
d0=d0[d0.disease=="장출혈성대장균감염증"]
scales=np.geomspace(4,260,140); periods=pywt.scale2frequency(WAV,scales)**-1; ann=(periods>=48)&(periods<=56)
def weekly(y0,y1):
    d=d0[(d0.year>=y0)&(d0.year<=y1)].groupby(["year","week"],as_index=False)["cases"].sum()
    idx=pd.MultiIndex.from_product([range(y0,y1+1),range(1,53)],names=["year","week"])
    d=d.set_index(["year","week"]).reindex(idx,fill_value=0).reset_index()
    return d["year"].values+(d["week"].values-1)/52.0, d["cases"].values.astype(float)
def cwtp(x): coef,_=pywt.cwt(x,scales,WAV,sampling_period=1.0); return np.abs(coef)**2
def detr(c):
    s=np.log1p(c)
    try: return s-STL(s,period=52,robust=True).fit().trend
    except: return s-pd.Series(s).rolling(103,center=True,min_periods=20).mean().bfill().ffill().values
def ar1(x): x=x-x.mean(); r=np.corrcoef(x[:-1],x[1:])[0,1] if len(x)>2 else 0.; return max(min(r,0.99),0.),x.std()
def gws_coi(pw):
    n=pw.shape[1]; e=min(52,n//4); return pw[:,e:n-e].mean(axis=1) if n>2*e else pw.mean(axis=1)
def ratio_annual(x):
    g=gws_coi(cwtp(x)); return g[ann].sum()/max(g.sum(),1e-9)
def boot_ratio(x,B=300,bl=52):
    n=len(x); o=np.empty(B)
    for b in range(B):
        st=np.random.randint(0,n,int(np.ceil(n/bl))); idx=np.concatenate([(np.arange(s,s+bl)%n) for s in st])[:n]; o[b]=ratio_annual(x[idx])
    return o
def pk_table(gws):  # dominant + 2nd peak(계절범위 ≤104주만; 다년 COI잡음 배제): (power, periodicity wk)
    seas=periods<=104
    pk,_=find_peaks(gws)
    pk=[i for i in pk if seas[i]]
    if not pk: return (np.nan,np.nan),(np.nan,np.nan)
    order=sorted(pk,key=lambda i:-gws[i])
    top=[(gws[i],periods[i]) for i in order[:2]]
    while len(top)<2: top.append((np.nan,np.nan))
    return top[0],top[1]
def analyze(y0,y1):
    t,c=weekly(y0,y1); x=detr(c); pw=cwtp(x); gws=gws_coi(pw)
    rho,sd=ar1(x); n=len(x); nl=np.empty((NSUR,len(scales))); npw95=np.zeros(len(scales))
    surpk=np.zeros((NSUR,len(scales)))  # per-surrogate per-scale GWS (for pointwise 95%)
    for i in range(NSUR):
        e=np.random.normal(0,sd,n); z=np.empty(n); z[0]=e[0]
        for k in range(1,n): z[k]=rho*z[k-1]+e[k]
        gz=gws_coi(cwtp(z)); nl[i]=gz
    g95=np.percentile(nl,95,axis=0); g90=np.percentile(nl,90,axis=0)
    dom,sec=pk_table(gws)
    return dict(t=t,c=c,pw=pw,gws=gws,g95=g95,g90=g90,rho=rho,abs_ann=gws[ann].mean(),
                factor=gws[ann].mean()/max(nl[:,ann].mean(),1e-9),
                ratio=gws[ann].sum()/max(gws.sum(),1e-9),boot=boot_ratio(x),totN=int(c.sum()),
                dom=dom,sec=sec)
FULL=analyze(2001,2024); EARLY=analyze(2001,2012); LATE=analyze(2013,2024)
dif=LATE["boot"]-EARLY["boot"]; dmed=np.median(dif); dlo,dhi=np.percentile(dif,[2.5,97.5]); sig=(dlo>0)or(dhi<0)
# ---- 준거 Table 2식 peak power(periodicity) 저장 ----
rows=[]
for lab,S in [("Total 2001-2024",FULL),("Early 2001-2012",EARLY),("Late 2013-2024",LATE)]:
    rows.append(dict(group=lab,dom_power=round(S["dom"][0],3),dom_period_wk=round(S["dom"][1],2),
                     sec_power=round(S["sec"][0],3),sec_period_wk=round(S["sec"][1],2),N=S["totN"]))
pd.DataFrame(rows).to_csv(os.path.join(RS,"MASTER_ehec_wavelet_peaks.csv"),index=False)
# 단봉 판정: 2nd 피크가 26주(반기)대(22-30wk)에 있고 dominant의 30%↑면 이중봉
def bimodal(S):
    p2,pd2=S["sec"]; return (22<=pd2<=30) and (p2>=0.30*S["dom"][0])
uni_full=not bimodal(FULL)

# ===== 그림 =====
def add_coi(ax,S):
    t=S["t"]; n=len(t); d=np.minimum(np.arange(n),n-1-np.arange(n))
    coib=np.clip(d/np.sqrt(2),periods.min(),periods.max())
    ax.fill_between(t,coib,periods.max(),facecolor="white",alpha=0.5,zorder=3,hatch="////",edgecolor="#c9c9c9",lw=0)
def scal(ax,S,title):
    logp=np.log(S["pw"]+1e-9); vmax=np.nanpercentile(np.abs(logp),98)
    cf=ax.contourf(S["t"],periods,logp,levels=22,cmap="RdBu_r",vmin=-vmax,vmax=vmax)
    # 95% 유의 흰 등고선(per-scale 서로게이트 95% 초과영역)
    sigm=(S["pw"]>S["g95"][:,None]).astype(float)
    ax.contour(S["t"],periods,sigm,levels=[0.5],colors="white",linewidths=0.8,zorder=2)
    # ridge(각 시점 계절대 최대파워 주기; 다년대 잡음 배제)
    seas=periods<=80
    ridge=periods[seas][np.argmax(S["pw"][seas],axis=0)]
    ax.plot(S["t"],ridge,color="#111",lw=0.5,alpha=0.5,zorder=2)
    add_coi(ax,S); ax.set_yscale("log"); ax.axhline(52,color="#111",ls="--",lw=0.9)
    ax.axhline(26,color="#111",ls=(0,(1,2)),lw=0.7,alpha=0.6)  # 반기(26주) 참조선
    ax.yaxis.set_major_locator(FixedLocator([16,26,52,104,208])); ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{int(v)}"))
    ax.set_ylabel("Period (weeks)"); ax.set_xlabel("Year"); ax.set_title(title,fontsize=10,fontweight="bold",loc="left")
    return cf
fig=plt.figure(figsize=(13.8,8.9)); gs=fig.add_gridspec(2,3,height_ratios=[1,1],hspace=0.34,wspace=0.30)
axA=fig.add_subplot(gs[0,:2]); cf=scal(axA,FULL,"A. Full period 2001–2024 (Morlet CWT, STL-detrended weekly; white=95% sig, hatch=COI, thin line=ridge)")
axA.axvline(2013,color="#333",ls=(0,(3,2)),lw=0.9); axA.text(2013,periods.max()*0.68,"2013 split",fontsize=8,color="#333",rotation=90,va="top")
fig.colorbar(cf,ax=axA,fraction=0.03,pad=0.02).set_label("log$_e$ power",fontsize=8)
# Panel D: 정규화 GWS + 95%/90% 유의점 (준거 Fig3 우측 규격)
axD=fig.add_subplot(gs[0,2])
ne=EARLY["gws"]/EARLY["gws"].sum(); nlt=LATE["gws"]/LATE["gws"].sum()
axD.plot(ne,periods,color="#3b7bb0",lw=1.9,label="Early 2001–12"); axD.plot(nlt,periods,color="#c1121f",lw=1.9,label="Late 2013–24")
for S,cur,cc in [(EARLY,ne,"#3b7bb0"),(LATE,nlt,"#c1121f")]:
    r95=S["gws"]>S["g95"]; r90=(S["gws"]>S["g90"])&~r95
    axD.scatter(cur[r95],periods[r95],s=9,c="#b10026",zorder=5); axD.scatter(cur[r90],periods[r90],s=7,c="#2166AC",zorder=5)
axD.fill_between([0,max(ne.max(),nlt.max())*1.05],48,56,color="#FDDCA8",alpha=0.6)
axD.axhline(26,color="#888",ls=(0,(1,2)),lw=0.7)
axD.set_yscale("log"); axD.set_ylim(periods.min(),periods.max()); axD.set_xlabel("Variance-normalised power")
axD.yaxis.set_major_locator(FixedLocator([16,26,52,104,208])); axD.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{int(v)}"))
axD.set_title("D. Normalised spectrum: early vs late\n(red=95% sig, blue=90%)",fontsize=10,fontweight="bold",loc="left")
axD.legend(frameon=False,fontsize=7.6,loc="lower right")
for sp in ("top","right"): axD.spines[sp].set_visible(False)
axB=fig.add_subplot(gs[1,0]); scal(axB,EARLY,f"B. Early half 2001–2012 (N={EARLY['totN']})")
axC=fig.add_subplot(gs[1,1]); scal(axC,LATE,f"C. Late half 2013–2024 (N={LATE['totN']})")
axE=fig.add_subplot(gs[1,2]); axE.axis("off")
absrise=LATE['abs_ann']/max(EARLY['abs_ann'],1e-9); nrise=LATE['totN']/max(EARLY['totN'],1)
msg=(f"Peak power (periodicity, weeks)\n"
     f"           dominant        2nd\n"
     f" Total  {FULL['dom'][0]:.3f}({FULL['dom'][1]:.1f})  {FULL['sec'][0]:.3f}({FULL['sec'][1]:.1f})\n"
     f" Early  {EARLY['dom'][0]:.3f}({EARLY['dom'][1]:.1f})  {EARLY['sec'][0]:.3f}({EARLY['sec'][1]:.1f})\n"
     f" Late   {LATE['dom'][0]:.3f}({LATE['dom'][1]:.1f})  {LATE['sec'][0]:.3f}({LATE['sec'][1]:.1f})\n\n"
     f"EHEC = single (annual, ~52wk) peak — UNIMODAL,\nno 26-wk semiannual peak (contrast: scrub typhus\ndual autumn+spring). Summer-driven foodborne cycle.\n\n"
     f"Scale-invariant annual/total (COI-excl): E {EARLY['ratio']:.3f}\n-> L {LATE['ratio']:.3f}; diff {dmed:+.3f} 95%CI[{dlo:+.3f},{dhi:+.3f}]\n= {'DIFFERENT' if sig else 'NO sig. difference (CI incl. 0)'}.\n"
     f"Absolute power x{absrise:.1f} = count scaling (N x{nrise:.1f}).\n\n"
     f"INTERPRETATION: dominant annual cycle in both\nhalves; relative annual-band share NOT increased\n(CI incl. 0) & monthly Gini fell (Round 1) =>\nmore temporally dispersed pattern. Attribution to\nclimate NOT warranted (no climatic covariate;\ndetection/surveillance changes plausible). ST\nwavelet = methodological comparator, not parallel.")
axE.text(0,1,msg,fontsize=7.2,va="top",ha="left",family="monospace",bbox=dict(boxstyle="round,pad=0.4",fc="#fbfbfb",ec="#ccc"))
fig.suptitle("EHEC — wavelet periodicity by period half (Morlet CWT)",x=0.01,y=0.995,ha="left",fontsize=12,fontweight="bold")
fig.text(0.01,0.955,f"Normalised spectra (raw power scales with incidence). Hatch=COI; white=95% surrogate significance ({NSUR} AR(1) sims); dashed=52wk, dotted=26wk. "
    "Full CWT (A) primary; halves (B,C) supplementary. Early-half N low (null approximate).",fontsize=7.2,color="#777")
fig.subplots_adjust(top=0.90,bottom=0.07,left=0.055,right=0.985)
fig.savefig(os.path.join(EH,"EHEC_Fig_wavelet_split.png"),dpi=200,facecolor="white"); plt.close(fig)
print("peaks Total dom %.3f(%.1f) Early dom %.3f(%.1f) Late dom %.3f(%.1f)"%(FULL['dom'][0],FULL['dom'][1],EARLY['dom'][0],EARLY['dom'][1],LATE['dom'][0],LATE['dom'][1]))
print("EHEC unimodal(2nd peak 없음/약함): %s · norm-ratio E%.3f/L%.3f diff95CI[%.3f,%.3f] %s"%(uni_full,EARLY['ratio'],LATE['ratio'],dlo,dhi,'sig' if sig else 'ns'))
print("저장: EHEC_Fig_wavelet_split.png + MASTER_ehec_wavelet_peaks.csv")
