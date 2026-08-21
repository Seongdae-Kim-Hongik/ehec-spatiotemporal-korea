#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 487 - Formal tests for the apparent early-versus-late difference in the wavelet spectrum.
#  The normalized spectrum looks different after 2014; the point is to test that impression rather than assert it.
#  전략(383 파이프라인 재사용: MFDS 주별, STL-디트렌드 log1p, Morlet cmor1.5-1.0):
#   (1) 연도별 분산정규화 연주기(48-56wk) 대역 점유율 시계열(COI 제외) 산출
#   (2) Pettitt 변화점 검정 — 변화 시점을 '사전지정 없이' 데이터가 고르게 함(2013 근방인지 확인)
#   (3) 순열검정 — 연도 라벨을 뒤섞어 전·후반 점유율 차이의 영분포 생성(블록=연 단위)
#   (4) ★count-matched 영모형 시뮬레이션(결정적): 계절 프로파일 고정 + 연도별 총 N만 관측대로 →
#       CWT 점유율 상승이 '계절성 불변'에서도 재현되는지. Pettitt·순열은 CWT 평활 자기상관으로 반보수적이므로
#       이 시뮬이 유일하게 공정한 검정. (5) 기존 블록부트스트랩 CI와 삼각 대조.
import os, numpy as np, pandas as pd
import pywt
from statsmodels.tsa.seasonal import STL
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); EH=os.path.join(ROOT,"EHEC"); RS=os.path.join(ROOT,"results_ST")
np.random.seed(42); WAV="cmor1.5-1.0"; NPERM=2000

d0=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["disease","year","week","cases"])
d0["cases"]=pd.to_numeric(d0["cases"],errors="coerce").fillna(0); d0=d0[(d0.week>=1)&(d0.week<=52)]
d0=d0[d0.disease=="장출혈성대장균감염증"]
Y0,Y1=2001,2024
idx=pd.MultiIndex.from_product([range(Y0,Y1+1),range(1,53)],names=["year","week"])
w=d0.groupby(["year","week"])["cases"].sum().reindex(idx,fill_value=0).reset_index()
x=np.log1p(w.cases.values.astype(float))
try: x=x-STL(x,period=52,robust=True).fit().trend
except Exception: x=x-pd.Series(x).rolling(103,center=True,min_periods=20).mean().bfill().ffill().values
years=w.year.values; nW=len(x)

scales=np.geomspace(4,260,140); periods=pywt.scale2frequency(WAV,scales)**-1
ann=(periods>=48)&(periods<=56)
coef,_=pywt.cwt(x,scales,WAV,sampling_period=1.0); P=np.abs(coef)**2   # (scale, time)
# COI: Morlet e-folding ≈ sqrt(2)*scale — 시계열 양끝에서 scale별 마스킹
coi_mask=np.ones_like(P,dtype=bool)
tt=np.arange(nW)
for si,s in enumerate(scales):
    ef=np.sqrt(2)*s
    coi_mask[si]=(tt>=ef)&(tt<=nW-1-ef)
Pm=np.where(coi_mask,P,np.nan)

# (1) 연도별 분산정규화 연주기 점유율 — ★주별 점유율 선산출 후 연평균 (COI 비대칭 마스킹 편향 제거:
#     합산비율 방식은 가장자리에서 연주기(넓은 COI)만 잘려 분자 과소 → 유효 주에서 s(t)=ann(t)/total(t) 먼저)
valid_t=np.isfinite(Pm[ann]).all(axis=0)          # 연주기 대역 전 스케일이 비마스킹인 주
sh_t=np.full(nW,np.nan)
tot_t=np.nansum(np.where(np.isfinite(Pm),Pm,0),axis=0)
ann_t=np.nansum(np.where(np.isfinite(Pm[ann]),Pm[ann],0),axis=0)
ok=valid_t&(tot_t>0)
sh_t[ok]=ann_t[ok]/tot_t[ok]
rows=[]
for y in range(Y0,Y1+1):
    sel=years==y
    cov=valid_t[sel].mean()
    sh=np.nanmean(sh_t[sel]) if cov>0 else np.nan
    rows.append(dict(year=y,share=sh,total=np.nansum(Pm[:,sel]),ann_coverage=cov))
S=pd.DataFrame(rows)
nfull=len(S)
S=S[S.ann_coverage>=0.5].reset_index(drop=True)   # ★엣지 아티팩트 제거: 연주기 대역이 절반 이상 COI에 잘린 해 제외
print(f"[COI 필터] 연주기 커버리지<50% 제외: {nfull-len(S)}개 연도 탈락 → 분석 {S.year.min()}-{S.year.max()} ({len(S)}년)")
S.to_csv(os.path.join(RS,"MASTER_ehec_wavelet_annshare.csv"),index=False)
print("=== 연도별 연주기 대역 점유율 (분산정규화, COI 제외) ===")
print(S.round(4).to_string(index=False))

# (2) Pettitt 변화점 검정
def pettitt(v):
    n=len(v); U=np.zeros(n)
    for t in range(n):
        U[t]=np.sum(np.sign(v[t+1:,None]-v[None,:t+1]))
    K=np.max(np.abs(U)); tau=int(np.argmax(np.abs(U)))
    p=2*np.exp(-6*K**2/(n**3+n**2))
    return tau,K,min(1.0,p)
tau,K,pP=pettitt(S.share.values)
print(f"\n[Pettitt 변화점] 최유력 변화점={S.year.iloc[tau]}년 K={K:.0f} p={pP:.3f} → {'유의한 변화점' if pP<0.05 else '변화점 근거 없음'}")

# (3) 순열검정: 연 블록 라벨 셔플 → 전(≤2012) vs 후(≥2013) 점유율 차이
obs=S[S.year>=2013].share.mean()-S[S.year<=2012].share.mean()
n1=(S.year<=2012).sum(); vals=S.share.values.copy()
null=[]
rng=np.random.default_rng(42)
for _ in range(NPERM):
    rng.shuffle(vals); null.append(vals[n1:].mean()-vals[:n1].mean())
null=np.array(null); pperm=np.mean(np.abs(null)>=abs(obs))
pp_txt=("<0.001" if pperm<0.0005 else f"{pperm:.3f}")
print(f"[순열검정(연 블록)] 후반-전반 점유율 차이={obs:+.4f}  p={pp_txt} (반보수적: CWT 평활 자기상관으로 교환가능성 위반)")
print("[기존 블록부트스트랩] 차이 +0.018, 95%CI [-0.066,+0.114] (383) — 삼각 대조용")

# (4) ★count-matched 영모형: 계절 프로파일 고정(전기간 풀링) + 연도별 N 관측대로 Poisson 재표집
prof=w.groupby("week")["cases"].sum(); prof=prof/prof.sum()
Nyr=w.groupby("year")["cases"].sum()
NSIM=2000
def share_series(xv):
    try: xd=xv-STL(xv,period=52,robust=True).fit().trend
    except Exception: xd=xv-pd.Series(xv).rolling(103,center=True,min_periods=20).mean().bfill().ffill().values
    Ps=np.abs(pywt.cwt(xd,scales,WAV,sampling_period=1.0)[0])**2
    Pms=np.where(coi_mask,Ps,np.nan)
    vt=np.isfinite(Pms[ann]).all(axis=0)
    tt2=np.nansum(np.where(np.isfinite(Pms),Pms,0),axis=0)
    at2=np.nansum(np.where(np.isfinite(Pms[ann]),Pms[ann],0),axis=0)
    st=np.full(nW,np.nan); ok2=vt&(tt2>0); st[ok2]=at2[ok2]/tt2[ok2]
    out={}
    for y in range(Y0,Y1+1):
        sel=years==y
        if vt[sel].mean()>=0.5: out[y]=np.nanmean(st[sel])
    return out
early=[y for y in S.year if y<=2012]; late=[y for y in S.year if y>=2013]
# 주간 과산포 추정(적률법): Var = mu + mu^2/k → NB 민감도용
mu_all=np.concatenate([prof.values*Nyr[y] for y in range(Y0,Y1+1)])
obs_all=w.cases.values.astype(float)
mask=mu_all>0.05
k_hat=np.nansum(mu_all[mask]**2)/max(np.nansum((obs_all[mask]-mu_all[mask])**2-mu_all[mask]),1e-9)
k_hat=max(k_hat,0.3)
print(f"[과산포 적률추정] NB size k={k_hat:.2f} (작을수록 과산포 큼)")
def one_sim(seed,family="pois"):
    r=np.random.default_rng(seed)
    lam=[prof.values*Nyr[y] for y in range(Y0,Y1+1)]
    if family=="pois":
        xs=np.concatenate([r.poisson(l) for l in lam]).astype(float)
    else:  # NB: gamma-Poisson 혼합
        xs=np.concatenate([r.poisson(r.gamma(k_hat,l/k_hat)) for l in lam]).astype(float)
    ss=share_series(np.log1p(xs))
    return np.nanmean([ss.get(y,np.nan) for y in late])-np.nanmean([ss.get(y,np.nan) for y in early])
simdiff=np.array([one_sim(i) for i in range(NSIM)])
NSIM_NB=1000
simdiff_nb=np.array([one_sim(10_000+i,"nb") for i in range(NSIM_NB)])
p_null=np.mean(simdiff>=obs)                              # 단측(사전지정 방향: 증가)
p_two=np.mean(np.abs(simdiff-np.nanmean(simdiff))>=abs(obs-np.nanmean(simdiff)))   # 양측(중심화)
print(f"\n[count-matched 영모형 {NSIM}회] 계절성 '불변'인데 사례수만 관측대로 → 기대 차이 {np.nanmean(simdiff):+.3f} "
      f"(95%범위 [{np.nanpercentile(simdiff,2.5):+.3f},{np.nanpercentile(simdiff,97.5):+.3f}])")
print(f"  관측 {obs:+.3f}: 단측 p={p_null:.3f} · 양측 p={p_two:.3f} (MC오차 ±{np.sqrt(p_null*(1-p_null)/NSIM):.3f})")
print(f"  → 상승분의 약 2/3(점추정 {np.nanmean(simdiff)/obs*100:.0f}%)는 count-scaling 기대치")
p_nb=np.mean(simdiff_nb>=obs)
print(f"[NB 과산포 민감도 {NSIM_NB}회, k={k_hat:.2f}] 기대 {np.nanmean(simdiff_nb):+.3f} "
      f"(95%범위 [{np.nanpercentile(simdiff_nb,2.5):+.3f},{np.nanpercentile(simdiff_nb,97.5):+.3f}]) · 단측 p={p_nb:.3f}")
pd.DataFrame({"sim_diff":simdiff}).to_csv(os.path.join(RS,"MASTER_ehec_wavelet_nullsim.csv"),index=False)


# (4) 그림: 점유율 시계열 + Pettitt U + 순열 영분포
fig,ax=plt.subplots(1,3,figsize=(15,4.3))
fig.suptitle("EHEC wavelet — statistical tests for the 'visible' early/late difference (annual-band share, variance-normalised, COI-excluded)",
             fontsize=12,fontweight="bold",x=0.02,ha="left")
a=ax[0]; a.plot(S.year,S.share,"o-",color="#2c3e50",lw=2,ms=5)
a.axvline(2013.5,ls="--",color="#c0392b",lw=1); a.text(2013.7,a.get_ylim()[1]*0.95,"2013|2014 split",fontsize=8,color="#c0392b")
a.set_title(f"A. Yearly annual-band share\n(Pettitt change-point: {S.year.iloc[tau]}, p={pP:.2f})",loc="left",fontsize=9)
a.set_xlabel("Year"); a.set_ylabel("Annual-band / total power")
b=ax[1]; b.hist(null,bins=40,color="#95a5a6",alpha=0.6,label="label permutation")
b.hist(simdiff,bins=40,color="#2980b9",alpha=0.55,label="count-matched null (constant seasonality)")
b.legend(frameon=False,fontsize=7.5)
b.axvline(obs,color="#c0392b",lw=2); b.set_title(f"B. Nulls: naive permutation vs count-matched\n(obs={obs:+.3f}; count-matched p={p_null:.3f})",loc="left",fontsize=9)
b.set_xlabel("late−early share difference")
c=ax[2]
c.plot(S.year,np.log10(S.total),"s-",color="#7f8c8d",ms=4,lw=1.5)
c.set_title("C. log10 total wavelet power by year\n(count-scaling driver)",loc="left",fontsize=9)
c.set_xlabel("Year"); c.set_ylabel("log10 total power")
for aa in ax: aa.spines["top"].set_visible(False); aa.spines["right"].set_visible(False)
plt.tight_layout(rect=[0,0,1,0.9])
out=os.path.join(EH,"EHEC_Fig_wavelet_changetest.png")
plt.savefig(out,dpi=150,bbox_inches="tight"); print("\n저장:",out)
