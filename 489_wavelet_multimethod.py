#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 489 - Seasonality change: nine complementary tests of the early-versus-late difference, reported side by side.
#  A single borderline p value is not a basis for a claim, so no method is privileged over the others.
#  기존(487): Pettitt / 연블록 순열 / count-matched Poisson·NB 영모형 / (383) 블록부트스트랩 / 조화회귀.
#  신규 3종:
#   [A] AR(1) 서로게이트 검정 — 디트렌드 계열의 lag-1 자기상관·분산을 보존한 적색잡음 1,000개로
#       '후반-전반 연주기 점유율 차이'의 영분포 생성 (Torrence&Compo 관례의 시계열 버전).
#       질문: "정상(stationary) 적색잡음만으로 이런 상승이 나오는가?"
#   [B] 전·후반 Global Wavelet Spectrum(연주기 대역 점유율) 차이 — 같은 AR(1) 서로게이트 영분포.
#   [C] STL 계절진폭의 연도별 추세 — 시간영역 교차검증(웨이블릿 무관): 연도별 계절성분 SD(log척도=상대진폭)
#       → Kendall·Pettitt. log 척도라 count-scaling 부분 통제, 잔여는 캐비앳.
#  산출: results_ST/MASTER_ehec_wavelet_multimethod.csv (방법×통계량×p×판정) + 콘솔 비교표.
#  ★사후(489b): 방법 C는 count-floor 아티팩트로 판명 — 계절불변+N관측 영모형에서 tau 평균 +0.934(관측 +0.942, p=0.533).
import os, numpy as np, pandas as pd
import pywt
from statsmodels.tsa.seasonal import STL
from scipy.stats import kendalltau
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); RS=os.path.join(ROOT,"results_ST")
np.random.seed(42); WAV="cmor1.5-1.0"; NSUR=1000

d0=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["disease","year","week","cases"])
d0["cases"]=pd.to_numeric(d0["cases"],errors="coerce").fillna(0); d0=d0[(d0.week>=1)&(d0.week<=52)]
d0=d0[d0.disease=="장출혈성대장균감염증"]
Y0,Y1=2001,2024
idx=pd.MultiIndex.from_product([range(Y0,Y1+1),range(1,53)],names=["year","week"])
w=d0.groupby(["year","week"])["cases"].sum().reindex(idx,fill_value=0).reset_index()
raw=np.log1p(w.cases.values.astype(float))
stl_full=STL(raw,period=52,robust=True).fit()
x=raw-stl_full.trend
years=w.year.values; nW=len(x)

scales=np.geomspace(4,260,140); periods=pywt.scale2frequency(WAV,scales)**-1
ann=(periods>=48)&(periods<=56)
tt=np.arange(nW)
coi_mask=np.ones((len(scales),nW),bool)
for si,s in enumerate(scales):
    ef=np.sqrt(2)*s; coi_mask[si]=(tt>=ef)&(tt<=nW-1-ef)

def yearly_share(xd):
    """487과 동일: 주별 점유율 선산출→연평균, 커버리지<50% 연도 제외. 반환 {year:share}"""
    P=np.abs(pywt.cwt(xd,scales,WAV,sampling_period=1.0)[0])**2
    Pm=np.where(coi_mask,P,np.nan)
    vt=np.isfinite(Pm[ann]).all(axis=0)
    tot=np.nansum(np.where(np.isfinite(Pm),Pm,0),axis=0)
    a=np.nansum(np.where(np.isfinite(Pm[ann]),Pm[ann],0),axis=0)
    st=np.full(nW,np.nan); ok=vt&(tot>0); st[ok]=a[ok]/tot[ok]
    out={}
    for y in range(Y0,Y1+1):
        sel=years==y
        if vt[sel].mean()>=0.5: out[y]=np.nanmean(st[sel])
    return out

obs_sh=yearly_share(x)
yrs=sorted(obs_sh); early=[y for y in yrs if y<=2012]; late=[y for y in yrs if y>=2013]
obs_diff=np.mean([obs_sh[y] for y in late])-np.mean([obs_sh[y] for y in early])
print(f"[관측] 후반-전반 연주기 점유율 차이 = {obs_diff:+.4f}  (분석연도 {yrs[0]}-{yrs[-1]}, n={len(yrs)})")

# ── [A]+[B] AR(1) 서로게이트 ──
rho=np.corrcoef(x[:-1],x[1:])[0,1]; sig=np.std(x)
print(f"[AR(1) 적합] lag-1 rho={rho:.3f}, sd={sig:.3f}")
def ar1_surrogate(r):
    e=r.normal(0,1,nW); s=np.zeros(nW); s[0]=e[0]
    for t in range(1,nW): s[t]=rho*s[t-1]+e[t]*np.sqrt(1-rho**2)
    return s*sig
def gws_share(xd,y_lo,y_hi):
    P=np.abs(pywt.cwt(xd,scales,WAV,sampling_period=1.0)[0])**2
    Pm=np.where(coi_mask,P,np.nan)
    sel=(years>=y_lo)&(years<=y_hi)
    g=np.nanmean(Pm[:,sel],axis=1)   # scale별 평균 = GWS
    return np.nansum(g[ann])/np.nansum(g)
obs_gws=gws_share(x,2013,2023)-gws_share(x,2002,2012)   # yearly_share 유효연도(2002-2023)와 동기화
nullA=[]; nullB=[]
rngA=np.random.default_rng(7)
for k in range(NSUR):
    ssur=ar1_surrogate(rngA)
    sh=yearly_share(ssur)
    e=[sh[y] for y in early if y in sh]; l=[sh[y] for y in late if y in sh]
    if e and l: nullA.append(np.mean(l)-np.mean(e))
    nullB.append(gws_share(ssur,2013,2023)-gws_share(ssur,2002,2012))
nullA=np.array(nullA); nullB=np.array(nullB)
pA=(np.sum(nullA>=obs_diff)+1)/(len(nullA)+1); pA2=(np.sum(np.abs(nullA-nullA.mean())>=abs(obs_diff-nullA.mean()))+1)/(len(nullA)+1)   # MC +1 보정
pB=(np.sum(nullB>=obs_gws)+1)/(len(nullB)+1); pB2=(np.sum(np.abs(nullB-nullB.mean())>=abs(obs_gws-nullB.mean()))+1)/(len(nullB)+1)
print(f"[A: AR(1) 서로게이트, 연점유율 차이] 관측 {obs_diff:+.4f} · 영 {nullA.mean():+.4f} · 단측 p={pA:.3f} · 양측 p={pA2:.3f}")
print(f"[B: AR(1) 서로게이트, GWS 점유율 차이] 관측 {obs_gws:+.4f} · 영 {nullB.mean():+.4f} · 단측 p={pB:.3f} · 양측 p={pB2:.3f}")

# ── [C] STL 계절진폭(연도별 SD, log척도) 추세 ──
seas=stl_full.seasonal
amp=[]
for y in range(Y0,Y1+1):
    sel=years==y; amp.append(dict(year=y,amp=np.std(seas[sel])))
A=pd.DataFrame(amp)
tC,pC=kendalltau(A.year,A.amp)
def pettitt(v):
    n=len(v); U=np.zeros(n)
    for t in range(n): U[t]=np.sum(np.sign(v[t+1:,None]-v[None,:t+1]))
    K=np.max(np.abs(U)); tau=int(np.argmax(np.abs(U)))
    return tau,min(1.0,2*np.exp(-6*K**2/(n**3+n**2)))
tauC,pPetC=pettitt(A.amp.values)
print(f"[C: STL 계절진폭(log, 연 SD) 추세] Kendall tau={tC:+.3f} p={pC:.3f} · Pettitt 변화점={A.year.iloc[tauC]} p={pPetC:.3f}")
A.to_csv(os.path.join(RS,"MASTER_ehec_stl_amp.csv"),index=False)

# ── 종합 비교표 ──
rows=[
 dict(method="Pettitt change-point (yearly share)",stat="CP=2014",p="0.003",note="반보수적: CWT 평활 자기상관"),
 dict(method="Year-label permutation",stat=f"diff={obs_diff:+.3f}",p="0.001",note="반보수적: 동일 사유"),
 dict(method="Count-matched null (Poisson, 2000)",stat="expected +0.057",p="one-sided 0.035 / two-sided 0.064",note="사례수 효과 통제 — 기준 검정"),
 dict(method="Count-matched null (NB k=2.08, 1000)",stat="expected +0.048",p="one-sided 0.028",note="과산포 민감도"),
 dict(method="Block bootstrap of half-spectra (383)",stat="bootstrap MEDIAN of half-share diff +0.026 (obs +0.111)",p="95% CI [-0.060,+0.124]",note="참조행: scale-invariant, 9방법 비교와 통계량 상이"),
 dict(method="Harmonic amplitude x period interaction",stat="LR df=2",p="0.175",note="시간영역, 진폭 관점"),
 dict(method="[NEW A] AR(1) surrogate (yearly share)",stat=f"obs {obs_diff:+.3f} vs null {nullA.mean():+.3f}",p=f"one-sided {pA:.3f} / two-sided {pA2:.3f}",note="정상 적색잡음 영모형 — 사례수 효과 미통제"),
 dict(method="[NEW B] AR(1) surrogate (GWS share)",stat=f"obs {obs_gws:+.3f} vs null {nullB.mean():+.3f}",p=f"one-sided {pB:.3f} / two-sided {pB2:.3f}",note="반기 GWS 관점"),
 dict(method="[NEW C-naive] STL seasonal amplitude trend",stat=f"tau={tC:+.3f}; CP {A.year.iloc[tauC]}",p=f"Kendall {pC:.3f} naive",note="사례수 미통제 — count-floor(489b로 판명)"),
 dict(method="[NEW C-fair] STL amplitude vs count-matched null",stat="null tau=+0.934 vs obs +0.942",p="one-sided 0.533 (489b)",note="계절불변+N관측 영모형 — naive 추세 완전 재현"),
]
M=pd.DataFrame(rows); M.to_csv(os.path.join(RS,"MASTER_ehec_wavelet_multimethod.csv"),index=False)
print("\n===== 방법 × 결과 비교표 =====")
for _,r in M.iterrows(): print(f"  {r['method']:<44} {r['stat']:<34} p: {r['p']:<34} {r['note']}")
print("\n저장: MASTER_ehec_wavelet_multimethod.csv, MASTER_ehec_stl_amp.csv")
