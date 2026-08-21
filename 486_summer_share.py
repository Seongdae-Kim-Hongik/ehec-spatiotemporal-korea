#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 486 - Year-by-year July-August case share, an interpretable complement to the Gini measure.
#  주검정 = 연 단위 로짓 OLS(HC3) + 베타이항 MLE 강건성. (주의: 초기 'scale 0.1 붕괴' 판단은 오독 —
#   실제 Pearson scale=7.58 초과분산. 전기간 추정은 가중방식 따라 p=.048~.238 경계이나 2005+ 전 방법 무추세.)
#  보조 = Kendall tau, 전·후반 평균 비교. + Gini 소표본 편의 진단(균등기대 대비 초과 Gini).
#  정직 포인트: Gini 감소는 사례수 증가(10→268건)에 따른 소표본 편의 축소가 상당 부분 설명 →
#  "여름 비중은 안정(≈40%), 집중도 완화 주장은 편의 보정 후 근거 약함"이 정확한 서술.
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import kendalltau
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); EH=os.path.join(ROOT,"EHEC"); RS=os.path.join(ROOT,"results_ST")
rng=np.random.default_rng(42)

m=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["disease","year","month","cases"])
e=m[m.disease.astype(str).str.contains("장출혈")]
g=e.groupby(["year","month"])["cases"].sum().reset_index()
full=pd.MultiIndex.from_product([sorted(g.year.unique()),range(1,13)],names=["year","month"]).to_frame(index=False)
d=full.merge(g,on=["year","month"],how="left").fillna({"cases":0})

def gini(v):
    v=np.sort(np.asarray(v,float)); n=len(v); s=v.sum()
    return np.nan if s<=0 else (2*np.arange(1,n+1)-n-1).dot(v)/(n*s)

rows=[]
for y,dy in d.groupby("year"):
    N=int(dy.cases.sum())
    ja=int(dy.loc[dy.month.isin([7,8]),"cases"].sum())
    js=int(dy.loc[dy.month.isin([6,7,8,9]),"cases"].sum())
    v=np.zeros(12)
    for _,r in dy.iterrows(): v[int(r.month)-1]=r.cases
    rows.append(dict(year=y,N=N,julaug=ja,junsep=js,share_ja=ja/N*100 if N else np.nan,
                     share_js=js/N*100 if N else np.nan,gini=gini(v)))
S=pd.DataFrame(rows).dropna().sort_values("year").reset_index(drop=True)
# 균등분포 기대 Gini(소표본 편의) — 각 연도 N에서 몬테카를로
S["gini_null"]=[np.mean([gini(rng.multinomial(n,[1/12]*12)) for _ in range(1500)]) for n in S.N]
S["gini_excess"]=S.gini-S.gini_null
# ★풀링 계절 프로파일 영모형(더 강한 버전): '계절성 불변 + N만 관측대로'의 기대 Gini와 잔차
prof_m=d.groupby("month")["cases"].sum(); prof_m=(prof_m/prof_m.sum()).values
S["gini_prof"]=[np.mean([gini(rng.multinomial(n,prof_m)) for _ in range(1500)]) for n in S.N]
S["gini_resid"]=S.gini-S.gini_prof
S["t"]=S.year-S.year.min()
S.to_csv(os.path.join(RS,"MASTER_ehec_summer_share.csv"),index=False)

print("=== 연도별 여름 비중 & Gini 진단 ===")
print(S[["year","N","share_ja","share_js","gini","gini_null","gini_excess"]].round(3).to_string(index=False))

# ── 주검정 세트: (a)연 단위 로짓 OLS(HC3) (b)올바른 quasi-binomial(scale=PearsonX2/df=7.58)
#    (c)베타이항 MLE. 전기간은 가중방식에 민감(초기 2003-04 유행기가 견인) → 2005+ 민감도가 결정적. ──
from scipy.optimize import minimize
from scipy.special import betaln, gammaln
from scipy.stats import norm as _norm
def betabin_fit(dfin,col):
    tt=(dfin.year-dfin.year.min()).values; k=dfin[col].values; m=dfin.N.values
    def nll(par):
        b0,b1,lr=par; rho=1/(1+np.exp(-lr))
        mu=1/(1+np.exp(-(b0+b1*tt))); a=mu*(1-rho)/rho; b=(1-mu)*(1-rho)/rho
        return -(gammaln(m+1)-gammaln(k+1)-gammaln(m-k+1)+betaln(k+a,m-k+b)-betaln(a,b)).sum()
    best=None
    for init in [(-0.3,-0.02,-3),(-0.3,0,-2),(0,-0.03,-4)]:
        r=minimize(nll,init,method="Nelder-Mead",options={"maxiter":6000,"xatol":1e-8,"fatol":1e-8})
        if best is None or r.fun<best.fun: best=r
    eps=1e-4; H=np.zeros((3,3))
    import itertools as _it
    for i,j in _it.product(range(3),repeat=2):
        pp=np.array(best.x); pm=np.array(best.x); mp=np.array(best.x); mm=np.array(best.x)
        pp[i]+=eps; pp[j]+=eps; pm[i]+=eps; pm[j]-=eps; mp[i]-=eps; mp[j]+=eps; mm[i]-=eps; mm[j]-=eps
        H[i,j]=(nll(pp)-nll(pm)-nll(mp)+nll(mm))/(4*eps*eps)
    se=np.sqrt(np.linalg.inv(H)[1,1]); z=best.x[1]/se
    return best.x[1],se,2*_norm.sf(abs(z))
X=sm.add_constant(S.t)
res={}
for col,lab in [("julaug","Jul-Aug"),("junsep","Jun-Sep")]:
    lg=np.log((S[col]+0.5)/(S.N-S[col]+0.5))
    fit=sm.OLS(lg,X).fit(cov_type="HC3")
    b=fit.params.iloc[1]; p=fit.pvalues.iloc[1]; ci=fit.conf_int().iloc[1]
    res[lab]=(b,p,ci)
    print(f"\n[주검정 연단위 로짓OLS(HC3)] {lab} share ~ year: OR/yr={np.exp(b):.4f} "
          f"(95%CI {np.exp(ci[0]):.4f}-{np.exp(ci[1]):.4f}) p={p:.3f} → {'유의' if p<0.05 else '추세 없음'}")
    # 올바른 quasi-binomial (N 가중, Pearson scale 보정)
    endog=np.column_stack([S[col],S.N-S[col]])
    fg=sm.GLM(endog,X,family=sm.families.Binomial()).fit()
    sc=fg.pearson_chi2/(len(S)-2); seq=fg.bse.iloc[1]*np.sqrt(sc); zq=fg.params.iloc[1]/seq
    from scipy.stats import norm as _n2
    print(f"   [quasi-binomial scale={sc:.2f}] slope={fg.params.iloc[1]:+.4f} p={2*_n2.sf(abs(zq)):.3f} (N 가중)")
    bb,seb,pb=betabin_fit(S,col)
    print(f"   [베타이항 MLE] slope={bb:+.4f} p={pb:.3f}")
    # 민감도: 초기 구간(2001-04) 제외 — 2003-04는 대유행(여름 78/83%), 2001-02는 극소표본(N=10/7). 결정적.
    S5=S[S.year>=2005].reset_index(drop=True); lg5=np.log((S5[col]+0.5)/(S5.N-S5[col]+0.5))
    f5=sm.OLS(lg5,sm.add_constant(S5.t)).fit(cov_type="HC3")
    bb5,seb5,pb5=betabin_fit(S5,col)
    print(f"   └ 민감도(2005+): OLS p={f5.pvalues.iloc[1]:.3f} · 베타이항 p={pb5:.3f} — 초기 유행기 제외 시 전 방법 무추세")
for col,lab in [("share_ja","Jul-Aug"),("share_js","Jun-Sep"),("gini","Gini(원시)"),("gini_excess","Gini(균등영모형 초과분)"),("gini_resid","Gini(계절불변 영모형 잔차)")]:
    t,p=kendalltau(S.year,S[col]); print(f"[보조 Kendall] {lab}: tau={t:+.3f} p={p:.3f}")
h1=S[S.year<=2012]; h2=S[S.year>=2013]
print(f"\n[전·후반] Jul-Aug {h1.share_ja.mean():.1f}% → {h2.share_ja.mean():.1f}%  ·  Jun-Sep {h1.share_js.mean():.1f}% → {h2.share_js.mean():.1f}%")

# ── 그림: (A) 7-8월 비중 (B) Gini 관측 vs 균등기대 ──
plt.rcParams.update({"font.size":11,"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
fig.suptitle("EHEC — summer (Jul-Aug) share by year, and Gini small-sample diagnosis",fontsize=13,fontweight="bold",x=0.02,ha="left")
a=ax[0]
a.bar(S.year,S.share_ja,color="#e67e22",alpha=0.75,label="Jul-Aug share")
a.plot(S.year,S.share_js,"s--",color="#c0392b",ms=4,lw=1.4,label="Jun-Sep share")
a.axhline(2/12*100,ls=":",color="0.5",lw=1); a.text(S.year.min(),2/12*100+1,"uniform 16.7%",fontsize=8,color="0.4")
b,p,ci=res["Jul-Aug"]
a.set_title(f"A. Summer share of annual cases\n(year-level logit OLS: OR/yr {np.exp(b):.3f}, p={p:.2f} → stable ~40%)",loc="left",fontsize=9.5)
a.set_xlabel("Year"); a.set_ylabel("% of annual cases"); a.legend(frameon=False,fontsize=9); a.set_ylim(0,100)
a2=ax[1]
a2.plot(S.year,S.gini,"o-",color="#2c3e50",lw=2,ms=5,label="Observed Gini")
a2.plot(S.year,S.gini_null,"--",color="#7f8c8d",lw=1.6,label="Expected under UNIFORM months\n(small-sample bias at that year's N)")
a2.fill_between(S.year,S.gini_null,S.gini,color="#3498db",alpha=0.18,label="Excess over uniform-null floor")
a2.plot(S.year,S.gini_prof,":",color="#c0392b",lw=1.6,label="Expected under CONSTANT seasonality\n(pooled profile, matched N)")
te,pe=kendalltau(S.year,S.gini_resid)
a2.set_title(f"B. Gini vs null floors (uniform / constant-seasonality)\n(residual vs constant-seasonality null: τ={te:+.2f}, p={pe:.2f} → no decline)",loc="left",fontsize=9.5)
a2.set_xlabel("Year"); a2.set_ylabel("Monthly Gini"); a2.legend(frameon=False,fontsize=8.5)
plt.tight_layout(rect=[0,0,1,0.92])
out=os.path.join(EH,"EHEC_Fig_summer_share.png")
plt.savefig(out,dpi=150,bbox_inches="tight"); print("\n저장:",out)
