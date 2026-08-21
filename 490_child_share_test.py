#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 490 - Formal test of the declining 0-9 year case share, with a demographic decomposition.
#  Null model: the share can fall purely because the child population is shrinking. Formalized as:
#   기대비중_y = (R_pooled × P_young,y) / (R_pooled × P_young,y + P_old,y)
#   (R_pooled = 전 기간 풀링 소아/비소아 발생률 비 — '상대위험 불변 + 인구구조만 변화' 영모형)
#  검정: (a)관측 비중의 베타이항 추세(486 재사용) + 창 민감도(2001+/2006+/2008+)
#        (b)관측-기대 로짓 잔차의 Kendall 추세 — 무추세면 "비중 감소 = 인구구조 효과" 입증.
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kendalltau, norm
from scipy.optimize import minimize
from scipy.special import betaln, gammaln
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); C=os.path.join(ROOT,"data_cache")
EH=os.path.join(ROOT,"EHEC"); RS=os.path.join(ROOT,"results_ST")

ac=pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv"))
pop=pd.read_csv(os.path.join(C,"ST_agepop.csv")).rename(columns={"age10":"age_group"})
cs=ac.groupby(["year","age_group"])["cases"].sum().reset_index()
cs["young"]=cs.age_group=="0~9세"
Y=cs.groupby(["year","young"])["cases"].sum().unstack(fill_value=0).rename(columns={True:"cy",False:"co"}).reset_index()
Y["N"]=Y.cy+Y.co; Y["share"]=Y.cy/Y.N*100
pp=pop.copy(); pp["young"]=pp.age_group=="0~9세"
P=pp.groupby(["year","young"])["pop"].sum().unstack().rename(columns={True:"Py",False:"Po"}).reset_index()
D=Y.merge(P,on="year",how="left")   # 인구는 2008+만 존재

# ── (a) 관측 비중 추세: 베타이항 + 창 민감도 ──
def betabin_fit(df):
    tt=(df.year-df.year.min()).values; k=df.cy.values; m=df.N.values
    def nll(par):
        b0,b1,lr=par; rho=1/(1+np.exp(-lr))
        mu=1/(1+np.exp(-(b0+b1*tt))); a=mu*(1-rho)/rho; b=(1-mu)*(1-rho)/rho
        return -(gammaln(m+1)-gammaln(k+1)-gammaln(m-k+1)+betaln(k+a,m-k+b)-betaln(a,b)).sum()
    best=None
    for init in [(0.2,-0.02,-3),(0,-0.01,-2),(0.3,0,-4)]:
        r=minimize(nll,init,method="Nelder-Mead",options={"maxiter":6000,"xatol":1e-8,"fatol":1e-8})
        if best is None or r.fun<best.fun: best=r
    eps=1e-4; H=np.zeros((3,3))
    import itertools as it
    for i,j in it.product(range(3),repeat=2):
        v=[np.array(best.x) for _ in range(4)]
        v[0][i]+=eps; v[0][j]+=eps; v[1][i]+=eps; v[1][j]-=eps; v[2][i]-=eps; v[2][j]+=eps; v[3][i]-=eps; v[3][j]-=eps
        H[i,j]=(nll(v[0])-nll(v[1])-nll(v[2])+nll(v[3]))/(4*eps*eps)
    se=np.sqrt(np.linalg.inv(H)[1,1]); z=best.x[1]/se
    return best.x[1],se,2*norm.sf(abs(z))
print("=== (a) 0-9세 비중 추세 — 베타이항 & Kendall, 창 민감도 ===")
for lo in (2001,2006,2008):
    sub=D[D.year>=lo]
    b,se,pbb=betabin_fit(sub)
    tk,pk=kendalltau(sub.year,sub.share)
    print(f"  {lo}-2024: 베타이항 slope={b:+.4f} p={pbb:.3f} · Kendall tau={tk:+.3f} p={pk:.3f}")

# ── (b) 인구구조 영모형: R_pooled 고정 + 관측 인구 → 기대 비중 (2008+) ──
D8=D[D.year>=2008].reset_index(drop=True)
ry=D8.cy.sum()/D8.Py.sum(); ro=D8.co.sum()/D8.Po.sum(); R=ry/ro
D8["exp_share"]=R*D8.Py/(R*D8.Py+D8.Po)*100
_po=(D8.cy+0.5)/(D8.N+1.0); _pe=(D8.exp_share/100).clip(1e-6,1-1e-6)   # Haldane 0.5 보정 + 클리핑(0/100% 방어)
D8["resid_logit"]=np.log(_po/(1-_po))-np.log(_pe/(1-_pe))
tr,pr=kendalltau(D8.year,D8.resid_logit)
te,pe=kendalltau(D8.year,D8.exp_share)
to,po=kendalltau(D8.year,D8.share)
print(f"\n=== (b) 인구구조 분해 (2008-2024, R_pooled={R:.1f}) ===")
print(D8[["year","N","share","exp_share","resid_logit"]].round(3).to_string(index=False))
print(f"  관측 비중 추세:  tau={to:+.3f} p={po:.3f}")
print(f"  기대 비중 추세(인구구조만): tau={te:+.3f} p={pe:.3f}  ← 인구만으로 생기는 감소")
print(f"  잔차(관측-기대) 추세: tau={tr:+.3f} p={pr:.3f}  ← 무추세면 '비중 감소=인구 효과' 입증")
D8.to_csv(os.path.join(RS,"MASTER_ehec_childshare_decomp.csv"),index=False)

# ── 그림 ──
plt.rcParams.update({"font.size":11,"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(12.5,4.6))
fig.suptitle("EHEC — is the falling child (0–9y) case share real risk change or demography?",fontsize=12.5,fontweight="bold",x=0.02,ha="left")
a=ax[0]
a.plot(D.year,D.share,"o-",color="#2c3e50",lw=2,ms=5,label="Observed child share")
a.plot(D8.year,D8.exp_share,"--",color="#c0392b",lw=1.8,label="Expected from demography alone\n(constant pooled rate ratio, observed populations)")
a.set_xlabel("Year"); a.set_ylabel("Child (0–9y) share of cases, %"); a.legend(frameon=False,fontsize=8.5)
a.set_title(f"A. Observed vs demography-only expectation\n(expected-share trend τ={te:+.2f}; observed τ={to:+.2f})",loc="left",fontsize=9.5)
b_=ax[1]
b_.axhline(0,ls=":",color="0.5")
b_.plot(D8.year,D8.resid_logit,"s-",color="#e67e22",lw=2,ms=5)
b_.set_xlabel("Year"); b_.set_ylabel("Residual (logit observed − logit expected)")
b_.set_title(f"B. Residual after removing demography\n(Kendall τ={tr:+.2f}, p={pr:.2f} → no residual trend)",loc="left",fontsize=9.5)
plt.tight_layout(rect=[0,0,1,0.9])
out=os.path.join(EH,"EHEC_Fig_childshare_decomp.png")
plt.savefig(out,dpi=150,bbox_inches="tight"); print("\n저장:",out)
