# -*- coding: utf-8 -*-
# 33_persistence_jkms_100k.py — (2026-10-21, /333 반영) 32_persistence_jkms.py와 같은 방식(SIR>1 연수 대 인구 보존 Poisson 영모형, 한쪽 P=(초과 수+1)/(B+1), BH FDR)에서
#  영모형 반복만 2,000 → 100,000회로 늘려 기존 E(540 규칙)와 v2 E 모두 다시 계산. 2,000회에서는 P 최솟값이 1/2001이라 BH q가 반복 수 해상도에 묶였음.
#  원고 시드 20260913 + 안정성 시드 4개. 같은 규칙(FDR<0.05)으로 두 E의 집합을 비교. 쓰기는 _v2_refit_261021/JKMS/ 에만.
#  사용: python3 33_persistence_jkms_100k.py [B=100000]
import pandas as pd, numpy as np, os, sys, warnings, time
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
GATE=os.environ.get("EHEC_GATE", os.path.join(os.environ.get("EHEC_ROOT","."),"gate")); sys.path.insert(0,GATE); import sigungu_key as sk
ROOT=os.environ.get("EHEC_ROOT","."); C=os.path.join(ROOT,"data_cache"); OUT=os.path.join(ROOT,"_v2_rerun_261022"); IN=os.path.join(ROOT,"_v2_refit_261021")
YY=list(range(2011,2025)); B=int(sys.argv[1]) if len(sys.argv)>1 else 100000; CH=5000
SEEDS=[20260913,1,2,3,4]
KEY=pd.read_csv(os.path.join(GATE,"sigungu_key.csv"),dtype={"sigungu_cd":str}); NM=dict(zip(KEY.sigungu_cd,KEY.canonical_region)); CODES=KEY.sigungu_cd.tolist()
def coded(df):
    d=df.copy(); d["region"]=d.region.astype(str).str.replace(r"\s+","",regex=True); d=sk.attach_code(d,"region",None); d["sigungu_cd"]=d.sigungu_cd.astype(str); return d
cs=coded(pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv")).groupby(["region","year"],as_index=False).cases.sum())
E0=coded(pd.read_csv(os.path.join(C,"ST_Estd_panel.csv")))
grid=pd.MultiIndex.from_product([CODES,YY],names=["sigungu_cd","year"]).to_frame(index=False)
D=grid.merge(cs[["sigungu_cd","year","cases"]],how="left").merge(E0[["sigungu_cd","year","E_std"]].rename(columns={"E_std":"E_old"}),how="left")
D["cases"]=D.cases.fillna(0); mn=D.E_old[D.E_old>0].min(); D["E_old"]=D.E_old.where(D.E_old>0).fillna(mn)   # 540 규칙 재현(32와 같음)
v2=pd.read_csv(os.path.join(IN,"inputs","E_cases_EHEC_v2.csv"),dtype={"sigungu_cd":str}); D=D.merge(v2[["sigungu_cd","year","E_std"]].rename(columns={"E_std":"E_v2"}),how="left"); assert D.E_v2.notna().all()
O=D.pivot(index="sigungu_cd",columns="year",values="cases").reindex(CODES)
def persist(Ecol,seed):
    E=D.pivot(index="sigungu_cd",columns="year",values=Ecol).reindex(CODES)
    nat=O.sum()/E.sum(); obs=((O>=1)&(O/E>1)).sum(axis=1).values
    rng=np.random.default_rng(seed); mu=(E*nat).values; Ev=E.values
    ge=np.zeros(len(CODES),dtype=np.int64); tot=np.zeros(len(CODES))
    done=0
    while done<B:   # 32와 같은 분포: 칸마다 독립 Poisson(E×전국 연도 비율). 메모리 때문에 CH회씩 나눠 뽑음
        b=min(CH,B-done); sims=rng.poisson(mu[None,:,:],size=(b,)+mu.shape)
        cnt=((sims>=1)&(sims/Ev[None,:,:]>1)).sum(axis=2)
        ge+=(cnt>=obs[None,:]).sum(axis=0); tot+=cnt.sum(axis=0); done+=b
    p=(ge+1)/(B+1); q=multipletests(p,method="fdr_bh")[1]
    return pd.DataFrame({"E":Ecol,"seed":seed,"B":B,"sigungu_cd":CODES,"region":[NM[c] for c in CODES],"obs_years":obs,"exp_null":(tot/B).round(2),"n_ge":ge,"P":p,"FDR":q})
rows=[]; t0=time.time()
for seed in SEEDS:
    for ec in ["E_old","E_v2"]:
        r=persist(ec,seed); rows.append(r)
        s=r[r.FDR<0.05].sort_values("P")
        print(f"[seed {seed} · {ec} · B={B}] FDR<0.05 {len(s)}곳: "+", ".join(f"{a}({o}/{e:.1f}, P {p:.5f}, q {q:.4f})" for a,o,e,p,q in zip(s.region,s.obs_years,s.exp_null,s.P,s.FDR))+f"  [{time.time()-t0:.0f}s]",flush=True)
R=pd.concat(rows)
R.to_csv(os.path.join(OUT,"JKMS",f"JKMS_persistence_B{B//1000}k_allseeds.csv"),index=False,encoding="utf-8-sig")
# 원고 시드 결과(S23 형식)
m=R[R.seed==20260913].copy(); m=m.sort_values(["E","P","obs_years"],ascending=[True,True,False])
m.to_csv(os.path.join(OUT,"JKMS",f"JKMS_persistence_S23_B{B//1000}k_seed20260913.csv"),index=False,encoding="utf-8-sig")
# 시드 안정성: 시드별 FDR<0.05 비율, P·FDR 범위
st=R.groupby(["E","sigungu_cd","region"]).agg(sig_share=("FDR",lambda x:(x<0.05).mean()),obs=("obs_years","first"),exp=("exp_null","mean"),P_min=("P","min"),P_max=("P","max"),FDR_min=("FDR","min"),FDR_max=("FDR","max")).reset_index()
st=st[st.FDR_min<0.2].sort_values(["E","FDR_max"])
st.to_csv(os.path.join(OUT,"JKMS",f"JKMS_persistence_seed_stability_B{B//1000}k.csv"),index=False,encoding="utf-8-sig")
print(st.to_string(index=False))
