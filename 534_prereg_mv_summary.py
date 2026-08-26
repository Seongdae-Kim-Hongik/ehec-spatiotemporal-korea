# -*- coding: utf-8 -*-
# 534 — 사전지정 멀티버스 합산·판정. 복합키(seed×spec). 판정 규칙은 PREREGISTERED_SPECIFICATION에 고정된 그대로.
import pandas as pd, numpy as np, glob, re, json, os, sys
RS="results_ST"; C="data_cache"
def load(pat,idcol):
    out=[]
    for f in sorted(glob.glob(os.path.join(RS,pat))):
        seed=int(re.search(r"seed(\d+)",f).group(1)); d=pd.read_csv(f); d["seed"]=seed; d["uid"]=seed*10000+d[idcol]; out.append(d)
    return pd.concat(out,ignore_index=True)
R=load("PREREG_MV_results_seed1[01][0-9].csv","spec"); S=load("PREREG_MV_specs_seed1[01][0-9].csv","id")
n0=len(R)
# 같은 (사양,변수)가 두 번 기록된 경우 = 같은 모형의 반복 실행(DIC 동일, IRR 소수4자리 반올림 차이). 첫 기록 채택.
# 원인: 동일 seed 워커가 이중 기동되어 같은 사양을 두 번 적합(변수집합 동일, INLA 수렴점 미세 차이). 사양당 첫 적합(첫 DIC) 채택.
first_dic=R.groupby("uid").DIC.transform("first"); R=R[R.DIC==first_dic]
R=R.drop_duplicates(subset=["uid","term"],keep="first").reset_index(drop=True)
assert R.groupby("uid").DIC.nunique().max()==1
print(f"이중 적합 정리: {n0}→{len(R)}행 · 사양당 첫 적합 채택 · 최종 사양 {R.uid.nunique():,}")
DES=json.load(open(os.path.join(C,"PREREG_DESIGN.json"))); M=pd.read_csv(os.path.join(C,"ST_prereg_panel_META.csv")).set_index("key")
dom={v:d for d,x in DES["domains"].items() for v in x["vars"]}; dirn={v:x["direction"][v] for d,x in DES["domains"].items() for v in x["vars"]}
dom.update({"rum_TMid":"D1 reservoir","rum_THigh":"D1 reservoir","kchs_pc1":"axis","sex_ratio":"forced"}); dirn.update({"rum_TMid":"+","rum_THigh":"+","kchs_pc1":"±","sex_ratio":"±"})
R["log"]=np.log(R.IRR)
g=R.groupby("term").agg(n_spec=("uid","nunique"),n_cred=("cred","sum"),median_IRR=("IRR","median"),q25=("IRR",lambda x:x.quantile(.25)),q75=("IRR",lambda x:x.quantile(.75)),pct_above1=("IRR",lambda x:(x>1).mean()*100)).reset_index()
g["cred_pct"]=g.n_cred/g.n_spec*100
# 방향 일관: credible 사양들 중 같은 방향 비율
cc=R[R.cred].groupby("term").IRR.apply(lambda x: max((x>1).mean(),(x<1).mean())*100).rename("cred_dir_consistency")
g=g.merge(cc,on="term",how="left")
def verdict(r):
    if r.n_spec<20: return "insufficient (n<20)"
    if r.cred_pct>=50 and (r.cred_dir_consistency or 0)>=90: return "robust"
    if r.cred_pct>=50: return "robust, direction unstable"
    if r.cred_pct>=20: return "suggestive"
    return "null"
g["verdict"]=g.apply(verdict,axis=1); g["domain"]=g.term.map(dom); g["hypothesis"]=g.term.map(dirn)
g["observed_dir"]=np.where(g.median_IRR>1,"+","-"); g["dir_match"]=g.apply(lambda r: "" if r.hypothesis in("±","0") else ("✓" if r.hypothesis==r.observed_dir else "✗ opposite"),axis=1)
g=g.sort_values(["cred_pct","n_spec"],ascending=False)
g.to_csv(os.path.join(RS,"PREREG_MV_summary.csv"),index=False,encoding="utf-8-sig")
print(f"사양 {R.uid.nunique():,} · 변수 {len(g)} · 판정: {g.verdict.value_counts().to_dict()}\n")
print(f"{'변수':40s} {'도메인':14s} {'n':>5s} {'cred%':>6s} {'IRR중앙':>8s} {'IQR':>14s} {'방향일관':>8s}  가설 판정")
for _,r in g.iterrows():
    if r.cred_pct>=20 or r.term in ("soap_sanitizer_use","beef_cattle_density","rum_THigh","fiscal_independence","grdp_per_capita","basic_livelihood_recipient_rate","drinking_water_groundwater_use","total_living_use_groundwater_volume","bars_pubs_per_1000_pop"):
        print(f"{r.term[:40]:40s} {str(r.domain)[:14]:14s} {r.n_spec:5d} {r.cred_pct:6.1f} {r.median_IRR:8.3f} {r.q25:6.3f}-{r.q75:6.3f} {r.cred_dir_consistency if r.cred_dir_consistency==r.cred_dir_consistency else 0:8.0f}  {r.hypothesis:2s} {r.dir_match:10s} {r.verdict}")
# 음성대조 보고
print("\n── 음성대조 ──")
for v in DES["negative_controls"]:
    r=g[g.term==v]
    if len(r): r=r.iloc[0]; print(f"  {v:36s} cred {r.cred_pct:5.1f}%  IRR {r.median_IRR:.3f}  {'✓ null' if r.cred_pct<20 else '✗ 교란 의심'}")
print("\n저장: results_ST/PREREG_MV_summary.csv")
