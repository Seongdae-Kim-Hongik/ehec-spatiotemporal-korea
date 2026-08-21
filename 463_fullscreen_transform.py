#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 463 — EHEC 전수 음성스크린 v2: 프로덕션 공간코드(fbd-korea-tongcode)의 변환메뉴를 그대로 적용.
#  이전 461은 per-100k raw + Spearman(횡단면·변환 없음)이었음.  본 스크립트는:
#   (1) 2011-2023 패널(region×year) — "2011부터 있는 변수만"(의료시설 2022-24 제외).
#   (2) 각 지표 → 5형태 {raw, log1p, binary, T3(tertile), Q4(quartile)} 생성
#       (tongcode 로직: %변수는 log 제외; zero>20%면 binary=有無·T3=0/저/고).
#   (3) 각 형태 단변량 NB GLM(cases ~ z(form) + offset(log E_std)) → 최적형태(최소 p) 선택.
#   (4) best-form 표준화 log-IRR(±)로 랭킹 → 그림/CSV.  누수(발생/장출혈성/감염증) 제거.
#  ★ 목적: "변환을 제대로 했나?"에 답 — 모든 변수에 최선의 변환기회를 줘도(=더 공격적 낚시)
#     생물학적 타당 결정요인이 안 뜨면 negative-control 결론이 더 강해진다.
import os, numpy as np, pandas as pd, warnings
from multiprocessing import Pool
import statsmodels.api as sm
from statsmodels.genmod.families import Poisson
warnings.filterwarnings("ignore")
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); C=os.path.join(ROOT,"data_cache"); EH=os.path.join(ROOT,"EHEC"); R=os.path.join(ROOT,"results_ST")
Y0,Y1=2011,2023

# ===== 패널 조립 =====
case=pd.read_csv(os.path.join(C,"ST_EHEC_agecase.csv")).groupby(["region","year"])["cases"].sum().reset_index()
Est=pd.read_csv(os.path.join(C,"ST_Estd_panel.csv"))
pop=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["region","year","population"]).drop_duplicates(["region","year"])
H=pd.read_csv(os.path.join(BASE,"국민건강결과_최종.csv"))
idc=[c for c in H.columns if c in ("year","sido","sigungu","region")]
ind0=[c for c in H.columns if c not in idc]
# 누수 제거: 감염병 발생/co-surveillance + EHEC 자신
ind=[c for c in ind0 if ("발생" not in c) and ("장출혈성" not in c) and ("감염증" not in c) and ("식중독" not in c)]
print(f"지표 {len(ind0)} → 누수제거 후 {len(ind)}  |  패널 {Y0}-{Y1} (의료시설 2022-24 제외: 2011부터 없음)")

# 카운트형 → per-100k;  이미 비율/율/지수/평균 등은 그대로(rate).
# ★/333 수정: 구 휴리스틱의 "세"(→세대수·인구수_20-44세)·"점수"(→점포수) 오매칭 제거(428개 카운트 오분류).
def is_pct(nm): return nm.endswith(("율","률")) or any(t in nm for t in ("백분율","비율","%"))
def is_rate(nm): return is_pct(nm) or any(t in nm for t in ("십만","10만","천명당","인구당","기대여명","지수","평균","표준화율","조율"))
base=(pd.merge(case,Est,on=["region","year"],how="outer")
      .merge(pop,on=["region","year"],how="left"))
base=base[base.year.between(Y0,Y1)].copy()
base["cases"]=base["cases"].fillna(0)
base=base[(base.E_std>0)&(base.population>0)].copy()
H=H[H.year.between(Y0,Y1)].copy()
# ★지역명 정규화 — 미적용 시 인천 미추홀구/남구가 조인에서 탈락(1개 시군구 손실)
for _df in (base,H):
    _df["region"]=_df["region"].astype(str).str.replace(r"\s+","",regex=True)
    _df["region"]=_df["region"].replace({"인천시미추홀구":"인천시남구"})
d=base.merge(H[idc+ind],on=["region","year"],how="inner")
d=d[np.isfinite(d.E_std)&(d.E_std>0)].reset_index(drop=True)
logE=np.log(d["E_std"].values); y=d["cases"].values.astype(float)
pv=d["population"].values
print(f"모델링 프레임: {len(d)} region-year cells  (지역 {d.region.nunique()} × 연도 {d.year.nunique()})")

# ===== 5형태 생성 (tongcode 로직) =====
# ★/333 수정: 결측(NaN) 반드시 보존 — binary NaN→0, digitize NaN→최상위 분위(3/4) 오코딩 버그 제거.
def _digit(xr,brk):  # finite 위치에만 digitize, 결측은 NaN 유지
    m=np.isfinite(xr); o=np.full_like(xr,np.nan,dtype=float)
    o[m]=np.digitize(xr[m],brk[1:-1],right=True).astype(float)+1; return o
def make_forms(nm,x):
    x=np.asarray(x,float); out={}
    xr = x/pv*1e5 if not is_rate(nm) else x    # 카운트→per-100k
    if np.nansum(np.isfinite(xr))<50 or np.nanstd(xr[np.isfinite(xr)])==0: return out
    fin=xr[np.isfinite(xr)]
    # ★/333: 퇴화/희소 변수 가드 — 비제로<3% 또는 distinct<10이면 스크린 무의미(이상치 지배 허구값 방지, tongcode MIN_OBS/coverage 취지).
    if (np.mean(fin>0)<0.03) or (len(np.unique(np.round(fin,6)))<10): return out
    out["raw"]=xr
    zp=np.mean(fin==0)*100; hz=zp>20; pt=is_pct(nm)
    if not pt:
        lv=np.log1p(np.clip(xr,0,None))  # np.log1p(nan)=nan → 결측 보존
        if np.nanstd(lv[np.isfinite(lv)])>0: out["log1p"]=lv
    b=np.full_like(xr,np.nan,dtype=float); m=np.isfinite(xr)
    thr=0.0 if hz else np.nanmedian(xr); b[m]=(xr[m]>thr).astype(float); out["binary"]=b
    if hz:
        nz=fin[fin>0]
        if len(nz)>10:
            mn=np.median(nz); out["T3"]=np.where(~np.isfinite(xr),np.nan,np.where(xr==0,1,np.where(xr<=mn,2,3)))
    else:
        q=np.nanquantile(xr,[1/3,2/3]); brk=np.unique(np.r_[-np.inf,q,np.inf])
        if len(brk)>=3: out["T3"]=_digit(xr,brk)
    q4=np.nanquantile(xr,[.25,.5,.75]); b4=np.unique(np.r_[-np.inf,q4,np.inf])
    if len(b4)>=3: out["Q4"]=_digit(xr,b4)
    return out

def zst(v):
    v=np.asarray(v,float); m=np.isfinite(v)
    if m.sum()<50: return None
    mu,sd=np.nanmean(v[m]),np.nanstd(v[m])
    if sd==0: return None
    z=np.full_like(v,np.nan); z[m]=(v[m]-mu)/sd; return z

def univ(z):
    # ★/333 수정: alpha=1 고정 NB는 점추정까지 alpha에 의존 → Poisson(계수 alpha-불변)로 랭킹.
    # p는 공간효과 미보정으로 무효(nominal) — 형태선택·랭킹 보조용일 뿐 유의판정 아님.
    m=np.isfinite(z)&np.isfinite(y)&np.isfinite(logE)
    if m.sum()<int(0.5*len(z)): return None
    X=sm.add_constant(z[m])
    try:
        r=sm.GLM(y[m],X,family=Poisson(),offset=logE[m]).fit()
        return dict(coef=r.params[1],p=r.pvalues[1],n=int(m.sum()))
    except Exception: return None

def work(nm):
    x=d[nm].values
    rate=is_rate(nm); forms=make_forms(nm,x)
    best=None
    for fn,fv in forms.items():
        z=zst(fv)
        if z is None: continue
        rr=univ(z)
        if rr is None: continue
        if (best is None) or (rr["p"]<best["p"]):
            best=dict(form=fn,coef=rr["coef"],p=rr["p"],n=rr["n"])
    if best is None: return None
    return dict(var=nm,is_rate=rate,coverage=round(best["n"]/len(x),3),**best)

if __name__=="__main__":
    with Pool(10) as pl:
        res=[r for r in pl.map(work,ind) if r is not None]
    M=pd.DataFrame(res).rename(columns={"p":"p_nominal"})
    # ★/333: p는 공간효과 미보정+best-of-5 형태선택으로 무효 → 'nominal'만. BH도 변수내 선택편향 미보정이라 참고용.
    from statsmodels.stats.multitest import multipletests
    M=M.sort_values("p_nominal").reset_index(drop=True)
    M["q_BH_nominal"]=multipletests(M["p_nominal"],method="fdr_bh")[1]
    M["logIRR_perSD"]=M["coef"]; M["abseff"]=M["coef"].abs()
    M=M.sort_values("abseff",ascending=False).reset_index(drop=True)
    M=M[["var","form","is_rate","coverage","n","logIRR_perSD","p_nominal","q_BH_nominal","abseff"]]
    M.to_csv(os.path.join(R,"MASTER_ehec_fullscreen_v2.csv"),index=False)
    print(f"\n=== 스크린 완료: {len(M)} 지표 (각 best-of-5형태, Poisson offset=logE_std) ===")
    print("형태 분포:",M["form"].value_counts().to_dict(),"| per-100k 변환(is_rate=False):",int((~M['is_rate']).sum()))
    print("\n상위 20 (|표준화 log-IRR|; p는 nominal·유의판정 아님):")
    for _,r in M.head(20).iterrows():
        print(f"  {r['logIRR_perSD']:+.3f} [{r['form']:>6}] p_nom={r['p_nominal']:.1e}  {r['var'][:46]}")

    # ===== 그림 (461 스타일 유지 + 변환명 주석) =====
    import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    KF=FontProperties(fname="/System/Library/Fonts/AppleSDGothicNeo.ttc") if os.path.exists("/System/Library/Fonts/AppleSDGothicNeo.ttc") else None
    mpl.rcParams.update({"font.size":9,"axes.edgecolor":"#444","axes.linewidth":0.8})
    top=M.head(15).iloc[::-1]
    fig,ax=plt.subplots(figsize=(11.6,6.4))
    cols=["#c1666b" if v>0 else "#6b8fc1" for v in top["logIRR_perSD"]]
    ax.barh(range(len(top)),top["logIRR_perSD"],color=cols,edgecolor="none")
    ax.set_yticks(range(len(top)))
    # ── y-라벨 영문화(하우스스타일: 영문 문서) ──
    def en_label(v):
        s=str(v)
        rep=[("가해","Assault"),("사망자수","deaths"),("표준화사망률","std. mortality"),
             ("1인가구율","Single-person household rate"),("내분비영양및대사질환","Endocrine/metabolic disease"),
             ("특수진료실","Special-care"),("응급실병실수","emergency-room beds"),("응급실병상수","emergency-room beds"),
             ("물리치료실병상수","physical-therapy beds"),("노후주택","Old housing"),("단독주택","detached"),
             ("달리분류되지않은증상징후","Ill-defined symptoms/signs"),("전염병","Notifiable-disease"),("진료실인원","clinic patients"),
             ("B형간염(HepB)예방접종률","HepB vaccination rate"),("결핵(BCG)예방접종률","BCG vaccination rate"),
             ("구강보건","Oral health"),("치면세마","dental scaling"),("수혜자수","recipients"),("실적","performance"),
             ("무면허교통사고","Unlicensed traffic accident"),("부상자수","injuries"),("보건소인력","Public-health-centre staff"),
             ("노인장기요양시설의사","LTC-facility physicians"),("노인장기요양","LTC"),("기관수","facilities"),
             ("위암결과","Gastric-cancer screening"),("간암결과","Liver-cancer screening"),("폐암결과","Lung-cancer screening"),
             ("자궁경부암결과","Cervical-cancer screening"),("대장내시경","colonoscopy"),("추적검사요망","follow-up advised"),
             ("결절선결절","nodular finding"),("기존암환자","prior-cancer patient"),("일상생활수행능력","ADL"),
             ("대중교통이용","public-transport use"),("유방의악성신생물","breast malignancy"),
             ("만1세","age 1"),("만2세","age 2"),("만3세","age 3"),("만6세","age 6"),("66세","age 66"),
             ("45-64세가구","(45–64y)"),("65세이상가구","(65+y)"),("전체","total"),("남자","male"),("여자","female"),
             ("경증","mild"),("기타","other")]
        for k,e in rep: s=s.replace(k,e)
        s=s.replace("_"," ").strip()
        return s if not any('가'<=c<='힣' for c in s) else str(v)  # 남은 한글 있으면 원문 유지(안전)
    labs=[en_label(v) for v in top["var"]]
    kor_left=sum(1 for l in labs if any('가'<=c<='힣' for c in l))
    if kor_left: print(f"⚠ 미번역 y-라벨 {kor_left}개(원문 유지)")
    ax.set_yticklabels(labs,fontsize=8)
    for i,(v,fm) in enumerate(zip(top["logIRR_perSD"],top["form"])):  # ★/333: 유의성 별표 제거(공간효과 미보정)
        ha="left" if v>0 else "right"; off=0.004 if v>0 else -0.004
        ax.text(v+off,i,f"{v:+.2f} [{fm}]",va="center",ha=ha,fontsize=7,color="#333")
    ax.axvline(0,color="#333",lw=0.8)
    ax.set_xlabel("Best-transform standardised log-IRR per SD (unadjusted Poisson, offset=log E_std; 2011–2023 panel; NOT spatially adjusted — ranking only)",fontsize=8.3)
    ax.set_title(f"Exploratory high-dimensional negative screen for omitted health/healthcare determinants ({len(M):,} indicators, 2011–2023 panel)\n"
                 "each indicator given its best of 5 transforms {raw, log, binary, tertile, quartile};\n"
                 "top signals are NON-SPECIFIC access / screening / urbanicity / population proxies, not EHEC-specific exposures",
                 fontsize=8.5,loc="left")
    mx=top["logIRR_perSD"].abs().max()*1.35; ax.set_xlim(-mx,mx)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    fig.text(0.01,0.008,"Ranking only — magnitudes and nominal significance are inflated by unmodelled spatial autocorrelation and best-of-5 (winner's-curse) selection, not causal.\n"
             "No EHEC-specific driver (livestock, agriculture, water, food handling) rises above these non-specific proxies. Definitive test = spatially-adjusted BYM2 (Table 1). Full method in Fig. S16 caption.",
             fontsize=6.2,color="#888")
    plt.tight_layout(rect=[0,0.03,1,1])
    plt.savefig(os.path.join(EH,"EHEC_Fig_fullscreen.png"),dpi=300,facecolor="white"); plt.close()
    print("\n저장: MASTER_ehec_fullscreen_v2.csv + EHEC_Fig_fullscreen.png (변환메뉴·2011+ 패널)")
