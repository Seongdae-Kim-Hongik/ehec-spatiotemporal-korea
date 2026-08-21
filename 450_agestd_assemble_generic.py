#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 450 — 연령표준화 자료조립(일반화). 인자: 질병한글명 영문명
#  EHEC 400을 파라미터화. KDCA 연령사례 + E_std(질병별 전국연령율×시군구인구) + 직접표준화 추세.
#  산출 data_cache/: ST_<eng>_agecase.csv · ST_<eng>_Estd.csv · ST_<eng>_natrate.csv · ST_<eng>_trend.csv
import pandas as pd, numpy as np, re, os, sys
DZ=sys.argv[1]; ENG=sys.argv[2]
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads")
ROOT=os.path.join(GD,"FBD_SpatioTemporal"); C=os.path.join(ROOT,"data_cache")
KDCA=os.path.join(GD,"건상/성별_종합 (3).csv"); FS=os.path.join(BASE,"식중독최종.csv")
AGE10=["0~9세","10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70세 이상"]
fs=pd.read_csv(FS,usecols=["sido","sigungu","region","year","population"])
REG229=sorted(fs.region.dropna().unique()); REGSET=set(REG229)
nat_totpop=fs.groupby(["year","region"],as_index=False)["population"].max().groupby("year")["population"].sum()
def p5to10(a):
    m={"0~4세":"0~9세","5~9세":"0~9세","10~14세":"10~19세","15~19세":"10~19세","20~24세":"20~29세","25~29세":"20~29세","30~34세":"30~39세","35~39세":"30~39세","40~44세":"40~49세","45~49세":"40~49세","50~54세":"50~59세","55~59세":"50~59세","60~64세":"60~69세","65~69세":"60~69세"}
    if a in m: return m[a]
    if a in ("70~74세","75~79세","80~84세","85~89세","90~94세","95~99세","100세 이상"): return "70세 이상"
    return None
# 인구(공유): 시군구×10세×연 (2008-2024) — 이미 400에서 만들었으면 재사용
POPF=os.path.join(C,"ST_agepop.csv")
if os.path.exists(POPF):
    pop=pd.read_csv(POPF)
else:
    SIDO_P={"서울특별시":"서울시","부산광역시":"부산시","대구광역시":"대구시","인천광역시":"인천시","광주광역시":"광주시","대전광역시":"대전시","울산광역시":"울산시","세종특별자치시":"세종시","경기도":"경기도","강원도":"강원도","강원특별자치도":"강원도","충청북도":"충청북도","충청남도":"충청남도","전라북도":"전라북도","전북특별자치도":"전라북도","전라남도":"전라남도","경상북도":"경상북도","경상남도":"경상남도","제주특별자치도":"제주도"}
    pp=pd.read_csv(os.path.join(GD,"인구데이터_통합.csv"),dtype={"인구수":str}).drop_duplicates(subset=["행정구역","날짜","구분","연령대"])
    pp=pp[(pp["구분"]=="계")&(pp["연령대"].map(lambda a:p5to10(a) is not None))]
    pp=pp[~pp["행정구역"].str.contains(r"\(3611",na=False)]
    def preg(hj):
        n=re.sub(r"\s*\(\d+\)\s*$","",str(hj)).strip();t=n.split();s=SIDO_P.get(t[0])
        if s is None:return None
        if s=="세종시":return "세종시세종시"
        r=s+"".join(t[1:]);return r if r in REGSET else None
    pp["region"]=pp["행정구역"].map(preg);pp=pp[pp.region.notna()].copy()
    pp["age10"]=pp["연령대"].map(p5to10);pp["인구수"]=pd.to_numeric(pp["인구수"].str.replace(",",""),errors="coerce")
    pp["yr"]=pp["날짜"].astype(str).str.extract(r"(\d{4})").astype(int);pp["mo"]=pp["날짜"].astype(str).str.extract(r"\d{4}[.\-]?(\d{2})").astype(float)
    pp=pp.groupby(["region","yr"],group_keys=False).apply(lambda g: g[g.mo==12] if (g.mo==12).any() else g[g.mo==g.mo.max()])
    pop=pp.groupby(["region","yr","age10"],as_index=False)["인구수"].sum().rename(columns={"yr":"year","인구수":"pop"})
    pop=pop[pop.year.between(2008,2024)];pop.to_csv(POPF,index=False)

# KDCA 연령사례(질병)
SIDO_K={"강원":"강원도","경기":"경기도","경남":"경상남도","경북":"경상북도","광주":"광주시","대구":"대구시","대전":"대전시","부산":"부산시","서울":"서울시","세종":"세종시","울산":"울산시","인천":"인천시","전남":"전라남도","전북":"전라북도","제주":"제주도","충남":"충청남도","충북":"충청북도"}
kd=pd.read_csv(KDCA,encoding="utf-8-sig");kd=kd[(kd.disease==DZ)&(kd.sex=="계")].copy()
def kreg(r):
    s=SIDO_K.get(str(r.sido).strip());sg=str(r.sigungu).strip()
    if s is None:return None
    if s=="세종시":return "세종시세종시"
    reg=s+sg
    if reg in REGSET:return reg
    reg2=reg.replace("미추홀구","남구")
    if reg2 in REGSET:return reg2
    if " " in sg:
        r3=s+sg.split()[0]
        if r3 in REGSET:return r3
    return None
kd["region"]=kd.apply(kreg,axis=1)
case=kd[kd.age_group.isin(AGE10)].dropna(subset=["region"])
case=case[case.region.isin(REGSET)].groupby(["region","year","age_group"],as_index=False)["cases"].sum()
case=case[case.year.between(2001,2024)];case.to_csv(os.path.join(C,f"ST_{ENG}_agecase.csv"),index=False)
n1124=int(case[case.year.between(2011,2024)].cases.sum())

# 2015 표준인구(공유) + 질병별 전국율(2011-24) + E_std
std2015=pop[pop.year==2015].groupby("age10")["pop"].sum().reindex(AGE10)
c_sp=case[case.year.between(2011,2024)].groupby("age_group")["cases"].sum().reindex(AGE10).fillna(0)
p_sp=pop[pop.year.between(2011,2024)].groupby("age10")["pop"].sum().reindex(AGE10)
natrate=(c_sp/p_sp);natrate.name="rate";natrate.to_csv(os.path.join(C,f"ST_{ENG}_natrate.csv"))
pw=pop.pivot_table(index=["region","year"],columns="age10",values="pop",fill_value=0).reindex(columns=AGE10,fill_value=0)
Estd=(pw*natrate.values).sum(axis=1).reset_index();Estd.columns=["region","year","E_std"]
Estd.to_csv(os.path.join(C,f"ST_{ENG}_Estd.csv"),index=False)

# 직접표준화 전국 시계열
tot_case=case.groupby("year")["cases"].sum();rows=[]
for y in range(2001,2025):
    O=int(tot_case.get(y,0));Ncr=nat_totpop.get(y,np.nan)
    crude=O/Ncr*1e5 if Ncr==Ncr and Ncr>0 else np.nan;astd=asl=asu=np.nan
    if 2008<=y<=2024:
        cy=case[case.year==y].groupby("age_group")["cases"].sum().reindex(AGE10,fill_value=0)
        py=pop[pop.year==y].groupby("age10")["pop"].sum().reindex(AGE10,fill_value=0)
        ry=(cy/py).values;w=std2015.values;W=w.sum()
        astd=np.nansum(ry*w)/W*1e5;var=np.nansum((w**2)*(cy.values)/(py.values**2));se=np.sqrt(var)/W*1e5
        asl=astd-1.96*se;asu=astd+1.96*se
    rows.append(dict(year=y,cases=O,crude=crude,agestd=astd,agestd_lo=asl,agestd_hi=asu))
pd.DataFrame(rows).to_csv(os.path.join(C,f"ST_{ENG}_trend.csv"),index=False)
Estd_sum=Estd[Estd.year.between(2011,2024)].E_std.sum()
print(f"[{ENG}] 2011-24 사례 {n1124} · E_std합 {Estd_sum:.0f} (SIR≈{n1124/max(Estd_sum,1):.3f}) · region {case[case.year.between(2011,2024)].region.nunique()}/229")
