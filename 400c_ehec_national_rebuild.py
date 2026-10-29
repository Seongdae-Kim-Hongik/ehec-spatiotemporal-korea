# 400c — EHEC 국가 시계열·연령 집계를 감염병포털 전국 완전집계로 재구축 (2026-09-09). 시군구 파일(ST_EHEC_agecase.csv)·E_std는 유지.
import pandas as pd,numpy as np,os
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads"); C=os.environ.get("FBD_CACHE",os.path.join(GD,"FBD_SpatioTemporal","data_cache"))   # (2026-10-21) 검증용 출력 위치 지정; 인구 = ST_agepop(v2, 코드 결합)
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
AGE10=["0~9세","10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70세 이상"]
n=pd.read_csv(os.path.join(GD,"FBD_SpatioTemporal","data_cache","KDCA_portal_national_age10_by_year_2001_2024.csv")); n=n[n.disease=="ehec"].copy(); n["age_group"]=n.age.str.replace("70세이상","70세 이상",regex=False)
tot=n[n.age=="계"].set_index("year").cases; case=n[n.age_group.isin(AGE10)][["year","age_group","cases"]].copy(); case.insert(0,"region","전국")
case.to_csv(os.path.join(C,"ST_EHEC_agecase_national.csv"),index=False)
fs=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["region","year","population"]); natpop=fs.groupby(["year","region"],as_index=False).population.max().groupby("year").population.sum()
pop=pd.read_csv(os.path.join(C,"ST_agepop.csv")); std=pd.read_csv(os.path.join(GD,"FBD_SpatioTemporal","data_cache","ST_stdpop2015.csv")).set_index("age10")["인구수"].reindex(AGE10)
oldp=os.path.join(C,"ST_nat_agestd_trend.csv"); old=pd.read_csv(oldp if os.path.exists(oldp) else os.path.join(GD,"FBD_SpatioTemporal","data_cache","ST_nat_agestd_trend.csv"))
if not os.path.exists(oldp.replace(".csv","_districtfile_backup.csv")): old.to_csv(oldp.replace(".csv","_districtfile_backup.csv"),index=False)
rows=[]
for y in range(2001,2025):
    O=int(tot[y]); crude=O/natpop[y]*1e5; astd=asl=asu=np.nan
    if y>=2008:
        cy=case[case.year==y].groupby("age_group").cases.sum().reindex(AGE10,fill_value=0); py=pop[pop.year==y].groupby("age10")["pop"].sum().reindex(AGE10)
        ry=(cy/py).values; w=std.values; W=w.sum(); astd=np.nansum(ry*w)/W*1e5; se=np.sqrt(np.nansum((w**2)*cy.values/(py.values**2)))/W*1e5; asl,asu=astd-1.96*se,astd+1.96*se
    rows.append(dict(year=y,cases=O,crude=crude,agestd=astd,agestd_lo=asl,agestd_hi=asu))
new=pd.DataFrame(rows); new.to_csv(oldp,index=False); print("ST_nat_agestd_trend.csv rebuilt: N",int(new.cases.sum()),"| 2008 agestd %.3f crude %.3f | 2024 agestd %.3f crude %.3f | fold agestd %.2f crude %.2f"%(new.agestd[7],new.crude[7],new.agestd[23],new.crude[23],new.agestd[23]/new.agestd[7],new.crude[23]/new.crude[7]))
