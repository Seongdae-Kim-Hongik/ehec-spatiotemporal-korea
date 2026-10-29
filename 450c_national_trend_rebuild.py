# 450c — 국가 시계열(ST_<ENG>_trend.csv) 을 감염병포털 전국 연령별 완전집계로 재구축 (2026-09-09).
#  이유: 시군구×연령 교차표(ST_<ENG>_agecase)는 주소·연령 미상 사례 누락(이질 −497, 장티푸스 −90). 시군구 모형(E_std·natrate)은 건드리지 않음.
#  사용: python3 450c_national_trend_rebuild.py Shigellosis shigellosis | Typhoid typhoid
import pandas as pd,numpy as np,os,sys
ENG,dis=sys.argv[1],sys.argv[2]
GD=os.path.join(os.environ.get("EHEC_ROOT","."),"G.Downloads"); C=os.environ.get("FBD_CACHE",os.path.join(GD,"FBD_SpatioTemporal","data_cache")); C0=os.path.join(GD,"FBD_SpatioTemporal","data_cache")   # (2026-10-21) 검증용 출력 위치; 인구 = ST_agepop(v2)
BASE=os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP")
AGE10=["0~9세","10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70세 이상"]
n=pd.read_csv(os.path.join(C0,"KDCA_portal_national_age10_by_year_2001_2024.csv")); n=n[n.disease==dis].copy(); n["age_group"]=n.age.str.replace("70세이상","70세 이상",regex=False)
tot=n[n.age=="계"].set_index("year").cases; case=n[n.age_group.isin(AGE10)]
fs=pd.read_csv(os.path.join(BASE,"식중독최종.csv"),usecols=["region","year","population"]); nat_totpop=fs.groupby(["year","region"],as_index=False)["population"].max().groupby("year")["population"].sum()
pop=pd.read_csv(os.path.join(C,"ST_agepop.csv")); std2015=pd.read_csv(os.path.join(C0,"ST_stdpop2015.csv")).set_index("age10")["인구수"].reindex(AGE10)
old=pd.read_csv(os.path.join(C,f"ST_{ENG}_trend.csv")).set_index("year")
rows=[]
for y in range(2001,2025):
    O=int(tot.get(y,0)); Ncr=nat_totpop.get(y,np.nan); crude=O/Ncr*1e5 if Ncr==Ncr and Ncr>0 else np.nan; astd=asl=asu=np.nan
    if 2008<=y<=2024:
        cy=case[case.year==y].groupby("age_group")["cases"].sum().reindex(AGE10,fill_value=0); py=pop[pop.year==y].groupby("age10")["pop"].sum().reindex(AGE10,fill_value=0)
        ry=(cy/py).values; w=std2015.values; W=w.sum(); astd=np.nansum(ry*w)/W*1e5; var=np.nansum((w**2)*(cy.values)/(py.values**2)); se=np.sqrt(var)/W*1e5; asl=astd-1.96*se; asu=astd+1.96*se
    rows.append(dict(year=y,cases=O,crude=crude,agestd=astd,agestd_lo=asl,agestd_hi=asu))
new=pd.DataFrame(rows); _bk=os.path.join(C0,f"ST_{ENG}_trend_districtfile_backup.csv")
if not os.path.exists(_bk): old.to_csv(_bk)   # (2026-10-21) 백업은 처음 한 번만(다시 돌리면 원 백업을 덮어쓰던 결함 수정)
new.to_csv(os.path.join(C,f"ST_{ENG}_trend.csv"),index=False)
cmp=new.set_index("year")[["cases","crude","agestd"]].join(old[["cases","crude","agestd"]],rsuffix="_old"); cmp["dcases"]=cmp.cases-cmp.cases_old
print(f"[{ENG}] rebuilt ST_{ENG}_trend.csv | total {int(new.cases.sum())} (old {int(old.cases.sum())})"); print(cmp[cmp.dcases!=0].round(3).to_string())
