#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 545 — 농림어업총조사(2010·2015·2020) 시군구 한우·젖소 사육농가 → 229 프레임 패널.
#  목적: 교수 규칙("2011년부터 있는 자료") 충족 검정용. 2011–2024 = 총조사 3시점 선형보간(2021+는 2020 유지).
import pandas as pd, numpy as np, os, re, sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import v2_fixes as v2
# ★2026-10-21 v2: 밀도 분모 = final.shp 기하 면적(코드 결합). 옛 판은 그때의 ST_prereg_panel area(연도별 하수도 면적 평균, 2017+ 이중 합산)를 써서
#   548이 넣은 한육우·젖소 농가 밀도 × shp 면적이 정수(농가 수)가 아니었음. 출력 위치 = FBD_CACHE, 총조사 원자료는 data_cache에서 읽음.
CIN="data_cache"; C=v2.cache_dir()
SIDO={"11":"서울시","21":"부산시","22":"대구시","23":"인천시","24":"광주시","25":"대전시","26":"울산시","29":"세종시","31":"경기도","32":"강원도","33":"충청북도","34":"충청남도","35":"전라북도","36":"전라남도","37":"경상북도","38":"경상남도","39":"제주특별자치도"}
ALIAS={"경상북도군위군":"대구시군위군","제주특별자치도서귀포시":"제주도서귀포시","제주특별자치도제주시":"제주도제주시","인천시미추홀구":"인천시남구","충청북도청원군":"충청북도청주시","충청남도연기군":"세종시세종시","세종시세종특별자치시":"세종시세종시","경상남도마산시":"경상남도창원시","경상남도진해시":"경상남도창원시","경기도여주군":"경기도여주시","충청남도당진군":"충청남도당진시","인천시남구":"인천시남구"}
regs=sorted(pd.read_csv(os.path.join(C,"ST_prereg_panel.csv"),usecols=["region"]).region.unique())
area=v2.area_by_region(regs)   # ★2026-10-21 shp 면적
CITY_GU=re.compile(r".+시(중구|동구|남구|서구|북구|.{1,3}구)$")  # 도 지역 시 산하 일반구 걸러내기용은 코드 계층으로 대체
rows=[]
for y in (2010,2015,2020):
    d=pd.read_csv(os.path.join(CIN,f"census{y}_utf8.csv"),skiprows=1,dtype={"C행정구역별":str})
    d["code"]=d["C행정구역별"].str.strip("'"); d["name"]=d["행정구역별"].str.strip()
    nat=pd.to_numeric(d.loc[d["code"]=="00","한우:농가 (가구)"],errors="coerce").iloc[0]
    d=d[~d["name"].isin(["동부","읍부","면부","전국"])]
    sg=d[d["code"].str.len()==5].copy()
    # 도(31~39) 산하 '일반구'(시 산하 구) 제거: 같은 시도 안에서 앞선 '~시' 행 다음에 오는 '~구' 행.
    keep=[]
    for _,r in sg.iterrows():
        s2=r["code"][:2]
        if s2 in ("11","21","22","23","24","25","26","29"): keep.append(True)          # 광역: 자치구·군 유지
        else: keep.append(not r["name"].endswith("구"))                                  # 도: 일반구 제외(시·군만)
    sg=sg[pd.Series(keep,index=sg.index)]
    sg["key"]=(sg["code"].str[:2].map(SIDO)+sg["name"]).map(lambda k:ALIAS.get(k,k))
    for c in ("한우:농가 (가구)","젖소암컷:농가 (가구)"): sg[c]=pd.to_numeric(sg[c],errors="coerce")
    g=sg.groupby("key").agg(beef_farms=("한우:농가 (가구)",lambda x:x.sum(min_count=1)),dairy_farms=("젖소암컷:농가 (가구)",lambda x:x.sum(min_count=1))).reset_index()
    assert abs(g[g.key.isin(regs)].beef_farms.sum()-nat)<1e-6, f"{y} 전국합 불일치"   # 이중집계·누락 가드
    g["year"]=y; rows.append(g)
X=pd.concat(rows,ignore_index=True)
miss=[r for r in regs if r not in set(X.key)]
extra=sorted(set(X.key)-set(regs))
print("프레임 미커버:",miss); print("잉여 키(무시):",extra[:12],"…" if len(extra)>12 else "")
X=X[X.key.isin(regs)]
# 밀도(호/km²) → 2011–2024 보간
out=[]
for r in regs:
    a=area[r]; sub=X[X.key==r].set_index("year").reindex([2010,2015,2020])
    sub=sub.fillna(0)  # 센서스는 제로셀을 빈칸으로 생략(전국합 일치로 검증됨) → 빈셀=해당 축종 농가 0. 매칭 실패는 위 miss/전국합 assert가 잡음.
    s=pd.Series(sub.beef_farms.values,index=[2010,2015,2020]).reindex(range(2010,2025)).interpolate().ffill()
    sd=pd.Series(sub.dairy_farms.values,index=[2010,2015,2020]).reindex(range(2010,2025)).interpolate().ffill()
    for y in range(2011,2025): out.append({"region":r,"year":y,"census_beef_farms_km2":s[y]/a,"census_dairy_farms_km2":sd[y]/a})
P=pd.DataFrame(out); P.to_csv(os.path.join(C,"ST_census_livestock.csv"),index=False,encoding="utf-8-sig")
# 검증: 기존 시간불변 농가밀도(2019-23 평균)와의 상관
old=pd.read_csv(os.path.join(C,"ST_prereg_panel.csv"),usecols=["region","year","beef_cattle_farms"])
m=P[P.year==2020].merge(old[old.year==2020],on="region")
print(f"229 커버 완료 · 2020 총조사 vs 기존(가축동향 2019-23) 상관 r={m.census_beef_farms_km2.corr(m.beef_cattle_farms):.3f}")
print(P.groupby("year").census_beef_farms_km2.mean().round(3).to_dict())
