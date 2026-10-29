#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 548 — 본모형 교체: ST_prereg_panel의 beef_cattle_farms·dairy_farms를 농림어업총조사 시변 밀도(호/km²)로 교체.
#  (교수 규칙 '2011년부터 있는 자료' 충족. 구 가축동향 값은 *_trend2019 컬럼으로 보존 → 민감도용.)
import pandas as pd, os, shutil
# (2026-10-21) 출력 위치 = FBD_CACHE(기본 data_cache; 검증용)
C=os.environ.get("FBD_CACHE","data_cache"); P=os.path.join(C,"ST_prereg_panel.csv")
bak=os.path.join(C,"ST_prereg_panel_bak_trendcattle.csv")
if not os.path.exists(bak): shutil.copy(P,bak); print("백업:",bak)
d=pd.read_csv(P); cz=pd.read_csv(os.path.join(C,"ST_census_livestock.csv"))
if "beef_cattle_farms_trend2019" not in d.columns:
    d=d.rename(columns={"beef_cattle_farms":"beef_cattle_farms_trend2019","dairy_farms":"dairy_farms_trend2019"})
else:
    d=d.drop(columns=["beef_cattle_farms","dairy_farms"],errors="ignore")
d=d.merge(cz.rename(columns={"census_beef_farms_km2":"beef_cattle_farms","census_dairy_farms_km2":"dairy_farms"}),on=["region","year"],how="left")
assert d.beef_cattle_farms.notna().all() and d.dairy_farms.notna().all()
d.to_csv(P,index=False,encoding="utf-8-sig")
# META 갱신
M=os.path.join(C,"ST_prereg_panel_META.csv"); m=pd.read_csv(M,keep_default_na=False)
for k,lab in (("beef_cattle_farms","Beef-cattle farms"),("dairy_farms","Dairy farms")):
    m.loc[m.key==k,"file"]="농림어업총조사 2010/2015/2020 (KOSIS DT_1AG501/1AG15501/1AG20501)"
    m.loc[m.key==k,"normalization"]="per km2 (census farm counts, linear interpolation 2011-2020, held at 2020 after)"
m.to_csv(M,index=False,encoding="utf-8-sig")
print("패널 교체 완료. 상관(신·구, 2020):",round(d[d.year==2020].beef_cattle_farms.corr(d[d.year==2020].beef_cattle_farms_trend2019),3))
