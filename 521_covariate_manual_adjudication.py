# -*- coding: utf-8 -*-
# 521 — 전수 공변량 수동 판정. 468개(1차 CANDIDATE+REVIEW)를 EHEC 병인학 기준으로 한 줄씩 판정.
#  출력: COVARIATE_PREREG_POOL.csv  — 도메인·방향가설·정규화·제외사유. 이 단계에서 어떤 IRR도 계산하지 않는다.
#  원칙: (1) 제외는 사유 필수 (2) 방향가설은 IRR 보기 전에 기록 (3) 같은 지표의 총량/비율/분위 중 하나만 (4) 결과변수·누출 변수 금지
import pandas as pd, re, os
C="data_cache"; Q=pd.read_csv(os.path.join(C,"COVARIATE_TRIAGE_v1.csv"))
S=Q[Q.status.isin(["CANDIDATE","REVIEW"])].copy()

# ── 포함 목록: (컬럼 패턴, 도메인, 영문 라벨, 방향가설, 정규화, 근거) ──
# 방향: + 위험, - 보호, ± 양방향/검출, 0 음성대조
INC=[
 # D1 반추동물 리저버 — 가축두수는 별도 파일(가축두수_전처리)에서 이미 패널에 있음(d_cattle/d_dairy/d_swine/d_poultry). 여기선 토지만.
 (r"^목장용지$","D1 reservoir","Pasture land area","+","per km2 of district","ruminant grazing land; manure-to-environment pathway"),
 # D2 농경지·토지이용 — 분뇨→농경지→농산물/관개수
 (r"^전$","D2 agri-land","Upland field (dry-field) area","+","share of district area","manure application; produce contamination"),
 (r"^답$","D2 agri-land","Paddy field area","±","share of district area","irrigation water; weaker manure link than upland"),
 (r"^과수원$","D2 agri-land","Orchard area","+","share of district area","manure-amended horticulture; fresh produce"),
 (r"^임야$","D2 agri-land","Forest area","-","share of district area","low human-livestock density; dilution"),
 (r"^하천$","D2 agri-land","River/stream land area","+","share of district area","surface-water transport of runoff"),
 (r"^유지$","D2 agri-land","Reservoir/pond land area","+","share of district area","standing water; irrigation source"),
 (r"^구거$","D2 agri-land","Irrigation ditch area","+","share of district area","irrigation channel density"),
 (r"^양어장$","D2 agri-land","Fish-farm land area","0","share of district area","negative control: aquaculture, not EHEC pathway"),
 (r"^대$","D2 agri-land","Residential building-site area","-","share of district area","urbanicity proxy (context, not exposure)"),
 (r"^공장용지$","D2 agri-land","Factory land area","0","share of district area","negative control: industrial land"),
 # D3 물·위생
 (r"^상수도보급률$","D3 water-sanit","Piped water supply coverage","-","%","municipal treated water reduces groundwater reliance"),
 (r"^하수도보급률$","D3 water-sanit","Sewerage coverage","-","%","sewered excreta removal"),
 (r"^공공하수처리구역인구보급률\(%\)$","D3 water-sanit","Public sewage-treatment coverage","-","%","treated sewage; duplicate-check vs 하수도보급률"),
 (r"^고도처리인구보급률\(%\)$","D3 water-sanit","Advanced wastewater treatment coverage","-","%","nutrient/pathogen removal"),
 (r"^하수처리구역외_정화조인구$","D3 water-sanit","Septic-tank population","+","% of district pop","on-site sanitation; groundwater contamination"),
 (r"^하수처리구역외_미처리인구$","D3 water-sanit","Untreated-sewage population","+","% of district pop","direct environmental release (low coverage 0.15; sensitivity only)"),
 (r"^하수도설치율\(%\)$","D3 water-sanit","Sewer installation rate","-","%","infrastructure; check collinearity with coverage"),
 (r"^총이용량$","D3 water-sanit","Drinking-water groundwater use","+","m3 per 1000 pop","untreated groundwater as drinking source; canonical EHEC waterborne route"),
 (r"^가정용$","D3 water-sanit","Household drinking groundwater use","+","m3 per 1000 pop","household-level untreated water"),
 (r"^간이상수도용$","D3 water-sanit","Simple-waterworks groundwater use","+","m3 per 1000 pop","small rural systems, minimal treatment"),
 (r"^학교용$","D3 water-sanit","School groundwater use","+","m3 per 1000 pop","child exposure route; EHEC is pediatric"),
 (r"^농업·생활겸용$","D3 water-sanit","Mixed agricultural-domestic groundwater","+","m3 per 1000 pop","irrigation-domestic cross-use"),
 (r"^가정용_이용량$","D3 water-sanit","Household (living-use) groundwater volume","+","m3 per 1000 pop","alternate series for household use; pick one after overlap check"),
 (r"^총 계_이용량$","D3 water-sanit","Total living-use groundwater volume","+","m3 per 1000 pop","alternate total series"),
 (r"^적합률\(%\)$","D3 water-sanit","Groundwater quality pass rate","-","%","inverse of contamination; replaces 부적합 count"),
 (r"^수밀검사·CCTV조사량\(m\)$","D3 water-sanit","Sewer inspection length","-","m per km2","maintenance intensity; leakage control"),
 (r"^개·보수관로_전체보수\(m\)_합류식$","D3 water-sanit","Combined-sewer repair length","±","m per km2","combined sewers overflow in summer rain"),
 # D4 개인위생·식품취급 — 손씻기 4종은 손씻기_전처리.csv 별도(이미 복구된 c_hw_*). 여기선 외식·식품.
 (r"^인구천명당주점업수$","D4 hygiene-food","Bars/pubs per 1000 pop","+","per 1000 pop","out-of-home eating; density form (not count)"),
 (r"^인구천명당패스트푸드점수$","D4 hygiene-food","Fast-food outlets per 1000 pop","+","per 1000 pop","ground-beef exposure; canonical EHEC food vehicle"),
 (r"^영양표시독해율_표준화율$","D4 hygiene-food","Nutrition-label literacy","-","% (age-std)","food-safety attentiveness"),
 (r"^아침결식예방인구율_표준화율$","D4 hygiene-food","Breakfast-eating rate","±","% (age-std)","home-prepared food proxy"),
 (r"^유아천명당보육시설수$","D4 hygiene-food","Childcare facilities per 1000 infants","+","per 1000 infants","daycare person-to-person transmission; pediatric EHEC"),
 (r"^유치원원아수$","D4 hygiene-food","Kindergarten enrolment","+","% of pop 3-5","congregate child settings"),
 (r"^초등학교학생수$","D4 hygiene-food","Primary-school enrolment","+","% of pop 6-11","school meals; congregate exposure"),
 (r"^학급당학생수_유치원$","D4 hygiene-food","Kindergarten class size","+","pupils per class","crowding in child settings"),
 (r"^학급당학생수_초등학교$","D4 hygiene-food","Primary class size","+","pupils per class","crowding"),
 # D5 SES·인구
 (r"^재정자립도$","D5 SES","Fiscal independence","-","%","local fiscal capacity; sanitation investment"),
 (r"^재정자주도$","D5 SES","Fiscal autonomy","-","%","alternate fiscal measure; pick one"),
 (r"^일반회계중보건세출비중$","D5 SES","Health share of general budget","-","%","public-health spending intensity"),
 (r"^일반회계중사회복지예산비중$","D5 SES","Welfare share of general budget","±","%","deprivation response vs deprivation level"),
 (r"^고용률$","D5 SES","Employment rate","-","%","economic level (coverage 0.37; sensitivity)"),
 (r"^실업률$","D5 SES","Unemployment rate","+","%","deprivation (coverage 0.37; sensitivity)"),
 (r"^인구천명당사업체수$","D5 SES","Businesses per 1000 pop","-","per 1000","economic activity"),
 (r"^인구천명당종사자수$","D5 SES","Workers per 1000 pop","-","per 1000","economic activity"),
 (r"^전문대학및대학교수$","D5 SES","Tertiary institutions","-","per 100k pop","education level proxy (supervisor hypothesis)"),
 (r"^인구천명당사설학원수$","D5 SES","Private academies per 1000 pop","-","per 1000","education investment; urban affluence"),
 (r"^고령인구비율$","D5 SES","Older-adult share","±","%","rural aging; excluded from main (offset absorbs age) but kept for urbanicity axis"),
 (r"^도시지역인구비율$","D5 SES","Urban population share","±","%","urbanicity axis (context)"),
 (r"^인구밀도$","D5 SES","Population density","±","per km2","contact rate vs detection"),
 (r"^성비$","D5 SES","Sex ratio","±","males per 100 females","forced demographic control"),
 (r"^인구천명당외국인수$","D5 SES","Foreign residents per 1000 pop","±","per 1000","dietary/hygiene heterogeneity; migrant labor in livestock"),
 (r"^순이동인구$","D5 SES","Net migration","±","per 1000 pop","population turnover"),
 (r"^주민등록세대수$","D5 SES","Households","±","persons per household (derive)","household size from pop/households"),
 (r"^소멸위험지수$","D5 SES","Local extinction-risk index","+","index","rural decline composite"),
 (r"^지가변동률$","D5 SES","Land-price change","-","%","development pressure"),
 (r"^주민1인당생활폐기물배출량$","D5 SES","Household waste per capita","+","kg per person","consumption/sanitation load"),
 (r"^합계출산율$","D5 SES","Total fertility rate","+","children per woman","young-child population (check vs offset)"),
 # D6 건강행태·예방 (사전지정 R7 개정에 따라 포함)
 (r"^연간인플루엔자예방접종률_표준화율$","D6 behaviour","Influenza vaccination coverage","-","% (age-std)","preventive-behaviour vs detection (bidirectional test)"),
 (r"^영유아건강검진율_전체$","D6 behaviour","Infant check-up rate","±","%","caregiver behaviour vs pediatric diagnosis opportunity"),
 (r"^걷기실천율_표준화율$","D6 behaviour","Walking practice","0","% (age-std)","negative control: general behaviour"),
 (r"^현재흡연율_표준화율$","D6 behaviour","Current smoking","+","% (age-std)","inverse health-behaviour indicator"),
 (r"^월간음주율_표준화율$","D6 behaviour","Monthly drinking","+","% (age-std)","eating-out / social dining exposure"),
 (r"^비만율_자가보고_표준화율$","D6 behaviour","Obesity (self-reported)","±","% (age-std)","diet pattern proxy"),
 (r"^건강생활실천율_표준화율$","D6 behaviour","Healthy-lifestyle practice","-","% (age-std)","composite preventive behaviour"),
 (r"^주관적건강수준인지율_표준화율$","D6 behaviour","Good self-rated health","±","% (age-std)","health status vs health-seeking"),
 (r"^연간보건기관이용률_표준화율$","D6 behaviour","Public-health-center use","±","% (age-std)","health-seeking; detection"),
 # D7 검출·의료접근 (통제/민감도 전용 — 결정요인 아님)
 (r"^인구천명당의료기관종사의사수$","D7 detection","Physicians per 1000 pop","±","per 1000","ascertainment proxy"),
 (r"^인구천명당의료기관병상수$","D7 detection","Hospital beds per 1000 pop","±","per 1000","ascertainment"),
 (r"^내과전문의$","D7 detection","Internal-medicine specialists","±","per 100k pop","GI presentation pathway"),
 (r"^진단검사의학과전문의$","D7 detection","Laboratory-medicine specialists","±","per 100k pop","culture/toxin testing capacity; most direct EHEC ascertainment proxy"),
 (r"^병리과전문의$","D7 detection","Pathology specialists","±","per 100k pop","lab capacity"),
 (r"^예방의학과전문의$","D7 detection","Preventive-medicine specialists","±","per 100k pop","surveillance capacity"),
 (r"^응급의학과전문의$","D7 detection","Emergency-medicine specialists","±","per 100k pop","acute presentation"),
 (r"^요양기관수_의원$","D7 detection","Clinics","±","per 100k pop","primary-care access"),
 (r"^요양기관수_병원$","D7 detection","Hospitals","±","per 100k pop","secondary access"),
 (r"^요양기관수_종합병원$","D7 detection","General hospitals","±","per 100k pop","tertiary access; HUS referral"),
 (r"^요양기관수_보건소$","D7 detection","Public health centers","±","per 100k pop","public surveillance node"),
 (r"^의원_내과$","D7 detection","Internal-medicine clinics","±","per 100k pop","GI primary care"),
 (r"^격리병실_병상수$","D7 detection","Isolation beds","±","per 100k pop","infectious-disease capacity"),
 (r"^전염병_진료실인원$","D7 detection","Infectious-disease outpatients","±","per 1000 pop","health-seeking for infection (leakage check: exclude if EHEC-inclusive)"),
 (r"^관내진료실인원_외래$","D7 detection","In-district outpatient volume","±","per 1000 pop","local health-seeking"),
 (r"^관외진료실인원_외래$","D7 detection","Out-of-district outpatient volume","±","per 1000 pop","referral outflow; rural detection leakage"),
 (r"^의사수$","D7 detection","Physicians (count)","±","per 100k pop","alternate to per-1000 series; pick one"),
 (r"^간호사수$","D7 detection","Nurses","±","per 100k pop","health workforce"),
 (r"^약사수$","D7 detection","Pharmacists","±","per 100k pop","OTC self-treatment → under-ascertainment"),
 (r"^요양기관수_약국$","D7 detection","Pharmacies","±","per 100k pop","self-medication access"),
]

# ── 명시적 제외(사유) — 패턴 ──
EXC=[
 (r"사망|결과_|진료비|_대상인원|_수검인원|검진판정|암검진|암_진료|잠혈|용종|이상소견|유질환자|일반질환의심|정상A","outcome/screening result or cost: not an exposure"),
 (r"^연도$|^분기$|^value$|^값$|^합계$|^적합$|^부적합$|총면적|총인구|남자인구|여자인구|노인인구수|행정구역면적|총 계_개소수|개소수\(총합\)","structural/denominator column, not a covariate"),
 (r"_개소수$|_시설수$","facility count: volume series retained instead"),
 (r"온천|광천|염전|사적지|유원지|묘지|종교용지|주차장|철도|창고|체육용지|학교용지|수도용지|잡종지|제방|주유소|도로|보행|운전|자동차|교통","land use / transport with no plausible EHEC pathway"),
 (r"구강|치과|치주|불소|스케일링|치위생|한의|한방|물리치료|작업치료|재활|정신|우울|스트레스|사회복지사|노인장기요양|노인여가|방문건강|모자보건|산부인과|조산원|출생|출산|이혼|혼인","unrelated clinical/social domain"),
 (r"성형|안과|이비인후|피부|비뇨|정형|신경|흉부|마취|영상의학|방사선|핵의학|가정의학|직업환경|외과전문의|결핵","specialty with no EHEC diagnostic role"),
 (r"고혈압|당뇨|간암|간질환|대장암|위암|유방암|자궁|폐암|관절|혈당|혈압|비만인지|체중조절|흡연자|간접흡연|흡연시작|남자현재|여자현재|신체활동|격렬|중등도|평생음주|연간음주","chronic-disease / behaviour sub-series: canonical series retained"),
 (r"영유아건강검진율_\d|영유아건강검진_|영유아검진","infant check-up sub-series: overall rate retained"),
 (r"전문대학 및 대학교|유치원 수|초등학교 수|초등학교 학생|초등학교교원|유치원교원|유치원수$|교원1인당|학급당학생수_(중|고|전체)","education sub-series: enrolment/class-size retained; supervisor GDP/education axis via 전문대학및대학교수"),
 (r"경제활동참가율|청년고용률|국민기초생활보장수급자","coverage < 0.10: unusable"),
 (r"요양기관수_(보건지소|보건진료소|보건의료원|상급종합|요양병원|치과|한방|한의원)|의원_(?!내과)","facility sub-type without EHEC role"),
 (r"의료기관진료|관내진료비|관외진료비|관내진료실인원_(입원|전체)|관외진료실인원_(입원|전체)|전염병_진료비|병상수_계|병실수|일반입원실|중환자실|정신과","utilisation series: outpatient volume retained"),
 (r"양전자|자기공명|초음파|체외충격|유방촬영","imaging equipment: no EHEC role"),
 (r"인구증가율|자연증가율|시도간|총전입|총전출|일인당자동차|인구천명당스포츠|일반공공행정|사고및중독|안전벨트","no plausible pathway"),
 (r"받이_|토실|개·보수관로_전체보수\(m\)_분류식|하수처리구역내|하수처리구역외_(계|오수)","sewer sub-component: canonical measures retained"),
 (r"공동주택용|민방위용|일반용|기타|상수도용|간이상수도용_이용량|학교용_이용량|농업·생활겸용_이용량|총계_","groundwater sub-series: drinking-water series retained"),
 (r"연간당뇨|고혈압관리|당뇨병관리|고혈압진단|당뇨병진단","chronic-disease management"),
]

rows=[]
for _,r in S.iterrows():
    c=r.col; hit=None
    for pat,dom,lab,dirn,norm,why in INC:
        if re.search(pat,c): hit=("INCLUDE",dom,lab,dirn,norm,why); break
    if hit is None:
        for pat,why in EXC:
            if re.search(pat,c): hit=("EXCLUDE","","","","",why); break
    if hit is None: hit=("REVIEW","","","","","not matched by include/exclude rules — needs a decision")
    rows.append(dict(file=r.file,col=c,n_years=r.n_years_obs,coverage=r.coverage,status=hit[0],domain=hit[1],label=hit[2],
                     direction=hit[3],normalization=hit[4],rationale=hit[5]))
P=pd.DataFrame(rows)
# 1차 triage의 EXCLUDE도 합쳐 전수 기록
E1=Q[Q.status=="EXCLUDE"][["file","col","n_years_obs","coverage","reason"]].rename(columns={"n_years_obs":"n_years","reason":"rationale"})
E1["status"]="EXCLUDE"; E1[["domain","label","direction","normalization"]]=""
ALL=pd.concat([P,E1[P.columns]],ignore_index=True)
ALL.to_csv(os.path.join(C,"COVARIATE_PREREG_POOL.csv"),index=False,encoding="utf-8-sig")
print(f"전수 {len(ALL):,}개 판정 완료 → data_cache/COVARIATE_PREREG_POOL.csv\n")
print(ALL.status.value_counts().to_string(),"\n")
I=ALL[ALL.status=="INCLUDE"]
print("── INCLUDE 도메인별 ──"); print(I.domain.value_counts().sort_index().to_string(),"\n")
print("── 방향가설 분포 ──"); print(I.direction.value_counts().to_string(),"\n")
rv=ALL[ALL.status=="REVIEW"]
print(f"── REVIEW 잔여 {len(rv)}개 ──")
for c in rv.col: print("   ",c)
