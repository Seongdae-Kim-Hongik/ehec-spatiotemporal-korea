# 공유 셋업 — ①~④ 시공간 후속분석 공통 (그래프·인접·prior·패널·RICH·Z빌더)
#  ★2026-10-21 v2(이전 판 _st_setup.R.bak_pre_v2, _v2_common.R): build_Z 평균 채움 금지(요청한 변수에 결측이 있으면 멈춤)·Z를 코드로 shp 라벨에 맞춤,
#   panel(ST_panel_15x229x14.rds)에 sigungu_cd 결합(이름이 다르면 멈춤). E·사례는 이 셋업을 쓰는 스크립트에서 v2_E_rc()·v2_cases_rc()·v2_finish()로(411·411b·451 수정).
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(Matrix)})
inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
source(file.path(ROOT,"_v2_common.R"))
EN<-c("노로바이러스"="Norovirus","A형간염"="HepatitisA","E형간염"="HepatitisE","황색포도상구균"="S_aureus","바실러스세레우스"="B_cereus","클로스트리디움퍼프린젠스"="C_perfringens","비브리오패혈증"="V_vulnificus","장염비브리오"="V_parahaemolyticus","장티푸스"="Typhoid","세균성이질"="Shigellosis","살모넬라"="Salmonella","병원성대장균"="Pathogenic_Ecoli","장출혈성대장균감염증"="EHEC","캠필로박터제주니"="Campylobacter","원충"="Kudoa")
KLAB<-c(a_water_supply="상수도보급",a_sewer="하수도보급",a_septic="정화조인구%",a_advtreat="고도하수처리",a_wastewater="폐수발생량",a_orgload="유기물부하",b_gw_use="지하수이용",b_gw_sites="지하수개소",b_gw_home="가정지하수",b_gw_simplewater="간이상수도",b_gw_unfit="지하수부적합%",c_hw_out="외출후손씻기",c_hw_meal="식전손씻기",c_hw_soap="비누사용",c_hw_toilet="화장실후손씻기",c_flu_vax="독감백신",d_cattle="한육우밀도",d_dairy="젖소밀도",d_swine="돼지밀도",d_poultry="가금밀도",d_paddy="논면적",d_field="밭면적",d_pasture="목장면적",e_urban_pop="도시인구%",e_rural_pop="농촌인구%",e_netmig="순이동",e_built="대지",e_road="도로",e_forest="임야",e_pop_dens="인구밀도",f_shellfish="패류생산",f_crustacean="갑각류생산",f_fish="어류생산",f_aqua="양식생산",g_elderly="고령%",g_solo_elder="독거노인%",g_sexratio="성비",g_fiscal="재정자립",g_basiclivelihood="기초생활%",g_singlehh="1인가구",g_infant="영유아(0-4)%",h_docs="의사수",h_pub="주점수",h_fastfood="패스트푸드",h_selfhealth="자가건강",temp="평균기온",sst="표층수온",precip="강수량")
# 영어 질병명(그림 제목용)
EN_DZ<-c("노로바이러스"="Norovirus","A형간염"="Hepatitis A","E형간염"="Hepatitis E","황색포도상구균"="Staphylococcus aureus","바실러스세레우스"="Bacillus cereus","클로스트리디움퍼프린젠스"="Clostridium perfringens","비브리오패혈증"="Vibrio vulnificus sepsis","장염비브리오"="Vibrio parahaemolyticus","장티푸스"="Typhoid fever","세균성이질"="Shigellosis","살모넬라"="Salmonellosis","병원성대장균"="Pathogenic E. coli","장출혈성대장균감염증"="EHEC","캠필로박터제주니"="Campylobacter","원충"="Kudoa")
endz<-function(k){v<-unname(EN_DZ[k]);ifelse(is.na(v),k,v)}
# 영어 공변량 라벨(그림/표 영어 표기용) — KLAB과 동일 키
ELAB<-c(a_water_supply="Piped-water coverage",a_sewer="Sewerage coverage",a_septic="Septic-tank population %",a_advtreat="Advanced sewage treatment",a_wastewater="Wastewater generation",a_orgload="Organic load",b_gw_use="Groundwater use",b_gw_sites="Groundwater sites",b_gw_home="Household groundwater",b_gw_simplewater="Village water supply",b_gw_unfit="Unfit groundwater %",c_hw_out="Handwashing after outing",c_hw_meal="Handwashing before meals",c_hw_soap="Soap use",c_hw_toilet="Handwashing after toilet",c_flu_vax="Influenza vaccination",d_cattle="Beef-cattle density",d_dairy="Dairy-cattle density",d_swine="Swine density",d_poultry="Poultry density",d_paddy="Paddy-field area",d_field="Upland-field area",d_pasture="Pasture area",e_urban_pop="Urban population %",e_rural_pop="Rural population %",e_netmig="Net migration",e_built="Built-up land",e_road="Road area",e_forest="Forest area",e_pop_dens="Population density",f_shellfish="Shellfish production",f_crustacean="Crustacean production",f_fish="Fish production",f_aqua="Aquaculture production",g_elderly="Elderly %",g_solo_elder="Solo-elderly %",g_sexratio="Sex ratio",g_fiscal="Fiscal independence",g_basiclivelihood="Basic-livelihood recipients %",g_singlehh="Single-person households",g_infant="Infants (0-4) %",h_docs="Physicians",h_pub="Bars/pubs",h_fastfood="Fast-food outlets",h_selfhealth="Self-rated health",temp="Mean temperature",sst="Sea-surface temperature",precip="Precipitation")
WATER<-c("a_water_supply","a_sewer","a_septic","a_advtreat","a_wastewater","a_orgload","b_gw_use","b_gw_sites","b_gw_home","b_gw_simplewater","b_gw_unfit","c_hw_out","c_hw_meal","c_hw_toilet","c_hw_soap","e_rural_pop","g_solo_elder","g_basiclivelihood","h_fastfood","h_pub")
TOXIN<-c("h_fastfood","h_pub","e_pop_dens","e_urban_pop","temp","c_hw_meal","c_hw_soap")
RICH<-list(
 "노로바이러스"=c("g_infant","a_sewer","a_septic","a_wastewater","a_orgload","b_gw_use","b_gw_sites","b_gw_home","b_gw_unfit","f_shellfish","f_aqua","f_crustacean","c_hw_meal","c_hw_toilet","c_hw_soap","e_urban_pop","e_pop_dens","h_pub","h_fastfood"),
 "A형간염"=WATER,"E형간염"=WATER,"장티푸스"=WATER,"세균성이질"=WATER,
 "장출혈성대장균감염증"=c("d_cattle","d_dairy","d_pasture","d_field","d_paddy","b_gw_use","b_gw_home","b_gw_unfit","g_infant","e_rural_pop","e_forest","a_septic","h_fastfood","h_pub"),
 "살모넬라"=c("d_poultry","d_swine","d_cattle","d_dairy","h_fastfood","h_pub","e_urban_pop","e_pop_dens","temp"),
 "캠필로박터제주니"=c("d_poultry","d_cattle","d_swine","d_dairy","b_gw_use","b_gw_home","e_rural_pop","temp","h_fastfood","h_pub"),
 "비브리오패혈증"=c("sst","temp","f_shellfish","f_aqua","f_crustacean","f_fish","g_elderly","g_solo_elder","e_rural_pop"),
 "장염비브리오"=c("sst","temp","f_shellfish","f_aqua","f_crustacean","f_fish","e_rural_pop"),
 "원충"=c("f_aqua","f_fish","f_shellfish","f_crustacean","sst","temp"),
 "황색포도상구균"=TOXIN,"바실러스세레우스"=TOXIN,"클로스트리디움퍼프린젠스"=TOXIN,
 "병원성대장균"=c("a_sewer","a_septic","b_gw_use","b_gw_home","c_hw_meal","c_hw_soap","e_urban_pop","e_pop_dens","h_fastfood","temp"))
ALLDZ<-names(RICH)
## 그래프
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
## 2026-10-21 재발 방지: final.dbf 인천 region 라벨(남구 23030 = 미추홀구 약 25 km², 남동구 23050 약 57 km²)이 면적과 맞는지 단언(2026-10-21 dbf 두 칸 교정, 백업 final.dbf.bak_pre_shpfix)
local({a<-as.numeric(sf::st_area(shp))/1e6;stopifnot(nrow(shp)==229,!anyDuplicated(shp$region),a[shp$region=="인천시남구"]<35,a[shp$region=="인천시남동구"]>45)})
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));NR<-nrow(shp);regs<-shp$region
build_nb<-function(type="queen",k=NULL){
  if(is.null(k)){nb<-poly2nb(shp,queen=(type=="queen"),snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
    for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
    cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
    nb<-make.sym.nb(nb)} else {nb<-make.sym.nb(knn2nb(knearneigh(coords,k=k)))}
  nb}
nb<-build_nb("queen");g<-tempfile(fileext=".graph");nb2INLA(g,nb);lw<-nb2listw(nb,style="W",zero.policy=TRUE)
Wmat<-nb2mat(nb,style="B",zero.policy=TRUE);Ra<-inla.scale.model(Diagonal(x=rowSums(Wmat))-Matrix(Wmat,sparse=TRUE),constr=list(A=matrix(1,1,NR),e=0))
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
panel<-readRDS(file.path(CACHE,"ST_panel_15x229x14.rds"))|>as.data.frame();panel$sigungu_cd<-v2_code(panel$region)   # ★v2
## 공변량 Z 빌더 (창 지정)
CV<-readRDS(file.path(CACHE,"ST_covariates_panel.rds"));cl<-readRDS(file.path(CACHE,"climate_covariates.rds"))
meltm<-function(M,nm){d<-as.data.frame(M);d$region<-rownames(M);pivot_longer(d,-region,names_to="year",values_to=nm)%>%mutate(year=as.integer(year))}
clim<-meltm(cl$temp,"temp")%>%left_join(meltm(cl$precip,"precip"),by=c("region","year"))%>%left_join(meltm(cl$sst,"sst"),by=c("region","year"))
basecov<-data.frame(region=CV$region,year=CV$year,CV$raw,check.names=FALSE)%>%left_join(clim,by=c("region","year"))
skewf<-function(x){x<-x[is.finite(x)];if(sd(x)==0)return(0);e1071::skewness(x)}
build_Z<-function(vars,y0,y1){b<-basecov[basecov$year>=y0&basecov$year<=y1,]
  Z<-sapply(vars,function(v)v2_zb(b[[v]],v))   # ★v2: 평균 채움 금지(결측이면 멈춤)
  Z<-as.data.frame(Z);Z$region<-b$region;Z$year<-b$year;v2_relabel(Z,"ST_covariates_panel",keep_code=FALSE)}   # ★v2: 코드로 shp 라벨
getsel<-function(){a<-commandArgs(trailingOnly=TRUE);if(length(a)>=1&&nchar(a[1])>0)ALLDZ[as.integer(strsplit(a[1],",")[[1]])] else ALLDZ}
cat("_st_setup 로드 완료 (NR=",NR,")\n")
