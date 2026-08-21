# ============================================================================
# 공간교란(spatial confounding) 민감도 — 제한적 공간회귀 RSR
#  /333 공통 HIGH 지적: 공간구조 공변량이 BYM2 ICAR 필드와 공유구조 → 고정효과 편향.
#  Reich-Hodges-Zadnik 2006 / Hughes-Haran 2013: 공간필드를 X에 직교(X'u=0)시켜
#  교란 제거한 IRR 회복. 표준 ICAR vs RSR 계수 비교(280 주분석 대비 민감도).
#  공간전용 모형(연도합산·지역평균 공변량, ICAR)으로 apples-to-apples 비교.
#  산출: results_ST/MASTER_RSR_sensitivity.csv (질병×공변량: IRR_std vs IRR_rsr)
# ============================================================================
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(Matrix)})
set.seed(20260707); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".")
CACHE<-file.path(ROOT,"data_cache")
EN<-c("노로바이러스"="Norovirus","A형간염"="HepatitisA","E형간염"="HepatitisE","황색포도상구균"="S_aureus","바실러스세레우스"="B_cereus","클로스트리디움퍼프린젠스"="C_perfringens","비브리오패혈증"="V_vulnificus","장염비브리오"="V_parahaemolyticus","장티푸스"="Typhoid","세균성이질"="Shigellosis","살모넬라"="Salmonella","병원성대장균"="Pathogenic_Ecoli","장출혈성대장균감염증"="EHEC","캠필로박터제주니"="Campylobacter","원충"="Kudoa")
KLAB<-c(a_water_supply="상수도보급",a_sewer="하수도보급",a_septic="정화조인구%",a_advtreat="고도하수처리",a_wastewater="폐수발생량",a_orgload="유기물부하",b_gw_use="지하수이용",b_gw_sites="지하수개소",b_gw_home="가정지하수",b_gw_simplewater="간이상수도",b_gw_unfit="지하수부적합%",c_hw_out="외출후손씻기",c_hw_meal="식전손씻기",c_hw_soap="비누사용",c_hw_toilet="화장실후손씻기",c_flu_vax="독감백신",d_cattle="한육우밀도",d_dairy="젖소밀도",d_swine="돼지밀도",d_poultry="가금밀도",d_paddy="논면적",d_field="밭면적",d_pasture="목장면적",e_urban_pop="도시인구%",e_rural_pop="농촌인구%",e_netmig="순이동",e_built="대지",e_road="도로",e_forest="임야",e_pop_dens="인구밀도",f_shellfish="패류생산",f_crustacean="갑각류생산",f_fish="어류생산",f_aqua="양식생산",g_elderly="고령%",g_solo_elder="독거노인%",g_sexratio="성비",g_fiscal="재정자립",g_basiclivelihood="기초생활%",g_singlehh="1인가구",g_infant="영유아(0-4)%",h_docs="의사수",h_pub="주점수",h_fastfood="패스트푸드",h_selfhealth="자가건강",temp="평균기온",sst="표층수온",precip="강수량")
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
 "황색포도상구균"=TOXIN,"바실러스세레우스"=TOXIN,"클로스트리디움퍼프린젠스"=TOXIN)

## 강화패널(native+기후) → 지역평균 공변량(연도합산 공간모형) + a priori 변환
CV<-readRDS(file.path(CACHE,"ST_covariates_panel.rds"));cl<-readRDS(file.path(CACHE,"climate_covariates.rds"))
meltm<-function(M,nm){d<-as.data.frame(M);d$region<-rownames(M);pivot_longer(d,-region,names_to="year",values_to=nm)%>%mutate(year=as.integer(year))}
clim<-meltm(cl$temp,"temp")%>%left_join(meltm(cl$precip,"precip"),by=c("region","year"))%>%left_join(meltm(cl$sst,"sst"),by=c("region","year"))
base0<-data.frame(region=CV$region,year=CV$year,CV$raw,check.names=FALSE)%>%left_join(clim,by=c("region","year"))
allv<-unique(unlist(RICH));skewf<-function(x){x<-x[is.finite(x)];if(sd(x)==0)return(0);e1071::skewness(x)}

## 그래프
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));NR<-nrow(shp);regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
panel<-readRDS(file.path(CACHE,"ST_panel_15x229x14.rds"))|>as.data.frame()

allrows<-list()
for(dz in names(RICH)){eng<-EN[dz];OUT<-file.path(ROOT,eng);cov<-RICH[[dz]]
  pd<-panel[panel$disease==dz,];fy<-unique(pd$fy);yy<-(2011:2024)[(2011:2024)>=fy]
  # 지역평균 공변량(분석기간), 지역합산 cases/E
  bb<-base0[base0$year%in%yy,]
  Xreg<-sapply(cov,function(v){ag<-tapply(bb[[v]],bb$region,function(x)mean(x,na.rm=TRUE))[regs];ag[!is.finite(ag)]<-mean(ag,na.rm=TRUE)
    xv<-if(min(ag,na.rm=TRUE)>=0&&skewf(ag)>1)log1p(pmax(ag,0))else ag;z<-scale(xv)[,1];z[!is.finite(z)]<-0;z})
  agg<-pd[pd$year%in%yy,]%>%group_by(region)%>%summarise(cases=sum(cases,na.rm=TRUE),pop=mean(pop,na.rm=TRUE),.groups="drop")
  d1<-data.frame(region=regs)%>%left_join(agg,by="region");d1$cases[is.na(d1$cases)]<-0
  d1$pop[is.na(d1$pop)|d1$pop<=0]<-median(d1$pop[d1$pop>0],na.rm=TRUE)
  d1$E<-pmax(d1$pop*sum(d1$cases)/sum(d1$pop),1);d1$area<-1:NR
  d1<-cbind(d1,as.data.frame(Xreg));rhs<-paste(sprintf("`%s`",cov),collapse="+")
  # 표준 ICAR
  f_std<-as.formula(paste0("cases~",rhs,"+f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)"))
  m_std<-tryCatch(inla(f_std,family="nbinomial",data=d1,E=d1$E,control.fixed=CF,control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
  # RSR: 공간필드를 X의 열공간에 직교 (Hughes-Haran). 공선성 대비 QR 정규직교 기저로 full-rank 제약.
  Xc<-as.matrix(Xreg)  # n×k (이미 표준화)
  qrX<-qr(Xc);rk<-qrX$rank;Qb<-qr.Q(qrX)[,seq_len(rk),drop=FALSE]  # n×rk 정규직교 기저
  A_rsr<-t(Qb)  # rk×n (full row rank → segfault 방지)
  f_rsr<-as.formula(paste0("cases~",rhs,"+f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,extraconstr=list(A=A_rsr,e=rep(0,rk)),hyper=pcp)"))
  m_rsr<-tryCatch(inla(f_rsr,family="nbinomial",data=d1,E=d1$E,control.fixed=CF,control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
  if(is.null(m_std)||is.null(m_rsr)){cat("ERR",eng,"\n");next}
  fs<-m_std$summary.fixed[cov,,drop=FALSE];fr<-m_rsr$summary.fixed[cov,,drop=FALSE]
  r<-data.frame(disease=dz,eng=eng,공변량=ifelse(cov%in%names(KLAB),KLAB[cov],cov),cov=cov,
    IRR_std=round(exp(fs$mean),3),lo_std=round(exp(fs$`0.025quant`),3),hi_std=round(exp(fs$`0.975quant`),3),cred_std=(fs$`0.025quant`>0|fs$`0.975quant`<0),
    IRR_rsr=round(exp(fr$mean),3),lo_rsr=round(exp(fr$`0.025quant`),3),hi_rsr=round(exp(fr$`0.975quant`),3),cred_rsr=(fr$`0.025quant`>0|fr$`0.975quant`<0))
  r$shift_pct<-round(100*(r$IRR_rsr-r$IRR_std)/r$IRR_std,1)
  write.csv(r,file.path(OUT,paste0(eng,"_RSR_sensitivity.csv")),row.names=FALSE);allrows[[eng]]<-r
  cat(sprintf("[%s] std↔rsr 완료 (부호전환 %d개)\n",dz,sum(r$cred_std!=r$cred_rsr)))}
M<-bind_rows(allrows);write.csv(M,file.path(ROOT,"results_ST/MASTER_RSR_sensitivity.csv"),row.names=FALSE)
cat("\n=== RSR 민감도: credible 상태 변경(공간교란 노출) ===\n")
ch<-M[M$cred_std!=M$cred_rsr,];if(nrow(ch)>0)print(ch[,c("disease","공변량","IRR_std","cred_std","IRR_rsr","cred_rsr","shift_pct")],row.names=FALSE) else cat("없음 — 공간교란에 강건\n")
cat("\n저장: MASTER_RSR_sensitivity.csv\n")
