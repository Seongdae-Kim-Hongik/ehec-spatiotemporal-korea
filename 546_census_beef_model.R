# 546 — 교수 규칙 검정 (type I 상호작용 = f(id, iid) 시공간 셀 무작위효과): 한우·젖소 농가밀도를 농림어업총조사(2010/2015/2020, 2011+ 커버)로 교체한 최종 20변수 모형.
#  (A) 시변 버전(보간 패널) (B) 시간불변 버전(3시점 평균). 나머지 18변수·구조·프라이어는 540과 동일.
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-file.path(ROOT,"results_ST")
V2_CACHE<-CACHE;source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지, KDCA형 청주 2011-2013 마스킹(CJ_MASK=0이면 끔). 이전 판 .bak_pre_v2
cont<-c("beef_cattle_farms","dairy_farms","reservoir_pond_land_area","upland_field_dry_field_area","forest_area","piped_water_supply_coverage","sewerage_coverage","septic_tank_population","drinking_water_groundwater_use","fast_food_outlets_per_1000_pop","childcare_facilities_per_1000_infants","soap_sanitizer_use","fiscal_independence","net_migration","older_adults_living_alone","sex_ratio","current_smoking","influenza_vaccination_coverage","kchs_pc1","physicians_per_1000_pop")
d0<-read.csv(file.path(CACHE,"ST_prereg_panel.csv"),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
cz<-read.csv(file.path(CACHE,"ST_census_livestock.csv"),fileEncoding="UTF-8-BOM");cz$region<-str_replace_all(cz$region,"\\s+","")
d0<-d0%>%left_join(cz,by=c("region","year"))
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
zs<-function(x){x<-v2_nofill(x,deparse(substitute(x)));xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;z}
# 548 패치 이후 패널의 beef/dairy는 이미 총조사 값 → baseline(모니터링)은 보존 컬럼 사용
if("beef_cattle_farms_trend2019" %in% names(d0)){d0$beef_cattle_farms<-d0$beef_cattle_farms_trend2019;d0$dairy_farms<-d0$dairy_farms_trend2019}
Z<-d0[,c("region","year")];for(v in cont)Z[[v]]<-zs(d0[[v]])
Z$census_beef_tv<-zs(d0$census_beef_farms_km2);Z$census_dairy_tv<-zs(d0$census_dairy_farms_km2)
mB<-d0%>%group_by(region)%>%summarise(b=mean(census_beef_farms_km2),dd=mean(census_dairy_farms_km2));d1<-d0%>%left_join(mB,by="region")
Z$census_beef_ti<-zs(d1$b);Z$census_dairy_ti<-zs(d1$dd)
# 연도별 표준화(시간 추세 완전 제거; 교차검토 민감도)
zy<-function(x,yr){out<-rep(NA_real_,length(x));for(y in unique(yr)){i<-yr==y;xi<-x[i];xi[is.na(xi)]<-mean(xi,na.rm=TRUE);xv<-if(min(xi,na.rm=TRUE)>=0&&skewf(xi)>1)log1p(pmax(xi,0))else xi;z<-scale(xv)[,1];z[!is.finite(z)]<-0;out[i]<-z};out}
Z$census_beef_yw<-zy(d0$census_beef_farms_km2,d0$year);Z$census_dairy_yw<-zy(d0$census_dairy_farms_km2,d0$year)
cat(sprintf("상관(표준화 후): 기존 vs 총조사 시변 %.3f · 시간불변 %.3f\n",cor(Z$beef_cattle_farms,Z$census_beef_tv),cor(Z$beef_cattle_farms,Z$census_beef_ti)))
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")];yy<-2011:2024
DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
stopifnot(!anyNA(DD$census_beef_tv),!any(duplicated(paste(DD$region,DD$year))))
DD<-v2_finish(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD)
fitmod<-function(cv){form<-as.formula(paste0("cases~",paste(sprintf("`%s`",cv),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(form,family="nbinomial",data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
rest<-setdiff(cont,c("beef_cattle_farms","dairy_farms"))
runs<-list(baseline_trend=c("beef_cattle_farms","dairy_farms",rest),census_timevarying=c("census_beef_tv","census_dairy_tv",rest),census_timeinvariant=c("census_beef_ti","census_dairy_ti",rest),census_yearwise_z=c("census_beef_yw","census_dairy_yw",rest))
out<-NULL
for(nm in names(runs)){m<-fitmod(runs[[nm]]);fx<-m$summary.fixed
  for(v in runs[[nm]][1:2]){r<-fx[v,];pp<-1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]])
    out<-rbind(out,data.frame(model=nm,cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(pp,3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),DIC=round(m$dic$dic,1)))}
  cat(sprintf("[%s] DIC %.1f\n",nm,m$dic$dic));for(v in runs[[nm]][1:2]){r<-fx[v,];cat(sprintf("   %-18s IRR %.3f (%.3f-%.3f)%s\n",v,exp(r$mean),exp(r$`0.025quant`),exp(r$`0.975quant`),ifelse(r$`0.025quant`>0|r$`0.975quant`<0," *","")))}}
write.csv(out,file.path(RS,"MASTER_census_beef_check.csv"),row.names=FALSE);cat("저장: MASTER_census_beef_check.csv\n")
