# 540s — (2026-10-22 v2c 재실행) 540 원고 모드의 사전 지정 민감도 사본. 공유 540은 그대로 두고 이 사본만 SENS 환경변수로 바꿈.
#  SENS=nopc1(조사지표 PC1 제외) · jejudummy(제주시·서귀포시 지시변수 표준화 추가) · cj0(CJ_MASK=0과 함께 = 청주 2011-2013 관측 0). 출력 = EHEC_RS. 음성대조는 생략.
# 540 — 최종 주모형 (사전지정 풀 멀티버스 v2 이후). 사전지정 도메인 대표(각 도메인 상위 credible 1~2개) + 도농축 + 성비 + 검출통제(의사/천명).
#  산출: Table 2(MASTER_final_table2.csv) · 지도용 regionSIR/exceedance · 적합통계 · 음성대조
#  ★2026-10-21 v2(이전 판 .bak_pre_v2, _v2_common.R): E·사례·공변량 sigungu_cd 결합, E 최솟값 채움·사례 0 채움 삭제(panel_gate), 공변량 평균 채움 금지(결측 있으면 멈춤),
#   KDCA형 청주 2011-2013 사례 NA 마스킹이 표준(CJ_MASK=0이면 관측 0 사용 = 민감도). 패널은 v2 정본(shp 면적·v2 인구·명시 채움, data_cache).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(car)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-Sys.getenv("EHEC_RS",unset=file.path(ROOT,"results_ST"))
source(file.path(ROOT,"_v2_common.R"))
# ── 사전지정 주모형: 도메인별 대표 (멀티버스 credible% 상위, 병인 경로 대표성) ──
DOM<-list(
 "1. Ruminant reservoir"=c("beef_cattle_farms","dairy_farms"),
 "2. Agriculture and land use"=c("reservoir_pond_land_area","upland_field_dry_field_area","forest_area"),
 "3. Water and sanitation"=c("piped_water_supply_coverage","sewerage_coverage","septic_tank_population","drinking_water_groundwater_use"),
 "4. Hygiene, food, and child settings"=c("fast_food_outlets_per_1000_pop","childcare_facilities_per_1000_infants","soap_sanitizer_use"),
 "5. Socioeconomic and demographic"=c("fiscal_independence","net_migration","older_adults_living_alone","sex_ratio"),
 "6. Health behavior"=c("current_smoking","influenza_vaccination_coverage"),
 "7. Context and ascertainment"=c("kchs_pc1","physicians_per_1000_pop"))
ENG<-c(beef_cattle_farms="Beef-cattle farm density",dairy_farms="Dairy farm density",reservoir_pond_land_area="Reservoir and pond area",upland_field_dry_field_area="Upland-field area",forest_area="Forest area",
 piped_water_supply_coverage="Piped water supply coverage",sewerage_coverage="Sewerage coverage",septic_tank_population="Septic-tank population",drinking_water_groundwater_use="Drinking-water groundwater use",
 fast_food_outlets_per_1000_pop="Fast-food and eating-out outlets",childcare_facilities_per_1000_infants="Childcare facilities",soap_sanitizer_use="Handwashing with soap",
 fiscal_independence="Fiscal independence",net_migration="Net migration",older_adults_living_alone="Older adults living alone",sex_ratio="Sex ratio",
 current_smoking="Current smoking",influenza_vaccination_coverage="Influenza vaccination coverage",kchs_pc1="Rural-urban axis (survey PC1)",physicians_per_1000_pop="Physicians per 1000 population")
NEG<-c("swine_density","poultry_density","walking_practice")
if(nzchar(Sys.getenv("EHEC_HM"))){   # ★2026-09-30 문헌 기반 풀: 700_hm_setcompare 최선 세트 + sets json 강제통제·음성대조로 DOM 교체
  HMD<-file.path(ROOT,"문헌가설_260930","EHEC");J<-jsonlite::fromJSON(file.path(HMD,"sets_EHEC.json"));BS<-sub("best_set=","",grep("^best_set=",readLines(file.path(HMD,"HM_best_EHEC.txt")),value=TRUE))
  HT<-read.csv(file.path(HMD,"hypothesis_table_EHEC.csv"));bv<-unlist(J$sets[[BS]]);dm<-HT$domain[match(bv,HT$variable_key)]
  DOM<-c(split(bv,dm),list("Forced controls"=J$forced));NEG<-J$negative_controls
  ME<-read.csv(file.path(CACHE,"HM_panel_META.csv"),fileEncoding="UTF-8-BOM");ENG<-setNames(ME$label,ME$key);ENG["kchs_pc1"]<-"Rural-urban axis (survey PC1)"
  cat("[HM] EHEC 문헌 최선 세트:",BS,"\n")}
cont<-unlist(DOM,use.names=FALSE)
SENS<-Sys.getenv("SENS");stopifnot(SENS%in%c("nopc1","jejudummy","cj0"))
if(SENS=="nopc1"){stopifnot("kchs_pc1"%in%cont);DOM<-lapply(DOM,function(x)setdiff(x,"kchs_pc1"));DOM<-DOM[lengths(DOM)>0];cont<-unlist(DOM,use.names=FALSE)}
if(SENS=="jejudummy"){DOM[["8. Jeju indicator"]]<-"jeju_dummy";ENG["jeju_dummy"]<-"Jeju-si/Seogwipo-si indicator";cont<-unlist(DOM,use.names=FALSE)}
if(SENS=="cj0")stopifnot(identical(Sys.getenv("CJ_MASK"),"0"))
d0<-read.csv(file.path(CACHE,Sys.getenv("EHEC_PANEL",unset="ST_prereg_panel.csv")),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
if(SENS=="jejudummy"){jd<-as.numeric(d0$region%in%c("제주도제주시","제주도서귀포시"));d0$jeju_dummy<-jd}
Z<-d0[,c("region","year")];for(v in c(cont,NEG))Z[[v]]<-if(v=="jeju_dummy")(d0[[v]]-mean(d0[[v]]))/sd(d0[[v]]) else v2_zb(d0[[v]],v)   # ★v2: 평균 채움 금지(결측이면 멈춤)
# VIF
vd<-Z[,cont];vd$.y<-rnorm(nrow(vd));lmf<-lm(as.formula(paste(".y~",paste(sprintf("`%s`",cont),collapse="+"))),data=vd);v<-car::vif(lmf);cat(sprintf("VIF 최대 %.2f (%s)\n",max(v),names(which.max(v))))
if(nzchar(Sys.getenv("EHEC_HM"))) write.csv(data.frame(cov=names(v),vif=round(as.numeric(v),2),forced=names(v)%in%DOM[["Forced controls"]]),file.path(Sys.getenv("EHEC_RS"),"EHEC_main_vif.csv"),row.names=FALSE)   # 2026-10-06
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
## 2026-10-21 재발 방지: final.dbf 인천 region 라벨(남구 23030 = 미추홀구 약 25 km², 남동구 23050 약 57 km²)이 면적과 맞는지 단언(2026-10-21 dbf 두 칸 교정, 백업 final.dbf.bak_pre_shpfix)
local({a<-as.numeric(sf::st_area(shp))/1e6;stopifnot(nrow(shp)==229,!anyDuplicated(shp$region),a[shp$region=="인천시남구"]<35,a[shp$region=="인천시남동구"]>45)})
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")];yy<-2011:2024   # ★v2: 코드 결합·게이트
DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
DD<-v2_finish(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD)   # ★v2: 채움 없이 단언 + 청주 마스킹. ★2026-10-21 v2b: v2 수정 때 time·id 대입이 주석 안으로 들어가 적합이 멈추던 것을 바로잡음(이전 판 .bak_pre_v2b)
fitmod<-function(cv){form<-as.formula(paste0("cases~",paste(sprintf("`%s`",cv),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(form,family=Sys.getenv("EHEC_FAM",unset="nbinomial"),data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
m<-fitmod(cont);fx<-m$summary.fixed
res<-do.call(rbind,lapply(names(DOM),function(dm)do.call(rbind,lapply(DOM[[dm]],function(v){r<-fx[v,];pp<-1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]])
  data.frame(domain=dm,cov=v,label=ENG[v],IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(pp,3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}))))
write.csv(res,file.path(RS,"MASTER_final_table2.csv"),row.names=FALSE)
cat(sprintf("\n===== 최종 주모형 %d변수 · DIC %.1f · WAIC %.1f · phi %.3f · EPV %.0f =====\n",length(cont),m$dic$dic,m$waic$waic,m$summary.hyperpar["Phi for area","mean"],sum(DD$cases,na.rm=TRUE)/length(cont)))
for(i in 1:nrow(res))cat(sprintf("  %-36s IRR %.3f (%.3f-%.3f)  P>1=%.3f%s\n",res$label[i],res$IRR[i],res$lo[i],res$hi[i],res$P_gt1[i],ifelse(res$cred[i]," *","")))
writeLines(c(sprintf("phi=%.3f",m$summary.hyperpar["Phi for area","mean"]),sprintf("DIC=%.1f",m$dic$dic),sprintf("WAIC=%.1f",m$waic$waic),sprintf("nvar=%d",length(cont)),sprintf("EPV=%.0f",sum(DD$cases,na.rm=TRUE)/length(cont))),file.path(RS,"final_fit.txt"))
# 지도용
uu<-m$summary.random$area[1:length(regs),];rr<-data.frame(region=regs,RR=round(exp(uu$mean),3),RR_lo=round(exp(uu$`0.025quant`),3),RR_hi=round(exp(uu$`0.975quant`),3));rr$hi_risk<-as.integer(rr$RR_lo>1)
write.csv(rr,file.path(RS,"MASTER_agestd_regionSIR.csv"),row.names=FALSE)
pep<-sapply(m$marginals.random$area[1:length(regs)],function(mg)1-INLA::inla.pmarginal(0,mg));write.csv(data.frame(region=regs,pep=round(pep,4)),file.path(RS,"MASTER_agestd_exceedance.csv"),row.names=FALSE)
cat(sprintf("\n지도용: 고위험 %d개 구 · 저장 완료\n",sum(rr$hi_risk)))
