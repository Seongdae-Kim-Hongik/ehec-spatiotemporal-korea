# 540cj — 청주 2011-2014 제외 민감도(DROP_CHEONGJU=1). 저자 지시 2026-09-28
# ★2026-09-29: 민감도 실행이 주모형 산출(적합통계·음성대조·지도용 RR·초과확률)을 같은 이름으로 덮어쓰던 결함 수정 → _noCheongju 접미사
# 540 — 최종 주모형 (사전지정 풀 멀티버스 v2 이후). 사전지정 도메인 대표(각 도메인 상위 credible 1~2개) + 조사지표 종합축 + 성비 + 검출통제(의사/천명).
#  산출: Table 2(MASTER_final_table2.csv) · 지도용 regionSIR/exceedance · 적합통계 · 음성대조
#  ★2026-10-21 v2(이전 판 .bak_pre_v2, _v2_common.R): E·사례·공변량 sigungu_cd 결합, E 최솟값 채움·사례 0 채움 삭제(panel_gate), 공변량 평균 채움 금지(결측 있으면 멈춤),
#   KDCA형 청주 2011-2013 사례 NA 마스킹이 표준(CJ_MASK=0이면 관측 0 사용 = 민감도). 패널은 v2 정본(shp 면적·v2 인구·명시 채움, data_cache).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(car)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-file.path(ROOT,"results_ST")
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
 current_smoking="Current smoking",influenza_vaccination_coverage="Influenza vaccination coverage",kchs_pc1="KCHS survey composite (PC1)",physicians_per_1000_pop="Physicians per 1000 population")
NEG<-c("swine_density","poultry_density","walking_practice")
cont<-unlist(DOM,use.names=FALSE)
d0<-read.csv(file.path(CACHE,"ST_prereg_panel.csv"),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")];for(v in c(cont,NEG))Z[[v]]<-v2_zb(d0[[v]],v)   # ★v2: 평균 채움 금지(결측이면 멈춤)
# VIF
vd<-Z[,cont];vd$.y<-rnorm(nrow(vd));lmf<-lm(as.formula(paste(".y~",paste(sprintf("`%s`",cont),collapse="+"))),data=vd);v<-car::vif(lmf);cat(sprintf("VIF 최대 %.2f (%s)\n",max(v),names(which.max(v))))
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
## 2026-10-21 재발 방지: final.dbf 인천 region 라벨(남구 23030 = 미추홀구 약 25 km², 남동구 23050 약 57 km²)이 면적과 맞는지 단언(2026-10-21 dbf 두 칸 교정, 백업 final.dbf.bak_pre_shpfix)
local({a<-as.numeric(sf::st_area(shp))/1e6;stopifnot(nrow(shp)==229,!anyDuplicated(shp$region),a[shp$region=="인천시남구"]<35,a[shp$region=="인천시남동구"]>45)})
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")];yy<-2011:2024   # ★v2: 코드 결합·게이트
DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
DD<-v2_finish(DD,kdca=TRUE)   # ★v2: 채움 없이 단언 + 청주 2011-2013 마스킹(표준). 아래 DROP_CHEONGJU=1은 2011-2014 행 제거 민감도(v2에서는 청주 기준이 반대 — 관측 0 민감도는 CJ_MASK=0)
# ★청주 민감도(저자 지시 2026-09-28): KDCA 시군구×연령 교차표에 청주가 2014 이전 부재 → 2011-2013 지구-연을 제외하고 1회 적합
CJ_DROP<-Sys.getenv("DROP_CHEONGJU")=="1"
if(CJ_DROP){n0<-nrow(DD);DD<-DD[!(DD$region=="충청북도청주시"&DD$year%in%2011:2014),];cat(sprintf("[민감도] 청주 2011-2014 제외: %d→%d행\n",n0,nrow(DD)))}
DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD)
fitmod<-function(cv){form<-as.formula(paste0("cases~",paste(sprintf("`%s`",cv),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(form,family="nbinomial",data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
m<-fitmod(cont);fx<-m$summary.fixed
res<-do.call(rbind,lapply(names(DOM),function(dm)do.call(rbind,lapply(DOM[[dm]],function(v){r<-fx[v,];pp<-1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]])
  data.frame(domain=dm,cov=v,label=ENG[v],IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(pp,3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}))))
write.csv(res,file.path(RS,if(CJ_DROP)"MASTER_final_table2_noCheongju.csv" else "MASTER_final_table2.csv"),row.names=FALSE)
cat(sprintf("\n===== 최종 주모형 %d변수 · DIC %.1f · WAIC %.1f · phi %.3f · EPV %.0f =====\n",length(cont),m$dic$dic,m$waic$waic,m$summary.hyperpar["Phi for area","mean"],sum(DD$cases,na.rm=TRUE)/length(cont)))
for(i in 1:nrow(res))cat(sprintf("  %-36s IRR %.3f (%.3f-%.3f)  P>1=%.3f%s\n",res$label[i],res$IRR[i],res$lo[i],res$hi[i],res$P_gt1[i],ifelse(res$cred[i]," *","")))
writeLines(c(sprintf("phi=%.3f",m$summary.hyperpar["Phi for area","mean"]),sprintf("DIC=%.1f",m$dic$dic),sprintf("WAIC=%.1f",m$waic$waic),sprintf("nvar=%d",length(cont)),sprintf("EPV=%.0f",sum(DD$cases,na.rm=TRUE)/length(cont))),file.path(RS,if(CJ_DROP)"final_fit_noCheongju.txt" else "final_fit.txt"))
# 음성대조
cat("\n[음성대조 개별 추가]\n");nc<-do.call(rbind,lapply(NEG,function(v){mm<-fitmod(c(cont,v));r<-mm$summary.fixed[v,];data.frame(cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}))
for(i in 1:nrow(nc))cat(sprintf("  %-18s IRR %.3f (%.3f-%.3f)%s\n",nc$cov[i],nc$IRR[i],nc$lo[i],nc$hi[i],ifelse(nc$cred[i]," * ← 교란 의심","")));write.csv(nc,file.path(RS,if(CJ_DROP)"MASTER_final_negcontrol_noCheongju.csv" else "MASTER_final_negcontrol.csv"),row.names=FALSE)
# 지도용
uu<-m$summary.random$area[1:length(regs),];rr<-data.frame(region=regs,RR=round(exp(uu$mean),3),RR_lo=round(exp(uu$`0.025quant`),3),RR_hi=round(exp(uu$`0.975quant`),3));rr$hi_risk<-as.integer(rr$RR_lo>1)
write.csv(rr,file.path(RS,if(CJ_DROP)"MASTER_agestd_regionSIR_noCheongju.csv" else "MASTER_agestd_regionSIR.csv"),row.names=FALSE)
pep<-sapply(m$marginals.random$area[1:length(regs)],function(mg)1-INLA::inla.pmarginal(0,mg));write.csv(data.frame(region=regs,pep=round(pep,4)),file.path(RS,if(CJ_DROP)"MASTER_agestd_exceedance_noCheongju.csv" else "MASTER_agestd_exceedance.csv"),row.names=FALSE)
cat(sprintf("\n지도용: 고위험 %d개 구 · 저장 완료\n",sum(rr$hi_risk)))
