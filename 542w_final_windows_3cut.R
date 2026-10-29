# 542w — EHEC 코로나 정렬 3창(2011-15/16-19/20-24) 순수공간 BYM2 — 542 블록 ① 대체 (2026-09-09 교수 지시). 산출 파일명은 542와 동일(MASTER_final_4yr_split/fitstat).
# 542 — 최종 주모형(20변수) 기준 기간 분석: ① 4년창 순수공간 ② 코로나 전후 분할 + 상호작용 ③ 주요 3변수 × 시기 상호작용
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-Sys.getenv("EHEC_RS",unset=file.path(ROOT,"results_ST"))
V2_CACHE<-CACHE;source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지, KDCA형 청주 2011-2013 마스킹(CJ_MASK=0이면 끔). 이전 판 .bak_pre_v2
T2<-read.csv(file.path(RS,"MASTER_final_table2.csv"));cont<-T2$cov;ENG<-setNames(T2$label,T2$cov)
d0<-read.csv(file.path(CACHE,Sys.getenv("EHEC_PANEL",unset="ST_prereg_panel.csv")),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")];for(v in cont){x<-d0[[v]];x<-v2_nofill(x,v);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")]
mk<-function(yy){DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
  DD<-v2_finish(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD);DD}
grab<-function(m,v)data.frame(cov=v,label=ENG[v],IRR=round(exp(m$summary.fixed[v,"mean"]),3),lo=round(exp(m$summary.fixed[v,"0.025quant"]),3),hi=round(exp(m$summary.fixed[v,"0.975quant"]),3),cred=(m$summary.fixed[v,"0.025quant"]>0|m$summary.fixed[v,"0.975quant"]<0))
KEY<-c("beef_cattle_farms","fiscal_independence","net_migration","current_smoking");if(nzchar(Sys.getenv("EHEC_HM")))KEY<-setdiff(cont,c("kchs_pc1","urban_population_share","physicians_per_1000","sex_ratio"))   # 2026-09-30 문헌 풀: 최종 세트 가설 변수+재정자주도
if(!nzchar(Sys.getenv("EHEC_HM")))v2_key_check(KEY,file.path(RS,"MASTER_final_table2.csv"),RS,"542w_final_windows_3cut.R")   # ★2026-10-21 v2b(이전 판 .bak_pre_v2b): 원고 모드 KEY 4개는 원고의 고정 목록이라 그대로 둠. 이번 540 결과의 1 배제 판정이 원고(9/3 Table 2)와 달라진 변수는 results_ST/MASTER_final_KEY_check.csv와 콘솔 경고로 표시(KEY 자체는 바꾸지 않음)
# ── ① 4년창 순수공간 ──
cat("═══ ① 4년창 (순수 공간 BYM2) ═══\n");W<-list("2011-2015"=2011:2015,"2016-2019"=2016:2019,"2020-2024"=2020:2024);split<-list();fs<-list()
for(w in names(W)){DD<-mk(W[[w]]);form<-as.formula(paste0("cases~",paste(sprintf("`%s`",cont),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)"))
  m<-inla(form,family=Sys.getenv("EHEC_FAM",unset="nbinomial"),data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE),control.inla=list(control.vb=list(emergency=50)))
  r<-do.call(rbind,lapply(cont,function(v)grab(m,v)));r$period<-w;r$N<-sum(DD$cases,na.rm=TRUE);split[[w]]<-r;fs[[w]]<-data.frame(period=w,N=sum(DD$cases,na.rm=TRUE),DIC=round(m$dic$dic,1),phi=round(m$summary.hyperpar["Phi for area","mean"],3))
  cat(sprintf("  [%s] N=%d · credible: %s\n",w,sum(DD$cases,na.rm=TRUE),paste(r$label[r$cred],collapse=", ")))}
write.csv(do.call(rbind,split),file.path(RS,"MASTER_final_4yr_split.csv"),row.names=FALSE);write.csv(do.call(rbind,fs),file.path(RS,"MASTER_final_4yr_fitstat.csv"),row.names=FALSE)
cat('done 3cut\n')
