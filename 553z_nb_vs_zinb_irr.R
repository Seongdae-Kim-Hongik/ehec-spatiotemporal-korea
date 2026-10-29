# 553z — 최종 20변수 시공간 모형 NB vs ZINB: 핵심 IRR·기대 영셀·적합도 (2026-09-14 외부 검토 대응). 헤더 = 553 1-22행.
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-Sys.getenv("EHEC_RS",unset=file.path(ROOT,"results_ST"))
V2_CACHE<-CACHE;source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지, KDCA형 청주 2011-2013 마스킹(CJ_MASK=0이면 끔). 이전 판 .bak_pre_v2
T2<-read.csv(file.path(RS,"MASTER_final_table2.csv"));cont<-T2$cov
d0<-read.csv(file.path(CACHE,Sys.getenv("EHEC_PANEL",unset="ST_prereg_panel.csv")),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")];for(v in cont){x<-d0[[v]];x<-v2_nofill(x,v);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")];yy<-2011:2024
DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
DD<-v2_finish(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD)
FX<-paste(sprintf("`%s`",cont),collapse="+")
fitf<-function(re,fam="nbinomial"){form<-as.formula(paste0("cases~",FX,if(nchar(re))paste0("+",re) else ""))
  inla(form,family=fam,data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
SP<-"f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)";TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(id,model='iid',hyper=pcp)"
full<-paste(SP,TI,IT,sep="+")
nbm<-fitf(full,"nbinomial");zim<-fitf(full,"zeroinflatednbinomial1")
fx<-function(m)data.frame(cov=cont,IRR=round(exp(m$summary.fixed[cont,"mean"]),3),lo=round(exp(m$summary.fixed[cont,"0.025quant"]),3),hi=round(exp(m$summary.fixed[cont,"0.975quant"]),3),cred=(m$summary.fixed[cont,"0.025quant"]>0|m$summary.fixed[cont,"0.975quant"]<0))
R<-merge(fx(nbm),fx(zim),by="cov",suffixes=c("_NB","_ZINB"));R$label<-T2$label[match(R$cov,T2$cov)];R<-R[match(cont,R$cov),]
hp<-function(m,pat){h<-m$summary.hyperpar;i<-grep(pat,rownames(h),ignore.case=TRUE);if(length(i))h$mean[i[1]] else NA}
zexp<-function(m,zi){mu<-DD$E_std*exp(m$summary.linear.predictor$mean);s<-hp(m,"size");p0<-(s/(s+mu))^s;if(zi){pz<-hp(m,"zero-prob");p0<-pz+(1-pz)*p0};sum(p0[!is.na(DD$cases)])}   # 서술 집계는 우도에 들어간 칸만(청주 NA 제외)
zz<-data.frame(model=c("NB","ZINB"),cells=sum(!is.na(DD$cases)),obs_zero=sum(DD$cases==0,na.rm=TRUE),exp_zero=round(c(zexp(nbm,FALSE),zexp(zim,TRUE)),1),
  size=round(c(hp(nbm,"size"),hp(zim,"size")),3),zero_prob=c(NA,round(hp(zim,"zero-prob"),3)),WAIC=round(c(nbm$waic$waic,zim$waic$waic),1),DIC=round(c(nbm$dic$dic,zim$dic$dic),1))
write.csv(R,file.path(RS,"MASTER_supp_nb_vs_zinb_irr.csv"),row.names=FALSE);write.csv(zz,file.path(RS,"MASTER_supp_nb_vs_zinb_fit.csv"),row.names=FALSE)
print(R);print(zz);cat("done 553z\n")
