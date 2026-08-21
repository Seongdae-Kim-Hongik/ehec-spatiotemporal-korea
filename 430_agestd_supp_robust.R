# 430 — EHEC 연령표준화 Supplementary 모형강건성: M1-M6 사다리 · 8-graph · NB vs ZINB
#  전부 offset=E_std, 결과변수=KDCA 신고사례, 13공변량(410과 동일). 산출 results_ST/MASTER_agestd_supp_*.csv
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260710); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
KLAB<-c(a_septic="정화조인구%",b_gw_use="지하수이용",b_gw_home="가정지하수",b_gw_unfit="지하수부적합%",d_paddy="논면적",d_field="밭면적",d_pasture="목장면적",e_rural_pop="농촌인구%",e_forest="임야",h_fastfood="패스트푸드",h_pub="주점수",g_fiscal="재정자립",g_basiclivelihood="기초생활%")
cov<-c("a_septic","b_gw_use","b_gw_home","b_gw_unfit","d_pasture","d_field","d_paddy","e_rural_pop","e_forest","h_fastfood","h_pub","g_fiscal","g_basiclivelihood")
# ---- Z ----
CV<-readRDS(file.path(CACHE,"ST_covariates_panel.rds"));cl<-readRDS(file.path(CACHE,"climate_covariates.rds"))
meltm<-function(M,nm){d<-as.data.frame(M);d$region<-rownames(M);pivot_longer(d,-region,names_to="year",values_to=nm)%>%mutate(year=as.integer(year))}
clim<-meltm(cl$temp,"temp")%>%left_join(meltm(cl$precip,"precip"),by=c("region","year"))%>%left_join(meltm(cl$sst,"sst"),by=c("region","year"))
base<-data.frame(region=CV$region,year=CV$year,CV$raw,check.names=FALSE)%>%left_join(clim,by=c("region","year"))
skewf<-function(x){x<-x[is.finite(x)];if(sd(x)==0)return(0);e1071::skewness(x)}
Z<-sapply(cov,function(v){x<-base[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;z});Z<-as.data.frame(Z);Z$region<-base$region;Z$year<-base$year
# ---- shp/graphs ----
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region;NR<-nrow(shp)
fixnb<-function(nb){k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))};cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};make.sym.nb(nb)}
mkgraph<-function(nb){g<-tempfile(fileext=".graph");nb2INLA(g,nb);g}
nbQ<-fixnb(poly2nb(shp,queen=TRUE,snap=0.01)); gQ<-mkgraph(nbQ)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
# ---- data ----
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"));yy<-2011:2024
d0<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
d0$cases[is.na(d0$cases)]<-0;d0$E_std[is.na(d0$E_std)|d0$E_std<=0]<-min(d0$E_std[d0$E_std>0],na.rm=TRUE)
d0$area<-match(d0$region,regs);d0$time<-match(d0$year,yy);d0$id<-1:nrow(d0)
rhs<-paste(sprintf("`%s`",cov),collapse="+")
runinla<-function(form,fam="nbinomial"){inla(form,family=fam,data=d0,E=d0$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}

# ===== M1-M6 사다리 =====
L<-list(
 M1=paste0("cases~",rhs,"+f(area,model='iid',hyper=pcp)"),
 M2=paste0("cases~",rhs,"+f(area,model='bym2',graph=gQ,scale.model=TRUE,constr=TRUE,hyper=pc)"),
 M3=paste0("cases~",rhs,"+f(time,model='rw1',scale.model=TRUE,hyper=pcp)"),
 M4=paste0("cases~",rhs,"+f(area,model='bym2',graph=gQ,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='iid',hyper=pcp)"),
 M5=paste0("cases~",rhs,"+f(area,model='bym2',graph=gQ,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)"),
 M6=paste0("cases~",rhs,"+f(area,model='bym2',graph=gQ,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
desc<-c(M1="Covariates + iid area",M2="+ BYM2 spatial",M3="+ RW1 temporal (no spatial)",M4="BYM2 + iid time",M5="BYM2 + RW1",M6="BYM2 + RW1 + Type I interaction")
lad<-lapply(names(L),function(k){m<-runinla(as.formula(L[[k]]));data.frame(model=k,structure=desc[k],DIC=round(m$dic$dic,1),WAIC=round(m$waic$waic,1))})
lad<-bind_rows(lad);lad$dDIC<-round(lad$DIC-min(lad$DIC),1)
write.csv(lad,file.path(ROOT,"results_ST/MASTER_agestd_supp_ladder.csv"),row.names=FALSE)
cat("[M1-M6]\n");print(lad,row.names=FALSE)

# ===== NB vs ZINB =====
f6<-as.formula(L$M6)
mnb<-runinla(f6,"nbinomial");mzi<-tryCatch(runinla(f6,"zeroinflatednbinomial1"),error=function(e)NULL)
cmp<-data.frame(family=c("NB","ZINB"),WAIC=c(round(mnb$waic$waic,1),if(is.null(mzi))NA else round(mzi$waic$waic,1)),DIC=c(round(mnb$dic$dic,1),if(is.null(mzi))NA else round(mzi$dic$dic,1)))
cmp$zero_pct<-round(mean(d0$cases==0)*100,1)
write.csv(cmp,file.path(ROOT,"results_ST/MASTER_agestd_supp_zinb.csv"),row.names=FALSE)
cat("\n[NB vs ZINB] zero%=",cmp$zero_pct[1],"\n");print(cmp,row.names=FALSE)

# ===== 8-graph 강건성 =====
graphs<-list(Queen=nbQ,Rook=fixnb(poly2nb(shp,queen=FALSE,snap=0.01)))
for(k in c(2,3,4,5,6)) graphs[[paste0("knn",k)]]<-make.sym.nb(knn2nb(knearneigh(coords,k=k)))
res8<-list()
for(gn in names(graphs)){gg<-mkgraph(graphs[[gn]])
  fm<-as.formula(paste0("cases~",rhs,"+f(area,model='bym2',graph=gg,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  m<-runinla(fm);fx<-m$summary.fixed[cov,]
  res8[[gn]]<-data.frame(graph=gn,cov=cov,credible=(fx$`0.025quant`>0|fx$`0.975quant`<0))}
r8<-bind_rows(res8)%>%group_by(cov)%>%summarise(n_credible=sum(credible),.groups="drop")%>%mutate(공변량=ifelse(cov%in%names(KLAB),KLAB[cov],cov))%>%arrange(desc(n_credible))
write.csv(r8,file.path(ROOT,"results_ST/MASTER_agestd_supp_8graph.csv"),row.names=FALSE)
cat("\n[8-graph] 신뢰횟수(/8)\n");print(as.data.frame(r8),row.names=FALSE)
cat("\n===== 430 완료 =====\n")
