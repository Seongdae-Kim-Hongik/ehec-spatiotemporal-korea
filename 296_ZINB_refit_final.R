# ============================================================================
# 296 (ⓐ) — 세균성이질·EHEC ZINB 최종 재적합 (202 결과 갱신, 제약·초과확률·검증 유지)
#  292서 ZINB 우세 확인된 두 질병만. Type I-IV(generic0 Kronecker+extraconstr) ZINB WAIC 선택.
#  산출: results_ST/ZINB_final_refit.csv + EXC(고위험) + CPO/PIT.
# ============================================================================
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(Matrix)})
set.seed(20260707); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
DZ<-c("세균성이질","장출혈성대장균감염증");FAM<-"zeroinflatednbinomial1"
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));NR<-nrow(shp);regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
Wmat<-nb2mat(nb,style="B",zero.policy=TRUE);Ra<-inla.scale.model(Diagonal(x=rowSums(Wmat))-Matrix(Wmat,sparse=TRUE),constr=list(A=matrix(1,1,NR),e=0))
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
panel<-readRDS(file.path(CACHE,"ST_panel_15x229x14.rds"))|>as.data.frame()
rows<-list()
for(d in DZ){pd<-panel[panel$disease==d,];fy<-unique(pd$fy);yy<-(2011:2024)[(2011:2024)>=fy];nyl<-length(yy)
  g2<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(pd[,c("region","year","cases","pop")],by=c("region","year"))
  g2$cases[is.na(g2$cases)]<-0;g2$pop[is.na(g2$pop)|g2$pop<=0]<-median(g2$pop,na.rm=TRUE)
  g2$E<-pmax(g2$pop*sum(g2$cases)/sum(g2$pop),1);g2$area<-match(g2$region,regs);g2$time<-match(g2$year,yy);g2$id<-1:nrow(g2)
  I_<-NR;T_<-nyl;NIT<-I_*T_;g2$st<-(g2$time-1)*I_+g2$area
  Rt<-inla.scale.model(crossprod(diff(Diagonal(T_))),constr=list(A=matrix(1,1,T_),e=0))
  Q_II<-kronecker(Rt,Diagonal(I_));Q_III<-kronecker(Diagonal(T_),Ra);Q_IV<-kronecker(Rt,Ra)
  A_t<-matrix(0,I_,NIT);for(a in 1:I_)A_t[a,((1:T_)-1)*I_+a]<-1;A_s<-matrix(0,T_,NIT);for(t in 1:T_)A_s[t,(t-1)*I_+(1:I_)]<-1
  ec_II<-list(A=A_t,e=rep(0,I_));ec_III<-list(A=A_s,e=rep(0,T_));ec_IV<-list(A=rbind(A_t,A_s),e=rep(0,I_+T_))
  base<-"cases~1+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)"
  forms<-list(I=paste0(base,"+f(id,model='iid',hyper=pcp)"),
    II=paste0(base,"+f(st,model='generic0',Cmatrix=Q_II,rankdef=",I_,",constr=FALSE,hyper=pcp,extraconstr=ec_II)"),
    III=paste0(base,"+f(st,model='generic0',Cmatrix=Q_III,rankdef=",T_,",constr=FALSE,hyper=pcp,extraconstr=ec_III)"),
    IV=paste0(base,"+f(st,model='generic0',Cmatrix=Q_IV,rankdef=",I_+T_-1,",constr=FALSE,hyper=pcp,extraconstr=ec_IV)"))
  best<-NULL;bw<-Inf;bt<-NA;waics<-c()
  for(tp in names(forms)){f<-tryCatch(inla(as.formula(forms[[tp]]),family=FAM,data=g2,E=g2$E,
    control.compute=list(waic=TRUE,cpo=TRUE,return.marginals.predictor=TRUE),control.predictor=list(compute=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
    w<-if(is.null(f))NA else f$waic$waic;waics[tp]<-round(w,1);if(!is.null(f)&&is.finite(w)&&w<bw){bw<-w;best<-f;bt<-tp}}
  if(is.null(best)){cat("FAIL",d,"\n");next}
  mf<-best$marginals.fitted.values;exc<-if(!is.null(mf))sapply(seq_len(nrow(g2)),function(r)1-INLA::inla.pmarginal(1,mf[[r]]))else NA
  rows[[d]]<-data.frame(disease=d,family="ZINB",best_type=bt,waic_I=waics["I"],waic_II=waics["II"],waic_III=waics["III"],waic_IV=waics["IV"],
    n_high=sum(exc>0.95,na.rm=TRUE),cpo_fail=sum(best$cpo$failure>0,na.rm=TRUE))
  cat(sprintf("[%s ZINB] best=%s n_high=%d cpo_fail=%d\n",d,bt,sum(exc>0.95,na.rm=TRUE),sum(best$cpo$failure>0,na.rm=TRUE)))}
M<-bind_rows(rows);write.csv(M,file.path(ROOT,"results_ST/ZINB_final_refit.csv"),row.names=FALSE)
cat("\n=== ZINB 최종 재적합 ===\n");print(as.data.frame(M),row.names=FALSE)
