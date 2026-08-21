# ============================================================================
# 292 — 희소질병 ZINB vs NB 비교 (CPO 대량실패 대응; /333 지적)
#  EHEC·비브리오패혈증(202서 CPO실패 3196/1271) + 대조 3개. 시공간 BYM2+RW1(Type I),
#  NB vs zeroinflatednbinomial1을 WAIC/CPO실패로 비교. 산출: results_ST/ZINB_vs_NB.csv
# ============================================================================
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA)})
set.seed(20260707); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
DZ<-c("장출혈성대장균감염증","비브리오패혈증","세균성이질","원충","E형간염")  # 희소·CPO실패 대상+대조
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));NR<-nrow(shp);regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
panel<-readRDS(file.path(CACHE,"ST_panel_15x229x14.rds"))|>as.data.frame()
rows<-list()
for(d in DZ){pd<-panel[panel$disease==d,];fy<-unique(pd$fy);yy<-(2011:2024)[(2011:2024)>=fy]
  g2<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(pd[,c("region","year","cases","pop")],by=c("region","year"))
  g2$cases[is.na(g2$cases)]<-0;g2$pop[is.na(g2$pop)|g2$pop<=0]<-median(g2$pop,na.rm=TRUE)
  g2$E<-pmax(g2$pop*sum(g2$cases)/sum(g2$pop),1);g2$area<-match(g2$region,regs);g2$time<-match(g2$year,yy);g2$id<-1:nrow(g2)
  zero<-mean(g2$cases==0)
  frm<-cases~1+f(area,model="bym2",graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model="rw1",scale.model=TRUE,hyper=pcp)+f(id,model="iid",hyper=pcp)
  fit<-function(fam)tryCatch(inla(frm,family=fam,data=g2,E=g2$E,control.compute=list(waic=TRUE,cpo=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
  mnb<-fit("nbinomial");mzi<-fit("zeroinflatednbinomial1")
  gv<-function(m)if(is.null(m))c(NA,NA)else c(round(m$waic$waic,1),sum(m$cpo$failure>0,na.rm=TRUE))
  a<-gv(mnb);b<-gv(mzi)
  rows[[d]]<-data.frame(disease=d,zero_pct=round(100*zero,1),waic_NB=a[1],cpofail_NB=a[2],waic_ZINB=b[1],cpofail_ZINB=b[2],
    선호=ifelse(is.na(b[1]),"NB",ifelse(b[1]<a[1]-2,"ZINB",ifelse(a[1]<b[1]-2,"NB","무차이"))))
  cat(sprintf("[%s] 0비율%.0f%% NB(waic%.0f,실패%d) ZINB(waic%.0f,실패%d)\n",d,100*zero,a[1],a[2],b[1],b[2]))}
M<-bind_rows(rows);write.csv(M,file.path(ROOT,"results_ST/ZINB_vs_NB.csv"),row.names=FALSE)
cat("\n=== ZINB vs NB ===\n");print(as.data.frame(M),row.names=FALSE)
