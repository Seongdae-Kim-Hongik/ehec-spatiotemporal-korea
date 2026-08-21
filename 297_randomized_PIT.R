# ============================================================================
# 297 (ⓑ) — Randomized PIT (Czado 2009) — 희소 이산자료 교정 재검증
#  비무작위 PIT는 count가 대부분 0일 때 인공적으로 1에 쏠림 → 오해. 무작위화로 교정.
#  rPIT = Fu*(y-1) + v*[Fu*(y)-Fu*(y-1)], v~U(0,1). 잘 교정=평평.
#  산출: results_ST/Fig_randomizedPIT_15diseases.png
# ============================================================================
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(ggplot2)})
EN_DZ<-c("노로바이러스"="Norovirus","A형간염"="Hepatitis A","E형간염"="Hepatitis E","황색포도상구균"="Staphylococcus aureus","바실러스세레우스"="Bacillus cereus","클로스트리디움퍼프린젠스"="Clostridium perfringens","비브리오패혈증"="Vibrio vulnificus sepsis","장염비브리오"="Vibrio parahaemolyticus","장티푸스"="Typhoid fever","세균성이질"="Shigellosis","살모넬라"="Salmonellosis","병원성대장균"="Pathogenic E. coli","장출혈성대장균감염증"="EHEC","캠필로박터제주니"="Campylobacter","원충"="Kudoa");endz<-function(k){v<-unname(EN_DZ[k]);ifelse(is.na(v),k,v)}
set.seed(20260707); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));NR<-nrow(shp);regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
panel<-readRDS(file.path(CACHE,"ST_panel_15x229x14.rds"))|>as.data.frame()
DZ<-unique(panel$disease);pits<-list()
for(d in DZ){pd<-panel[panel$disease==d,];fy<-unique(pd$fy);yy<-(2011:2024)[(2011:2024)>=fy]
  g2<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(pd[,c("region","year","cases","pop")],by=c("region","year"))
  g2$cases[is.na(g2$cases)]<-0;g2$pop[is.na(g2$pop)|g2$pop<=0]<-median(g2$pop,na.rm=TRUE)
  g2$E<-pmax(g2$pop*sum(g2$cases)/sum(g2$pop),1);g2$area<-match(g2$region,regs);g2$time<-match(g2$year,yy);g2$id<-1:nrow(g2)
  frm<-cases~1+f(area,model="bym2",graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model="rw1",scale.model=TRUE,hyper=pcp)+f(id,model="iid",hyper=pcp)
  m<-tryCatch(inla(frm,family="nbinomial",data=g2,E=g2$E,control.compute=list(config=TRUE,dic=FALSE),control.predictor=list(compute=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
  if(is.null(m)){cat("ERR",d,"\n");next}
  # 사후 fitted 평균률 → NB 파라미터로 rPIT (size=사후 과산포)
  mu<-m$summary.fitted.values$mean[1:nrow(g2)]*g2$E
  size<-tryCatch(m$summary.hyperpar["size for the nbinomial observations (1/overdispersion)","mean"],error=function(e)10)
  if(!is.finite(size))size<-10
  y<-g2$cases;Fy<-pnbinom(y,size=size,mu=mu);Fy1<-pnbinom(y-1,size=size,mu=mu)
  v<-runif(length(y));rpit<-Fy1+v*(Fy-Fy1);rpit<-pmin(pmax(rpit,0),1)
  pits[[d]]<-data.frame(disease=d,rpit=rpit);cat(sprintf("[%s] rPIT n=%d size=%.1f\n",d,length(rpit),size))}
P<-bind_rows(pits);P$disease<-factor(P$disease,levels=DZ,labels=endz(DZ))
gg<-ggplot(P,aes(rpit))+geom_histogram(bins=12,fill="#238B45",color="white",linewidth=0.2,boundary=0)+
  facet_wrap(~disease,scales="free_y",ncol=4)+labs(title="Randomized PIT (Czado 2009) — calibration for sparse discrete counts",subtitle="Removes the discreteness artefact of non-randomised PIT. Flat = well calibrated / skew = misfit. NB-BYM2-RW1.",x="Randomized PIT",y="Frequency")+
  theme_bw(base_size=9)+theme(strip.text=element_text(size=7.5),plot.title=element_text(face="bold"),text=element_text(family="Helvetica"))
ggsave(file.path(ROOT,"results_ST/Fig_randomizedPIT_15diseases.png"),gg,width=11,height=8,dpi=200,bg="white")
cat("저장: results_ST/Fig_randomizedPIT_15diseases.png\n")
