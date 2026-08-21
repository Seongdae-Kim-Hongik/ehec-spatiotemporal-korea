# 415 — EHEC 연령표준화 초과확률 지도 P(spatial RR>1) (보충용, 레퍼런스 JMIR 방식)
#  410 M6 재적합 후 area 사후 주변분포에서 P(u_i>0)=P(RR>1) 추출. fig_style bin_pep/PAL_PEP5.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA)})
set.seed(20260711); inla.setOption(num.threads="1:1")
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-file.path(GD,"FBD_SpatioTemporal");CACHE<-file.path(ROOT,"data_cache");EH<-file.path(ROOT,"EHEC")
source(file.path(GD,"EHP_Spatio/fig_style.R"))
cov<-c("a_septic","b_gw_use","b_gw_home","b_gw_unfit","d_pasture","d_field","d_paddy","e_rural_pop","e_forest","h_fastfood","h_pub","g_fiscal","g_basiclivelihood")
CV<-readRDS(file.path(CACHE,"ST_covariates_panel.rds"));cl<-readRDS(file.path(CACHE,"climate_covariates.rds"))
meltm<-function(M,nm){d<-as.data.frame(M);d$region<-rownames(M);pivot_longer(d,-region,names_to="year",values_to=nm)%>%mutate(year=as.integer(year))}
clim<-meltm(cl$temp,"temp")%>%left_join(meltm(cl$precip,"precip"),by=c("region","year"))%>%left_join(meltm(cl$sst,"sst"),by=c("region","year"))
base<-data.frame(region=CV$region,year=CV$year,CV$raw,check.names=FALSE)%>%left_join(clim,by=c("region","year"))
skewf<-function(x){x<-x[is.finite(x)];if(sd(x)==0)return(0);e1071::skewness(x)}
Z<-sapply(cov,function(v){x<-base[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;z});Z<-as.data.frame(Z);Z$region<-base$region;Z$year<-base$year
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region;NR<-nrow(shp)
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"));yy<-2011:2024
d0<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
d0$cases[is.na(d0$cases)]<-0;d0$E_std[is.na(d0$E_std)|d0$E_std<=0]<-min(d0$E_std[d0$E_std>0],na.rm=TRUE)
d0$area<-match(d0$region,regs);d0$time<-match(d0$year,yy);d0$id<-1:nrow(d0)
rhs<-paste(sprintf("`%s`",cov),collapse="+")
form<-as.formula(paste0("cases~",rhs,"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
m<-inla(form,family="nbinomial",data=d0,E=d0$E_std,control.fixed=CF,control.compute=list(config=TRUE),control.inla=list(control.vb=list(emergency=50)))
# BYM2 area: 첫 NR개 = 총 공간효과(u+v). P(효과>0)=P(RR>1)
pep<-sapply(1:NR,function(i) 1-inla.pmarginal(0,m$marginals.random$area[[i]]))
out<-data.frame(region=regs,pep=round(pep,4))
write.csv(out,file.path(ROOT,"results_ST/MASTER_agestd_exceedance.csv"),row.names=FALSE)
cat(sprintf("초과확률: P>0.95(hot) %d개 · P<0.05(cold) %d개\n",sum(pep>0.95),sum(pep<0.05)))
# ---- 지도 ----
sp<-left_join(shp[,c("region","geometry")],out,by="region");sp$cls<-bin_pep(sp$pep)
PLV<-c("<0.05 (cold)","0.05–0.5","0.5–0.8","0.8–0.95","≥0.95 (hot)")
gg<-ggplot(sp)+geom_sf(aes(fill=cls),color="grey78",linewidth=0.03)+lvl_layer(sp,PLV)+layer_sido(BASE)+
  scale_fill_manual(values=setNames(PAL_PEP5,PLV),breaks=PLV,limits=PLV,na.value="grey90",name="P(spatial RR > 1)",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+
  labs(title="EHEC posterior exceedance probability P(spatial RR > 1)",
       subtitle="Age-standardised BYM2 structured spatial effect. ≥0.95 = credible high-risk; <0.05 = credible low-risk. Korea 2000/UTM-K.")+
  theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
AGG<-if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else NULL
for(p in c(file.path(ROOT,"results_ST/Fig_agestd_exceedance.png"),file.path(EH,"EHEC_Fig_agestd_exceedance.png")))
  ggsave(p,gg,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)
cat("저장: Fig_agestd_exceedance (results_ST + EHEC)\n")
