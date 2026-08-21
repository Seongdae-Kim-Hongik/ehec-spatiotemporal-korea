# 451 — 연령표준화 모델+S6+지도+초과확률 (일반화). 인자: 질병한글 영문
#  EHEC 410+411+412+415 통합·파라미터화. cov=RICH[dz]-연령구조+SES. offset=E_std.
args<-commandArgs(trailingOnly=TRUE); DZ<-args[1]; ENG<-args[2]
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
source(file.path(GD,"FBD_SpatioTemporal/_st_setup.R"))    # nb,g,regs,NR,coords,build_Z,RICH,pc,pcp,CF,KLAB
source(file.path(GD,"EHP_Spatio/fig_style.R"))            # PAL_DIV6,bin_div,DIV_LABS,guide_fill_box,north_ehp,scale_ehp,layer_sido,bin_pep,PAL_PEP5,theme_ehp_map,lvl_layer,ELAB
suppressPackageStartupMessages({library(dplyr);library(sf);library(ggplot2);library(stringr)})
set.seed(20260711)
R<-file.path(ROOT,"results_ST");C<-CACHE;EH<-file.path(ROOT,ENG);dir.create(EH,showWarnings=FALSE)
AGG<-if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else NULL
# ---- 공변량: RICH-연령구조+SES ----
WATER_covs<-WATER
base_cov<-if(!is.null(RICH[[DZ]])) RICH[[DZ]] else WATER_covs
cov<-setdiff(base_cov,c("g_infant","g_elderly","g_solo_elder"))
cov<-unique(c(cov,"g_fiscal","g_basiclivelihood"))
lab<-function(cv) ifelse(cv%in%names(ELAB),ELAB[cv],ifelse(cv%in%names(KLAB),KLAB[cv],cv))
Z<-build_Z(cov,2011,2024)
# VIF
Xv<-Z[,cov];vif<-sapply(cov,function(v){r2<-tryCatch(summary(lm(Xv[[v]]~.,data=Xv[,setdiff(cov,v)]))$r.squared,error=function(e)NA);1/(1-r2)})
# ---- data ----
case<-read.csv(file.path(C,paste0("ST_",ENG,"_agecase.csv")))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(C,paste0("ST_",ENG,"_Estd.csv")));yy<-2011:2024
d0<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
d0$cases[is.na(d0$cases)]<-0;d0$E_std[is.na(d0$E_std)|d0$E_std<=0]<-min(d0$E_std[d0$E_std>0],na.rm=TRUE)
d0$area<-match(d0$region,regs);d0$time<-match(d0$year,yy);d0$id<-1:nrow(d0)
rhs<-paste(sprintf("`%s`",cov),collapse="+")
fitm<-function(covs,cfg=FALSE){f<-as.formula(paste0("cases~",paste(sprintf("`%s`",covs),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(f,family="nbinomial",data=d0,E=d0$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,config=cfg),control.inla=list(control.vb=list(emergency=50)))}
m<-fitm(cov,cfg=TRUE)
fx<-m$summary.fixed[cov,,drop=FALSE]
drv<-data.frame(cov=cov,label=lab(cov),IRR=round(exp(fx$mean),3),lo=round(exp(fx$`0.025quant`),3),hi=round(exp(fx$`0.975quant`),3),credible=(fx$`0.025quant`>0|fx$`0.975quant`<0),VIF=round(vif[cov],2))
drv<-drv[order(-drv$IRR),];write.csv(drv,file.path(R,paste0("MASTER_agestd_",ENG,"_drivers.csv")),row.names=FALSE)
uu<-m$summary.random$area[1:NR,];rr<-data.frame(region=regs,RR=round(exp(uu$mean),3),RR_lo=round(exp(uu$`0.025quant`),3),RR_hi=round(exp(uu$`0.975quant`),3))
write.csv(rr,file.path(R,paste0("MASTER_agestd_",ENG,"_regionSIR.csv")),row.names=FALSE)
pep<-sapply(1:NR,function(i) 1-inla.pmarginal(0,m$marginals.random$area[[i]]))
pepdf<-data.frame(region=regs,pep=round(pep,4));write.csv(pepdf,file.path(R,paste0("MASTER_agestd_",ENG,"_exceedance.csv")),row.names=FALSE)
phi<-m$summary.hyperpar["Phi for area","mean"];dic<-m$dic$dic;waic<-m$waic$waic
# children 민감도
CH<-read.csv(file.path(C,"ST_children_pct.csv"))
d0<-d0%>%left_join(CH[,c("region","year","children_0_4")],by=c("region","year"));d0$children<-as.numeric(scale(d0$children_0_4));d0$children[!is.finite(d0$children)]<-0
ms<-fitm(c(cov,"children"));fs<-ms$summary.fixed["children",]
sens<-data.frame(IRR=round(exp(fs$mean),3),lo=round(exp(fs$`0.025quant`),3),hi=round(exp(fs$`0.975quant`),3),credible=(fs$`0.025quant`>0|fs$`0.975quant`<0))
write.csv(sens,file.path(R,paste0("MASTER_agestd_",ENG,"_sens.csv")),row.names=FALSE)
# ---- EHSA + Moran (EB-SIR) ----
grid<-d0[,c("region","year","cases","E_std")];M<-matrix(0,NR,length(yy),dimnames=list(regs,yy))
for(ti in seq_along(yy)){py<-grid[grid$year==yy[ti],];idx<-match(regs,py$region);cs<-py$cases[idx];Ee<-py$E_std[idx];cs[is.na(cs)]<-0;Ee[is.na(Ee)|Ee<=0]<-min(Ee[Ee>0],na.rm=TRUE)
  eb<-tryCatch(EBest(cs,Ee)$estmm,error=function(e)cs/Ee);eb[!is.finite(eb)]<-0;M[,ti]<-eb}
nbs<-include.self(nb);lws<-nb2listw(nbs,style="W",zero.policy=TRUE);lw<-nb2listw(nb,style="W",zero.policy=TRUE)
mmc<-moran.mc(rowMeans(M),lw,nsim=999,zero.policy=TRUE)
classify<-function(zv,yrs){n<-length(zv);hot<-zv>1.96;cold<-zv< -1.96;fh<-mean(hot);fc<-mean(cold);mk<-suppressWarnings(tryCatch(cor.test(yrs,zv,method="kendall"),error=function(e)NULL));tau<-if(is.null(mk))NA else unname(mk$estimate);p<-if(is.null(mk))NA else mk$p.value;up<-(!is.na(p)&&p<0.05&&tau>0);dn<-(!is.na(p)&&p<0.05&&tau<0);lh<-hot[n];lc<-cold[n]
  if(any(hot)&&any(cold))return("Oscillating");if(lh){if(!any(hot[1:(n-1)]))return("New hot");if(fh>=0.9&&up)return("Intensifying hot");if(fh>=0.9)return("Persistent hot");if(dn)return("Diminishing hot");return("Sporadic/Consecutive hot")};if(!lh&&any(hot))return("Historical hot");if(lc){if(!any(cold[1:(n-1)]))return("New cold");if(fc>=0.9)return("Persistent cold");return("Sporadic cold")};if(!lc&&any(cold))return("Historical cold");"No pattern"}
Z2<-apply(M,2,function(x){x[!is.finite(x)]<-0;if(sd(x)==0)return(rep(0,length(x)));as.numeric(localG(x,lws,zero.policy=TRUE))});rownames(Z2)<-regs
cats<-sapply(1:NR,function(i)classify(Z2[i,],as.numeric(yy)))
write.csv(data.frame(region=regs,EHSA=cats),file.path(R,paste0("MASTER_agestd_",ENG,"_EHSAcat.csv")),row.names=FALSE)
ehsum<-data.frame(cat=names(table(cats)),n=as.integer(table(cats)));write.csv(ehsum,file.path(R,paste0("MASTER_agestd_",ENG,"_EHSAsum.csv")),row.names=FALSE)
# ---- 5년 구간 ----
PER<-list("2011-2015"=2011:2015,"2016-2020"=2016:2020,"2021-2024"=2021:2024);pr<-list()
for(pn in names(PER)){ys<-PER[[pn]];gg<-grid[grid$year%in%ys,]%>%group_by(region)%>%summarise(cases=sum(cases),E_std=sum(E_std),.groups="drop")
  Zw<-build_Z(cov,min(ys),max(ys))%>%group_by(region)%>%summarise(across(all_of(cov),mean),.groups="drop");d<-gg%>%left_join(Zw,by="region");d$area<-match(d$region,regs)
  fp<-as.formula(paste0("cases~",rhs,"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)"))
  mp<-tryCatch(inla(fp,family="nbinomial",data=d,E=d$E_std,control.fixed=CF,control.compute=list(dic=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
  if(!is.null(mp)){f2<-mp$summary.fixed[cov,,drop=FALSE];pr[[pn]]<-data.frame(period=pn,cov=cov,label=lab(cov),IRR=round(exp(f2$mean),3),lo=round(exp(f2$`0.025quant`),3),hi=round(exp(f2$`0.975quant`),3),credible=(f2$`0.025quant`>0|f2$`0.975quant`<0))}}
write.csv(bind_rows(pr),file.path(R,paste0("MASTER_agestd_",ENG,"_period.csv")),row.names=FALSE)
cat(sprintf("[%s] cov%d · credible %d · phi %.3f · Moran %.3f · EHSA hist %d/spor %d/new %d · children %.3f(%s)\n",ENG,length(cov),sum(drv$credible),phi,mmc$statistic,sum(cats=="Historical hot"),sum(cats=="Sporadic/Consecutive hot"),sum(cats=="New hot"),sens$IRR,ifelse(sens$credible,"cred","ns")))
# ================= 지도 3종 =================
sido_df<-suppressWarnings(shp%>%mutate(sido=substr(as.character(sigungu_cd),1,2))%>%group_by(sido)%>%summarise(.groups="drop")%>%st_set_crs(5179))
lay<-geom_sf(data=sido_df,fill=NA,color="grey55",linewidth=0.15,inherit.aes=FALSE)
enm<-endz(DZ)
# A) SIR
D<-rr;D$cls<-factor(bin_div(D$RR),levels=DIV_LABS);sp<-left_join(shp[,c("region","geometry")],D,by="region")
gA<-ggplot(sp)+geom_sf(aes(fill=cls),color="grey78",linewidth=0.03)+lvl_layer(sp,DIV_LABS)+lay+
  scale_fill_manual(values=setNames(PAL_DIV6,DIV_LABS),name="Age-standardised\nspatial RR",na.value="grey90",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+labs(title=paste0(enm," age-standardised spatial relative risk"),subtitle="Indirect standardisation (E_std offset). Colour = posterior RR; deep red = highest. Korea 2000/UTM-K.")+
  theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
ggsave(file.path(EH,paste0(ENG,"_Fig2_agestd_SIR.png")),gA,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)
# B) EHSA
lev<-c("New hot","Intensifying hot","Persistent hot","Sporadic/Consecutive hot","Diminishing hot","Historical hot","Oscillating","New cold","Persistent cold","Sporadic cold","Historical cold","No pattern")
pal<-c("New hot"="#B10026","Intensifying hot"="#E31A1C","Persistent hot"="#FC4E2A","Sporadic/Consecutive hot"="#FD8D3C","Diminishing hot"="#FEB24C","Historical hot"="#FED976","Oscillating"="#984EA3","New cold"="#08519C","Persistent cold"="#3182BD","Sporadic cold"="#6BAED6","Historical cold"="#BDD7E7","No pattern"="grey88")
spE<-left_join(shp[,c("region","geometry")],data.frame(region=regs,cls=factor(cats,levels=lev)),by="region")
gB<-ggplot(spE)+geom_sf(aes(fill=cls),color="grey78",linewidth=0.03)+lvl_layer(spE,lev)+lay+
  scale_fill_manual(values=pal[lev],breaks=lev,limits=lev,na.value="grey92",name="EHSA class",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+labs(title=paste0(enm," emerging hot spot analysis — age-standardised SIR (2011–2024)"),subtitle="Getis-Ord Gi* (EB-smoothed annual SIR) + Mann-Kendall. Korea 2000/UTM-K.")+
  theme_ehp_map(11)+guide_fill_box()+theme(text=element_text(family="Helvetica"),legend.position="right")
ggsave(file.path(EH,paste0(ENG,"_Fig3_agestd_EHSA.png")),gB,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)
# C) 초과확률
spX<-left_join(shp[,c("region","geometry")],pepdf,by="region");spX$cls<-bin_pep(spX$pep)
PLV<-c("<0.05 (cold)","0.05–0.5","0.5–0.8","0.8–0.95","≥0.95 (hot)")
gC<-ggplot(spX)+geom_sf(aes(fill=cls),color="grey78",linewidth=0.03)+lvl_layer(spX,PLV)+lay+
  scale_fill_manual(values=setNames(PAL_PEP5,PLV),breaks=PLV,limits=PLV,na.value="grey90",name="P(spatial RR > 1)",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+labs(title=paste0(enm," posterior exceedance probability P(spatial RR > 1)"),subtitle="Age-standardised BYM2 structured spatial effect. >=0.95 credible high-risk; <0.05 credible low-risk. Korea 2000/UTM-K.")+
  theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
ggsave(file.path(EH,paste0(ENG,"_Fig_agestd_exceedance.png")),gC,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)
# fit 요약
writeLines(sprintf("phi=%.3f DIC=%.1f WAIC=%.1f Moran=%.3f(p=%.3f) hi=%d lo=%d",phi,dic,waic,mmc$statistic,mmc$p.value,sum(rr$RR_lo>1),sum(rr$RR_hi<1)),file.path(R,paste0("agestd_",ENG,"_fit.txt")))
cat("완료:",ENG,"\n")
