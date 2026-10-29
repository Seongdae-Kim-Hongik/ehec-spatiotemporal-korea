# 595b_ehec_f46_redraw.R — EHEC 지속성 지도(F46)만 큰 범례로 재작도(메인 Fig 3C, A/B 지도와 글자 크기 맞춤). 입력: EHEC/series/EHEC_T34_persistence.csv
# 542 — 최종 주모형(20변수) 기준 기간 분석: ① 4년창 순수공간 ② 코로나 전후 분할 + 상호작용 ③ 주요 3변수 × 시기 상호작용
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-file.path(ROOT,"results_ST")
V2_CACHE<-CACHE;source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지 (기술 분석이라 청주 마스킹 없음 — 관측 0 사용). 이전 판 .bak_pre_v2
T2<-read.csv(file.path(RS,"MASTER_final_table2.csv"));cont<-T2$cov;ENG<-setNames(T2$label,T2$cov)
d0<-read.csv(file.path(CACHE,"ST_prereg_panel.csv"),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")];for(v in cont){x<-d0[[v]];x<-v2_nofill(x,v);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")]
mk<-function(yy){DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
  DD<-v2_finish(DD,kdca=FALSE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD);DD}
grab<-function(m,v)data.frame(cov=v,label=ENG[v],IRR=round(exp(m$summary.fixed[v,"mean"]),3),lo=round(exp(m$summary.fixed[v,"0.025quant"]),3),hi=round(exp(m$summary.fixed[v,"0.975quant"]),3),cred=(m$summary.fixed[v,"0.025quant"]>0|m$summary.fixed[v,"0.975quant"]<0))
KEY<-c("beef_cattle_farms","fiscal_independence","net_migration","current_smoking")
suppressPackageStartupMessages({library(ggplot2);library(ggspatial)});source(file.path(file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads"),"EHP_Spatio/fig_style.R"));AGG<-ragg::agg_png
SH<-file.path(ROOT,"EHEC","series");Pst<-read.csv(file.path(SH,"EHEC_T34_persistence.csv"),fileEncoding="UTF-8")
sido_df<-suppressWarnings(shp%>%mutate(sido=substr(as.character(sigungu_cd),1,2))%>%group_by(sido)%>%summarise(.groups="drop")%>%st_set_crs(5179))
Pst$cls<-cut(Pst$years_SIR_gt1,breaks=c(-1,0,2,4,6,9,14),labels=c("0","1-2","3-4","5-6","7-9","10-14"))
spP<-left_join(shp[,c("region","geometry")],Pst,by="region")
g46<-ggplot(spP)+geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+lvl_layer(spP,levels(Pst$cls))+geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.32,inherit.aes=FALSE)+
  scale_fill_manual(values=setNames(c("white",PAL_SEQ5),levels(Pst$cls)),breaks=levels(Pst$cls),limits=levels(Pst$cls),na.value="grey90",name="Years with >= 1 case\nand SIR > 1 (of 14)",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+theme_ehp_map(16)+guide_fill_box(rev=TRUE)+
  theme(text=element_text(family="Helvetica"),legend.position="right",legend.title=element_text(size=17,face="bold"),legend.text=element_text(size=16),legend.key.size=unit(0.75,"cm"))
d7<-spP[spP$fdr_years_SIR_gt1<0.05,];if(nrow(d7))g46<-g46+geom_sf(data=d7,aes(color="More years than expected\n(FDR < 0.05)"),fill=NA,linewidth=0.7,inherit.aes=FALSE)+scale_color_manual(name=NULL,values=c("More years than expected\n(FDR < 0.05)"="black"))+guides(color=guide_legend(override.aes=list(fill=NA,linewidth=0.7),order=2))
ggsave(file.path(SH,"EHEC_F46_map_persistence.png"),g46,width=11,height=11.2,dpi=300,bg="white",device=AGG);cat("F46 redrawn\n")
