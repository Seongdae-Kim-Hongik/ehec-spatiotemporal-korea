# 412b — JAMA용 지도 변형 (412 미러): baked-in 제목·부제 제거 + 미국식 철자 범례.
#  플레이북 §6-C(제목은 캡션이 담당) + JAMA(그림 내 편집문구 금지). 출력: EHEC/JAMA_mapA_RR.png, JAMA_mapB_EHSA.png
suppressPackageStartupMessages({library(dplyr);library(sf);library(stringr);library(ggplot2)})
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-file.path(GD,"FBD_SpatioTemporal");EH<-file.path(ROOT,"EHEC")
source(file.path(GD,"EHP_Spatio/fig_style.R"))
AGG<-if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else NULL
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
sido_df<-suppressWarnings(shp%>%mutate(sido=substr(as.character(sigungu_cd),1,2))%>%group_by(sido)%>%summarise(.groups="drop")%>%st_set_crs(5179))
lay_sido<-geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.32,inherit.aes=FALSE)

# (A) 연령표준화 공간 RR
D<-read.csv(file.path(ROOT,"results_ST/MASTER_agestd_regionSIR.csv"))
D$cls<-factor(bin_div(D$RR),levels=DIV_LABS); D$hi_risk<-D$RR_lo>1
sp<-left_join(shp[,c("region","geometry")],D,by="region")
gA<-ggplot(sp)+
  geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+
  lvl_layer(sp,DIV_LABS)+
  lay_sido+
  scale_fill_manual(values=setNames(PAL_DIV6,DIV_LABS),name="Age-standardized\nspatial RR",na.value="grey90",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+
  theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
ggsave(file.path(EH,"JAMA_mapA_RR.png"),gA,width=8.4,height=8.6,dpi=300,bg="white",device=AGG)

# (B) EHSA 유형
E<-read.csv(file.path(ROOT,"results_ST/EHSA_agestd_cat.csv"),fileEncoding="UTF-8")
lev<-c("New hot","Intensifying hot","Persistent hot","Sporadic/Consecutive hot","Diminishing hot","Historical hot","Oscillating","New cold","Persistent cold","Sporadic cold","Historical cold","No pattern")
pal<-c("New hot"="#B10026","Intensifying hot"="#E31A1C","Persistent hot"="#FC4E2A","Sporadic/Consecutive hot"="#FD8D3C","Diminishing hot"="#FEB24C","Historical hot"="#FED976","Oscillating"="#984EA3","New cold"="#08519C","Persistent cold"="#3182BD","Sporadic cold"="#6BAED6","Historical cold"="#BDD7E7","No pattern"="grey88")
E$cls<-factor(E$EHSA,levels=lev)
spE<-left_join(shp[,c("region","geometry")],E,by="region")
gB<-ggplot(spE)+
  geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+
  lvl_layer(spE,lev)+
  lay_sido+
  scale_fill_manual(values=pal[lev],breaks=lev,limits=lev,na.value="grey92",name="Emerging hot spot\nclass",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+
  theme_ehp_map(11)+guide_fill_box()+theme(text=element_text(family="Helvetica"),legend.position="right")
ggsave(file.path(EH,"JAMA_mapB_EHSA.png"),gB,width=8.4,height=8.6,dpi=300,bg="white",device=AGG)
cat(sprintf("저장: JAMA_mapA_RR + JAMA_mapB_EHSA. 고위험 %d개\n",sum(D$hi_risk)))
