# 412c — Fig. 3B 재작도: 표준(Esri) EHSA 분류(411b), 존재하는 범주만 범례 (2026-09-14). 헤더 = 412b 1-30행.
# 412b — JAMA용 지도 변형 (412 미러): baked-in 제목·부제 제거 + 미국식 철자 범례.
#  플레이북 §6-C(제목은 캡션이 담당) + JAMA(그림 내 편집문구 금지). 출력: EHEC/JAMA_mapA_RR.png, JAMA_mapB_EHSA.png
suppressPackageStartupMessages({library(dplyr);library(sf);library(stringr);library(ggplot2)})
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-file.path(GD,"FBD_SpatioTemporal");EH<-file.path(ROOT,"EHEC")
source(file.path(GD,"EHP_Spatio/fig_style.R"))
# 412b 전용 오버라이드: N 화살표 3배 (fig_style 공용은 불변)
north_big<-function(loc="tr") ggspatial::annotation_north_arrow(location=loc,which_north="true",
  height=unit(2.55,"cm"),width=unit(1.8,"cm"),
  style=ggspatial::north_arrow_orienteering(text_size=21,fill=c("white","grey20"),line_col="grey20"),
  pad_x=unit(0.25,"cm"),pad_y=unit(0.25,"cm"))
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
  north_big("tr")+scale_ehp("bl")+
  theme_ehp_map(20)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right",legend.text=element_text(size=17),legend.title=element_text(size=19,face="bold"),legend.key.size=unit(0.85,"cm"))
ggsave(file.path(EH,"JAMA_mapA_RR.png"),gA,width=11,height=11.2,dpi=300,bg="white",device=AGG)

E<-read.csv(file.path(ROOT,"results_ST/EHSA_agestd_cat_esri.csv"),fileEncoding="UTF-8")
lev_all<-c("New hot","Consecutive hot","Intensifying hot","Persistent hot","Diminishing hot","Sporadic hot","Oscillating hot","Historical hot",
           "New cold","Consecutive cold","Intensifying cold","Persistent cold","Diminishing cold","Sporadic cold","Oscillating cold","Historical cold","No pattern")
pal<-c("New hot"="#B10026","Consecutive hot"="#E31A1C","Intensifying hot"="#FC4E2A","Persistent hot"="#FD6A3C","Diminishing hot"="#FEB24C","Sporadic hot"="#FD8D3C","Oscillating hot"="#984EA3","Historical hot"="#FED976",
       "New cold"="#08519C","Consecutive cold"="#2171B5","Intensifying cold"="#4292C6","Persistent cold"="#3182BD","Diminishing cold"="#9ECAE1","Sporadic cold"="#6BAED6","Oscillating cold"="#807DBA","Historical cold"="#BDD7E7","No pattern"="grey88")
lev<-lev_all[lev_all %in% unique(E$EHSA_esri)]
E$cls<-factor(E$EHSA_esri,levels=lev)
spE<-left_join(shp[,c("region","geometry")],E,by="region")
gB<-ggplot(spE)+geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+lvl_layer(spE,lev)+lay_sido+
  scale_fill_manual(values=pal[lev],breaks=lev,limits=lev,na.value="grey92",name="Emerging hot spot\nclass",drop=FALSE)+
  north_big("tr")+scale_ehp("bl")+
  theme_ehp_map(20)+guide_fill_box()+theme(text=element_text(family="Helvetica"),legend.position="right",legend.text=element_text(size=17),legend.title=element_text(size=19,face="bold"),legend.key.size=unit(0.85,"cm"))
ggsave(file.path(EH,"JAMA_mapB_EHSA.png"),gB,width=11,height=11.2,dpi=300,bg="white",device=AGG)
cat("saved JAMA_mapB_EHSA (Esri):",paste(sprintf("%s %d",names(table(E$cls)),table(E$cls)),collapse=", "),"\n")
