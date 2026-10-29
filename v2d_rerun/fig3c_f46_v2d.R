# fig3c_f46_v2d.R — (2026-10-29 /333 7번) 메인 Fig 3C 지속성 지도 재작도: 범례·방위표·축척을 A·B(412b)와 같은 크기로, ">= 1" → "≥ 1".
#  입력 EHEC/series/EHEC_T34_persistence.csv(595 산출). 출력 JKMS_EHEC_v2d_261029/fig_src/EHEC_F46_map_persistence_v2d.png (공유 F46은 그대로 둠).
suppressPackageStartupMessages({library(dplyr);library(sf);library(stringr);library(ggplot2);library(ggspatial)})
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-file.path(GD,"FBD_SpatioTemporal");OUT<-file.path(ROOT,"JKMS_FBD","JKMS_EHEC_v2d_261029","fig_src");dir.create(OUT,showWarnings=FALSE)
source(file.path(GD,"EHP_Spatio/fig_style.R"));AGG<-ragg::agg_png
north_big<-function(loc="tr") ggspatial::annotation_north_arrow(location=loc,which_north="true",height=unit(2.55,"cm"),width=unit(1.8,"cm"),
  style=ggspatial::north_arrow_orienteering(text_size=21,fill=c("white","grey20"),line_col="grey20"),pad_x=unit(0.25,"cm"),pad_y=unit(0.25,"cm"))
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
sido_df<-suppressWarnings(shp%>%mutate(sido=substr(as.character(sigungu_cd),1,2))%>%group_by(sido)%>%summarise(.groups="drop")%>%st_set_crs(5179))
Pst<-read.csv(file.path(ROOT,"EHEC","series","EHEC_T34_persistence.csv"),fileEncoding="UTF-8");Pst$region<-str_replace_all(Pst$region,"\\s+","")
stopifnot(nrow(Pst)==229,all(Pst$region%in%shp$region))
Pst$cls<-cut(Pst$years_SIR_gt1,breaks=c(-1,0,2,4,6,9,14),labels=c("0","1-2","3-4","5-6","7-9","10-14"))
spP<-left_join(shp[,c("region","geometry")],Pst,by="region")
g46<-ggplot(spP)+geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+lvl_layer(spP,levels(Pst$cls))+geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.32,inherit.aes=FALSE)+
  scale_fill_manual(values=setNames(c("white",PAL_SEQ5),levels(Pst$cls)),breaks=levels(Pst$cls),limits=levels(Pst$cls),na.value="grey90",name="Years with ≥ 1 case\nand SIR > 1 (of 14)",drop=FALSE)+
  north_big("tr")+scale_ehp("bl")+theme_ehp_map(20)+guide_fill_box(rev=TRUE)+
  theme(text=element_text(family="Helvetica"),legend.position="right",legend.title=element_text(size=19,face="bold"),legend.text=element_text(size=17),legend.key.size=unit(0.85,"cm"))
d7<-spP[spP$fdr_years_SIR_gt1<0.05,];cat("FDR<0.05:",nrow(d7),"\n")
if(nrow(d7))g46<-g46+geom_sf(data=d7,aes(color="More years than expected\n(FDR < 0.05)"),fill=NA,linewidth=0.7,inherit.aes=FALSE)+scale_color_manual(name=NULL,values=c("More years than expected\n(FDR < 0.05)"="black"))+guides(color=guide_legend(override.aes=list(fill=NA,linewidth=0.7),order=2))
ggsave(file.path(OUT,"EHEC_F46_map_persistence_v2d.png"),g46,width=12.2,height=11.2,dpi=300,bg="white",device=AGG);cat("F46 v2d redrawn\n")
