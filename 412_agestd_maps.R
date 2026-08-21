# 412 — EHEC 연령표준화 지도 (fig_style.R 저널표준 적용, 304/333 미러)
#  (A) 공간 SIR RR 지도: 이산 PAL_DIV6+bin_div, 박스범례(고위험 상단), 고위험 dissolve 외곽선
#  (B) EHSA 유형 지도: 전 레벨 lvl_layer 강제 채움, 박스범례
#  공통: layer_sido·north_ehp·scale_ehp·theme_ehp_map·EPSG5179·ragg(en-dash)
suppressPackageStartupMessages({library(dplyr);library(sf);library(stringr);library(ggplot2)})
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-file.path(GD,"FBD_SpatioTemporal");EH<-file.path(ROOT,"EHEC")
source(file.path(GD,"EHP_Spatio/fig_style.R"))
AGG<-if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else NULL
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
# 위계: 시군구=얇은 회색 실선(보이게), 시도=그보다 진한 중간회색(검은 볼드는 아님)
sido_df<-suppressWarnings(shp%>%mutate(sido=substr(as.character(sigungu_cd),1,2))%>%group_by(sido)%>%summarise(.groups="drop")%>%st_set_crs(5179))
lay_sido<-geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.32,inherit.aes=FALSE)

# ================= (A) 공간 SIR RR 지도 =================
D<-read.csv(file.path(ROOT,"results_ST/MASTER_agestd_regionSIR.csv"))
D$cls<-factor(bin_div(D$RR),levels=DIV_LABS); D$hi_risk<-D$RR_lo>1
sp<-left_join(shp[,c("region","geometry")],D,by="region")
gA<-ggplot(sp)+
  geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+
  lvl_layer(sp,DIV_LABS)+
  lay_sido+
  scale_fill_manual(values=setNames(PAL_DIV6,DIV_LABS),name="Age-standardised\nspatial RR",na.value="grey90",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+
  labs(title="EHEC age-standardised spatial relative risk",
       subtitle="Indirect standardisation (E_std offset). Colour = posterior relative risk. Korea 2000 / UTM-K (EPSG:5179).")+
  theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
for(p in c(file.path(ROOT,"results_ST/Fig2_agestd_spatialSIR.png"),file.path(EH,"EHEC_Fig2_agestd_SIR.png")))
  ggsave(p,gA,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)

# ================= (B) EHSA 유형 지도 =================
E<-read.csv(file.path(ROOT,"results_ST/EHSA_agestd_cat.csv"),fileEncoding="UTF-8")
lev<-c("New hot","Intensifying hot","Persistent hot","Sporadic/Consecutive hot","Diminishing hot","Historical hot","Oscillating","New cold","Persistent cold","Sporadic cold","Historical cold","No pattern")
pal<-c("New hot"="#B10026","Intensifying hot"="#E31A1C","Persistent hot"="#FC4E2A","Sporadic/Consecutive hot"="#FD8D3C","Diminishing hot"="#FEB24C","Historical hot"="#FED976","Oscillating"="#984EA3","New cold"="#08519C","Persistent cold"="#3182BD","Sporadic cold"="#6BAED6","Historical cold"="#BDD7E7","No pattern"="grey88")
E$cls<-factor(E$EHSA,levels=lev)
spE<-left_join(shp[,c("region","geometry")],E,by="region")
gB<-ggplot(spE)+
  geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+
  lvl_layer(spE,lev)+
  lay_sido+
  scale_fill_manual(values=pal[lev],breaks=lev,limits=lev,na.value="grey92",name="EHSA class",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+
  labs(title="EHEC emerging hot spot analysis — age-standardised SIR (2011–2024)",
       subtitle="Getis-Ord Gi* (EB-smoothed annual SIR) + Mann-Kendall. Historical-hot dominant = clustering diffusing. Korea 2000 / UTM-K.")+
  theme_ehp_map(11)+guide_fill_box()+theme(text=element_text(family="Helvetica"),legend.position="right")
for(p in c(file.path(ROOT,"results_ST/Fig3_agestd_EHSA.png"),file.path(EH,"EHEC_Fig3_agestd_EHSA.png")))
  ggsave(p,gB,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)

# ================= (C) 초과확률 지도 (CSV에서, 재적합 불필요) =================
xf<-file.path(ROOT,"results_ST/MASTER_agestd_exceedance.csv")
if(file.exists(xf)){
  X<-read.csv(xf);spX<-left_join(shp[,c("region","geometry")],X,by="region");spX$cls<-bin_pep(spX$pep)
  PLV<-c("<0.05 (cold)","0.05–0.5","0.5–0.8","0.8–0.95","≥0.95 (hot)")
  gC<-ggplot(spX)+geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+lvl_layer(spX,PLV)+lay_sido+
    scale_fill_manual(values=setNames(PAL_PEP5,PLV),breaks=PLV,limits=PLV,na.value="grey90",name="P(spatial RR > 1)",drop=FALSE)+
    north_ehp("tr")+scale_ehp("bl")+
    labs(title="EHEC posterior exceedance probability P(spatial RR > 1)",
         subtitle="Age-standardised BYM2 structured spatial effect. >=0.95 = credible high-risk; <0.05 = credible low-risk. Korea 2000/UTM-K.")+
    theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
  for(p in c(file.path(ROOT,"results_ST/Fig_agestd_exceedance.png"),file.path(EH,"EHEC_Fig_agestd_exceedance.png")))
    ggsave(p,gC,width=8.4,height=9.2,dpi=300,bg="white",device=AGG)
  cat("저장: Fig_agestd_exceedance\n")
}
cat(sprintf("저장: Fig2_agestd_spatialSIR·Fig3_agestd_EHSA (results_ST + EHEC). 고위험 %d개\n",sum(D$hi_risk)))
