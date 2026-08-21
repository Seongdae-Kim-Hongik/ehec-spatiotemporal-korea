# 499 — 코로나 전후 공간 RR 2패널 지도 (498 [A] 산출 regionRR 사용). 메일 보고용.
suppressPackageStartupMessages({library(dplyr);library(sf);library(stringr);library(ggplot2);library(patchwork)})
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".")
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
g<-read.csv(file.path(ROOT,"results_ST","MASTER_ehec_covid_regionRR.csv"))
m<-shp%>%left_join(g,by="region")
lim<-max(abs(log(c(m$RR_pre,m$RR_post))),na.rm=TRUE)
mk<-function(col,hicol,ttl,n,ndist){
  ggplot(m)+geom_sf(aes(fill=log(.data[[col]])),color="grey55",linewidth=0.08)+
    geom_sf(data=m[m[[hicol]]==1,],fill=NA,color="black",linewidth=0.35)+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-lim,lim),
      breaks=log(c(0.5,1,2)),labels=c("0.5","1","2"),name="Spatial RR")+
    labs(title=ttl,subtitle=sprintf("N = %s cases; %d districts with credibly elevated RR (black outline)",format(n,big.mark=","),ndist))+
    theme_void(base_family="Helvetica")+
    theme(plot.title=element_text(face="bold",size=13),plot.subtitle=element_text(size=9,color="grey30"),
          legend.position=c(0.88,0.28),legend.title=element_text(size=9),legend.text=element_text(size=8))
}
pA<-mk("RR_pre","hi_pre","A  Pre-COVID (2011-2019)",880,sum(g$hi_pre))
pB<-mk("RR_post","hi_post","B  COVID/post (2020-2024)",1136,sum(g$hi_post))
gg<-pA+pB
ggsave(file.path(ROOT,"EHEC","EHEC_Fig_covid_split_RR.png"),gg,width=13,height=6.2,dpi=300,bg="white")
cat("저장: EHEC/EHEC_Fig_covid_split_RR.png\n")
