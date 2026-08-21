# 전국 월별 시계열 STL 분해 (보충): 추세+계절+잔차. 밀도 충분한 질병만. robust=TRUE (대유행→잔차).
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(dplyr)})
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-file.path(GD,"FBD_SpatioTemporal")
EN<-c("노로바이러스"="Norovirus","A형간염"="HepatitisA","E형간염"="HepatitisE","황색포도상구균"="S_aureus","바실러스세레우스"="B_cereus","클로스트리디움퍼프린젠스"="C_perfringens","비브리오패혈증"="V_vulnificus","장염비브리오"="V_parahaemolyticus","장티푸스"="Typhoid","세균성이질"="Shigellosis","살모넬라"="Salmonella","병원성대장균"="Pathogenic_Ecoli","장출혈성대장균감염증"="EHEC","캠필로박터제주니"="Campylobacter","원충"="Kudoa")
raw<-fread(file.path(BASE,"식중독최종.csv"),select=c("disease","year","month","cases"))[!is.na(month)&month>=1&month<=12]
summ<-data.frame()
for(dz in names(EN)){eng<-EN[dz]
  d<-raw[disease==dz,.(cs=sum(cases,na.rm=TRUE)),by=.(year,month)][order(year,month)]
  if(nrow(d)<24){cat("skip(짧음):",dz,"\n");next}
  y0<-min(d$year);y1<-max(d$year)
  full<-CJ(year=y0:y1,month=1:12);d<-merge(full,d,by=c("year","month"),all.x=TRUE);d$cs[is.na(d$cs)]<-0
  d<-d[order(year,month)]
  npts<-nrow(d);zero<-mean(d$cs==0);mn<-mean(d$cs)
  # 밀도 기준: 5년(60월)+ , 0월<55%, 월평균>=8
  ok<-(npts>=60 & zero<0.55 & mn>=8)
  summ<-rbind(summ,data.frame(disease=dz,eng=eng,y0=y0,y1=y1,npts=npts,zero_pct=round(zero*100),mean_mo=round(mn,1),STL=ok))
  if(!ok){cat("skip(희소):",dz,sprintf("(0월%.0f%%, 월평균%.1f)\n",zero*100,mn));next}
  ts0<-ts(d$cs,start=c(y0,1),frequency=12)
  fit<-stl(ts0,s.window="periodic",robust=TRUE)
  cp<-as.data.frame(fit$time.series);cp$obs<-as.numeric(ts0)
  cp$t<-seq(y0,by=1/12,length.out=nrow(cp))
  long<-rbind(
    data.frame(t=cp$t,val=cp$obs, comp="관측(Observed)"),
    data.frame(t=cp$t,val=cp$trend,comp="추세(Trend)"),
    data.frame(t=cp$t,val=cp$seasonal,comp="계절(Seasonal)"),
    data.frame(t=cp$t,val=cp$remainder,comp="잔차(Remainder)"))
  long$comp<-factor(long$comp,levels=c("관측(Observed)","추세(Trend)","계절(Seasonal)","잔차(Remainder)"))
  gg<-ggplot(long,aes(t,val))+geom_line(color="#2166AC",linewidth=0.4)+
    geom_vline(xintercept=2020,linetype="dashed",color="#B2182B",linewidth=0.3)+
    geom_hline(data=subset(long,comp%in%c("계절(Seasonal)","잔차(Remainder)")),aes(yintercept=0),color="grey70",linewidth=0.2)+
    facet_wrap(~comp,ncol=1,scales="free_y")+
    labs(title=paste0(dz," — 전국 월별 시계열 STL 분해 (",y0,"–",y1,")"),
         subtitle="robust STL(계절 주기 고정). 빨강파선=COVID(2020). 대유행은 잔차로 분리. 월별 신고건수.",
         x="연도",y="월별 건수")+
    theme_bw(base_size=11)+theme(text=element_text(family="Helvetica"),
      strip.text=element_text(face="bold"),plot.title=element_text(face="bold"),
      plot.subtitle=element_text(size=8,color="grey30"),panel.grid.minor=element_blank())
  ggsave(file.path(ROOT,eng,paste0(eng,"_FigS_STL.png")),gg,width=8,height=8.5,dpi=200,bg="white")
  cat("STL 저장:",dz,"\n")}
write.csv(summ,file.path(ROOT,"results_ST/STL_feasibility.csv"),row.names=FALSE)
cat("\n=== STL 적합성 요약 ===\n");print(summ[,c("disease","y0","y1","zero_pct","mean_mo","STL")])
