# ④ M1–M6 모형선택 사다리 (저자 서식): 비공간→ICAR→BYM→+상호작용→+RW1→full, DIC/WAIC
source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_st_setup.R"))
DZ<-getsel();rows<-list()
for(dz in DZ){pd<-panel[panel$disease==dz,];fy<-unique(pd$fy);yy<-(2011:2024)[(2011:2024)>=fy]
  d0<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(pd[,c("region","year","cases","pop")],by=c("region","year"))
  d0$cases[is.na(d0$cases)]<-0;d0$pop[is.na(d0$pop)|d0$pop<=0]<-median(d0$pop,na.rm=TRUE)
  d0$E<-pmax(d0$pop*sum(d0$cases)/sum(d0$pop),1);d0$area<-match(d0$region,regs);d0$time<-match(d0$year,yy);d0$id<-1:nrow(d0)
  forms<-list(
   M1="cases~1",
   M2="cases~1+f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)",
   M3="cases~1+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)",
   M4="cases~1+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(id,model='iid',hyper=pcp)",
   M5="cases~1+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)",
   M6="cases~1+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)")
  dic<-c();wa<-c()
  for(nm in names(forms)){f<-tryCatch(inla(as.formula(forms[[nm]]),family="nbinomial",data=d0,E=d0$E,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
    dic[nm]<-if(is.null(f))NA else round(f$dic$dic,1);wa[nm]<-if(is.null(f))NA else round(f$waic$waic,1)}
  best<-names(which.min(dic))
  rows[[dz]]<-data.frame(disease=dz,eng=EN[dz],DIC_M1=dic["M1"],DIC_M2=dic["M2"],DIC_M3=dic["M3"],DIC_M4=dic["M4"],DIC_M5=dic["M5"],DIC_M6=dic["M6"],WAIC_M6=wa["M6"],best_DIC=best,dDIC_M1_M6=round(dic["M1"]-dic["M6"],1))
  cat(sprintf("[%s] M1..M6 DIC: %s best=%s\n",dz,paste(dic,collapse="/"),best))}
M<-bind_rows(rows);tag<-paste(match(DZ,ALLDZ),collapse="_")
write.csv(M,file.path(ROOT,"results_ST",paste0("M1M6_",tag,".csv")),row.names=FALSE);cat("저장 M1M6_",tag,".csv\n",sep="")
