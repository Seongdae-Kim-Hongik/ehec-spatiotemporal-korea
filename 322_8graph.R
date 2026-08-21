# ③ 8-graph 인접 민감도 (저자 서식): credible 드라이버가 8개 이웃정의서 유지되나 → regional vs local
source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_st_setup.R"))
DZ<-getsel()
GRAPHS<-list()
GRAPHS[["Queen"]]<-build_nb("queen");GRAPHS[["Rook"]]<-build_nb("rook")
for(k in 2:7)GRAPHS[[paste0("knn",k)]]<-build_nb("knn",k=k)
gfiles<-lapply(GRAPHS,function(nbx){f<-tempfile(fileext=".graph");nb2INLA(f,nbx);f})
rows<-list()
for(dz in DZ){cov<-RICH[[dz]];pd<-panel[panel$disease==dz,];fy<-unique(pd$fy);yy<-(2011:2024)[(2011:2024)>=fy]
  Z<-build_Z(cov,2011,2024)
  d0<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(pd[,c("region","year","cases","pop")],by=c("region","year"))%>%left_join(Z,by=c("region","year"))
  d0$cases[is.na(d0$cases)]<-0;d0$pop[is.na(d0$pop)|d0$pop<=0]<-median(d0$pop[d0$pop>0],na.rm=TRUE)
  d0$E<-pmax(d0$pop*sum(d0$cases)/sum(d0$pop),1);d0$area<-match(d0$region,regs);d0$time<-match(d0$year,yy)
  rhs<-paste(sprintf("`%s`",cov),collapse="+")
  credmat<-matrix(FALSE,length(cov),length(GRAPHS),dimnames=list(cov,names(GRAPHS)))
  for(gi in seq_along(GRAPHS)){gg<-gfiles[[gi]]
    form<-as.formula(paste0("cases~",rhs,"+f(area,model='bym2',graph=gg,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)"))
    m<-tryCatch(inla(form,family="nbinomial",data=d0,E=d0$E,control.fixed=CF,control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
    if(is.null(m))next;fx<-m$summary.fixed[cov,,drop=FALSE];credmat[,gi]<-(fx$`0.025quant`>0|fx$`0.975quant`<0)}
  ng<-rowSums(credmat)
  for(j in seq_along(cov))if(ng[j]>0)rows[[paste(dz,cov[j])]]<-data.frame(disease=dz,eng=EN[dz],공변량=ifelse(cov[j]%in%names(KLAB),KLAB[cov[j]],cov[j]),cov=cov[j],
    n_graphs_credible=ng[j],classification=ifelse(ng[j]==8,"Regional (8/8)",ifelse(ng[j]>=6,paste0("Regional (",ng[j],"/8)"),ifelse(credmat[j,"Queen"]|credmat[j,"Rook"],"Local (contiguity-only)","Graph-sensitive"))))
  cat(sprintf("[%s] 8-graph 완료 (regional8/8: %d)\n",dz,sum(ng==8)))}
M<-bind_rows(rows);tag<-paste(match(DZ,ALLDZ),collapse="_")
write.csv(M,file.path(ROOT,"results_ST",paste0("graph8_",tag,".csv")),row.names=FALSE);cat("저장 graph8_",tag,".csv\n",sep="")
