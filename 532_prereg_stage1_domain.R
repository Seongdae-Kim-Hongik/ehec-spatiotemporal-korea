# 532 — 사전지정 1단계: 도메인별 BYM2 시공간 모형. 각 도메인 변수 + 도농축(kchs_pc1) + 성비(forced).
#  목적: 도메인 안에서 어떤 변수가 신호를 갖는지 — 2단계 결합 모형의 대표 선정 근거(사전규칙: credible 또는 |log IRR| 상위 3).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(jsonlite);library(car)})
set.seed(20260822); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA",unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT",unset="."); CACHE<-file.path(ROOT,"data_cache"); RS<-file.path(ROOT,"results_ST")
DES<-fromJSON(file.path(CACHE,"PREREG_DESIGN.json"))
d0<-read.csv(file.path(CACHE,"ST_prereg_panel.csv"),fileEncoding="UTF-8-BOM"); d0$region<-str_replace_all(d0$region,"\\s+","")
allv<-unique(c(unlist(lapply(DES$domains,function(x)x$vars)),"kchs_pc1","sex_ratio"))
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in allv){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
# 그래프
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"));yy<-2011:2024
fitmod<-function(cont){
  dd<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
  dd$cases[is.na(dd$cases)]<-0;dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
  dd$area<-match(dd$region,regs);dd$time<-match(dd$year,yy);dd$id<-1:nrow(dd)
  form<-as.formula(paste0("cases~",paste(sprintf("`%s`",cont),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  tryCatch(inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e){cat("  [실패]",conditionMessage(e),"\n");NULL})}
grab<-function(m,vs,tag,dom){fx<-m$summary.fixed;do.call(rbind,lapply(vs,function(v)if(v%in%rownames(fx)){r<-fx[v,];data.frame(model=tag,domain=dom,term=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),DIC=round(m$dic$dic,1))}else NULL))}
# VIF 사전정리 함수(도메인 내 >10 반복제거)
vifprune<-function(vs){keep<-vs;repeat{if(length(keep)<2)break;vd<-Z[,keep];vd$.y<-rnorm(nrow(vd));v<-tryCatch(car::vif(lm(.y~.,data=vd)),error=function(e)NULL);if(is.null(v)||max(v)<10)break;drop<-names(which.max(v));cat("   VIF>10 제거:",drop,sprintf("(%.1f)\n",max(v)));keep<-setdiff(keep,drop)};keep}

R<-list(); cat("═══ 1단계: 도메인별 모형 (+ kchs_pc1 + sex_ratio) ═══\n")
for(dom in names(DES$domains)){
  vs<-setdiff(DES$domains[[dom]]$vars,c("sex_ratio","kchs_pc1")); vs<-vifprune(vs)
  cont<-unique(c(vs,"kchs_pc1","sex_ratio")); cat(sprintf("\n[%s] %d변수\n",dom,length(vs)))
  m<-fitmod(cont); if(is.null(m))next
  r<-grab(m,cont,paste0("S1:",dom),dom); R[[dom]]<-r
  for(i in 1:nrow(r))cat(sprintf("  %-40s IRR %.3f (%.3f-%.3f)%s\n",r$term[i],r$IRR[i],r$lo[i],r$hi[i],ifelse(r$cred[i]," *","")))
  cat(sprintf("  DIC %.1f\n",r$DIC[1]))}
out<-do.call(rbind,R); write.csv(out,file.path(RS,"PREREG_S1_domain.csv"),row.names=FALSE)
# 2단계 대표 선정(사전규칙): 도메인별 credible 전부 + 비credible 중 |log IRR| 상위 2 (최대 5/도메인); kchs_pc1·sex_ratio 제외
rep<-out%>%filter(!term%in%c("kchs_pc1","sex_ratio"))%>%group_by(domain)%>%arrange(desc(cred),desc(abs(log(IRR))))%>%mutate(rk=row_number())%>%filter(cred|rk<=2)%>%slice_head(n=5)%>%ungroup()
write.csv(rep,file.path(RS,"PREREG_S1_representatives.csv"),row.names=FALSE)
cat("\n═══ 2단계 대표 변수 (사전규칙) ═══\n");for(dm in unique(rep$domain))cat(sprintf("  %-16s %s\n",dm,paste(rep$term[rep$domain==dm],collapse=", ")))
cat("\n저장: PREREG_S1_domain.csv · PREREG_S1_representatives.csv\n")
