# 514 — c_hw_soap 부호 반전 추적. 원시 r(SIR)=-0.083 → 조정 IRR 1.529 의 전환점을 찾는다.
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260821); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
cont_cov<-c("d_pasture","d_field","d_paddy","e_forest","a_septic","a_sewer","b_gw_use","b_gw_unfit","h_fastfood","h_pub","h_selfhealth","g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder")
d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv")); d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
d0$ruminant<-d0$d_cattle+d0$d_dairy; qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in unique(c(cont_cov,"c_hw_soap"))){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE)
  xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1) log1p(pmax(x,0)) else x
  z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
Z$rum_T<-factor(d0$rum_T,levels=c("Low","Mid","High"))
cat("── c_hw_soap 과 주모형 공변량의 상관 (|r|>0.35) ──\n")
for(v in cont_cov){r<-cor(Z$c_hw_soap,Z[[v]]);if(abs(r)>0.35)cat(sprintf("   %-20s r=%+.3f\n",v,r))}
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"));yy<-2011:2024
fitmod<-function(cont,useT=TRUE){
  dd<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
  dd$cases[is.na(dd$cases)]<-0;dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
  dd$area<-match(dd$region,regs);dd$time<-match(dd$year,yy);dd$id<-1:nrow(dd)
  dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High"))
  tm<-c(if(useT)"rum_T",if(length(cont))sprintf("`%s`",cont))
  form<-as.formula(paste0("cases~",paste(tm,collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
sh<-function(m,tag){fx<-m$summary.fixed;r<-fx["c_hw_soap",]
  cat(sprintf("  %-34s IRR %.3f (%.3f-%.3f)%s  DIC %.1f\n",tag,exp(r$mean),exp(r$`0.025quant`),exp(r$`0.975quant`),
    ifelse(r$`0.025quant`>0|r$`0.975quant`<0," ★",""),m$dic$dic))
  data.frame(step=tag,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),
    cred=(r$`0.025quant`>0|r$`0.975quant`<0),DIC=round(m$dic$dic,1))}
R<-list()
cat("\n═══ 단계적 축적: 어디서 부호가 뒤집히는가 ═══\n")
R[[1]]<-sh(fitmod("c_hw_soap",useT=FALSE),"① 단독 (BYM2 + soap)")
R[[2]]<-sh(fitmod("c_hw_soap",useT=TRUE), "② + 반추동물 tertile")
GRP<-list("③ + 농업토지(D2)"=c("d_pasture","d_field","d_paddy","e_forest"),
          "④ + 물위생(D3)"=c("a_septic","a_sewer","b_gw_use","b_gw_unfit"),
          "⑤ + 위생식품(D4)"=c("h_fastfood","h_pub","h_selfhealth"),
          "⑥ + SES(D5) = 주모형"=c("g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder"))
acc<-c()
for(nm in names(GRP)){acc<-c(acc,GRP[[nm]]);R[[length(R)+1]]<-sh(fitmod(c("c_hw_soap",acc),useT=TRUE),nm)}
cat("\n═══ 단일 변수 제거 검사 (주모형에서 하나씩 빼기) ═══\n")
for(v in cont_cov){m<-fitmod(c("c_hw_soap",setdiff(cont_cov,v)),useT=TRUE);fx<-m$summary.fixed;r<-fx["c_hw_soap",]
  irr<-exp(r$mean); if(irr<1.30) cat(sprintf("  %-20s 제거 시 → IRR %.3f (%.3f-%.3f)%s  ← 크게 감쇠\n",v,irr,exp(r$`0.025quant`),exp(r$`0.975quant`),ifelse(r$`0.025quant`>0|r$`0.975quant`<0," ★","")))}
out<-do.call(rbind,R); write.csv(out,file.path(ROOT,"results_ST","MASTER_ehec_soap_signflip.csv"),row.names=FALSE)
cat("\n저장: results_ST/MASTER_ehec_soap_signflip.csv\n")
