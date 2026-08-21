# 508 — 검출(detection) 대리변수 통제 후 반추동물 저장소 신호 검정
#  근거: 멀티버스에서 유일한 '강건' 변수가 자기보고 위생지표(c_hw_soap, 78%, IRR 1.269, 반직관 방향)였고,
#        같은 도메인 손씻기 3종은 credible 0% → 검출/감시 강도 신호로 판단.
#        반추동물은 credible 29-33%이나 방향 일관 97-100% → 검출 통제 시 신호가 드러나는지 사전지정 검정.
#  ★ 사후 선택 아님: PREREGISTERED §2에서 h_docs·h_selfhealth·c_flu_vax를 '검출 대리변수, 음성대조 전용'으로 이미 분류했다.
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260821); inla.setOption(num.threads="2:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
DET<-c("h_docs","h_selfhealth","c_flu_vax","c_hw_soap")     # 검출/감시 강도 대리
EXPO<-c("d_cattle","d_dairy","d_pasture")                    # 저장소(관심 노출)
CTRL<-c("d_field","b_gw_use","g_fiscal","e_pop_dens")        # 도메인 대표 통제(멀티버스 상위)
d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv"));d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in unique(c(DET,EXPO,CTRL))){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE)
  xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1) log1p(pmax(x,0)) else x
  z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
d0$ruminant<-d0$d_cattle+d0$d_dairy
qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
Z$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m)
 D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,]
 a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)))
pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"))
dd<-expand.grid(region=regs,year=2011:2024,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%
    left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
dd$cases[is.na(dd$cases)]<-0;dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
for(v in unique(c(DET,EXPO,CTRL))) dd[[v]][is.na(dd[[v]])]<-0
dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High"));dd$rum_T[is.na(dd$rum_T)]<-"Low"
dd$area<-match(dd$region,regs);dd$yr<-dd$year-2010;dd$area_time<-seq_len(nrow(dd))
fit<-function(vars,useT,tag){
  terms<-c(if(useT)"rum_T" else NULL,sprintf("`%s`",vars))
  f<-as.formula(paste0("cases~",paste(terms,collapse="+"),
    "+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)",
    "+f(yr,model='rw1',scale.model=TRUE,constr=TRUE,hyper=pcp)+f(area_time,model='iid',hyper=pcp)"))
  m<-inla(f,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),
          control.inla=list(control.vb=list(emergency=50)))
  fx<-m$summary.fixed
  out<-do.call(rbind,lapply(rownames(fx),function(nm){r<-fx[nm,]
    data.frame(model=tag,term=nm,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),
      hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),DIC=round(m$dic$dic,1),stringsAsFactors=FALSE)}))
  cat(sprintf("\n[%s] DIC=%.1f\n",tag,m$dic$dic))
  for(i in seq_len(nrow(out))) if(out$term[i]!="(Intercept)")
    cat(sprintf("   %-22s IRR %.3f (%.3f-%.3f)%s\n",out$term[i],out$IRR[i],out$lo[i],out$hi[i],ifelse(out$cred[i],"  *","")))
  out}
res<-rbind(
 fit(c(EXPO,CTRL),FALSE,"M0_검출통제없음_연속"),
 fit(c(EXPO,CTRL,DET),FALSE,"M1_검출통제_연속"),
 fit(c(setdiff(EXPO,c("d_cattle","d_dairy")),CTRL),TRUE,"M2_검출통제없음_3분위"),
 fit(c(setdiff(EXPO,c("d_cattle","d_dairy")),CTRL,DET),TRUE,"M3_검출통제_3분위"))
write.csv(res,file.path(ROOT,"results_ST","MASTER_detection_adjusted.csv"),row.names=FALSE)
cat("\n저장: MASTER_detection_adjusted.csv\n")
