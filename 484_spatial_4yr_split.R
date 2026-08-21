# 484 - Calendar-window split: 4-year pure-spatial BYM2 models (no temporal terms).
#  통합(2011-2024)은 시기별 위험요인 변화를 평균내 희석 → 4년 구간별로 순수 공간모형.
#  구간: 2011-2014 / 2015-2018 / 2019-2022 / 2023-2024. 각 구간 O=사례합, offset=E_std합.
#  공변량=Enriched A(주모형과 동일 16개); per-SD 표준화는 전체기간 기준(구간간 IRR 비교 가능)→구간평균.
#  시간항(RW1)·상호작용 제거 = 순수 공간 BYM2만. 각 구간 credible 요인 + 지역RR 지도용 산출.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260805); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
PERIODS<-list("2011-2014"=c(2011,2014),"2015-2018"=c(2015,2018),"2019-2022"=c(2019,2022),"2023-2024"=c(2023,2024))
DOM<-list("1. Livestock / ruminant reservoir"=c("rum_T","d_pasture"),
 "2. Agriculture & land use"=c("d_field","d_paddy","e_forest"),
 "3. Water & sanitation"=c("a_septic","a_sewer","b_gw_use","b_gw_unfit"),
 "4. Hygiene & food handling"=c("h_fastfood","h_pub","h_selfhealth"),
 "5. Socioeconomic & demographic"=c("g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder"))
ENG<-c(rum_T="Ruminant density (cattle+dairy), tertile",d_pasture="Pasture area",d_field="Upland-field (cropland) area",e_forest="Forest area",
 a_septic="Septic-tank population",a_sewer="Public-sewerage coverage",b_gw_use="Groundwater use",b_gw_unfit="Unfit groundwater",
 h_fastfood="Fast-food / eating-out outlets",h_pub="Bars/pubs",d_paddy="Paddy-field area",h_selfhealth="Poor self-rated health",
 g_solo_elder="Elderly living alone",g_sexratio="Sex ratio",g_fiscal="Fiscal independence",g_basiclivelihood="Basic-livelihood recipients")
cont_cov<-c("d_pasture","d_field","d_paddy","e_forest","a_septic","a_sewer","b_gw_use","b_gw_unfit","h_fastfood","h_pub","h_selfhealth",
            "g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder")

# ── 데이터: 전체기간 기준 per-SD 표준화 (구간간 비교 가능) ──
d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv")); d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
d0$ruminant<-d0$d_cattle+d0$d_dairy; qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in cont_cov){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
Z$rum_T<-factor(d0$rum_T,levels=c("Low","Mid","High"))

# ── 그래프 (418 동일) ──
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"))

fit_period<-function(lo,hi){
  # 구간 집계: O=사례합, E=E_std합, 공변량=구간평균 z
  cc2<-case%>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(O=sum(cases),.groups="drop")
  ee <-Est %>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(E=sum(E_std),.groups="drop")
  zz <-Z   %>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(across(all_of(cont_cov),mean),.groups="drop")
  # 반추동물 tertile은 시간불변 → 첫 값
  rt <-Z%>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(rum_T=first(rum_T),.groups="drop")
  dd<-data.frame(region=regs)%>%left_join(cc2,by="region")%>%left_join(ee,by="region")%>%left_join(zz,by="region")%>%left_join(rt,by="region")
  dd$O[is.na(dd$O)]<-0; dd$E[is.na(dd$E)|dd$E<=0]<-min(dd$E[dd$E>0],na.rm=TRUE)
  for(v in cont_cov){dd[[v]][is.na(dd[[v]])]<-0}
  dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High")); dd$rum_T[is.na(dd$rum_T)]<-"Low"
  dd$area<-match(dd$region,regs)
  terms<-c("rum_T",sprintf("`%s`",cont_cov))
  form<-as.formula(paste0("O~",paste(terms,collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)"))
  m<-inla(form,family="nbinomial",data=dd,E=dd$E,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))
  list(m=m,N=sum(dd$O),dd=dd)
}
grabrow<-function(fx,cv){
  nm<-if(cv=="rum_T:Mid")"rum_TMid" else if(cv=="rum_T:High")"rum_THigh" else cv
  if(!nm%in%rownames(fx))return(NULL); r<-fx[nm,]
  data.frame(cov=cv,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}

allres<-data.frame(); fitstat<-data.frame(); regRR<-data.frame(region=regs)
for(pn in names(PERIODS)){
  pr<-PERIODS[[pn]]; cat(sprintf("\n===== %s 적합 중 =====\n",pn))
  R1<-fit_period(pr[1],pr[2]); fx<-R1$m$summary.fixed
  out<-list()
  for(dm in names(DOM)) for(cv in DOM[[dm]]){
    if(cv=="rum_T"){for(lv in c("rum_T:Mid","rum_T:High")){rr<-grabrow(fx,lv);if(!is.null(rr)){rr$domain<-dm;rr$label<-paste0(ENG["rum_T"],ifelse(lv=="rum_T:Mid"," — Mid vs Low"," — High vs Low"));out[[length(out)+1]]<-rr}}}
    else{rr<-grabrow(fx,cv);if(!is.null(rr)){rr$domain<-dm;rr$label<-ENG[cv];out[[length(out)+1]]<-rr}}}
  rp<-do.call(rbind,out); rp$period<-pn; rp$N<-R1$N; allres<-rbind(allres,rp)
  phi<-R1$m$summary.hyperpar["Phi for area","mean"]
  fitstat<-rbind(fitstat,data.frame(period=pn,N=R1$N,DIC=round(R1$m$dic$dic,1),WAIC=round(R1$m$waic$waic,1),phi=round(phi,3)))
  uu<-R1$m$summary.random$area[1:length(regs),]
  regRR[[paste0("RR_",pn)]]<-round(exp(uu$mean),3)
  regRR[[paste0("hi_",pn)]]<-as.integer(uu$`0.025quant`>0)
  cat(sprintf("  N=%d DIC=%.1f phi=%.3f · credible: %s\n",R1$N,R1$m$dic$dic,phi,paste(rp$label[rp$cred],collapse=" | ")))
}
write.csv(allres[,c("period","N","domain","cov","label","IRR","lo","hi","cred")],file.path(ROOT,"results_ST","MASTER_ehec_4yr_split.csv"),row.names=FALSE)
write.csv(fitstat,file.path(ROOT,"results_ST","MASTER_ehec_4yr_fitstat.csv"),row.names=FALSE)
write.csv(regRR,file.path(ROOT,"results_ST","MASTER_ehec_4yr_regionRR.csv"),row.names=FALSE)
cat("\n===== 구간별 credible 요인 요약 =====\n")
for(pn in names(PERIODS)){s<-allres[allres$period==pn&allres$cred,];cat(sprintf("[%s] N=%d: %s\n",pn,s$N[1],if(nrow(s))paste(sprintf("%s %.3f(%.3f-%.3f)",s$label,s$IRR,s$lo,s$hi),collapse=" | ")else "없음"))}
cat("\n저장: MASTER_ehec_4yr_split.csv · _fitstat.csv · _regionRR.csv\n")
