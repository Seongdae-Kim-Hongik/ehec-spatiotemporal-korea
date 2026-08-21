# 498 - Pre- versus post-COVID-19 spatial split (2011-2019 vs 2020-2024); the pooled model is retained.
#  Period-specific models are descriptive; the formal test of effect change is the covariate-by-period interaction.
#  단서: "EHEC는 코비드 나인틴 때 영향을 거의 안 받았었죠" → 차이 없을 가능성도 사전 예상.
#  구성: [A] 구간별 순수 공간 BYM2(484 미러, 구간만 2개) + [B] 통합모형 + 공변량×COVID 상호작용 정식검정(485 미러).
#  ※ 7/20 교훈: 구간 비교만으로 시기변화를 판정하지 않는다(검정력 손실 아티팩트). [B]가 확정 검정.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260820); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
PERIODS<-list("2011-2019 (pre-COVID)"=c(2011,2019),"2020-2024 (COVID/post)"=c(2020,2024))
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

# ── 데이터: 전체기간 기준 per-SD 표준화(구간간 IRR 비교 가능) ──
d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv")); d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
d0$ruminant<-d0$d_cattle+d0$d_dairy; qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in cont_cov){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
Z$rum_T<-factor(d0$rum_T,levels=c("Low","Mid","High"))

# ── 그래프(418/484 동일) ──
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)))
pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"))

# ══════════ [A] 구간별 순수 공간 BYM2 ══════════
fit_period<-function(lo,hi){
  cc2<-case%>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(O=sum(cases),.groups="drop")
  ee <-Est %>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(E=sum(E_std),.groups="drop")
  zz <-Z   %>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(across(all_of(cont_cov),mean),.groups="drop")
  rt <-Z   %>%filter(year>=lo,year<=hi)%>%group_by(region)%>%summarise(rum_T=first(rum_T),.groups="drop")
  dd<-data.frame(region=regs)%>%left_join(cc2,by="region")%>%left_join(ee,by="region")%>%left_join(zz,by="region")%>%left_join(rt,by="region")
  dd$O[is.na(dd$O)]<-0; dd$E[is.na(dd$E)|dd$E<=0]<-min(dd$E[dd$E>0],na.rm=TRUE)
  for(v in cont_cov){dd[[v]][is.na(dd[[v]])]<-0}
  dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High")); dd$rum_T[is.na(dd$rum_T)]<-"Low"
  dd$area<-match(dd$region,regs)
  terms<-c("rum_T",sprintf("`%s`",cont_cov))
  form<-as.formula(paste0("O~",paste(terms,collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)"))
  m<-inla(form,family="nbinomial",data=dd,E=dd$E,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))
  list(m=m,N=sum(dd$O),nz=sum(dd$O>0),dd=dd)
}
grabrow<-function(fx,cv){
  nm<-if(cv=="rum_T:Mid")"rum_TMid" else if(cv=="rum_T:High")"rum_THigh" else cv
  if(!nm%in%rownames(fx))return(NULL); r<-fx[nm,]
  data.frame(cov=cv,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}

allres<-data.frame(); fitstat<-data.frame(); regRR<-data.frame(region=regs)
for(pn in names(PERIODS)){
  pr<-PERIODS[[pn]]; cat(sprintf("\n===== [A] %s 적합 중 =====\n",pn))
  R1<-fit_period(pr[1],pr[2]); fx<-R1$m$summary.fixed
  out<-list()
  for(dm in names(DOM)) for(cv in DOM[[dm]]){
    if(cv=="rum_T"){for(lv in c("rum_T:Mid","rum_T:High")){rr<-grabrow(fx,lv);if(!is.null(rr)){rr$domain<-dm;rr$label<-paste0(ENG["rum_T"],ifelse(lv=="rum_T:Mid"," — Mid vs Low"," — High vs Low"));out[[length(out)+1]]<-rr}}}
    else{rr<-grabrow(fx,cv);if(!is.null(rr)){rr$domain<-dm;rr$label<-ENG[cv];out[[length(out)+1]]<-rr}}}
  rp<-do.call(rbind,out); rp$period<-pn; rp$N<-R1$N; allres<-rbind(allres,rp)
  phi<-R1$m$summary.hyperpar["Phi for area","mean"]
  fitstat<-rbind(fitstat,data.frame(period=pn,N=R1$N,districts_nonzero=R1$nz,DIC=round(R1$m$dic$dic,1),WAIC=round(R1$m$waic$waic,1),phi=round(phi,3)))
  uu<-R1$m$summary.random$area[1:length(regs),]
  tag<-ifelse(grepl("pre",pn),"pre","post")
  regRR[[paste0("RR_",tag)]]<-round(exp(uu$mean),3)
  regRR[[paste0("hi_",tag)]]<-as.integer(uu$`0.025quant`>0)
  cat(sprintf("  N=%d (비제로 %d개 시군구) DIC=%.1f phi=%.3f\n  credible: %s\n",R1$N,R1$nz,R1$m$dic$dic,phi,
      if(any(rp$cred))paste(sprintf("%s %.3f(%.3f-%.3f)",rp$label[rp$cred],rp$IRR[rp$cred],rp$lo[rp$cred],rp$hi[rp$cred]),collapse=" | ") else "없음"))
}
write.csv(allres[,c("period","N","domain","cov","label","IRR","lo","hi","cred")],file.path(ROOT,"results_ST","MASTER_ehec_covid_split.csv"),row.names=FALSE)
write.csv(fitstat,file.path(ROOT,"results_ST","MASTER_ehec_covid_fitstat.csv"),row.names=FALSE)
write.csv(regRR,file.path(ROOT,"results_ST","MASTER_ehec_covid_regionRR.csv"),row.names=FALSE)

# ══════════ [B] 통합모형 + 공변량×COVID 상호작용 정식검정 ══════════
cat("\n\n===== [B] 통합모형 + COVID 상호작용 정식검정 =====\n")
dd<-expand.grid(region=regs,year=2011:2024,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
dd$cases[is.na(dd$cases)]<-0
dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
for(v in cont_cov){dd[[v]][is.na(dd[[v]])]<-0}
dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High")); dd$rum_T[is.na(dd$rum_T)]<-"Low"
dd$area<-match(dd$region,regs); dd$yr<-dd$year-2010; dd$area_time<-seq_len(nrow(dd))
dd$post<-as.integer(dd$year>=2020)   # COVID 이후 지시자
INT<-c("h_pub","g_fiscal","h_fastfood")
for(v in INT) dd[[paste0(v,"_x_post")]]<-dd[[v]]*dd$post
terms<-c("rum_T",sprintf("`%s`",cont_cov),"post",paste0(INT,"_x_post"))
form<-as.formula(paste0("cases~",paste(terms,collapse="+"),
  "+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)",
  "+f(yr,model='rw1',scale.model=TRUE,constr=TRUE,hyper=pcp)",
  "+f(area_time,model='iid',hyper=pcp)"))
m<-inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,
        control.compute=list(dic=TRUE,waic=TRUE,config=TRUE),control.inla=list(control.vb=list(emergency=50)))
fx<-m$summary.fixed; res<-data.frame()
show<-function(nm,lab){ if(!nm%in%rownames(fx))return(invisible(NULL)); r<-fx[nm,]
  cred<-(r$`0.025quant`>0|r$`0.975quant`<0)
  cat(sprintf("  %-26s IRR %.3f (%.3f-%.3f)%s\n",lab,exp(r$mean),exp(r$`0.025quant`),exp(r$`0.975quant`),ifelse(cred,"  *credible","")))
  res<<-rbind(res,data.frame(term=nm,label=lab,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=cred))}
cat(sprintf("[통합모형 DIC=%.1f (주모형 pooled 5335.0 대비)]\n",m$dic$dic))
cat("[주효과 = 2011-2019 기준]\n")
show("h_pub","Bars/pubs"); show("g_fiscal","Fiscal independence"); show("h_fastfood","Fast-food outlets")
cat("[상호작용 = COVID 이후(2020-) 추가 배율] ← 시기 변화의 정식검정\n")
show("h_pub_x_post","Bars/pubs x post-COVID"); show("g_fiscal_x_post","Fiscal x post-COVID"); show("h_fastfood_x_post","Fast-food x post-COVID")
cat("[COVID 시기 자체]\n"); show("post","Post-COVID period")
# 선형결합: COVID 이후 총효과 = 주효과 + 상호작용
lcs<-list(); for(v in INT){ l<-list(); l[[v]]<-1; l[[paste0(v,"_x_post")]]<-1; lcs[[paste0(v,"_post_total")]]<-l }
lc<-do.call(c,lapply(names(lcs),function(k){x<-do.call(inla.make.lincomb,lcs[[k]]);names(x)<-k;x}))
m2<-inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,lincomb=lc,
         control.compute=list(dic=TRUE),control.inla=list(control.vb=list(emergency=50)))
cat("\n[선형결합] COVID 이후 총효과\n")
for(k in names(lcs)){ r<-m2$summary.lincomb.derived[k,]
  cred<-(r$`0.025quant`>0|r$`0.975quant`<0)
  cat(sprintf("  %-26s IRR %.3f (%.3f-%.3f)%s\n",k,exp(r$mean),exp(r$`0.025quant`),exp(r$`0.975quant`),ifelse(cred,"  *credible","")))
  res<-rbind(res,data.frame(term=k,label=paste0(k," (post-COVID total)"),IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=cred))}
write.csv(res,file.path(ROOT,"results_ST","MASTER_ehec_covid_interaction.csv"),row.names=FALSE)
cat("\n저장: MASTER_ehec_covid_split.csv · _fitstat.csv · _regionRR.csv · _interaction.csv\n")
