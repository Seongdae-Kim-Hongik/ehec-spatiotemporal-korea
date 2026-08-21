# 512 — D6 건강행태 도메인 검정 + 도시성 축(KCHS PC1) 교란 검정.
#  Q1. 건강행태·예방행동 지표가 EHEC 위험과 연관되는가(개별 투입).
#  Q2. c_hw_soap의 멀티버스 78% 신호가 도시성 통제 후에도 남는가(=행태 신호 vs 도시성 대리).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(car)})
set.seed(20260821); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
cont_cov<-c("d_pasture","d_field","d_paddy","e_forest","a_septic","a_sewer","b_gw_use","b_gw_unfit","h_fastfood","h_pub","h_selfhealth","g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder")
D6<-c("c_flu_vax","c_infant_chk","c_nutrilabel","c_walking","c_smoking","c_drinking")
HW<-c("c_hw_soap","c_hw_toilet","c_hw_out","c_hw_meal")
ENG<-c(c_flu_vax="Influenza vaccination coverage",c_infant_chk="Infant health check-up rate",c_nutrilabel="Nutrition-label literacy",
 c_walking="Walking practice rate",c_smoking="Current smoking rate",c_drinking="Monthly drinking rate",
 c_hw_soap="Handwashing with soap",kchs_pc1="Rural-urban axis (survey PC1)")

d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv"))
d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
# D6 · PC1 병합(시간불변)
v6<-read.csv(file.path(CACHE,"D6_health_behaviour_values.csv"),fileEncoding="UTF-8-BOM")
p1<-read.csv(file.path(CACHE,"D6_kchs_pc1.csv"),fileEncoding="UTF-8-BOM")
fixr<-function(x){x<-str_replace_all(as.character(x),"\\s+","");ifelse(x=="인천시미추홀구","인천시남구",x)}
v6$X_rk<-fixr(v6[[1]]); p1$X_rk<-fixr(p1[[1]])
# 패널에 동명 컬럼(c_flu_vax 등)이 이미 있으면 조인이 .x/.y로 분리됨 → 게이트 통과본으로 교체
dupc<-intersect(c(D6,"kchs_pc1"),names(d0)); if(length(dupc)>0){cat("기존 컬럼 교체:",paste(dupc,collapse=","),"\n"); d0<-d0[,setdiff(names(d0),dupc)]}
d0<-d0%>%left_join(v6[,c("X_rk",D6)],by=c("region"="X_rk"))%>%left_join(p1[,c("X_rk","kchs_pc1")],by=c("region"="X_rk"))
cat(sprintf("병합 확인: D6 결측 %s · PC1 결측 %d\n",paste(sapply(D6,function(v)sum(is.na(d0[[v]]))),collapse="/"),sum(is.na(d0$kchs_pc1))))

d0$ruminant<-d0$d_cattle+d0$d_dairy
qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in unique(c(cont_cov,D6,HW,"kchs_pc1"))){ if(is.null(d0[[v]]))next
  x<-d0[[v]]; x[is.na(x)]<-mean(x,na.rm=TRUE)
  xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1) log1p(pmax(x,0)) else x
  z<-scale(xv)[,1]; z[!is.finite(z)]<-0; Z[[v]]<-z}
Z$rum_T<-factor(d0$rum_T,levels=c("Low","Mid","High"))

shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3))
for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"));yy<-2011:2024
fitmod<-function(cont){
  dd<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
  dd$cases[is.na(dd$cases)]<-0; dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
  dd$area<-match(dd$region,regs);dd$time<-match(dd$year,yy);dd$id<-1:nrow(dd)
  dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High"))
  form<-as.formula(paste0("cases~rum_T+",paste(sprintf("`%s`",cont),collapse="+"),
    "+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
grab<-function(m,vs,tag){fx<-m$summary.fixed;do.call(rbind,lapply(vs,function(v) if(v%in%rownames(fx)){r<-fx[v,]
  data.frame(model=tag,term=v,label=ifelse(is.na(ENG[v]),v,ENG[v]),IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),
   hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),DIC=round(m$dic$dic,1))} else NULL))}

R<-list()
cat("\n═══ Q1. D6 개별 투입 (기존 15변수 + 1) ═══\n")
for(v in D6){ m<-fitmod(c(cont_cov,v)); r<-grab(m,v,paste0("A:",v)); R[[length(R)+1]]<-r
  cat(sprintf("  %-14s IRR %.3f (%.3f-%.3f)%s  DIC %.1f\n",v,r$IRR,r$lo,r$hi,ifelse(r$cred," ★",""),r$DIC))}
cat("\n═══ Q2. c_hw_soap × 도시성 통제 ═══\n")
mB1<-fitmod(c(cont_cov,"c_hw_soap")); rB1<-grab(mB1,"c_hw_soap","B1:soap"); R[[length(R)+1]]<-rB1
cat(sprintf("  통제 없음      IRR %.3f (%.3f-%.3f)%s  DIC %.1f\n",rB1$IRR,rB1$lo,rB1$hi,ifelse(rB1$cred," ★",""),rB1$DIC))
mB2<-fitmod(c(cont_cov,"c_hw_soap","kchs_pc1")); rB2<-grab(mB2,c("c_hw_soap","kchs_pc1"),"B2:soap+PC1"); R[[length(R)+1]]<-rB2
for(i in 1:nrow(rB2)) cat(sprintf("  +도시성 통제   %-10s IRR %.3f (%.3f-%.3f)%s  DIC %.1f\n",rB2$term[i],rB2$IRR[i],rB2$lo[i],rB2$hi[i],ifelse(rB2$cred[i]," ★",""),rB2$DIC[i]))
cat("\n═══ Q3. D6 전체 + 도시성 (VIF>10 제거) ═══\n")
cand<-c(cont_cov,D6,"kchs_pc1"); vd<-Z[,cand]; vd$.y<-rnorm(nrow(vd)); keep<-cand
repeat{ lmf<-lm(as.formula(paste(".y~",paste(sprintf("`%s`",keep),collapse="+"))),data=vd)
  v<-tryCatch(car::vif(lmf),error=function(e)NULL); if(is.null(v)||max(v)<10)break
  drop<-names(which.max(v)); drop<-gsub("`","",drop); cat(sprintf("  VIF>10 제거: %s (%.1f)\n",drop,max(v))); keep<-setdiff(keep,drop)}
mC<-fitmod(keep); rC<-grab(mC,intersect(keep,c(D6,"kchs_pc1")),"C:full"); R[[length(R)+1]]<-rC
for(i in seq_len(nrow(rC))) cat(sprintf("  %-14s IRR %.3f (%.3f-%.3f)%s\n",rC$term[i],rC$IRR[i],rC$lo[i],rC$hi[i],ifelse(rC$cred[i]," ★","")))
cat(sprintf("  모형 DIC %.1f\n",rC$DIC[1]))
out<-do.call(rbind,R); write.csv(out,file.path(ROOT,"results_ST","MASTER_ehec_D6_urbanicity.csv"),row.names=FALSE)
cat("\n저장: results_ST/MASTER_ehec_D6_urbanicity.csv\n")
