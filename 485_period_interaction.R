# 485 — /333 3사 공통처방: 통합 시공간모형 + 공변량×시기 상호작용 정식검정.
#  4년분할(484)은 검정력 손실로 시기변화를 판정 못함 → 통합모형(전체 데이터·검정력 유지)에
#  '초기(2011-2014) vs 이후' 상호작용을 넣어 효과의 시기차이를 직접 검정한다.
#  대상: 임야(분할의 유일한 시기변화 후보) + 재정자립도·주점(헤드라인 확인용).
#  구조는 418 주모형과 동일(NB BYM2+RW1+iid, offset=E_std, Enriched A 16공변량).
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260805); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
cont_cov<-c("d_pasture","d_field","d_paddy","e_forest","a_septic","a_sewer","b_gw_use","b_gw_unfit","h_fastfood","h_pub","h_selfhealth",
            "g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder")

d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv")); d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
d0$ruminant<-d0$d_cattle+d0$d_dairy; qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in cont_cov){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
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

dd<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(Z,by=c("region","year"))
dd$cases[is.na(dd$cases)]<-0; dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
dd$area<-match(dd$region,regs);dd$time<-match(dd$year,yy);dd$id<-1:nrow(dd)
dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High"))
dd$early<-as.integer(dd$year<=2014)   # 초기(2011-2014) 지시자
# 상호작용 항 (per-SD z × early)
for(v in c("e_forest","g_fiscal","h_pub")) dd[[paste0(v,"_x_early")]]<-dd[[v]]*dd$early

terms<-c("rum_T",sprintf("`%s`",cont_cov),"early","e_forest_x_early","g_fiscal_x_early","h_pub_x_early")
form<-as.formula(paste0("cases~",paste(terms,collapse="+"),
  "+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
m<-inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))
fx<-m$summary.fixed
show<-function(nm,lab){if(!nm%in%rownames(fx))return(invisible());r<-fx[nm,]
  cat(sprintf("  %-24s IRR %.3f (%.3f-%.3f)%s\n",lab,exp(r$mean),exp(r$`0.025quant`),exp(r$`0.975quant`),
      ifelse(r$`0.025quant`>0|r$`0.975quant`<0," *","")))}
cat(sprintf("===== 통합모형 + 시기 상호작용 (DIC=%.1f, 주모형 5335.4 대비) =====\n",m$dic$dic))
cat("[주효과 (2015-2024 기준)]\n")
show("e_forest","Forest area");show("g_fiscal","Fiscal independence");show("h_pub","Bars/pubs")
cat("[상호작용 (초기 2011-2014에서의 추가 배율)]\n")
show("e_forest_x_early","Forest x early");show("g_fiscal_x_early","Fiscal x early");show("h_pub_x_early","Bars/pubs x early")
cat("[초기 시기 자체]\n"); show("early","Early period (11-14)")
# 초기 시기의 임야 총효과 = e_forest + e_forest_x_early (선형결합)
lc<-inla.make.lincomb(e_forest=1,e_forest_x_early=1); names(lc)<-"forest_early_total"
m2<-inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,lincomb=lc,control.inla=list(control.vb=list(emergency=50)))
r<-m2$summary.lincomb.derived["forest_early_total",]
cat(sprintf("\n[선형결합] 초기 임야 총효과: IRR %.3f (%.3f-%.3f)%s\n",exp(r$mean),exp(r$`0.025quant`),exp(r$`0.975quant`),
    ifelse(r$`0.025quant`>0|r$`0.975quant`<0," *","")))
res<-data.frame(term=rownames(fx),IRR=round(exp(fx$mean),3),lo=round(exp(fx$`0.025quant`),3),hi=round(exp(fx$`0.975quant`),3),
                cred=(fx$`0.025quant`>0|fx$`0.975quant`<0))
res<-rbind(res,data.frame(term="forest_early_total",IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0)))
write.csv(res,file.path(ROOT,"results_ST","MASTER_ehec_period_interaction.csv"),row.names=FALSE)
cat("\n저장: MASTER_ehec_period_interaction.csv\n")
