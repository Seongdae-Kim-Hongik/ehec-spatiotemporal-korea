# s3_child04_refit_v2d.R — (2026-10-29 /333 8번) 표 S3 재적합: 최종 20 공변량 + 0-4세 인구 비율(v2 코드 결합 패널 _internal/s3_refit/ST_children_pct_v2d.csv) 1회.
#  540s_sensitivity.R(v2c)와 같은 셋업(패널·E·그래프·사전·청주 마스킹)에 children_0_4만 더함. 출력 _internal/s3_refit/ (공유 results_ST는 건드리지 않음).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(car)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
H<-file.path(ROOT,"JKMS_FBD","JKMS_EHEC_v2d_261029");RS<-file.path(H,"_internal","s3_refit");dir.create(RS,showWarnings=FALSE,recursive=TRUE)
source(file.path(ROOT,"_v2_common.R"))
T2<-read.csv(file.path(ROOT,"results_ST","MASTER_final_table2.csv"));cont<-T2$cov;stopifnot(length(cont)==20);ENG<-setNames(T2$label,T2$cov)
d0<-read.csv(file.path(CACHE,"ST_prereg_panel.csv"),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
CH<-read.csv(file.path(RS,"ST_children_pct_v2d.csv"),fileEncoding="UTF-8-BOM",colClasses=c(sigungu_cd="character"))
d0$sigungu_cd<-v2_code(d0$region);d0<-left_join(d0,CH[,c("sigungu_cd","year","children_0_4")],by=c("sigungu_cd","year"))
cat(sprintf("패널 %d행, 연도 %s, children 결측 %d\n",nrow(d0),paste(range(d0$year),collapse="-"),sum(is.na(d0$children_0_4))))
cv<-c(cont,"children_0_4");ENG["children_0_4"]<-"Children aged 0-4 years, % of population"
Z<-d0[,c("region","year")];for(v in cv)Z[[v]]<-v2_zb(d0[[v]],v)
vd<-Z[,cv];vd$.y<-rnorm(nrow(vd));v<-car::vif(lm(as.formula(paste(".y~",paste(sprintf("`%s`",cv),collapse="+"))),data=vd));cat(sprintf("VIF 최대 %.2f (%s); children VIF %.2f\n",max(v),names(which.max(v)),v["`children_0_4`"]))
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
local({a<-as.numeric(sf::st_area(shp))/1e6;stopifnot(nrow(shp)==229,!anyDuplicated(shp$region),a[shp$region=="인천시남구"]<35,a[shp$region=="인천시남동구"]>45)})
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")];yy<-2011:2024
DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
DD<-v2_finish(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD);cat(sprintf("모형 칸 %d (사례 NA %d)\n",nrow(DD),sum(is.na(DD$cases))))
form<-as.formula(paste0("cases~",paste(sprintf("`%s`",cv),collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
m<-inla(form,family="nbinomial",data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))
fx<-m$summary.fixed
res<-do.call(rbind,lapply(cv,function(v){r<-fx[v,];pp<-1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]]);data.frame(cov=v,label=ENG[v],IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(pp,3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}))
write.csv(res,file.path(RS,"MASTER_final_table2.csv"),row.names=FALSE)
writeLines(c(sprintf("phi=%.3f",m$summary.hyperpar["Phi for area","mean"]),sprintf("DIC=%.1f",m$dic$dic),sprintf("WAIC=%.1f",m$waic$waic),sprintf("nvar=%d",length(cv)),sprintf("n=%d",sum(!is.na(DD$cases)))),file.path(RS,"final_fit.txt"))
cat(sprintf("\n===== S3 재적합 %d변수 · DIC %.1f · WAIC %.1f · phi %.3f =====\n",length(cv),m$dic$dic,m$waic$waic,m$summary.hyperpar["Phi for area","mean"]))
for(i in 1:nrow(res))cat(sprintf("  %-40s IRR %.3f (%.3f-%.3f)%s\n",res$label[i],res$IRR[i],res$lo[i],res$hi[i],ifelse(res$cred[i]," *","")))
