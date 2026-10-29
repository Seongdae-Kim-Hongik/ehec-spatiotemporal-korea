# _hm_setup.R — (2026-09-30) 문헌 기반 풀 재구성 공용 셋업. HM_TAG = EHEC|SHIG|TYPH|SALM|NORO.
#  사례·기대건수·그래프·prior는 각 질환 기존 셋업 그대로(EHEC는 540과 동일 조립), 공변량만 마스터 패널(HM_panel.csv, 563_hm 게이트)로 교체.
#  kchs_pc1(손씻기 종합지표·조사지표 종합축, 질환별 정의는 _fbd_labels.py)은 질환별 기존 값 유지(MFDS 계열=NORO/SALM, KDCA 계열=SHIG/TYPH, EHEC=ST_prereg_panel).
#  ★2026-10-21 v2(이전 판 _hm_setup.R.bak_pre_v2, _v2_common.R): E·사례·공변량·kchs_pc1 sigungu_cd 결합, E 최솟값 채움 삭제(panel_gate), 공변량 평균 채움 금지
#   (결측 열은 Z에서 빠지고 모형에 넣으면 멈춤), KDCA형(EHEC·SHIG·TYPH) 청주 2011-2013 사례 NA 마스킹(CJ_MASK=0이면 끔). kchs_pc1은 v2 정본(229곳; 손씻기 원자료 교정 뒤 재계산).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
TAG<-Sys.getenv("HM_TAG");stopifnot(TAG%in%c("EHEC","SHIG","TYPH","SALM","NORO"))
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
source(file.path(ROOT,"_v2_common.R"))
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
zb<-function(d0,vars){v2_Z(d0,vars,"HM_panel")}   # ★v2: 평균 채움 없음(결측 열은 Z에서 빠지고 이름 출력)
if(TAG=="EHEC"){
  set.seed(20260822);inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
  shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
  ## 2026-10-21 재발 방지: final.dbf 인천 region 라벨(남구 23030 = 미추홀구 약 25 km², 남동구 23050 약 57 km²)이 면적과 맞는지 단언(2026-10-21 dbf 두 칸 교정, 백업 final.dbf.bak_pre_shpfix)
  local({a<-as.numeric(sf::st_area(shp))/1e6;stopifnot(nrow(shp)==229,!anyDuplicated(shp$region),a[shp$region=="인천시남구"]<35,a[shp$region=="인천시남동구"]>45)})
  coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
  nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
  cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
  pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
  YY<-2011:2024
  DD<-v2_DD(regs,YY,"ST_Estd_panel.csv","ST_EHEC_agecase.csv",kdca=TRUE,label="EHEC")   # ★v2: 코드 결합·게이트, E 최솟값 채움 삭제
  DD<-v2_cj_mask(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,YY);DD$id<-1:nrow(DD)
  pe<-read.csv(file.path(CACHE,"ST_prereg_panel.csv"),fileEncoding="UTF-8-BOM");pe$region<-str_replace_all(pe$region,"\\s+","");PC1<-v2_relabel(pe[,c("region","year","kchs_pc1")],"EHEC kchs_pc1")
  SP<-"f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)";TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(id,model='iid',hyper=pcp)";FULLRE<-paste(SP,TI,IT,sep="+")
  MAINFAM<-"nbinomial";SH<-file.path(ROOT,"results_ST")
} else {
  src<-c(SHIG="_shig_setup.R",TYPH="_typh_setup.R",SALM="_salm_setup.R",NORO="_noro_setup.R")[[TAG]]
  source(file.path(ROOT,src))
  pcf<-c(SHIG="SHIG_kchs_pc1.csv",TYPH="TYPH_kchs_pc1.csv",SALM="SALM_kchs_pc1.csv",NORO="NORO_kchs_pc1.csv")[[TAG]]
  PC1<-read.csv(file.path(CACHE,pcf),colClasses=c(sigungu_cd="character"));PC1<-v2_relabel(PC1,pcf)   # ★v2: 코드(229곳)
  DD<-DD[,c("region","year","sigungu_cd","cases","cases_obs","cj","E_std","area","time","id")]
  MAINFAM<-"zeroinflatednbinomial1"
}
hm<-read.csv(file.path(CACHE,"HM_panel.csv"),fileEncoding="UTF-8-BOM");hm$region<-str_replace_all(hm$region,"\\s+","");hm$sigungu_cd<-v2_code(hm$region)
HMETA<-read.csv(file.path(CACHE,"HM_panel_META.csv"),fileEncoding="UTF-8-BOM")
HMV<-setdiff(names(hm),c("region","year","pop","area","sigungu_cd"))
hm<-hm%>%left_join(PC1[,intersect(c("sigungu_cd","year","kchs_pc1"),names(PC1))],by=intersect(names(PC1),c("sigungu_cd","year")))   # ★v2: 코드 결합
ZH<-zb(hm,c(HMV,"kchs_pc1"))
DD<-v2_join(DD,ZH,"HM 공변량")
stopifnot(!anyNA(DD[,intersect(c(HMV,"kchs_pc1"),names(DD))]),"kchs_pc1"%in%names(DD))
fitf<-function(cv=character(0),fam=MAINFAM,cpo=TRUE,data=DD){
  v2_need(cv,data)   # ★v2
  fx<-if(length(cv))paste(sprintf("`%s`",cv),collapse="+") else "1"
  form<-as.formula(paste0("cases~",fx,"+",FULLRE))
  inla(form,family=fam,data=data,E=data$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,cpo=cpo,config=FALSE),control.inla=list(control.vb=list(emergency=50)))}
irr<-function(m,v){r<-m$summary.fixed[v,];data.frame(cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]]),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}
safefit<-function(...)tryCatch(fitf(...),error=function(e){message("  ✗ fit failed: ",conditionMessage(e));NULL})
cat(sprintf("[hm_setup %s] %d cells · cases %d · zero %.1f%% · HM covariates %d · main family %s\n",TAG,nrow(DD),sum(DD$cases,na.rm=TRUE),100*mean(DD$cases==0,na.rm=TRUE),length(HMV),MAINFAM))
