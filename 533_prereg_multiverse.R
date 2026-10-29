# 533 — 사전지정 대규모 멀티버스. 84개 변수 전부 대상. 선정 단계 없음.
#  사양 = 변수조합 × 도농축 × 검출통제 × 반추동물형태 × 공간구조 × 인접그래프. 판정: ≥50% credible + 방향일관 = 강건.
#  실행: N_SPEC(환경변수) 개 사양을 층화 무작위로 뽑아 적합. 체크포인트 저장으로 중단 복구 가능.
#  ★2026-10-21 v2b(이전 판 .bak_pre_v2b): 풀 결측 변수는 _v2_promote_261021/15_세트풀_제외변수.csv의 (3) 결정대로 v2_pool_check가 빼고, 사양 생성도 남은 풀(allv) 안에서만 뽑음
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(jsonlite)})
inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
N_SPEC<-as.integer(Sys.getenv("N_SPEC",unset="50")); SEED<-as.integer(Sys.getenv("SEED",unset="20260822")); set.seed(SEED)
BASE<-file.path(Sys.getenv("EHEC_DATA",unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT",unset="."); CACHE<-file.path(ROOT,"data_cache"); RS<-Sys.getenv("EHEC_RS",unset=file.path(ROOT,"results_ST"))
V2_CACHE<-CACHE;source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지, KDCA형 청주 2011-2013 마스킹(CJ_MASK=0이면 끔). 이전 판 .bak_pre_v2
OUT<-file.path(RS,sprintf("PREREG_MV_results_seed%d.csv",SEED)); SPECF<-file.path(RS,sprintf("PREREG_MV_specs_seed%d.csv",SEED))
DES<-fromJSON(file.path(CACHE,Sys.getenv("EHEC_DESIGN",unset="PREREG_DESIGN.json")))
d0<-read.csv(file.path(CACHE,Sys.getenv("EHEC_PANEL",unset="ST_prereg_panel.csv")),fileEncoding="UTF-8-BOM"); d0$region<-str_replace_all(d0$region,"\\s+","")
DOMS<-setdiff(names(DES$domains),"D7 detection"); DET<-DES$domains[["D7 detection"]]$vars
allv<-unique(c(unlist(lapply(DES$domains,function(x)x$vars)),"kchs_pc1","sex_ratio",DES$forced))
allv<-v2_pool_check(d0,allv,"533_prereg_multiverse.R")   # ★v2: 풀 변수 결측이면 멈춤(V2_POOL_DROP_NA=1이면 빼고 기록)
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")]
for(v in allv){x<-d0[[v]];x<-v2_nofill(x,v);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
# 반추동물 3분위
rum<-d0$beef_cattle_density+ifelse(is.na(d0$dairy_cattle_density),0,d0$dairy_cattle_density); qt<-quantile(rum,c(1/3,2/3),na.rm=TRUE)
Z$rum_T<-factor(cut(rum,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE),levels=c("Low","Mid","High"))
# 그래프 3종
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
fixnb<-function(nb){k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
  cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};make.sym.nb(nb)}
mkg<-function(nb){g<-tempfile(fileext=".graph");nb2INLA(g,nb);g}
G<-list(queen=mkg(fixnb(poly2nb(shp,queen=TRUE,snap=0.01))),rook=mkg(fixnb(poly2nb(shp,queen=FALSE,snap=0.01))),knn5=mkg(make.sym.nb(knn2nb(knearneigh(coords,k=5)))))
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv")
Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")];yy<-2011:2024
DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
DD<-v2_finish(DD,kdca=TRUE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD)

# ── 사양 생성(층화 무작위): 각 변수가 고르게 등장하도록 도메인별 추출 ──
DETSETS<-list(none=character(0),docs="physicians_per_1000_pop",lab=c("laboratory_medicine_specialists","pathology_specialists"),full=DET)
DETSETS<-lapply(DETSETS,function(s)intersect(s,allv))
gen_spec<-function(i){
  nd<-sample(min(3,length(DOMS)):length(DOMS),1); doms<-sample(DOMS,nd)
  vars<-unlist(lapply(doms,function(d){v<-intersect(DES$domains[[d]]$vars,allv);sample(v,min(length(v),sample(1:3,1)))}))   # ★2026-10-21 v2b: 풀에서 뺀 변수(15_세트풀_제외변수.csv)는 사양에 넣지 않음(넣으면 적합 실패가 조용히 버려짐)
  list(id=i,vars=vars,pc1=sample(c(TRUE,FALSE),1),det=sample(names(DETSETS),1),rum=sample(c("cont","tertile"),1),
       st=sample(c("bym2","bym2_int"),1),graph=sample(names(G),1))}
specs<-lapply(seq_len(N_SPEC),gen_spec)
write.csv(do.call(rbind,lapply(specs,function(s)data.frame(id=s$id,n_vars=length(s$vars),vars=paste(s$vars,collapse="|"),pc1=s$pc1,det=s$det,rum=s$rum,st=s$st,graph=s$graph))),SPECF,row.names=FALSE)
fit1<-function(s){
  cont<-unique(c(s$vars,if(s$pc1)"kchs_pc1","sex_ratio",DETSETS[[s$det]],DES$forced))   # DES$forced: 2026-09-30 문헌 풀(도시인구비율·재정자주도 모든 사양 고정)
  useT<-s$rum=="tertile"&&any(grepl("cattle",cont)); if(useT)cont<-setdiff(cont,c("beef_cattle_density","dairy_cattle_density"))
  tm<-c(if(useT)"rum_T",sprintf("`%s`",cont))
  re<-"+f(area,model='bym2',graph=G[[s$graph]],scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)"
  if(s$st=="bym2_int")re<-paste0(re,"+f(id,model='iid',hyper=pcp)")
  form<-as.formula(paste0("cases~",paste(tm,collapse="+"),re))
  m<-tryCatch(inla(form,family=Sys.getenv("EHEC_FAM",unset="nbinomial"),data=DD,E=DD$E_std,control.fixed=CF,control.compute=list(dic=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e)NULL)
  if(is.null(m))return(NULL)
  fx<-m$summary.fixed; fx<-fx[rownames(fx)!="(Intercept)",]
  data.frame(spec=s$id,term=rownames(fx),IRR=round(exp(fx$mean),4),lo=round(exp(fx$`0.025quant`),4),hi=round(exp(fx$`0.975quant`),4),cred=(fx$`0.025quant`>0|fx$`0.975quant`<0),DIC=round(m$dic$dic,1),stringsAsFactors=FALSE)}
done<-if(file.exists(OUT))unique(read.csv(OUT)$spec)else integer(0); t0<-Sys.time(); ok<-0
for(s in specs){ if(s$id%in%done)next
  r<-fit1(s); if(!is.null(r)){write.table(r,OUT,sep=",",row.names=FALSE,col.names=!file.exists(OUT),append=file.exists(OUT));ok<-ok+1}
  if(s$id%%10==0)cat(sprintf("  %d/%d 적합 · 수렴 %d · %.1f분\n",s$id,N_SPEC,ok,as.numeric(difftime(Sys.time(),t0,units="mins"))))}
cat(sprintf("\n완료: 사양 %d · 수렴 %d · %.1f분 · %s\n",N_SPEC,ok,as.numeric(difftime(Sys.time(),t0,units="mins")),basename(OUT)))
