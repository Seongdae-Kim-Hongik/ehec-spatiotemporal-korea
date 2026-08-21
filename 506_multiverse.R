# 506 — EHEC 공간 결정요인 멀티버스 (사전지정: PREREGISTERED_SPECIFICATION_2026-08-21.md)
#  ★ 사후 선택 방지: 도메인·후보·판정규칙을 IRR 보기 전에 확정했고 이 스크립트는 그대로 집행만 한다.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260821); inla.setOption(num.threads="2:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")
N_SPEC<-as.integer(Sys.getenv("N_SPEC","300"))

DOM<-list(D1=c("d_cattle","d_dairy","d_pasture"),
          D2=c("d_field","d_paddy","e_forest"),
          D3=c("a_sewer","a_septic","a_water_supply","a_advtreat","b_gw_use","b_gw_unfit","b_gw_home","a_wastewater","a_orgload"),
          D4=c("c_hw_toilet","c_hw_meal","c_hw_out","c_hw_soap","h_fastfood","h_pub"),
          D5=c("g_fiscal","g_basiclivelihood","g_solo_elder","g_sexratio","g_foreign","e_urban_pop","e_rural_pop","e_pop_dens"))
ALLV<-unlist(DOM,use.names=FALSE)

d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv"))
d0<-d0[d0$year>=2011&d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
## §5-Z 가드 재확인
for(v in ALLV){ zr<-tapply(d0[[v]],d0$year,function(x) mean(x==0,na.rm=TRUE))
  if(any(zr>0.9,na.rm=TRUE)) stop(sprintf("[가드] %s 오염 잔존",v)) }
cat("§5-Z 가드 통과 ·",length(ALLV),"변수\n")

skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
mkZ<-function(logtr){Z<-d0[,c("region","year")]
  for(v in ALLV){x<-d0[[v]];x[is.na(x)]<-mean(x,na.rm=TRUE)
    xv<-if(logtr&&min(x,na.rm=TRUE)>=0&&skewf(x)>1) log1p(pmax(x,0)) else x
    z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z};Z}
Z_log<-mkZ(TRUE);Z_raw<-mkZ(FALSE)
d0$ruminant<-d0$d_cattle+d0$d_dairy
qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)

shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
mkgraph<-function(kind){
  nb<-if(kind=="queen") poly2nb(shp,queen=TRUE,snap=0.01) else knn2nb(knearneigh(coords,k=5))
  k3<-knn2nb(knearneigh(coords,k=3))
  for(i in which(card(nb)==0)) for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
  cc<-n.comp.nb(nb)
  while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m)
    D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf
    ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]]
    nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)}
  nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb);g}
G_queen<-mkgraph("queen");G_knn<-mkgraph("knn")
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)))
pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-read.csv(file.path(CACHE,"ST_EHEC_agecase.csv"))%>%group_by(region,year)%>%summarise(cases=sum(cases),.groups="drop")
Est<-read.csv(file.path(CACHE,"ST_Estd_panel.csv"))
popd<-d0%>%transmute(region,year,pop_off=pmax(e_pop_dens,1e-6))

## ── 사양 생성 (층화 무작위: 각 도메인 변수 균등 출현) ──
grid<-expand.grid(D1=DOM$D1,D2=DOM$D2,D3=DOM$D3,D4=DOM$D4,D5=DOM$D5,stringsAsFactors=FALSE)
idx<-sample(seq_len(nrow(grid)),min(N_SPEC,nrow(grid)))
specs<-grid[idx,]
specs$rum_form<-rep(c("cont","tertile"),length.out=nrow(specs))
specs$transform<-rep(c("log1p","raw"),length.out=nrow(specs))
specs$structure<-rep(c("full","rw1","spatial"),length.out=nrow(specs))
specs$offset   <-rep(c("agestd","agestd","pop"),length.out=nrow(specs))
specs$graph    <-rep(c("queen","queen","queen","knn"),length.out=nrow(specs))
specs$spec_id<-seq_len(nrow(specs))
write.csv(specs,file.path(ROOT,"results_ST","MASTER_multiverse_specs.csv"),row.names=FALSE)
cat("사양",nrow(specs),"개 생성\n")

fit_one<-function(s){
  Z<-if(s$transform=="log1p") Z_log else Z_raw
  vars<-c(s$D1,s$D2,s$D3,s$D4,s$D5)
  dd<-expand.grid(region=regs,year=2011:2024,stringsAsFactors=FALSE)%>%
    left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%
    left_join(Z,by=c("region","year"))%>%left_join(popd,by=c("region","year"))
  dd$cases[is.na(dd$cases)]<-0
  dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
  dd$pop_off[is.na(dd$pop_off)|dd$pop_off<=0]<-min(dd$pop_off[dd$pop_off>0],na.rm=TRUE)
  for(v in vars) dd[[v]][is.na(dd[[v]])]<-0
  dd$rum_T<-d0$rum_T[match(paste(dd$region,dd$year),paste(d0$region,d0$year))]
  dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High")); dd$rum_T[is.na(dd$rum_T)]<-"Low"
  dd$area<-match(dd$region,regs);dd$yr<-dd$year-2010;dd$area_time<-seq_len(nrow(dd))
  useT<-(s$rum_form=="tertile")
  terms<-c(if(useT)"rum_T" else NULL, sprintf("`%s`",if(useT) setdiff(vars,c("d_cattle","d_dairy")) else vars))
  gph<-if(s$graph=="queen") G_queen else G_knn
  f<-paste0("cases~",paste(terms,collapse="+"),"+f(area,model='bym2',graph=gph,scale.model=TRUE,constr=TRUE,hyper=pc)")
  if(s$structure %in% c("full","rw1")) f<-paste0(f,"+f(yr,model='rw1',scale.model=TRUE,constr=TRUE,hyper=pcp)")
  if(s$structure=="full") f<-paste0(f,"+f(area_time,model='iid',hyper=pcp)")
  E<-if(s$offset=="agestd") dd$E_std else dd$pop_off
  m<-try(inla(as.formula(f),family="nbinomial",data=dd,E=E,control.fixed=CF,
        control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50))),silent=TRUE)
  if(inherits(m,"try-error")) return(NULL)
  fx<-m$summary.fixed
  res<-lapply(rownames(fx),function(nm){
    r<-fx[nm,];data.frame(spec_id=s$spec_id,term=nm,IRR=exp(r$mean),lo=exp(r$`0.025quant`),hi=exp(r$`0.975quant`),
      cred=(r$`0.025quant`>0|r$`0.975quant`<0),DIC=m$dic$dic,WAIC=m$waic$waic,stringsAsFactors=FALSE)})
  do.call(rbind,res)
}
out<-list();t0<-Sys.time()
for(i in seq_len(nrow(specs))){
  r<-fit_one(specs[i,]); if(!is.null(r)) out[[length(out)+1]]<-r
  if(i%%20==0) cat(sprintf("  %d/%d (%.1f분)\n",i,nrow(specs),as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
R<-do.call(rbind,out)
write.csv(R,file.path(ROOT,"results_ST","MASTER_multiverse_results.csv"),row.names=FALSE)
S<-R%>%filter(term!="(Intercept)")%>%group_by(term)%>%
  summarise(n_spec=n(),n_cred=sum(cred),cred_pct=round(100*mean(cred),1),
            median_IRR=round(median(IRR),3),pct_above1=round(100*mean(IRR>1),1),.groups="drop")%>%
  arrange(desc(cred_pct))
write.csv(S,file.path(ROOT,"results_ST","MASTER_multiverse_summary.csv"),row.names=FALSE)
cat("\n═══ 변수별 credible 비율 (사전지정 규칙: ≥50% 강건 / 20-50% 시사적 / <20% null) ═══\n")
print(as.data.frame(S),row.names=FALSE)
cat(sprintf("\n총 %d 사양 · %.1f분\n",length(out),as.numeric(difftime(Sys.time(),t0,units="mins"))))
