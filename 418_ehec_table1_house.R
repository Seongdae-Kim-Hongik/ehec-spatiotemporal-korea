# 418 — EHEC 하우스스타일 Table 1 (문헌기반 7도메인, /333 종합 확정 설계)
#  내러티브(Zotero FBD_04_EHEC공간 73편): 반추동물(소) 장내 리저버(Ward2025·Friesema2010·Mulder2020)
#   → 분뇨매개 농경지/농산물·관개수 오염(Sarnino2025·Berry2015·An2006) → 물/하수/지하수(수인성)
#   → 사람간·식품취급(손씻기·외식) → SES/도농 결핍(Cleary2021) → 소아 HUS 우세(Ylinen·Gould)
#   → 감시체계 검출편의(Graham2024·Gibbons2014).  소는 리저버, 돼지·가금류는 아님 → 음성대조.
#  방법(하우스=HAV/NORO): 단일 결합 BYM2, per-SD 표준화·log1p 변환, 반추동물 tertile 범주형(EID 선례),
#   forced confounder, VIF>10 반복제거(NORO), 8-graph는 후속. offset=E_std(연령표준화).
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071);library(car)})
set.seed(20260717); inla.setOption(num.threads="1:1")
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache")

# ── 라벨(도메인·영문·forced) ──
DOM<-list(   # ★ Enriched A (버전 A base + 논면적·주관건강·성비·독거노인; 재정·주점 credible 유지·DIC 최소)
 "1. Livestock / ruminant reservoir"=c("rum_T","d_pasture"),
 "2. Agriculture & land use"=c("d_field","d_paddy","e_forest"),
 "3. Water & sanitation"=c("a_septic","a_sewer","b_gw_use","b_gw_unfit"),
 "4. Hygiene & food handling"=c("h_fastfood","h_pub","h_selfhealth"),
 "5. Socioeconomic & demographic"=c("g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder"))
# 주: 소아 0-4(g_infant)는 offset(연령표준화 E_std)이 연령을 이미 흡수 → 이중조정이므로 주모형 제외.
#     소아 우세는 §3.7 연령별 발생률(Fig S14)로 서술.
ENG<-c(rum_T="Ruminant density (cattle+dairy), tertile", d_pasture="Pasture area",
 d_field="Upland-field (cropland) area", e_forest="Forest area",
 a_septic="Septic-tank population", a_sewer="Public-sewerage coverage", b_gw_use="Groundwater use", b_gw_unfit="Unfit groundwater",
 h_fastfood="Fast-food / eating-out outlets", h_pub="Bars/pubs",
 d_paddy="Paddy-field area", h_selfhealth="Poor self-rated health", g_solo_elder="Elderly living alone",
 g_sexratio="Sex ratio", g_elderly="Elderly population share", e_urban_pop="Urban population rate", g_fiscal="Fiscal independence", g_basiclivelihood="Basic-livelihood recipients",
 g_infant="Children aged 0-4 years", h_docs="Physicians (per capita)")
FORCED<-c("g_sexratio")   # Enriched A: 인구보정(성비)만 nuisance로 표기; 재정자립은 결정요인이라 forced 아님.
cont_cov<-c("d_pasture","d_field","d_paddy","e_forest","a_septic","a_sewer","b_gw_use","b_gw_unfit","h_fastfood","h_pub","h_selfhealth",
            "g_sexratio","g_fiscal","g_basiclivelihood","g_solo_elder")   # Enriched A (소아·도시화·고령·의사 제외: 헤드라인 흐림/이중조정)

# ── 데이터 ──
d0<-read.csv(file.path(CACHE,"ST_curated46_panel.csv"))
d0<-d0[d0$year>=2011 & d0$year<=2024,]
d0$region<-str_replace_all(as.character(d0$region),"\\s+","");d0$region<-ifelse(d0$region=="인천시미추홀구","인천시남구",d0$region)
# 반추동물 = 소+젖소, tertile(저/중/고)
d0$ruminant<-d0$d_cattle+d0$d_dairy
qt<-quantile(d0$ruminant,c(1/3,2/3),na.rm=TRUE)
d0$rum_T<-cut(d0$ruminant,breaks=unique(c(-Inf,qt,Inf)),labels=c("Low","Mid","High"),include.lowest=TRUE)
cat(sprintf("반추동물 tertile 경계: [%.1f, %.1f]; n(Low/Mid/High)=%s\n",qt[1],qt[2],paste(table(d0$rum_T),collapse="/")))

skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
# 연속형 → (skew>1이면 log1p) → per-SD z
Z<-d0[,c("region","year")]
## [502 수정] SES 탐색 모형(135행)이 cont_cov 밖 변수를 참조 → Z에 함께 표준화(주모형은 cont_cov만 사용, 영향 없음)
EXTRA_Z<-intersect(c("g_elderly","e_urban_pop","h_docs"),names(d0))
for(v in unique(c(cont_cov,EXTRA_Z))){x<-d0[[v]]; x[is.na(x)]<-mean(x,na.rm=TRUE)
  xv<-if(min(x,na.rm=TRUE)>=0 && skewf(x)>1) log1p(pmax(x,0)) else x
  z<-scale(xv)[,1]; z[!is.finite(z)]<-0; Z[[v]]<-z}
Z$rum_T<-factor(d0$rum_T,levels=c("Low","Mid","High"))

# ── VIF 사전 스크린 (>10 반복제거; forced·rum_T 보호) ──
vifdat<-Z[,c(cont_cov)]; vifdat$rum_num<-as.integer(Z$rum_T)  # VIF용 tertile은 정수 근사
vifdat$.y<-rnorm(nrow(vifdat))
keep<-c("rum_num",cont_cov); repeat{
  lmf<-lm(as.formula(paste(".y~",paste(keep,collapse="+"))),data=vifdat)
  v<-tryCatch(car::vif(lmf),error=function(e)NULL); if(is.null(v))break
  v<-v[setdiff(names(v),c(FORCED,"rum_num"))]  # forced·reservoir 제외하고 후보 중 최대
  if(length(v)==0||max(v)<10)break
  drop<-names(which.max(v)); cat(sprintf("VIF>10 제거: %s (VIF=%.1f)\n",drop,max(v))); keep<-setdiff(keep,drop)}
final_cont<-setdiff(keep,"rum_num")
cat("최종 연속형:",length(final_cont),"개 + rum_T tertile\n")

# ── 그래프 (410/414 동일) ──
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
  dd$cases[is.na(dd$cases)]<-0; dd$E_std[is.na(dd$E_std)|dd$E_std<=0]<-min(dd$E_std[dd$E_std>0],na.rm=TRUE)
  dd$area<-match(dd$region,regs);dd$time<-match(dd$year,yy);dd$id<-1:nrow(dd)
  dd$rum_T<-factor(dd$rum_T,levels=c("Low","Mid","High"))
  terms<-c(if(useT)"rum_T",sprintf("`%s`",cont))
  form<-as.formula(paste0("cases~",paste(terms,collapse="+"),"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(id,model='iid',hyper=pcp)"))
  inla(form,family="nbinomial",data=dd,E=dd$E_std,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}

m<-fitmod(final_cont,useT=TRUE)
fx<-m$summary.fixed
# 결과표 조립 (도메인 순서)
row1<-function(cv){
  nm<-if(cv=="rum_T:Mid")"rum_TMid" else if(cv=="rum_T:High")"rum_THigh" else cv
  if(!nm%in%rownames(fx))return(NULL)
  r<-fx[nm,]; data.frame(cov=cv,IRR=exp(r$mean),lo=exp(r$`0.025quant`),hi=exp(r$`0.975quant`),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}
out<-list()
for(dm in names(DOM)) for(cv in DOM[[dm]]){
  if(cv=="rum_T"){for(lv in c("rum_T:Mid","rum_T:High")){rr<-row1(lv);if(!is.null(rr)){rr$domain<-dm;rr$label<-paste0(ENG["rum_T"],ifelse(lv=="rum_T:Mid"," — Mid vs Low"," — High vs Low"));rr$forced<-FALSE;out[[length(out)+1]]<-rr}}}
  else if(cv%in%final_cont){rr<-row1(cv);if(!is.null(rr)){rr$domain<-dm;rr$label<-ENG[cv];rr$forced<-cv%in%FORCED;out[[length(out)+1]]<-rr}}
}
res<-do.call(rbind,out)
res$IRR<-round(res$IRR,3);res$lo<-round(res$lo,3);res$hi<-round(res$hi,3)
write.csv(res[,c("domain","cov","label","IRR","lo","hi","cred","forced")],file.path(ROOT,"results_ST","MASTER_ehec_table1_house.csv"),row.names=FALSE)
# ---- 지도용 사후 지역 RR + 초과확률 (17변수 주모형 잔차 공간효과) ----
uu<-m$summary.random$area[1:length(regs),]
rr<-data.frame(region=regs,RR=round(exp(uu$mean),3),RR_lo=round(exp(uu$`0.025quant`),3),RR_hi=round(exp(uu$`0.975quant`),3))
rr$hi_risk<-as.integer(rr$RR_lo>1)
write.csv(rr,file.path(ROOT,"results_ST","MASTER_agestd_regionSIR.csv"),row.names=FALSE)
pep<-sapply(m$marginals.random$area[1:length(regs)],function(mg) 1-INLA::inla.pmarginal(0,mg))
write.csv(data.frame(region=regs,pep=round(pep,4)),file.path(ROOT,"results_ST","MASTER_agestd_exceedance.csv"),row.names=FALSE)
cat(sprintf("지도용 갱신: regionSIR(고위험 %d개)·exceedance from 17변수 주모형\n",sum(rr$hi_risk)))
# ---- 적합통계(본문 §3.2용): phi/DIC/WAIC ----
phi_m<-m$summary.hyperpar["Phi for area","mean"]
writeLines(c(sprintf("phi for area mean=%.3f",phi_m),sprintf("DIC=%.1f",m$dic$dic),sprintf("WAIC=%.1f",m$waic$waic),
             sprintf("model=house 17-var (2011-2024, EPV=%.0f)",sum(case$cases[case$year>=2011])/(length(final_cont)+2))),
           file.path(ROOT,"results_ST","agestd_fit.txt"))
# (주: §3.2 Moran은 연령표준화 SIR O/E_std 기반=공변량모형 무관, 411에서 산출. 전모형 iid항 때문에
#  full-model 잔차-Moran은 well-defined 지표가 아니라 보고하지 않음.)
cat(sprintf("\n===== EHEC Table 1 결합 BYM2 (DIC=%.1f WAIC=%.1f, EPV=%.0f) =====\n",m$dic$dic,m$waic$waic,sum(case$cases[case$year>=2011])/(length(final_cont)+2)))
cur<-"";for(i in 1:nrow(res)){if(res$domain[i]!=cur){cur<-res$domain[i];cat("\n[",cur,"]\n")};cat(sprintf("  %-42s IRR %.3f (%.3f-%.3f)%s%s\n",res$label[i],res$IRR[i],res$lo[i],res$hi[i],ifelse(res$cred[i]," *",""),ifelse(res$forced[i]," ‡","")))}

# ── 음성대조: 돼지·가금류 개별 추가 (Friesema: 소는 리저버, 돼지·가금류 아님) ──
Z$d_swine<-scale(log1p(pmax(d0$d_swine,0)))[,1];Z$d_poultry<-scale(log1p(pmax(d0$d_poultry,0)))[,1]
Z$d_swine[!is.finite(Z$d_swine)]<-0;Z$d_poultry[!is.finite(Z$d_poultry)]<-0
nc<-data.frame()
for(v in c("d_swine","d_poultry")){mm<-fitmod(c(final_cont,v),useT=TRUE);r<-mm$summary.fixed[v,]
  nc<-rbind(nc,data.frame(nc_var=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0)))}
write.csv(nc,file.path(ROOT,"results_ST","MASTER_ehec_table1_negcontrol.csv"),row.names=FALSE)
cat("\n[음성대조 (개별추가)]\n");for(i in 1:nrow(nc))cat(sprintf("  %-10s IRR %.3f (%.3f-%.3f) cred=%s\n",nc$nc_var[i],nc$IRR[i],nc$lo[i],nc$hi[i],nc$cred[i]))

# ============ /333 종합 민감도 (반추동물 null 정당화·강건성) ============
Z$rum_cont<-scale(log1p(pmax(d0$ruminant,0)))[,1]; Z$rum_cont[!is.finite(Z$rum_cont)]<-0   # 연속 per-SD
Z$rum_ord<-as.integer(factor(d0$rum_T,levels=c("Low","Mid","High")))                        # 순서형 추세(1/2/3)
grab<-function(m,nm){r<-m$summary.fixed[nm,];data.frame(IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}
sens<-data.frame()
# [515 수정] 개별 민감도 모형이 수렴 실패해도 나머지를 저장하도록 보호(Z8: 산출물 보존 우선)
safefit<-function(...) tryCatch(fitmod(...),error=function(e){cat("  [건너뜀] 수렴 실패:",conditionMessage(e),"\n");NULL})
addrow<-function(m,nm,mlab,vlab){ if(is.null(m))return(invisible(NULL))
  r<-tryCatch(grab(m,nm),error=function(e)NULL); if(is.null(r))return(invisible(NULL))
  sens<<-rbind(sens,cbind(model=mlab,var=vlab,r)) }
# S1 연속 per-SD 반추동물 (tertile 대체)
m_c<-safefit(c(final_cont,"rum_cont"),useT=FALSE); addrow(m_c,"rum_cont","S1 ruminant continuous (per-SD, log1p)","ruminant")
# 순서형 추세(저→중→고 선형)
m_o<-safefit(c(final_cont,"rum_ord"),useT=FALSE); addrow(m_o,"rum_ord","S1 ruminant ordinal trend (per level)","ruminant")
# S3 매개 과보정 분리: 리저버 tertile + 최소교란(forced)만 = 총효과
m_p<-safefit(c("g_sexratio","g_elderly","e_urban_pop","g_fiscal","h_docs"),useT=TRUE)
for(lv in c("rum_TMid","rum_THigh")) addrow(m_p,lv,"S3 parsimonious (reservoir+forced only)",lv)
# (소아 0-4는 주모형에서 제외됐으므로 children-excluded 민감도는 불필요)
# S4 대안 인접행렬(knn k=5) 강건성 → 주요계수
nbk<-knn2nb(knearneigh(coords,k=5)); nbk<-make.sym.nb(nbk); g<-tempfile(fileext=".graph"); nb2INLA(g,nbk)  # g 재지정
m_k<-safefit(final_cont,useT=TRUE)
for(v in c("rum_THigh","h_pub","g_fiscal")) addrow(m_k,v,"S4 alt adjacency (kNN k=5)",v)
write.csv(sens,file.path(ROOT,"results_ST","MASTER_ehec_table1_sens.csv"),row.names=FALSE)
cat("\n===== /333 민감도 =====\n"); for(i in 1:nrow(sens))cat(sprintf("  [%-38s] %-10s IRR %.3f (%.3f-%.3f)%s\n",sens$model[i],sens$var[i],sens$IRR[i],sens$lo[i],sens$hi[i],ifelse(sens$cred[i]," *","")))

# GVIF (실제 더미 design matrix; gpt 지적)
vf<-Z[,c(final_cont)]; vf$rum_T<-factor(d0$rum_T,levels=c("Low","Mid","High")); vf$.y<-rnorm(nrow(vf))
gv<-car::vif(lm(as.formula(paste(".y~rum_T+",paste(final_cont,collapse="+"))),data=vf))
cat("\n[GVIF 최대 (더미 design)]:",round(max(gv[,ncol(gv)]^2,na.rm=TRUE),2)," (rum_T GVIF^(1/2df)=",round(gv["rum_T",3],2),")\n")
cat("\n저장: house.csv + negcontrol.csv + sens.csv\n")
