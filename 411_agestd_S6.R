# ★SUPERSEDED EHSA CLASSIFICATION (2026-09-14): 이 파일의 EHSA 분류는 Esri 비표준. 원고용 EHSA는 411b_ehsa_esri.R + 412c_mapB_esri.R 산출만 사용.
# 411 — EHEC 연령표준화 S6: SIR 기반 EHSA + Moran + 3구간 risk/preventive
#  결과변수=KDCA 신고 사례(O), offset=E_std. SIR=O/E_std, EB평활 → Gi*/Moran/EHSA.
#  구간=COVID정렬 3기(2011-15 / 2016-19 / 2020-24). 각 구간 top risk + top preventive.
#  산출 results_ST/: EHSA_agestd_cat.csv · EHSA_agestd_sum.csv · agestd_SIRpanel.csv
#                    agestd_moran.txt · MASTER_agestd_period_drivers.csv
source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_st_setup.R"))
source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지 (EHSA·Gi*는 결측을 못 다뤄 청주 마스킹 없음). 이전 판 .bak_pre_v2
set.seed(20260710)
MC<-v2_main_cov("장출혈성대장균감염증");cov<-MC$vars   # ★2026-10-21 v2b(이전 판 .bak_pre_v2b): 기간 동인은 옛 13공변량(ST_covariates_panel.rds, 결측 있어 v2에서 멈춤) 대신 원고 최종 20변수(540 DOM)·v2 패널
yy<-2011:2024
O<-v2_cases_rc("ST_EHEC_agecase.csv")
E<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")]
grid<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(O,by=c("region","year"))%>%left_join(E,by=c("region","year"))
grid<-v2_finish(grid,kdca=FALSE)

# ===== EB 평활 SIR 패널 (region×year) =====
M<-matrix(0,NR,length(yy),dimnames=list(regs,yy))
for(ti in seq_along(yy)){py<-grid[grid$year==yy[ti],];idx<-match(regs,py$region)
  cs<-py$cases[idx];Ee<-py$E_std[idx];cs[is.na(cs)]<-0;Ee[is.na(Ee)|Ee<=0]<-min(Ee[Ee>0],na.rm=TRUE)
  eb<-tryCatch(EBest(cs,Ee)$estmm,error=function(e)cs/Ee);eb[!is.finite(eb)]<-0;M[,ti]<-eb}    # 평활 SIR
sirp<-as.data.frame(M);sirp$region<-regs
write.csv(sirp,file.path(ROOT,"results_ST/agestd_SIRpanel.csv"),row.names=FALSE)

# ===== Moran's I (moran.mc 999) — 최종연도·기간평균 SIR =====
sink(file.path(ROOT,"results_ST/agestd_moran.txt"))
for(lbl in c("2024","period_mean")){
  x<-if(lbl=="2024")M[,ncol(M)] else rowMeans(M)
  mc<-moran.mc(x,lw,nsim=999,zero.policy=TRUE)
  cat(sprintf("[SIR %s] Moran I=%.3f  p=%.3f\n",lbl,mc$statistic,mc$p.value))
}
sink()

# ===== EHSA (Gi* 연별 + Mann-Kendall + ESRI 유형) =====
nbs<-include.self(nb);lws<-nb2listw(nbs,style="W",zero.policy=TRUE)
classify<-function(zv,yrs){n<-length(zv);hot<-zv>1.96;cold<-zv< -1.96;fh<-mean(hot);fc<-mean(cold)
  mk<-suppressWarnings(tryCatch(cor.test(yrs,zv,method="kendall"),error=function(e)NULL))
  tau<-if(is.null(mk))NA else unname(mk$estimate);p<-if(is.null(mk))NA else mk$p.value
  up<-(!is.na(p)&&p<0.05&&tau>0);dn<-(!is.na(p)&&p<0.05&&tau<0);last_hot<-hot[n];last_cold<-cold[n]
  if(any(hot)&&any(cold))return("Oscillating")
  if(last_hot){if(!any(hot[1:(n-1)]))return("New hot");if(fh>=0.9&&up)return("Intensifying hot");if(fh>=0.9)return("Persistent hot");if(dn)return("Diminishing hot");return("Sporadic/Consecutive hot")}
  if(!last_hot&&any(hot))return("Historical hot")
  if(last_cold){if(!any(cold[1:(n-1)]))return("New cold");if(fc>=0.9)return("Persistent cold");return("Sporadic cold")}
  if(!last_cold&&any(cold))return("Historical cold");"No pattern"}
Z<-apply(M,2,function(x){x[!is.finite(x)]<-0;if(sd(x)==0)return(rep(0,length(x)));as.numeric(localG(x,lws,zero.policy=TRUE))});rownames(Z)<-regs
pz<-2*pnorm(-abs(Z[,ncol(Z)]));padj<-p.adjust(pz,method="BH")
cats<-sapply(1:NR,function(i)classify(Z[i,],as.numeric(yy)))
cent<-t(sapply(seq_along(yy),function(ti){h<-which(Z[,ti]>1.96);if(length(h)==0)c(NA,NA)else colMeans(coords[h,,drop=FALSE])}))
ok<-which(!is.na(cent[,1]));shift<-if(length(ok)>=2)sqrt(sum((cent[ok[length(ok)],]-cent[ok[1],])^2))/1000 else NA
cat.df<-data.frame(region=regs,EHSA=cats,giz_final=round(Z[,ncol(Z)],2),p_BH=round(padj,3))
write.csv(cat.df,file.path(ROOT,"results_ST/EHSA_agestd_cat.csv"),row.names=FALSE)
sumdf<-as.data.frame(t(table(factor(cats))));sumdf$hotspot_shift_km<-round(shift,1)
write.csv(data.frame(cat=names(table(cats)),n=as.integer(table(cats))),file.path(ROOT,"results_ST/EHSA_agestd_sum.csv"),row.names=FALSE)

# ===== 5년 구간 risk/preventive (교수 Table 11 경계: 2011-15/2016-20/2021-24) =====
PER<-list("2011-2015"=2011:2015,"2016-2020"=2016:2020,"2021-2024"=2021:2024)
lab<-function(cv) ifelse(cv%in%names(MC$lab),MC$lab[cv],ifelse(cv%in%names(KLAB),KLAB[cv],cv))   # ★v2b
allper<-list()
for(pn in names(PER)){ys<-PER[[pn]]
  gg<-grid[grid$year%in%ys,]%>%mutate(.cj=sigungu_cd==V2_CJ_CODE&year%in%V2_CJ_YEARS&v2_cj_on())%>%group_by(region)%>%summarise(cases=if(any(.cj))NA_real_ else sum(cases),E_std=sum(E_std),.groups="drop")   # ★v2b: 모형 적합만 KDCA형 청주 2011-2013 마스킹(그 칸이 든 구간 합은 NA → INLA가 우도에서 뺌; CJ_MASK=0이면 끔). EHSA·Moran은 관측값 그대로
  Zw<-v2_build_Z_main("장출혈성대장균감염증",min(ys),max(ys),MC)%>%group_by(region)%>%summarise(across(all_of(cov),mean),.groups="drop")   # ★v2b: v2 패널, 창 안 표준화
  d<-gg%>%left_join(Zw,by="region");d$area<-match(d$region,regs)
  rhs<-paste(sprintf("`%s`",cov),collapse="+")
  form<-as.formula(paste0("cases~",rhs,"+f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)"))
  m<-inla(form,family="nbinomial",data=d,E=d$E_std,control.fixed=CF,control.compute=list(dic=TRUE),control.inla=list(control.vb=list(emergency=50)))
  fx<-m$summary.fixed[cov,,drop=FALSE]
  r<-data.frame(period=pn,cov=cov,공변량=lab(cov),IRR=round(exp(fx$mean),3),lo=round(exp(fx$`0.025quant`),3),hi=round(exp(fx$`0.975quant`),3),
    credible=(fx$`0.025quant`>0|fx$`0.975quant`<0))
  allper[[pn]]<-r
  # 구간별 최대크기 요인 1개씩(credible=* 표시); 교수 Table11 형식
  risk<-r[r$IRR>1,,drop=FALSE];risk<-risk[which.max(risk$IRR),]
  prev<-r[r$IRR<1,,drop=FALSE];prev<-prev[which.min(prev$IRR),]
  mk<-function(x) if(nrow(x))sprintf("%s %.2f(%.2f-%.2f)%s",x$공변량,x$IRR,x$lo,x$hi,ifelse(isTRUE(x$credible),"*","")) else "—"
  cat(sprintf("[%s] risk↑ %s | preventive↓ %s\n",pn,mk(risk),mk(prev)))}
write.csv(bind_rows(allper),file.path(ROOT,"results_ST/MASTER_agestd_period_drivers.csv"),row.names=FALSE)

cat("\n===== EHSA 유형 분포 =====\n");print(table(cats))
cat(sprintf("hotspot 이동 %.0f km\n",ifelse(is.na(shift),0,shift)))
cat("Moran:\n");cat(readLines(file.path(ROOT,"results_ST/agestd_moran.txt")),sep="\n")
cat("\n완료 → EHSA_agestd_*, agestd_SIRpanel/moran, MASTER_agestd_period_drivers\n")
