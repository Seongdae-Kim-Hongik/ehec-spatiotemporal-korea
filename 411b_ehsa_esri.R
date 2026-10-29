# 411b — EHEC EHSA 표준(Esri) 분류 재계산 + 연도별 Moran I (2026-09-14 외부 검토 대응)
#  411의 classify()는 "마지막 해 hot 아님 & 과거 1회라도 hot"을 Historical hot으로 분류(표준 아님) → Esri 정의로 재분류.
#  Esri: https://doc.esri.com/en/arcgis-pro/latest/tool-reference/space-time-pattern-mining/learnmoreemerging.html
#  입력·Gi* 계산은 411과 동일(EB 평활 연령표준화 SIR, include.self queen W, z>1.96). 산출 results_ST/: EHSA_agestd_cat_esri.csv, EHSA_agestd_sum_esri.csv, EHSA_old_vs_esri.csv, agestd_moran_yearly.csv
source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_st_setup.R"))
source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지 (EHSA·Gi*는 결측을 못 다뤄 청주 마스킹 없음). 이전 판 .bak_pre_v2
set.seed(20260914)
yy<-2011:2024
O<-v2_cases_rc("ST_EHEC_agecase.csv")
E<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")]
grid<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(O,by=c("region","year"))%>%left_join(E,by=c("region","year"))
grid<-v2_finish(grid,kdca=FALSE)
M<-matrix(0,NR,length(yy),dimnames=list(regs,yy))
for(ti in seq_along(yy)){py<-grid[grid$year==yy[ti],];idx<-match(regs,py$region)
  cs<-py$cases[idx];Ee<-py$E_std[idx];cs[is.na(cs)]<-0;Ee[is.na(Ee)|Ee<=0]<-min(Ee[Ee>0],na.rm=TRUE)
  eb<-tryCatch(EBest(cs,Ee)$estmm,error=function(e)cs/Ee);eb[!is.finite(eb)]<-0;M[,ti]<-eb}
nbs<-include.self(nb);lws<-nb2listw(nbs,style="W",zero.policy=TRUE)
Z<-apply(M,2,function(x){x[!is.finite(x)]<-0;if(sd(x)==0)return(rep(0,length(x)));as.numeric(localG(x,lws,zero.policy=TRUE))});rownames(Z)<-regs
esri<-function(zv,yrs){n<-length(zv);hot<-zv>1.96;cold<-zv< -1.96;nh<-sum(hot);nc<-sum(cold);thr<-ceiling(0.9*n)
  mk<-suppressWarnings(tryCatch(cor.test(yrs,zv,method="kendall"),error=function(e)NULL))
  tau<-if(is.null(mk))NA else unname(mk$estimate);p<-if(is.null(mk))NA else mk$p.value
  up<-isTRUE(p<0.05&&tau>0);dn<-isTRUE(p<0.05&&tau<0)
  runs<-function(b){r<-rle(b);list(last=if(tail(r$values,1))tail(r$lengths,1) else 0)}
  if(hot[n]){
    if(nh>=thr){if(up)return("Intensifying hot");if(dn)return("Diminishing hot");return("Persistent hot")}
    if(nh==1)return("New hot")
    if(runs(hot)$last>=2&&runs(hot)$last==nh)return("Consecutive hot")
    if(any(cold))return("Oscillating hot")
    return("Sporadic hot")}
  if(cold[n]){
    if(nc>=thr){if(dn)return("Intensifying cold");if(up)return("Diminishing cold");return("Persistent cold")}
    if(nc==1)return("New cold")
    if(runs(cold)$last>=2&&runs(cold)$last==nc)return("Consecutive cold")
    if(any(hot))return("Oscillating cold")
    return("Sporadic cold")}
  if(nh>=thr)return("Historical hot");if(nc>=thr)return("Historical cold")
  "No pattern"}
old<-read.csv(file.path(ROOT,"results_ST/EHSA_agestd_cat.csv"),fileEncoding="UTF-8")
cats<-sapply(1:NR,function(i)esri(Z[i,],as.numeric(yy)))
nhot<-rowSums(Z>1.96);ncold<-rowSums(Z< -1.96)
mkp<-t(sapply(1:NR,function(i){k<-suppressWarnings(tryCatch(cor.test(as.numeric(yy),Z[i,],method="kendall"),error=function(e)NULL));if(is.null(k))c(NA,NA) else c(unname(k$estimate),k$p.value)}))
cat.df<-data.frame(region=regs,EHSA_esri=cats,EHSA_old=old$EHSA[match(regs,old$region)],n_hot_years=nhot,n_cold_years=ncold,giz_final=round(Z[,ncol(Z)],2),mk_tau=round(mkp[,1],3),mk_p=round(mkp[,2],4))
write.csv(cat.df,file.path(ROOT,"results_ST/EHSA_agestd_cat_esri.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(data.frame(cat=names(table(cats)),n=as.integer(table(cats))),file.path(ROOT,"results_ST/EHSA_agestd_sum_esri.csv"),row.names=FALSE)
write.csv(as.data.frame.matrix(table(old=cat.df$EHSA_old,esri=cat.df$EHSA_esri)),file.path(ROOT,"results_ST/EHSA_old_vs_esri.csv"))
cat("== Esri 분류 ==\n");print(table(cats));cat("== old × esri ==\n");print(table(old=cat.df$EHSA_old,esri=cat.df$EHSA_esri))
cat("== hot 연수 분포(1회 이상 hot 구) ==\n");print(table(nhot[nhot>0]))
# 연도별 Moran I (EB SIR, 999 순열) + 추세
my<-do.call(rbind,lapply(seq_along(yy),function(ti){mc<-moran.mc(M[,ti],lw,nsim=999,zero.policy=TRUE);data.frame(year=yy[ti],moran_I=round(unname(mc$statistic),3),p=round(mc$p.value,3))}))
kt<-cor.test(my$year,my$moran_I,method="kendall");my$trend_tau<-round(unname(kt$estimate),3);my$trend_p<-round(kt$p.value,4)
pm<-moran.mc(rowMeans(M),lw,nsim=999,zero.policy=TRUE)
write.csv(my,file.path(ROOT,"results_ST/agestd_moran_yearly.csv"),row.names=FALSE)
cat(sprintf("== 연도별 Moran I: 범위 %.3f-%.3f, 추세 tau %.3f P %.4f; 기간평균 SIR I %.3f P %.3f ==\n",min(my$moran_I),max(my$moran_I),kt$estimate,kt$p.value,pm$statistic,pm$p.value));print(my)
