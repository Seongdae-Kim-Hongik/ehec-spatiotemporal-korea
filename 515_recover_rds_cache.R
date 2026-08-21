# 515 — rds 공변량 캐시 복구. 500-505의 복구는 CSV 패널에만 반영되어 rds는 오염 상태로 남아 있었다.
#  410/411/416/417/430은 CSV가 아니라 이 rds를 읽으므로, 부록 민감도 분석 전체가 오염 결과였다.
CACHE<-"data_cache"
src<-file.path(CACHE,"ST_covariates_panel.rds")
bak<-file.path(CACHE,"ST_covariates_panel_CONTAMINATED_BACKUP.rds")
if(!file.exists(bak)){ file.copy(src,bak); cat("백업 생성:",basename(bak),"\n") } else cat("백업 이미 존재\n")
CV<-readRDS(src); d<-read.csv(file.path(CACHE,"ST_curated46_panel.csv"))
nsp<-function(x) gsub("[[:space:]]+","",as.character(x))
stopifnot(identical(paste(nsp(CV$region),CV$year), paste(nsp(d$region),d$year)))   # 정렬 동일 확인
V<-colnames(CV$raw); stopifnot(all(V %in% names(d)))
old<-CV$raw
new<-as.matrix(d[,V]); rownames(new)<-rownames(old); colnames(new)<-V
# Ms = scale(raw) (원 정의와 동일)
Ms<-apply(new,2,function(x){ s<-sd(x,na.rm=TRUE); if(!is.finite(s)||s==0) rep(0,length(x)) else as.numeric(scale(x)) })
colnames(Ms)<-V
cover<-apply(new,2,function(x) mean(is.finite(x)&!is.na(x)))
CV$raw<-new; CV$Ms<-Ms; CV$cover<-cover
saveRDS(CV,src)
cat("\n── 변경된 변수 (중앙값 기준) ──\n")
ch<-0
for(v in V){ o<-median(old[,v],na.rm=TRUE); n<-median(new[,v],na.rm=TRUE)
  if(!isTRUE(all.equal(o,n,tolerance=1e-8))){ ch<-ch+1
    cat(sprintf("  %-20s %12.4f → %12.4f\n",v,o,n)) } }
cat(sprintf("\n총 %d/%d 변수 갱신 · rds 복구 완료\n",ch,length(V)))
# Z1 게이트 재검증
cat("\n── Z1 재검증 (연도별 0비율>90%) ──\n"); bad<-0
for(v in V){ zr<-tapply(new[,v],CV$year,function(z) mean(z==0,na.rm=TRUE))
  if(any(zr>0.9,na.rm=TRUE)){ bad<-bad+1; cat(sprintf("  ✗ %s (%d년)\n",v,sum(zr>0.9,na.rm=TRUE))) } }
cat(if(bad==0) "  통과 — 오염 변수 없음\n" else sprintf("  → %d개 잔존\n",bad))
