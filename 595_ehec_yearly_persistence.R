# 595_ehec_yearly_persistence.R — EHEC 지속성(595_shig 복제, 2026-09-10): F45/F45b 연도별 SIR 지도, F46 지속성 지도, T34/T34b/T34c.
# 595_shig_yearly_persistence.R — [EHEC] B+F(2026-09-04 교수 피드백): 연도별 SIR 소지도(공통 계급) + 클러스터 지속성.
#  B: F45 = 2011-2024 14패널 EB평활 연령표준화 SIR, 모든 해 동일 E_std 기준·공통 로그대칭 6계급(<0.5 … >=2), 2013 인천 표기.
#  F: T34 지속성 = 구별 (i) 원 SIR>1 연수, (ii) EB SIR 상위 5분위 연수, (iii) 연도간 순위 일치(Kendall W, EB SIR 순위);
#     F46 = 상위분위 연수 지도(+ 7년 이상 외곽선). T34b = 시도별 연도 발생률(/10만, PP 인구). 577 헬퍼 재사용. INLA 미실행.
#  EHEC_CASES=domestic 이면 국내발생 기준(setup에서 cases·E_std·SH·CAPF 전환).
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(ggplot2);library(ggspatial);library(jsonlite);library(tidyr)})
GD<-file.path(Sys.getenv("EHEC_ROOT", unset="."),"G.Downloads")
# ── EHEC 최소 셋업(542 헤더 재사용: shp·regs·g·mk) ──
# 542 — 최종 주모형(20변수) 기준 기간 분석: ① 4년창 순수공간 ② 코로나 전후 분할 + 상호작용 ③ 주요 3변수 × 시기 상호작용
suppressPackageStartupMessages({library(dplyr);library(sf);library(spdep);library(stringr);library(INLA);library(e1071)})
set.seed(20260822); inla.setOption(num.threads=Sys.getenv("INLA_THREADS","2:1"))
BASE<-file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP")
ROOT<-Sys.getenv("EHEC_ROOT", unset=".");CACHE<-file.path(ROOT,"data_cache");RS<-Sys.getenv("EHEC_RS",unset=file.path(ROOT,"results_ST"))
V2_CACHE<-CACHE;source(file.path(Sys.getenv("EHEC_ROOT", unset="."),"_v2_common.R"))   # ★2026-10-21 v2: E·사례 코드 결합·게이트, 채움 금지 (기술 분석이라 청주 마스킹 없음 — 관측 0 사용). 이전 판 .bak_pre_v2
T2<-read.csv(file.path(RS,"MASTER_final_table2.csv"));cont<-T2$cov;ENG<-setNames(T2$label,T2$cov)
d0<-read.csv(file.path(CACHE,Sys.getenv("EHEC_PANEL",unset="ST_prereg_panel.csv")),fileEncoding="UTF-8-BOM");d0$region<-str_replace_all(d0$region,"\\s+","")
skewf<-function(x){x<-x[is.finite(x)];if(length(x)<3||sd(x)==0)return(0);e1071::skewness(x)}
Z<-d0[,c("region","year")];for(v in cont){x<-d0[[v]];x<-v2_nofill(x,v);xv<-if(min(x,na.rm=TRUE)>=0&&skewf(x)>1)log1p(pmax(x,0))else x;z<-scale(xv)[,1];z[!is.finite(z)]<-0;Z[[v]]<-z}
shp<-st_read(file.path(BASE,"final.shp"),quiet=TRUE)%>%mutate(region=str_replace_all(as.character(region),"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))%>%st_set_crs(5179)
coords<-suppressWarnings(st_coordinates(st_centroid(st_geometry(shp))));regs<-shp$region
nb<-poly2nb(shp,queen=TRUE,snap=0.01);k3<-knn2nb(knearneigh(coords,k=3));for(i in which(card(nb)==0))for(j in k3[[i]]){nb[[i]]<-as.integer(unique(c(setdiff(nb[[i]],0L),j)));nb[[j]]<-as.integer(sort(unique(c(setdiff(nb[[j]],0L),i))))}
cc<-n.comp.nb(nb);while(cc$nc>1){m<-as.integer(names(which.max(table(cc$comp.id))));mn<-which(cc$comp.id==m);ot<-which(cc$comp.id!=m);D<-as.matrix(dist(coords));D[ot,ot]<-Inf;D[mn,mn]<-Inf;ix<-which(D[ot,mn,drop=FALSE]==min(D[ot,mn]),arr.ind=TRUE)[1,];a<-ot[ix[1]];b<-mn[ix[2]];nb[[a]]<-as.integer(sort(unique(c(setdiff(nb[[a]],0L),b))));nb[[b]]<-as.integer(sort(unique(c(setdiff(nb[[b]],0L),a))));cc<-n.comp.nb(nb)};nb<-make.sym.nb(nb);g<-tempfile(fileext=".graph");nb2INLA(g,nb)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
case<-v2_cases_rc("ST_EHEC_agecase.csv");Est<-v2_E_rc("ST_Estd_panel.csv")[,c("region","year","E_std")]
mk<-function(yy){DD<-expand.grid(region=regs,year=yy,stringsAsFactors=FALSE)%>%left_join(case,by=c("region","year"))%>%left_join(Est,by=c("region","year"))%>%left_join(v2_relabel(Z,"공변량",keep_code=FALSE),by=c("region","year"))
  DD<-v2_finish(DD,kdca=FALSE);DD$area<-match(DD$region,regs);DD$time<-match(DD$year,yy);DD$id<-1:nrow(DD);DD}
grab<-function(m,v)data.frame(cov=v,label=ENG[v],IRR=round(exp(m$summary.fixed[v,"mean"]),3),lo=round(exp(m$summary.fixed[v,"0.025quant"]),3),hi=round(exp(m$summary.fixed[v,"0.975quant"]),3),cred=(m$summary.fixed[v,"0.025quant"]>0|m$summary.fixed[v,"0.975quant"]<0))
KEY<-c("beef_cattle_farms","fiscal_independence","net_migration","current_smoking");if(nzchar(Sys.getenv("EHEC_HM")))KEY<-setdiff(cont,c("kchs_pc1","urban_population_share","physicians_per_1000","sex_ratio"))   # 2026-09-30 문헌 풀: 최종 세트 가설 변수+재정자주도
DD<-mk(2011:2024);ROOT<-Sys.getenv("EHEC_ROOT", unset=".")
SH<-Sys.getenv("EHEC_SH",unset=file.path(ROOT,"EHEC","series"));dir.create(SH,showWarnings=FALSE,recursive=TRUE);CAPF<-Sys.getenv("EHEC_CAPF",unset=file.path(CACHE,"EHEC_SERIES_CAPTIONS.json"));DOMESTIC<-FALSE   # ★v2b: 시험 실행 때 공용 캡션 파일을 건드리지 않도록 캡션 경로 인자화(EHEC_CAPF)
suppressPackageStartupMessages({library(ggplot2);library(ggspatial);library(jsonlite);library(tidyr)});source(file.path(GD,"EHP_Spatio/fig_style.R"));AGG<-ragg::agg_png

PP<-readRDS(file.path(CACHE,"ST_panel_15x229x14.rds"))|>as.data.frame()|>filter(disease=="장출혈성대장균감염증")|>
  mutate(region=str_replace_all(region,"\\s+",""),region=if_else(region=="인천시미추홀구","인천시남구",region))|>select(region,year,pop)
DD<-left_join(DD,PP,by=c("region","year"));stopifnot(!any(is.na(DD$pop)))
sido_df<-suppressWarnings(shp%>%mutate(sido=substr(as.character(sigungu_cd),1,2))%>%group_by(sido)%>%summarise(.groups="drop")%>%st_set_crs(5179))
DIVL<-c("<0.5","0.5-0.67","0.67-1","1-1.5","1.5-2",">=2");bin_div_ascii<-function(x)cut(x,breaks=c(-Inf,0.5,0.67,1,1.5,2,Inf),labels=DIVL)
upd_caps<-function(new){f<-CAPF;cur<-if(file.exists(f))tryCatch(fromJSON(f,simplifyVector=FALSE),error=function(e)list()) else list();for(k in names(new))cur[[k]]<-new[[k]]
  writeLines(enc2utf8(as.character(toJSON(cur,auto_unbox=TRUE,pretty=TRUE))),f,useBytes=TRUE)}
savepng<-function(p,f,w=8.4,h=7.8)ggsave(file.path(SH,f),p,width=w,height=h,dpi=300,bg="white",device=AGG)
YRS<-2011:2024;LAB<-if(DOMESTIC)"domestic" else "notified"
# ---- 연도별 EB 평활 SIR (Marshall global, 연도 내) ----
Y<-lapply(YRS,function(y){s<-DD%>%filter(year==y)%>%arrange(match(region,regs));eb<-tryCatch(EBest(s$cases,s$E_std)$estmm,error=function(e)s$cases/s$E_std);eb[!is.finite(eb)]<-0
  data.frame(region=s$region,year=y,cases=s$cases,E=s$E_std,rawSIR=s$cases/s$E_std,ebSIR=eb)})%>%bind_rows()
Y$cls<-bin_div_ascii(Y$ebSIR);Y$panel<-factor(as.character(Y$year),levels=as.character(YRS))
spY<-left_join(shp[,c("region","geometry")],Y,by="region")
g45<-ggplot(spY)+geom_sf(aes(fill=cls),color="grey70",linewidth=0.06)+geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.18,inherit.aes=FALSE)+facet_wrap(~panel,nrow=3)+
  scale_fill_manual(values=setNames(PAL_DIV6,DIVL),breaks=DIVL,limits=DIVL,na.value="grey90",name="EB-smoothed\nage-standardized SIR",drop=FALSE)+
  theme_ehp_map(10)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),strip.text=element_text(face="bold",size=9),legend.position="right")
savepng(g45,"EHEC_F45_map_yearly_SIR.png",w=13,h=8.6)
# F45b: 연도 내 상대 SIR (그 해 전국 SIR = 1로 정규화) — 전국 수준 변화를 제거한 순수 지리 패턴
nat<-Y%>%group_by(year)%>%summarise(natSIR=sum(cases)/sum(E),.groups="drop");Yr<-Y%>%left_join(nat,by="year")%>%mutate(relSIR=ebSIR/natSIR,cls=bin_div_ascii(relSIR))
spYr<-left_join(shp[,c("region","geometry")],Yr,by="region")
g45b<-ggplot(spYr)+geom_sf(aes(fill=cls),color="grey70",linewidth=0.06)+geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.18,inherit.aes=FALSE)+facet_wrap(~panel,nrow=3)+
  scale_fill_manual(values=setNames(PAL_DIV6,DIVL),breaks=DIVL,limits=DIVL,na.value="grey90",name="EB-smoothed SIR\nrelative to national\nSIR of the same year",drop=FALSE)+
  theme_ehp_map(10)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),strip.text=element_text(face="bold",size=9),legend.position="right")
savepng(g45b,"EHEC_F45b_map_yearly_relSIR.png",w=13,h=8.6);write.csv(Yr,file.path(SH,"EHEC_F45_data_yearly_SIR.csv"),row.names=FALSE,fileEncoding="UTF-8")
# ---- 지속성 ----
# 지속성 지표: 희소자료(영셀 77-86%)에서 EB 평활 순위는 소규모 군의 0건 셀이 상향되는 인공물이 있어, 순위는 원 SIR(0건=최하위 동률), 상위분위는 '해당 연도에 1건 이상 & 원 SIR 상위 20%'로 정의.
# 상위분위 = 그 해 1건 이상 발생한 구들 중 원 SIR 상위 20% (전체 229구 기준 80분위는 영셀>80%인 해에 0으로 퇴화 → /333 지적 반영)
W<-Y%>%group_by(year)%>%mutate(rk=rank(-rawSIR,ties.method="average"),q80=if(any(cases>0))quantile(rawSIR[cases>0],0.8) else Inf,top5=cases>0&rawSIR>=q80,gt1=rawSIR>1&cases>0)%>%ungroup()
Pst<-W%>%group_by(region)%>%summarise(cases=sum(cases),E=sum(E),years_SIR_gt1=sum(gt1),years_top_quintile=sum(top5),mean_rank=mean(rk),.groups="drop")%>%mutate(SIR_2011_2024=cases/E)%>%arrange(desc(years_top_quintile),mean_rank)
# Kendall W (연도=평가자, 지역=대상), 2013 제외 민감도
kendallW_raw<-function(M){m<-nrow(M);n<-ncol(M);R<-colSums(M);12*sum((R-mean(R))^2)/(m^2*(n^3-n))}
# 동률 보정 W(2026-09-14 검토 반영): 분모에서 m*sum(t^3-t) 차감 — 영셀 동률이 많아 표준 보정 필요
kendallW<-function(M){m<-nrow(M);n<-ncol(M);R<-colSums(M);Tt<-sum(apply(M,1,function(r){t<-table(r);sum(t^3-t)}));12*sum((R-mean(R))^2)/(m^2*(n^3-n)-m*Tt)}
RM<-W%>%select(region,year,rk)%>%pivot_wider(names_from=region,values_from=rk)%>%arrange(year);M<-as.matrix(RM[,-1])
SETS<-list("2011-2024"=RM$year>=2011,"2011-2024 excl. 2020"=RM$year!=2020,"2011-2019"=RM$year<=2019,"2020-2024"=RM$year>=2020)
B<-as.integer(Sys.getenv("PERSIST_B",unset="100000"));set.seed(as.integer(Sys.getenv("PERSIST_SEED",unset="20260904")))   # ★2026-10-21 v2b(이전 판 .bak_pre_v2b): 영모형 반복 수 인자화, 기본 100,000(원고 S23·Fig 3C, 33_persistence_jkms_100k.py와 같은 수). 옛 판은 2,000 고정
# 영모형 2종: (i) 순열 = 연도 내 구 라벨 셔플(구 교환가능성; 동률 구조 보존), (ii) 파라메트릭 = cases_iy ~ Poisson(E_std_iy × natSIR_y) (인구규모·연도 수준 보존, 공간 이질성 없음)
#  ★v2b 계산 방식(분포·검정 통계량은 옛 판과 같음): (i) 순열해도 연도별 동률 구조가 그대로라 동률 보정항은 고정, 열 합만 다시 계산. (ii) Poisson 표본 한 번으로 4개 기간 집합의 W와 구별 'SIR>1 연수'를 함께 계산(구별 결과는 B×229 행렬 대신 누적).
natS<-Y%>%group_by(year)%>%summarise(natSIR=sum(cases)/sum(E),.groups="drop");Ey<-Y%>%select(region,year,E)%>%left_join(natS,by="year")
Em<-matrix(Ey$E*Ey$natSIR,nrow=length(YRS),byrow=TRUE)   # rows=years (Y is ordered region-within-year)
EmatR<-matrix(Ey$E,nrow=length(YRS),byrow=TRUE)   # 모의 환자수를 관측과 같은 SIR(cases/E_std)로 바꿔 순위화하기 위한 E (2026-09-14 버그 수정: 이전엔 환자수 자체를 순위화)
Emat<-EmatR   # rows=years, cols=regions(regs order)
tie_t<-function(v){l<-rle(sort(v))$lengths;sum(l^3-l)}
nR<-ncol(M);t0p<-Sys.time()
Wperm<-lapply(names(SETS),function(nm){Mi<-M[SETS[[nm]],,drop=FALSE];m<-nrow(Mi);den<-m^2*(nR^3-nR)-m*sum(apply(Mi,1,tie_t));out<-numeric(B)
  for(bb in seq_len(B)){R<-numeric(nR);for(i in seq_len(m))R<-R+Mi[i,sample.int(nR)];out[bb]<-12*sum((R-mean(R))^2)/den};out});names(Wperm)<-names(SETS)
cat(sprintf("[595] 순열 영모형 B=%d 완료 %.1f분\n",B,as.numeric(difftime(Sys.time(),t0p,units="mins"))))
obsG<-Pst$years_SIR_gt1[match(regs,Pst$region)];geG<-numeric(nR);totG<-numeric(nR);Wq<-matrix(NA_real_,length(SETS),B,dimnames=list(names(SETS),NULL))
selm<-sapply(SETS,identity);nY<-length(YRS);t0q<-Sys.time()
for(bb in seq_len(B)){cs<-matrix(rpois(length(Em),Em),nrow=nY)
  cnt<-colSums(cs>0&cs>Emat);geG<-geG+(cnt>=obsG);totG<-totG+cnt
  sr<-cs/EmatR;Rk<-matrix(0,nY,nR);Tr<-numeric(nY);for(i in seq_len(nY)){v<- -sr[i,];Rk[i,]<-rank(v,ties.method="average");Tr[i]<-tie_t(v)}
  for(k in seq_along(SETS)){sel<-selm[,k];m<-sum(sel);R<-colSums(Rk[sel,,drop=FALSE]);Wq[k,bb]<-12*sum((R-mean(R))^2)/(m^2*(nR^3-nR)-m*sum(Tr[sel]))}}
cat(sprintf("[595] Poisson 영모형 B=%d 완료 %.1f분\n",B,as.numeric(difftime(Sys.time(),t0q,units="mins"))))
KW<-bind_rows(lapply(names(SETS),function(nm){Mi<-M[SETS[[nm]],,drop=FALSE];Wv<-kendallW(Mi);Wp<-Wperm[[nm]];wq<-Wq[nm,]
  data.frame(set=nm,n_years=nrow(Mi),kendall_W=round(Wv,3),kendall_W_uncorrected=round(kendallW_raw(Mi),3),perm_p=round((sum(Wp>=Wv)+1)/(B+1),6),perm_null_97.5=round(quantile(Wp,0.975),3),perm_z=round((Wv-mean(Wp))/sd(Wp),1),
             poisson_p=round((sum(wq>=Wv)+1)/(B+1),6),poisson_null_97.5=round(quantile(wq,0.975),3),poisson_z=round((Wv-mean(wq))/sd(wq),1),B=B)}))
Wall<-KW$kendall_W[1];Wno13<-KW$kendall_W[2];Wpre<-KW$kendall_W[3];Wpost<-KW$kendall_W[4]
# 구별 지속성의 인구규모 보정: 영모형 (ii)에서 '1건 이상 & SIR>1'인 연수의 분포 → 관측 연수의 상단 P, BH-FDR
Pst$exp_years_SIR_gt1<-round(totG/B,2)[match(Pst$region,regs)]
Pst$p_years_SIR_gt1<-((geG+1)/(B+1))[match(Pst$region,regs)]
Pst$null_B<-B
Pst$fdr_years_SIR_gt1<-p.adjust(Pst$p_years_SIR_gt1,"BH")
Pst<-Pst%>%arrange(p_years_SIR_gt1,desc(years_SIR_gt1))
write.csv(Pst,file.path(SH,"EHEC_T34_persistence.csv"),row.names=FALSE,fileEncoding="UTF-8");write.csv(KW,file.path(SH,"EHEC_T34c_kendallW.csv"),row.names=FALSE)
# 시도별 연도 발생률
SD<-DD%>%mutate(sido=str_extract(region,"^(서울시|부산시|대구시|인천시|광주시|대전시|울산시|세종시|경기도|강원도|충청북도|충청남도|전라북도|전라남도|경상북도|경상남도|제주도)"))%>%
  group_by(sido,year)%>%summarise(cases=sum(cases),pop=sum(pop),.groups="drop")%>%mutate(rate=round(cases/pop*1e5,3))
T34b<-SD%>%select(sido,year,rate)%>%pivot_wider(names_from=year,values_from=rate)%>%left_join(SD%>%group_by(sido)%>%summarise(total_cases=sum(cases),rate_2011_2024=round(sum(cases)/sum(pop)*1e5,3),.groups="drop"),by="sido")%>%arrange(desc(rate_2011_2024))
write.csv(T34b,file.path(SH,"EHEC_T34b_sido_yearly_rate.csv"),row.names=FALSE,fileEncoding="UTF-8")
# F46 지속성 지도: '1건 이상 & SIR>1' 연수(관측) + Poisson 영모형 FDR<0.05 외곽선
Pst$cls<-cut(Pst$years_SIR_gt1,breaks=c(-1,0,2,4,6,9,14),labels=c("0","1-2","3-4","5-6","7-9","10-14"))
spP<-left_join(shp[,c("region","geometry")],Pst,by="region")
g46<-ggplot(spP)+geom_sf(aes(fill=cls),color="grey62",linewidth=0.14)+lvl_layer(spP,levels(Pst$cls))+geom_sf(data=sido_df,fill=NA,color="grey40",linewidth=0.32,inherit.aes=FALSE)+
  scale_fill_manual(values=setNames(c("white",PAL_SEQ5),levels(Pst$cls)),breaks=levels(Pst$cls),limits=levels(Pst$cls),na.value="grey90",name="Years with >= 1 case\nand SIR > 1 (of 14)",drop=FALSE)+
  north_ehp("tr")+scale_ehp("bl")+theme_ehp_map(11)+guide_fill_box(rev=TRUE)+theme(text=element_text(family="Helvetica"),legend.position="right")
d7<-spP[spP$fdr_years_SIR_gt1<0.05,];if(nrow(d7))g46<-g46+geom_sf(data=d7,aes(color="More years than expected (FDR < 0.05)"),fill=NA,linewidth=0.5,inherit.aes=FALSE)+scale_color_manual(name=NULL,values=c("More years than expected (FDR < 0.05)"="black"))+guides(color=guide_legend(override.aes=list(fill=NA,linewidth=0.5),order=2))
savepng(g46,"EHEC_F46_map_persistence.png")
top<-head(Pst,12);sig<-Pst[Pst$fdr_years_SIR_gt1<0.05,];cat(sprintf("[595] Kendall W %.3f (perm P=%.3g; perm 97.5%% %.3f; Poisson-null P=%.3g) excl2020 %.3f pre %.3f post %.3f · FDR<0.05 persistent: %d districts: %s\n",Wall,KW$perm_p[1],KW$perm_null_97.5[1],KW$poisson_p[1],Wno13,Wpre,Wpost,nrow(sig),paste(sprintf("%s %d/%.1f",sig$region,sig$years_SIR_gt1,sig$exp_years_SIR_gt1),collapse=", ")))
print(head(Pst,15));print(KW);print(T34b[,c("sido","total_cases","rate_2011_2024")])
sw<-sig%>%mutate(sido=str_extract(region,"^[가-힣]+?(시|도)"))%>%count(sido)%>%arrange(desc(n))
upd_caps(list(
 F45=list(file="EHEC_F45_map_yearly_SIR.png",title=sprintf("Yearly EB-smoothed age-standardized SIR of %s enterohemorrhagic E. coli infection by district, 2011-2024.",LAB),
   caption=sprintf("Empirical-Bayes-smoothed age-standardized standardized incidence ratio (SIR) of %s enterohemorrhagic E. coli infection by district for each year 2011-2024, Republic of Korea, on a common scale.",LAB),
   legend="Same index in every panel: SIR = observed / age-standardized expected cases (national age-specific reference rates applied to each district's age structure; identical E_std construction across years), Marshall global empirical-Bayes smoothing within each year, common log-symmetric classes around 1 (RdBu). 2020 = largest outbreak year (Ansaon outbreak year. Korea 2000 / Unified CS (EPSG:5179).",
   yoji=sprintf("연도별 동일 인덱스(연령표준화 SIR, 공통 계급). 인구규모 보정 후 기대보다 많은 해에 SIR>1인 구(FDR<0.05) %d개(%s); 연도간 순위 일치 Kendall W = %.3f (순열 P = %.4f, Poisson 영모형 P = %.2f).",nrow(sig),if(nrow(sig))paste(sprintf("%s %d",sw$sido,sw$n),collapse=", ") else "없음",Wall,KW$perm_p[1],KW$poisson_p[1])),
 F45b=list(file="EHEC_F45b_map_yearly_relSIR.png",title=sprintf("Yearly within-year relative SIR of %s enterohemorrhagic E. coli infection by district, 2011-2024.",LAB),
   caption=sprintf("EB-smoothed age-standardized SIR of %s enterohemorrhagic E. coli infection divided by the national SIR of the same year, by district and year, 2011-2024.",LAB),
   legend="Same index as F45 but rescaled so that the national level of each year equals 1; this removes the national decline (2020-2024) and shows the geographic pattern within each year on a common log-symmetric scale.",
   yoji=sprintf("연도 내 상대 SIR: 국가 수준 감소를 제거한 지리 패턴. 전국 SIR(절대) %s.",paste(sprintf("%d %.2f",nat$year,nat$natSIR),collapse=", "))),
 F46=list(file="EHEC_F46_map_persistence.png",title="Persistence of high enterohemorrhagic E. coli infection risk by district, 2011-2024.",
   caption=sprintf("Number of years (of 14) in which each district had at least one case and an age-standardized SIR above 1, %s enterohemorrhagic E. coli infection, 2011-2024; black outline = more such years than expected under a Poisson null that preserves district population size (BH-FDR < 0.05).",LAB),
   legend=paste0("Null: counts drawn as Poisson(E_std × national SIR of the year) for every district-year (",format(B,big.mark=",")," draws), so a small district that reaches SIR > 1 with a single case is not rewarded; P = (b + 1)/(B + 1), Benjamini-Hochberg across 229 districts. Companion tables: T34 (district list with expected years and FDR), T34c (Kendall W of yearly district ranks)."),
   yoji=sprintf("인구규모 보정 후 지속 고위험(FDR<0.05) %d개 구: %s. SIR>1 연수 최다: %s.",nrow(sig),if(nrow(sig))paste(sprintf("%s %d년(기대 %.1f)",sig$region,sig$years_SIR_gt1,sig$exp_years_SIR_gt1),collapse=", ") else "없음",paste(sprintf("%s %d",head(Pst%>%arrange(desc(years_SIR_gt1)),5)$region,head(Pst%>%arrange(desc(years_SIR_gt1)),5)$years_SIR_gt1),collapse=", "))),
 T34=list(file="EHEC_T34_persistence.csv",title="District-level persistence of high enterohemorrhagic E. coli infection risk, 2011-2024.",caption="Years with at least one case and raw SIR above 1 (observed vs expected under the population-preserving Poisson null, with permutation P and BH-FDR), years in the top quintile of raw SIR among districts with cases, mean yearly rank and cumulative SIR by district (sorted by P).",legend="Expected years and P from 2,000 Poisson(E_std × national SIR) draws; rank 1 = highest raw SIR in that year among 229 districts (zero-case districts tied at the bottom); top quintile = 80th percentile of raw SIR among districts with at least one case that year.",yoji=sprintf("Poisson 영모형 대비 SIR>1 연수 초과 상위: %s (FDR<0.05 %d개).",paste(sprintf("%s %d/%.1f년 P=%.3f",top$region[1:6],top$years_SIR_gt1[1:6],top$exp_years_SIR_gt1[1:6],top$p_years_SIR_gt1[1:6]),collapse="; "),nrow(sig))),
 T34b=list(file="EHEC_T34b_sido_yearly_rate.csv",title="Yearly enterohemorrhagic E. coli infection incidence per 100,000 by province, 2011-2024.",caption=sprintf("Crude %s enterohemorrhagic E. coli infection incidence per 100,000 by province and year, with the 2011-2024 cumulative rate.",LAB),legend="Population = district resident registration population summed by province.",yoji=sprintf("2011-24 누적 발생률 상위 시도: %s.",paste(sprintf("%s %.2f",T34b$sido[1:5],T34b$rate_2011_2024[1:5]),collapse=", "))),
 T34c=list(file="EHEC_T34c_kendallW.csv",title="Year-to-year concordance of district risk ranks (Kendall W).",caption="Kendall coefficient of concordance of yearly district ranks of raw age-standardized SIR (zero-case districts tied at the bottom), all years, excluding 2013, pre- and post-2020, against two null models.",legend=paste0("W is corrected for ties (zero-case districts share the bottom rank); the uncorrected value is also given. Observed and simulated ranks use the same SIR (cases / age-standardized expected count) and tie handling. Null (i) permutation: district labels shuffled within each year (",format(B,big.mark=",")," draws), i.e. districts exchangeable given the year; null (ii) parametric: counts drawn as Poisson(E_std × national SIR of the year), which preserves each district's population size and the yearly national level but has no spatial heterogeneity. P = (b + 1)/(B + 1)."),yoji=sprintf("Kendall W 전체 %.3f: 순열 영모형 z = %.1f (P = %.4f), Poisson(E×전국SIR) 영모형 z = %.1f (P = %.4f); 2011-19 z = %.1f/%.1f, 2020-24 z = %.1f/%.1f (순열/Poisson). 동률이 많아 W 절대값은 작고 집합 간 비교 불가, z로 판단.",Wall,KW$perm_z[1],KW$perm_p[1],KW$poisson_z[1],KW$poisson_p[1],KW$perm_z[3],KW$poisson_z[3],KW$perm_z[4],KW$poisson_z[4]))))
