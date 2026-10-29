#!/bin/zsh
# JKMS EHEC v2d 재실행 체인(2026-10-21 대장 결정 (가)). 워커 2개.
#  W2 = 다중우주 씨드 큐(처음부터). W1 = A단계(공유 수정과 무관한 단계) → 다중우주 합류, 사이사이 shared_fix_done.flag가 생기면 B단계(na.rm 공유 수정 대상 단계와 그 산출을 쓰는 그림).
ROOT="${EHEC_ROOT:-.}"
R="$ROOT/_v2_rerun_261022/JKMS"; L="$R/logs"; FLAG="$ROOT/_v2_promote_261021/shared_fix_done.flag"; cd "$ROOT"
export TMPDIR="${TMPDIR:-/tmp}/ehec_rtmp"; mkdir -p "$TMPDIR" "$L/mvclaim_v2d"
unset EHEC_HM EHEC_FAM EHEC_PANEL EHEC_RS EHEC_DESIGN EHEC_SH EHEC_CASES CJ_MASK V2_POOL_DROP_NA
export INLA_THREADS="2:1" PERSIST_B=100000 PERSIST_SEED=20260913 EHEC_CAPF="$R/EHEC_SERIES_CAPTIONS_JKMS.json"
say(){ echo "$(date '+%m-%d %H:%M:%S')  [v2d] $*" >> "$L/chain.log"; }
runR(){ local s=$1; shift; local t0=$SECONDS; Rscript "$s" "$@" > "$L/v2d_${s%.R}.log" 2>&1; local e=$?; say "$s $* exit=$e ($(( (SECONDS-t0)/60 ))분)"; }
runP(){ local s=$1; shift; local t0=$SECONDS; python3 "$s" "$@" > "$L/v2d_$(basename ${s%.py}).log" 2>&1; local e=$?; say "$(basename $s) exit=$e ($(( (SECONDS-t0)/60 ))분)"; }
SEEDS=(101 102 103 104 105 106 107 108 109 110 111 112)
mvone(){ local w=$1
  for sd in $SEEDS; do
    mkdir "$L/mvclaim_v2d/seed$sd" 2>/dev/null || continue
    local t0=$SECONDS; say "533 seed$sd 시작 ($w)"
    SEED=$sd N_SPEC=250 Rscript 533_prereg_multiverse.R > "$L/v2d_533_seed$sd.log" 2>&1; local e=$?
    say "533 seed$sd exit=$e ($w, $(( (SECONDS-t0)/60 ))분)"; return 0
  done; return 1; }
phaseB(){
  say "--- B단계 시작(shared_fix_done.flag 확인) ---"
  for s in 542_final_period_analyses.R 542w_final_windows_3cut.R 542c_covid_phi_ci.R 553_final_ladder_zinb.R 553z_nb_vs_zinb_irr.R 546b_monitor_beef_sens.R 499_covid_map.R; do runR $s; done
  runP 592_ehec_window_grid.py
  runP 544_det_windows_final.py
  say "--- B단계 끝 ---"; }
say "=== 시작(v2d 패널 md5 $(md5 -q data_cache/ST_prereg_panel.csv)) ==="
( while mvone W2; do :; done; say "MV worker W2 끝" ) &; P2=$!
# A단계
for s in 415_agestd_exceedance.R 540_final_main_model.R; do runR $s; done
for sv in nopc1 jejudummy cj0; do mkdir -p "$R/sens_$sv"; t0=$SECONDS
  if [[ $sv == cj0 ]]; then CJ_MASK=0 SENS=$sv EHEC_RS="$R/sens_$sv" Rscript "$R/540s_sensitivity.R" > "$L/v2d_540s_$sv.log" 2>&1; else SENS=$sv EHEC_RS="$R/sens_$sv" Rscript "$R/540s_sensitivity.R" > "$L/v2d_540s_$sv.log" 2>&1; fi
  say "540s SENS=$sv exit=$? ($(( (SECONDS-t0)/60 ))분)"; done
for s in 542x_ehec_window_interaction.R 595_ehec_yearly_persistence.R 595b_ehec_f46_redraw.R; do runR $s; done
runP "$R/33_persistence_jkms_100k_rerun.py"
runR 411_agestd_S6.R
runR 411b_ehsa_esri.R
runR 412c_mapB_esri.R
runR 451_agestd_model_generic.R 장출혈성대장균감염증 EHEC
runR 546_census_beef_model.R
say "=== A단계 끝, W1 다중우주 합류 ==="
BDONE=0
while true; do
  if [[ $BDONE == 0 && -f "$FLAG" ]]; then phaseB; BDONE=1; fi
  mvone W1 || break
done
if [[ $BDONE == 0 ]]; then say "B단계: shared_fix_done.flag 기다림"; until [[ -f "$FLAG" ]]; do sleep 60; done; phaseB; fi
wait $P2; say "=== 다중우주 끝 ==="
runP 534_prereg_mv_summary.py
runP 543_final_spec_curve.py
say "=== 전체 끝 ==="
