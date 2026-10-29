#!/bin/zsh
# mv_pool_v2d.sh — (2026-10-21 대장 속도 조정) 533 다중우주 씨드 101-112를 워커 6개로 돌리는 관리자. run_chain_v2d.sh의 다중우주·마무리 부분을 넘겨받음.
#  이미 도는 101·102(INLA 2:1, run_chain_v2d가 띄움)는 그대로 두고 PID로 지켜봄. 새 워커는 INLA_THREADS=1:1. 하나 끝나면 다음 씨드.
#  워커 PID는 logs/mv_workers_v2d.pids에 기록, 끝날 때·중단될 때 남은 워커를 정리. 12씨드 모두 '완료' 확인 뒤 534 → 543 → '=== 전체 끝 ==='(finalize_after_mv.sh가 이어받음).
ROOT="${EHEC_ROOT:-.}"
R="$ROOT/_v2_rerun_261022/JKMS"; L="$R/logs"; CL="$L/mvclaim_v2d"; PF="$L/mv_workers_v2d.pids"; cd "$ROOT"
export TMPDIR="${TMPDIR:-/tmp}/ehec_rtmp"
unset EHEC_HM EHEC_FAM EHEC_PANEL EHEC_RS EHEC_DESIGN EHEC_SH EHEC_CASES CJ_MASK V2_POOL_DROP_NA
say(){ echo "$(date '+%m-%d %H:%M:%S')  [v2d] $*" >> "$L/chain.log"; }
MAXW=6; SEEDS=(101 102 103 104 105 106 107 108 109 110 111 112)
typeset -A PID; for kv in ${=EXIST}; do PID[${kv%%:*}]=${kv##*:}; done   # EXIST="101:16781 102:14379"
cleanup(){ for sd in ${(k)PID}; do kill -0 ${PID[$sd]} 2>/dev/null && { pkill -TERM -P ${PID[$sd]} 2>/dev/null; kill ${PID[$sd]} 2>/dev/null; say "정리: seed$sd 워커 PID ${PID[$sd]} 종료"; }; done; }
trap 'cleanup; say "mv_pool 중단(신호)"; exit 1' INT TERM
done_seed(){ grep -q "^완료: 사양 250" "$L/v2d_533_seed$1.log" 2>/dev/null; }
launch(){ local sd=$1; ( SEED=$sd N_SPEC=250 INLA_THREADS=1:1 Rscript 533_prereg_multiverse.R >> "$L/v2d_533_seed$sd.log" 2>&1; echo "exit=$?" >> "$L/v2d_533_seed$sd.log" ) &
  local wp=$!; PID[$sd]=$wp; echo "$sd ${PID[$sd]} $(date '+%m-%d %H:%M:%S') INLA_THREADS=1:1" >> "$PF"; say "533 seed$sd 시작(mv_pool, PID ${PID[$sd]}, INLA 1:1)"; }
say "=== mv_pool 시작(워커 최대 $MAXW, 기존 ${(kv)PID}) ==="
for sd in ${(k)PID}; do echo "$sd ${PID[$sd]} 기존(run_chain_v2d, INLA 2:1)" >> "$PF"; done
TRIES=()
while true; do
  for sd in ${(k)PID}; do
    if ! kill -0 ${PID[$sd]} 2>/dev/null; then
      if done_seed $sd; then say "533 seed$sd 완료(mv_pool 확인)"; else say "533 seed$sd 워커 끝났으나 '완료' 없음 → 다시 대기열"; rmdir "$CL/seed$sd" 2>/dev/null; fi
      unset "PID[$sd]"
    fi
  done
  for sd in $SEEDS; do
    (( ${#PID} >= MAXW )) && break
    done_seed $sd && continue
    [[ -n ${PID[$sd]} ]] && continue
    n=${#${(M)TRIES:#$sd}}; (( n >= 3 )) && continue
    mkdir "$CL/seed$sd" 2>/dev/null || continue
    TRIES+=($sd); launch $sd
  done
  ndone=0; for sd in $SEEDS; do done_seed $sd && ((ndone++)); done
  (( ndone == ${#SEEDS} )) && break
  (( ${#PID} == 0 )) && { say "mv_pool: 남은 워커 없음, 완료 ${ndone}/12 — 중단"; exit 1; }
  sleep 60
done
cleanup
say "=== 다중우주 끝(12/12, mv_pool) ==="
t0=$SECONDS; python3 534_prereg_mv_summary.py > "$L/v2d_534_prereg_mv_summary.log" 2>&1; say "534_prereg_mv_summary.py exit=$? ($(( (SECONDS-t0)/60 ))분)"
t0=$SECONDS; python3 543_final_spec_curve.py > "$L/v2d_543_final_spec_curve.log" 2>&1; say "543_final_spec_curve.py exit=$? ($(( (SECONDS-t0)/60 ))분)"
say "=== 전체 끝 ==="
