#!/bin/zsh
# finalize_after_mv.sh — 다중우주(run_chain.sh의 534·543까지) 끝난 뒤 S15 계산 → 보충표·그림 갱신 → 새 버전 폴더 원고 묶음 재생성 → 00_결과요약.md 다중우주 절 채움.
ROOT="${EHEC_ROOT:-.}"; R="$ROOT/_v2_rerun_261022/JKMS"; L="$R/logs"; N="$R/JKMS_FBD_v2c"
say(){ echo "$(date '+%m-%d %H:%M:%S')  $*" >> "$L/chain.log"; }
until grep -q "=== 전체 끝 ===" "$L/chain.log"; do sleep 60; done
cd "$ROOT"
n=$(ls results_ST/PREREG_MV_results_seed1[01][0-9].csv | wc -l | tr -d ' '); [[ $n == 12 ]] || { say "finalize 중단: 씨드 결과 파일 $n개"; exit 1; }
python3 "$R/update_table2_v2c.py" > "$L/update_table2.log" 2>&1; say "update_table2_v2c exit=$?"
python3 "$R/s15_axis_robustness.py" > "$L/s15.log" 2>&1; say "s15_axis_robustness exit=$?"
python3 "$R/update_supp_tables_v2c.py" > "$L/update_supp_tables.log" 2>&1; say "update_supp_tables_v2c exit=$? ($(tail -1 $L/update_supp_tables.log | cut -c1-200))"
python3 "$R/update_supp_figures_v2c.py" > "$L/update_supp_figures.log" 2>&1; say "update_supp_figures_v2c exit=$? ($(tail -1 $L/update_supp_figures.log))"
cd "$N"; python3 gen_jkms_supp_ehec.py > "$L/gen_supp_v2c.log" 2>&1; say "gen_jkms_supp_ehec(v2c 사본) exit=$?"
python3 gen_jkms_fbd.py EHEC > "$L/gen_fbd_v2c.log" 2>&1; say "gen_jkms_fbd(v2c 사본) exit=$?"
python3 tools/split_manuscript.py EHEC > "$L/split_v2c.log" 2>&1; say "split_manuscript(v2c 사본) exit=$?"
cd "$ROOT"; python3 "$R/mv_section.py" > "$L/mv_section.log" 2>&1; say "mv_section exit=$?"
say "=== finalize 끝 ==="
