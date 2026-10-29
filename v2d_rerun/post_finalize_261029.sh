#!/bin/zsh
# post_finalize_261029.sh — finalize 끝난 뒤 00_결과요약.md 진행 상태 표를 '끝'으로 바꾸고 완료 표시 파일을 남김.
ROOT="${EHEC_ROOT:-.}"; J="$ROOT/_v2_rerun_261022/JKMS"; L="$J/logs"
until grep -q "=== finalize 끝 ===" "$L/chain.log"; do sleep 120; done
python3 - <<'PY'
import re, os
F=os.path.join(os.environ.get("EHEC_ROOT","."),"_v2_rerun_261022/JKMS/00_결과요약.md")
s=open(F,encoding="utf-8").read()
s=s.replace("| 다중우주 533 씨드 101-112 × 250 | 진행 중(워커 2개, 씨드당 약 3.3시간) |","| 다중우주 533 씨드 101-112 × 250 | 끝(10/21 중단 뒤 10/29 이어하기, 워커 6개 INLA 1:1, `logs/mv_workers_v2d.pids`) |")
s=s.replace("| 534·543 → finalize(S15·보충표·그림·Table 2·원고 묶음·이 문서 다중우주 절) | `run_chain_v2d.sh`와 `finalize_after_mv.sh`가 자동으로 함. 진행은 `logs/chain.log`의 [v2d] 줄 |","| 534·543 → finalize(S15·보충표·그림·Table 2·원고 묶음·이 문서 다중우주 절) | 끝(`logs/chain.log` '=== finalize 끝 ===') |")
s=s.replace("| 원고 조립(새 버전 폴더) | 1차 생성 끝(다중우주 칸 제외) |","| 원고 조립(새 버전 폴더) | 끝(다중우주 칸 포함, `JKMS_FBD_v2c/JKMS_EHEC/`) |")
s=s.replace("- **다중우주 3,000사양은 진행 중임**('진행 상태' 참고). 끝나면 자동 처리 스크립트가 이 문서의 다중우주 절을 채우고 보충표·그림·원고 묶음을 다시 만듦.","- **다중우주 3,000사양 결과는 아래 다중우주 절에 있음**(보충표 S12-S15·Fig S17·원고 묶음 갱신 완료).")
s=s.replace("보충그림은 S14(초과확률), S15(COVID 전후 RR), S16(창별 결정요인)을 새 그림으로 바꿨고, S17(사양 곡선)은 다중우주 뒤에 바꿈.","보충그림은 S14(초과확률), S15(COVID 전후 RR), S16(창별 결정요인), S17(사양 곡선)을 새 그림으로 바꿨음.")
open(F,"w",encoding="utf-8").write(s)
PY
echo "$(date '+%m-%d %H:%M:%S')  [v2d] post_finalize: 00_결과요약.md 진행 상태 갱신 exit=$?" >> "$L/chain.log"
date > "$J/ALL_DONE_261029.flag"; cp "$L/chain.log" "$J/ALL_DONE_261029.flag.chainlog"
