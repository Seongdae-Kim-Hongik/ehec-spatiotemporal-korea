# 시군구 결합 규칙 (2026-10-21)

시군구 자료는 이름으로 결합하지 않는다. 모든 자료는 먼저 `sigungu_cd`(5자리, 229 체계)로 바꾸고, 결합한 패널은 `panel_gate`를 통과한 뒤에만 모형에 넣는다.

계기는 이름 결합이 개명·통합 지역에서 조용히 빠지고, 빠진 칸을 전국 최솟값 등으로 조용히 채운 일이다. EHEC 인천 미추홀구 2019-2024 기대건수가 1.14 대신 0.0198로 들어가 '고위험'으로 나왔고, final.shp에서는 인천 23030·23050의 region 라벨이 서로 바뀌어 있었다(`FBD_SpatioTemporal/_shpfix_impact_261020/`).

## 네 가지 규칙

1. **이름 결합 금지.** 어떤 자료든 `to_code()` 또는 `attach_code()`로 코드화한 뒤 `sigungu_cd`(와 연도)로 결합한다. `region` 문자열끼리 `merge`/`left_join`/`match` 하지 않는다. INLA 그래프 순서(`DD$area`)도 코드로 맞춘다.
2. **결측 채움 금지.** `is.na(x) <- 최솟값`, `approx(rule=2)`, `fillna(median)`, 전체 평균 대입처럼 조용한 채움을 쓰지 않는다. 꼭 필요하면 `explicit_fill(df, col, how, reason, log_path)`로 채우고, 채운 칸 목록을 분석 폴더에 남기고 Methods에 적는다. 이 함수는 'zero'(구조적 0)와 'interp'(같은 지역 안 보간)만 제공한다.
3. **모든 패널은 `panel_gate` 통과 뒤 모형.** 인구·기대건수 E·결과 열이 229 × 분석 연도 모든 칸에 있고 유한하며 인구·E > 0이어야 한다. 실패하면 멈춘다. 경고(채움 의심, 급변, 통합 연도 단절)는 하나씩 원자료로 확인하고 확인 결과를 기록한다.
4. **shp는 `read_shp_checked()`로만 읽는다.** `st_read("final.shp")`를 직접 쓰지 않는다. 이 함수는 `sigungu_cd`로 `canonical_region`을 다시 붙이고 코드·이름·면적(±5%)을 단언하므로, dbf의 region 라벨이 틀려도 결과가 맞다.

### 함께 지킬 것

- 통합·개명이 있는 해의 합산은 **(코드, 연도)마다 같은 기준 시점의 행만** 더한다. 이름마다 '12월, 없으면 마지막 달'을 따로 고르면 옛 이름(6월)과 새 이름(12월)이 함께 더해져 인구가 두 배가 된다. 주민등록 시험에서 미추홀 2018, 세종 2012, 군위 2023, 청주 2014가 이렇게 부풀었다.
- 통합 전 단위(청원군·연기군·마산시·진해시·북제주군·남제주군)는 같은 코드로 붙는다. 그 해들은 `(sigungu_cd, year)`로 합산해야 하고, 합산하지 않으면 `panel_gate`가 중복으로 멈춘다.
- 일반구·출장소·보건소 구역 표기(예: 고양시 덕양구, 송탄출장소, 제주시 동부)는 기본으로 막힌다. 시 합계 행과 함께 더하면 두 번 셈해지기 때문이다. 시 합계 행이 없는 자료에서만 `subdistrict="city"`로 소속 시 코드를 붙인다.
- 시도 합계·전국 행은 `drop_kinds=c("시도","전국")`로 명시해서 뺀다. 그러면 뺀 행 수가 출력된다.
- 시도 없이 쓴 중복 이름(중구·동구·서구·남구·북구·강서구·고성군)과 '광주시' 단독 표기는 모호해서 실패한다. 원자료에서 시도를 붙여 다시 만든다.

## 파일

```
sigungu_key.csv            229 코드표: sigungu_cd, sido, sigungu, canonical_region(FBD 패널 region 표기), area_km2(final.shp 기하, EPSG:5179), sido_2024, shp_sigungu_nm, note
sigungu_alias.csv          별칭 2,329줄: alias, sigungu_cd, valid_from, valid_to, source_note, needs_review, kind, sido_part, sigungu_part (needs_review=1 78줄은 표 끝)
sigungu_key.py / .R        to_code, attach_code, name_kind, explicit_fill, read_shp_checked (두 언어 결과 동일, 2,979개 표기로 대조함)
panel_gate.py / .R         결합 뒤 모형 전 게이트
_build_sigungu_tables.py   두 표를 다시 만드는 생성기(--names _관측_시군구표기_261021.csv)
_harvest_sigungu_names.py  원자료·캐시에서 시군구 표기를 모으는 수집기(읽기만 함)
_관측_시군구표기_261021.csv 수집된 표기 3,398개와 나온 파일·연도
_시험_261021.py / .R       실데이터 시험. 결과 요약은 시험결과_261021.md
```

## 사용 예

R
```r
source("gate/panel_gate.R")  # sigungu_key.R도 함께 읽음
shp <- read_shp_checked()                                           # region = 코드 기준 canonical_region
E   <- attach_code(read.csv("ST_Estd_panel.csv"), "region", "year")  # 못 찾는 이름이 있으면 여기서 멈춤
DD  <- merge(expand.grid(sigungu_cd = shp$sigungu_cd, year = 2011:2024, stringsAsFactors = FALSE), E[, c("sigungu_cd", "year", "E_std")], all.x = TRUE)
panel_gate(DD, cols = "E_std", positive = "E_std", years = 2011:2024, code_col = "sigungu_cd")  # 빈칸 있으면 stop
DD$area <- match(DD$sigungu_cd, shp$sigungu_cd)                     # 그래프 순서도 코드로
```

Python
```python
import sys; sys.path.insert(0, "gate")
from sigungu_key import attach_code, read_shp_checked; from panel_gate import panel_gate
pop = attach_code(raw, "행정구역", "year", drop_kinds=["시도", "전국", "일반구"])        # 주민등록 원자료
pop = pop[pop.mo == pop.groupby(["sigungu_cd", "year"]).mo.transform("max")].groupby(["sigungu_cd", "year"], as_index=False)["인구수"].sum()
panel_gate(pop, cols=["인구수"], positive=["인구수"], years=range(2011, 2025), code_col="sigungu_cd")
```

명령줄
```
python3 sigungu_key.py "인천광역시 미추홀구" "경상북도 군위군" --year 2020
python3 panel_gate.py ST_Estd_panel.csv --name-col region --cols E_std --positive E_std --years 2011-2024
```

## to_code가 하는 일

1. NFC로 바꾸고, 괄호 속 숫자 코드(`(2817700000)`)와 공백·`-`·`_`를 지운다.
2. 별칭 원문과 똑같으면 그 코드를 쓴다.
3. 앞의 시도 표기(서울특별시/서울시/서울, 강원도/강원특별자치도/강원, 전라북도/전북특별자치도/전북 등)를 약칭으로 통일하고, 같은 시도가 두 번 적힌 경우(`경기도경기 여주군`)도 벗긴 뒤 (시도, 시군구)로 찾는다.
4. 그래도 없으면 끝의 시·군·구를 생략한 식약처형 표기(`경남 마산`, `강원 양구`)로 찾는다.
5. 시도 없는 이름은 전국에 하나뿐일 때만 받는다(`부산진구`, `제주시`).
6. 별칭의 유효기간 밖 연도면 경고한다(같은 코드로 붙이되 소급 표기인지 확인하라는 뜻). `strict_period=TRUE`면 멈춘다.
7. 하나라도 실패하면 못 찾은 이름과 이유를 모두 보여 주고 멈춘다.

## 별칭 표의 kind

| kind | 뜻 | 예 |
|---|---|---|
| 표준 | 코드표 이름의 표준 표기 | 인천시남구, 인천광역시 남구, 인천 남구 |
| 표기변형 | 원자료에서 실제로 나온 다른 표기 | 인천광역시 남구 (2817000000), 경기도 수원시 |
| 개명 | 이름이 바뀐 단위 | 인천 미추홀구(2018-07-01부터) |
| 시도변경 | 소속 시도가 바뀐 단위 | 경상북도 군위군(2023-06-30까지), 대구광역시 군위군 |
| 승격 | 군→시, 출장소→시·군 | 여주군, 당진군, 계룡출장소, 증평출장소 |
| 통합전 | 다른 단위와 합쳐진 옛 단위(합산 필요) | 청원군, 연기군, 마산시, 진해시 |
| 일반구 | 시의 하위 구·출장소·보건소 구역(기본 차단) | 고양시 덕양구, 송탄출장소, 제주시 동부 |
| 검토필요 | 자동 매칭 금지(needs_review=1) | 남구, 고성군, 거제군, 충북 보령, 인천광역시 제물포구 |

## 별칭 추가

`to_code`가 멈추면 원자료의 표기와 실제 단위를 확인한 뒤 `sigungu_alias.csv`에 한 줄을 더한다. 출처(파일·연도)를 `source_note`에 적고, 애매하면 `needs_review=1`로 둔다. 표를 통째로 다시 만들 때는 `_build_sigungu_tables.py`를 쓴다.

## panel_gate 판정

| 구분 | 항목 |
|---|---|
| 실패 | 코드화 실패, 229에 없는 코드, (코드, 연도) 중복, 지정 열의 행 없음·결측·비유한, 인구·E ≤ 0 |
| 경고 | 전국 최솟값·최댓값과 같은 값의 연속·반복, 양수 열의 같은 값 3년 이상 연속(상수 채움·rule=2 외삽 흔적) |
| 경고 | 그해 전국 중앙 배율 대비 1.3배 넘는 지역 변화(전국 발생률을 따라 매년 움직이는 E도 지역 고유 급변만 걸림) |
| 경고 | 통합 연도(청주 2014, 세종 2012, 창원 2010, 제주·서귀포 2006) 직전→그해 ±10% 넘는 단절 |
| 경고(선택) | `area_col`이 sigungu_key.csv 기하 면적과 ±5% 넘게 다름 |
