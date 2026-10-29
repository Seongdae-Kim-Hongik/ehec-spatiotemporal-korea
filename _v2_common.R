# _v2_common.R — (2026-10-21) FBD 시공간 셋업·모형 스크립트 공통 v2 규칙(논문_마스터자료/00_분석전_데이터게이트/README_시군구결합_규칙.md).
#  1) 기대건수 E·사례·공변량은 sigungu_cd로 결합. 표의 region 라벨은 코드에서 shp(read_shp_checked) 라벨로 다시 붙임 → 이름이 달라 조용히 빠지는 일 없음.
#  2) 결측 E를 전국 최솟값·이웃 해·같은 지역 보간으로 채우지 않음. E·사례 파일은 읽을 때 panel_gate(229 × 연도 모든 칸, E > 0)를 통과해야 함.
#  3) 공변량 표준화(v2_zb)는 평균 채움을 하지 않고 결측이 있으면 멈춤. 패널 전체를 표준화하는 셋업은 결측 있는 열을 Z에서 빼고(이름 출력),
#     그 변수를 모형에 넣으면 멈춤(v2_need). 채움이 필요하면 빌더에서 explicit_fill로 채우고 data_cache/_v2_logs에 기록.
#  4) KDCA형(EHEC·이질·장티푸스·파라티푸스 등 질병청 교차표) 청주 2011-2013은 사례를 NA로 마스킹(교차표에 옛 청주시·청원군 신고가 없음).
#     INLA는 이 칸을 우도에서 빼고 예측만 함. 원래 관측(0)은 DD$cases_obs에 남김. CJ_MASK=0이면 마스킹하지 않음(민감도 '청주 관측 0').
#  사용: source(file.path(ROOT,"_v2_common.R")) 뒤 v2_DD()·v2_Z()·v2_join()·v2_need() 또는 인라인 스크립트용 v2_cases_rc()·v2_E_rc()·v2_relabel()·v2_finish()
V2_ROOT <- Sys.getenv("EHEC_ROOT", unset=".")
if (!exists("V2_CACHE")) V2_CACHE <- file.path(V2_ROOT, "data_cache")
V2_GATE_DIR <- Sys.getenv("EHEC_GATE", unset=file.path(Sys.getenv("EHEC_ROOT", unset="."),"gate"))
if (!exists("to_code", mode = "function")) source(file.path(V2_GATE_DIR, "sigungu_key.R"), local = globalenv())   # panel_gate.R의 자기 위치 찾기가 중첩 source에서 틀리므로 먼저 읽음
if (!exists("panel_gate", mode = "function")) source(file.path(V2_GATE_DIR, "panel_gate.R"), local = globalenv())
assign(".SGK_DIR", V2_GATE_DIR, envir = globalenv())   # sigungu_key.R도 sys.frame(1)$ofile로 자기 위치를 찾아, 셋업 안에서 source하면 셋업 폴더를 가리킴 → 게이트 폴더로 고정
V2_CJ_CODE <- "33010"; V2_CJ_YEARS <- 2011:2013
if (!exists("V2_NA_VARS")) V2_NA_VARS <- character(0)
v2_cj_on <- function() !identical(Sys.getenv("CJ_MASK"), "0")
.v2e <- new.env()
.v2_cd <- function(x) formatC(as.integer(sub("\\.0$", "", trimws(as.character(x)))), width = 5, flag = "0")

v2_shp_lab <- function() {   # read_shp_checked 행 순서의 (코드, region 라벨) — 셋업의 shp$region(공백 제거·미추홀→남구)과 같은 표기
  if (is.null(.v2e$lab)) { s <- suppressWarnings(suppressMessages(read_shp_checked())); .v2e$lab <- data.frame(sigungu_cd = .v2_cd(s$sigungu_cd), region = gsub("\\s+", "", as.character(s$region)), stringsAsFactors = FALSE) }
  .v2e$lab
}
v2_codes_for_regs <- function(regs) {   # 셋업 그래프 순서(regs)의 sigungu_cd. 라벨이 shp 코드 라벨과 하나라도 다르면 멈춤
  L <- v2_shp_lab(); k <- match(regs, L$region)
  if (anyNA(k)) stop("[v2] 셋업 shp region과 read_shp_checked 라벨 불일치: ", paste(regs[is.na(k)], collapse = ", "), call. = FALSE)
  L$sigungu_cd[k]
}
v2_code <- function(nm) {   # region 표기 → sigungu_cd(못 찾으면 to_code가 멈춤)
  nm <- gsub("\\s+", "", as.character(nm)); u <- unique(nm)
  if (is.null(.v2e$cc)) .v2e$cc <- character(0)
  new <- setdiff(u, names(.v2e$cc)); if (length(new)) .v2e$cc <- c(.v2e$cc, setNames(to_code(new), new))
  unname(.v2e$cc[nm])
}
v2_relabel <- function(tab, what = "표", keep_code = TRUE) {   # 코드 기준으로 region을 shp 라벨로 다시 붙임(sigungu_cd 열이 있으면 그것을, 없으면 region을 코드화)
  tab <- as.data.frame(tab)
  cd <- if ("sigungu_cd" %in% names(tab)) .v2_cd(tab$sigungu_cd) else v2_code(tab$region)
  L <- v2_shp_lab(); k <- match(cd, L$sigungu_cd)
  if (anyNA(k)) stop(sprintf("[v2] %s: 229 그래프에 없는 코드 %s", what, paste(unique(cd[is.na(k)]), collapse = ",")), call. = FALSE)
  tab$region <- L$region[k]; tab$sigungu_cd <- cd
  if ("year" %in% names(tab) && anyDuplicated(tab[, c("sigungu_cd", "year")])) stop(sprintf("[v2] %s: (sigungu_cd, year) 중복", what), call. = FALSE)
  if (!keep_code) tab$sigungu_cd <- NULL
  tab
}
.v2_read <- function(file) { f <- if (file.exists(file)) file else file.path(V2_CACHE, file); read.csv(f, fileEncoding = "UTF-8-BOM") }
v2_E_rc <- function(file = "ST_Estd_panel.csv", years = 2011:2024) {   # E(v2 정본, sigungu_cd 포함) → region(shp 라벨)·year·E_std·sigungu_cd, 게이트 통과
  E <- .v2_read(file)
  if (!"sigungu_cd" %in% names(E)) stop(sprintf("[v2] %s에 sigungu_cd 없음 — v2 이전 판. 400/450/450b 빌더(v2)로 다시 만들 것", file), call. = FALSE)
  E <- E[E$year %in% years, c("sigungu_cd", "year", "E_std")]; E$sigungu_cd <- .v2_cd(E$sigungu_cd)
  suppressWarnings(panel_gate(E, cols = "E_std", positive = "E_std", years = years, code_col = "sigungu_cd", label = file, verbose = FALSE))
  v2_relabel(E, file)[, c("region", "year", "E_std", "sigungu_cd")]
}
v2_cases_rc <- function(file = "ST_EHEC_agecase.csv", years = 2011:2024) {   # 사례(v2 정본, 0 명시) → region·year·cases(연령 합)·sigungu_cd, 게이트 통과
  cs <- .v2_read(file)
  if (!"sigungu_cd" %in% names(cs)) stop(sprintf("[v2] %s에 sigungu_cd 없음 — v2 이전 판", file), call. = FALSE)
  cs$sigungu_cd <- .v2_cd(cs$sigungu_cd); cs <- cs[cs$year %in% years, ]
  a <- aggregate(cases ~ sigungu_cd + year, data = cs, FUN = sum)
  suppressWarnings(panel_gate(a, cols = "cases", years = years, code_col = "sigungu_cd", label = file, verbose = FALSE))
  v2_relabel(a, file)[, c("region", "year", "cases", "sigungu_cd")]
}
v2_cj_mask <- function(DD, kdca) {   # KDCA형 청주 2011-2013 사례 NA(표준). 원 관측은 cases_obs
  cj <- DD$sigungu_cd == V2_CJ_CODE & DD$year %in% V2_CJ_YEARS
  DD$cases_obs <- DD$cases; DD$cj <- as.integer((isTRUE(kdca) && v2_cj_on()) & cj)
  if (any(DD$cj == 1)) { DD$cases[DD$cj == 1] <- NA; message(sprintf("[v2] 청주 2011-2013 KDCA형 사례 NA 마스킹 %d칸(관측 %s; CJ_MASK=0이면 끔)", sum(DD$cj), paste(DD$cases_obs[DD$cj == 1], collapse = ","))) }
  else if (isTRUE(kdca) && !v2_cj_on()) message("[v2] CJ_MASK=0: 청주 2011-2013 관측 0을 그대로 씀(민감도)")
  DD
}
v2_DD <- function(regs, years, Efile, casefile, kdca, label = "") {   # 셋업용: shp 순서 격자(region이 빠르게 도는 expand.grid) + 코드 결합 E·사례 + 게이트(마스킹은 v2_cj_mask로 따로)
  DD <- expand.grid(region = regs, year = years, stringsAsFactors = FALSE)
  DD$sigungu_cd <- v2_codes_for_regs(regs)[match(DD$region, regs)]
  E <- v2_E_rc(Efile, years); cs <- v2_cases_rc(casefile, years); key <- paste(DD$sigungu_cd, DD$year)
  DD$cases <- cs$cases[match(key, paste(cs$sigungu_cd, cs$year))]; DD$E_std <- E$E_std[match(key, paste(E$sigungu_cd, E$year))]
  g <- suppressWarnings(panel_gate(DD, cols = c("E_std", "cases"), positive = "E_std", years = years, code_col = "sigungu_cd", label = label, verbose = FALSE))
  message(sprintf("[v2] %s E·사례 코드 결합 · panel_gate 통과(%d×%d칸, 경고 %d종: %s)", label, length(regs), length(years), length(g$warns), paste(names(g$warns), collapse = "; ")))
  DD
}
.v2_skew <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3 || sd(x) == 0) return(0); e1071::skewness(x) }
v2_zb <- function(x, v = "?") {   # 셋업 zb와 같은 변환(왜도 > 1이고 최솟값 ≥ 0이면 log1p → 표준화)에서 평균 채움만 뺀 것. 결측 있으면 멈춤
  if (anyNA(x)) stop(sprintf("[v2] 공변량 '%s' 결측 %d칸 — 평균 채움 금지. 빌더에서 explicit_fill(기록)로 채우거나 이 변수를 빼고 다시 실행", v, sum(is.na(x))), call. = FALSE)
  xv <- if (min(x) >= 0 && .v2_skew(x) > 1) log1p(pmax(x, 0)) else x; z <- scale(xv)[, 1]; z[!is.finite(z)] <- 0; z
}
v2_nofill <- function(x, v = "?") { if (anyNA(x)) stop(sprintf("[v2] 공변량 '%s' 결측 %d칸 — 평균 채움 금지(빌더에서 explicit_fill로 채울 것)", v, sum(is.na(x))), call. = FALSE); x }
v2_Z <- function(d0, vars, label = "패널") {   # 패널 전체 표준화: 결측 있는 열은 Z에서 빼고 이름을 남김(V2_NA_VARS) — 그 변수를 모형에 넣으면 v2_need가 멈춤
  Z <- d0[, c("region", "year")]; nav <- character(0)
  for (v in vars) { x <- d0[[v]]; if (is.null(x)) stop(sprintf("[v2] %s에 열 없음: %s", label, v), call. = FALSE); if (anyNA(x)) { nav <- c(nav, v); next }; Z[[v]] <- v2_zb(x, v) }
  if (length(nav)) message(sprintf("[v2] %s: 결측 있는 공변량 %d개는 Z에서 뺌(평균 채움 금지; 모형에 넣으면 멈춤): %s", label, length(nav), paste(nav, collapse = ", ")))
  V2_NA_VARS <<- unique(c(V2_NA_VARS, nav)); attr(Z, "na_vars") <- nav; Z
}
v2_join <- function(DD, tab, what = "공변량") {   # tab(region 표기·year)을 코드로 DD(sigungu_cd·year)에 붙임
  t2 <- v2_relabel(tab, what); t2$region <- NULL
  if (!"year" %in% names(t2)) return(dplyr::left_join(DD, t2, by = "sigungu_cd"))
  dplyr::left_join(DD, t2, by = c("sigungu_cd", "year"))
}
v2_need <- function(cv, data) {   # 모형에 넣을 공변량이 모두 있고 결측 없는지
  miss <- setdiff(cv, names(data))
  if (length(miss)) stop(sprintf("[v2] 모형 공변량 없음: %s%s", paste(miss, collapse = ", "), if (any(miss %in% V2_NA_VARS)) " — 결측이 있어 Z에서 뺀 변수. 빌더에서 explicit_fill로 채우거나 사양에서 뺄 것" else ""), call. = FALSE)
  bad <- cv[vapply(cv, function(v) anyNA(data[[v]]), TRUE)]
  if (length(bad)) stop(sprintf("[v2] 모형 공변량 결측: %s", paste(bad, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}
## ★2026-10-21 v2b(이전 판 _v2_common.R.bak_pre_v2b): 정지 지점 해결 규칙 (1) 원자료 sigungu_cd 결합 → (2) 구조적 결측은 외부 출처 대체 → (3) 그래도 남으면 세트·풀에서 빼고 기록.
##  (3)의 결정은 한 곳(_v2_promote_261021/15_세트풀_제외변수.csv)에만 적고, 스크립트는 scope_regex가 맞는 행을 읽어 뺌(평균 채움 금지, 결측이 실제로 있을 때만 뺌).
V2_EXCL_FILE <- file.path(V2_ROOT, "_v2_promote_261021", "15_세트풀_제외변수.csv")
v2_excl <- function(scope) {   # scope(예: "533_prereg_multiverse.R", "700:SHIG:REF_previous_main_260928")에 해당하는 (3) 제외 변수
  if (!file.exists(V2_EXCL_FILE)) return(character(0))
  x <- utils::read.csv(V2_EXCL_FILE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  unique(x$variable[vapply(x$scope_regex, function(r) grepl(r, scope, perl = TRUE), TRUE)])
}
v2_excl_log <- function(script, vars, data) {   # 실제로 뺀 변수 기록(_v2_logs/pool_dropped_na_vars.csv, 같은 형식)
  if (!length(vars)) return(invisible())
  lf <- file.path(V2_CACHE, "_v2_logs", "pool_dropped_na_vars.csv"); dir.create(dirname(lf), showWarnings = FALSE)
  utils::write.table(data.frame(time = format(Sys.time()), script = paste0(script, " [15_세트풀_제외변수.csv]"), variable = vars, na_cells = vapply(vars, function(v) if (is.null(data[[v]])) NA_integer_ else sum(is.na(data[[v]])), 1L)), lf, sep = ",", append = file.exists(lf), col.names = !file.exists(lf), row.names = FALSE, fileEncoding = "UTF-8")
}
v2_pool_check <- function(d0, vars, label = "후보 풀") {   # 다중우주·1단계처럼 풀 전체를 표준화하는 스크립트용: 결측 변수가 있으면 멈춤. V2_POOL_DROP_NA=1이면 빼고 기록
  nav <- vars[vapply(vars, function(v) !is.null(d0[[v]]) && anyNA(d0[[v]]), TRUE)]
  if (!length(nav)) { st <- intersect(v2_excl(label), vars); if (length(st)) message(sprintf("[v2b] %s: 제외 목록에 있으나 결측이 없어 유지: %s (15_세트풀_제외변수.csv 갱신 필요)", label, paste(st, collapse = ", "))); return(vars) }
  ex <- intersect(v2_excl(label), nav)   # ★v2b: (3) 결정된 변수는 기록하고 뺌
  if (length(ex)) {
    message(sprintf("[v2b] %s: 결측 변수 %d개를 (3) 결정대로 풀에서 뺌(15_세트풀_제외변수.csv): %s", label, length(ex), paste(sprintf("%s(%d)", ex, vapply(ex, function(v) sum(is.na(d0[[v]])), 1L)), collapse = ", ")))
    v2_excl_log(label, ex, d0); V2_NA_VARS <<- unique(c(V2_NA_VARS, ex)); vars <- setdiff(vars, ex); nav <- setdiff(nav, ex)
    if (!length(nav)) return(vars)
  }
  msg <- sprintf("[v2] %s에 결측 있는 변수 %d개: %s", label, length(nav), paste(sprintf("%s(%d)", nav, vapply(nav, function(v) sum(is.na(d0[[v]])), 1L)), collapse = ", "))
  if (!identical(Sys.getenv("V2_POOL_DROP_NA"), "1")) stop(msg, "\n  평균 채움 금지. 빌더에서 explicit_fill로 채우거나, 이 변수들을 풀에서 빼려면 V2_POOL_DROP_NA=1로 다시 실행(뺀 목록은 기록됨).", call. = FALSE)
  message(msg, " → V2_POOL_DROP_NA=1: 풀에서 뺌"); V2_NA_VARS <<- unique(c(V2_NA_VARS, nav))
  lf <- file.path(V2_CACHE, "_v2_logs", "pool_dropped_na_vars.csv"); dir.create(dirname(lf), showWarnings = FALSE)
  utils::write.table(data.frame(time = format(Sys.time()), script = label, variable = nav, na_cells = vapply(nav, function(v) sum(is.na(d0[[v]])), 1L)), lf, sep = ",", append = file.exists(lf), col.names = !file.exists(lf), row.names = FALSE, fileEncoding = "UTF-8")
  setdiff(vars, nav)
}
v2_finish <- function(DD, kdca = TRUE) {   # 인라인 스크립트용: 코드 결합한 E·사례가 모든 칸에 있는지 단언하고(채움 없음) 청주 마스킹
  if (!"sigungu_cd" %in% names(DD)) DD$sigungu_cd <- v2_codes_for_regs(DD$region)
  if (anyNA(DD$E_std) || any(DD$E_std <= 0)) stop(sprintf("[v2] E_std 결측·0 %d칸 — 최솟값 채움 금지(E 파일·결합 확인)", sum(is.na(DD$E_std) | DD$E_std <= 0, na.rm = TRUE)), call. = FALSE)
  if (anyNA(DD$cases)) stop(sprintf("[v2] 사례 결측 %d칸 — 0 채움 금지(사례 파일은 0을 명시해야 함)", sum(is.na(DD$cases))), call. = FALSE)
  if (anyDuplicated(DD[, c("sigungu_cd", "year")])) stop("[v2] (sigungu_cd, year) 중복", call. = FALSE)
  v2_cj_mask(DD, kdca)
}

v2_set_drop <- function(tag, sets, data_na_vars, script = "700_hm_setcompare.R") {   # 700 비교 세트: (3) 결정 변수 중 실제 결측(Z에서 빠진 변수)만 세트에서 뺌. 반환 list(sets, dropped data.frame)
  dr <- NULL
  for (sn in names(sets)) {
    ex <- intersect(intersect(v2_excl(sprintf("700:%s:%s", tag, sn)), unlist(sets[[sn]])), data_na_vars)
    if (length(ex)) { sets[[sn]] <- setdiff(unlist(sets[[sn]]), ex); dr <- rbind(dr, data.frame(set = sn, variable = ex)); message(sprintf("[v2b] %s %s: (3) 결정대로 세트에서 뺌 — %s", tag, sn, paste(ex, collapse = ", "))) }
  }
  list(sets = sets, dropped = dr)
}
v2_vif_set <- function(S, data, tag, outd = NULL, thr = 10) {   # 'VIF 10 초과 제외' 민감도 세트: 사전 규칙(주모형 가설 변수 중 VIF > thr를 뺌) 그대로, 대상만 현재 자료로 자동 계산
  vss <- S$rules$vif_sensitivity_set; if (is.null(vss) || !nzchar(vss)) return(S)
  main <- unlist(S$sets[[S$rules$main_set]]); cv <- unique(c(main, unlist(S$forced))); v2_need(cv, data)
  vd <- as.data.frame(data)[, cv, drop = FALSE]; set.seed(1); vd$.y <- stats::rnorm(nrow(vd))
  v <- car::vif(stats::lm(stats::as.formula(paste(".y~", paste(sprintf("`%s`", cv), collapse = "+"))), data = vd)); names(v) <- gsub("`", "", names(v))
  hi <- setdiff(names(v)[v > thr], unlist(S$forced)); old <- unlist(S$sets[[vss]]); new <- setdiff(main, hi)
  rec <- data.frame(set = vss, cov = names(v), vif = round(as.numeric(v), 2), forced = names(v) %in% unlist(S$forced), dropped = names(v) %in% hi, in_json_set = names(v) %in% old)
  if (!is.null(outd)) { dir.create(outd, showWarnings = FALSE, recursive = TRUE); utils::write.csv(rec, file.path(outd, sprintf("HM_vif_sensitivity_%s.csv", tag)), row.names = FALSE) }
  message(sprintf("[v2b] %s VIF 민감도 세트 %s = 주모형 가설 %d개 − VIF>%g %d개(%s) → %d개%s", tag, vss, length(main), thr, length(hi), if (length(hi)) paste(sprintf("%s %.1f", hi, v[hi]), collapse = ", ") else "없음", length(new),
                  if (setequal(new, old)) " (json 고정 목록과 같음)" else sprintf(" (json 고정 목록 %d개와 다름: json에서만 %s / 새로 %s)", length(old), paste(setdiff(old, new), collapse = ","), paste(setdiff(new, old), collapse = ","))))
  S$sets[[vss]] <- new; S$rules$vif_sensitivity_dropped <- hi; S
}
v2_jkms_main <- function() {   # JKMS 원고 주모형(540 원고 모드) 도메인·라벨을 540에서 그대로 읽음(목록을 한 곳에서만 관리)
  ex <- parse(file.path(V2_ROOT, "540_final_main_model.R"), keep.source = FALSE); e <- new.env()
  for (x in ex) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) && as.character(x[[2]]) %in% c("DOM", "ENG", "NEG")) eval(x, e)
  stopifnot(exists("DOM", e), exists("ENG", e)); list(DOM = e$DOM, vars = unlist(e$DOM, use.names = FALSE), ENG = e$ENG, NEG = if (exists("NEG", e)) e$NEG else character(0))
}
V2_DZ_TAG <- c("장출혈성대장균감염증" = "EHEC", "EHEC" = "EHEC", "세균성이질" = "SHIG", "장티푸스" = "TYPH", "노로바이러스" = "NORO", "살모넬라" = "SALM")
v2_main_cov <- function(dz) {   # 질환의 v2 주모형 공변량·패널(기간 동인 분석 411·415·451용). EHEC = JKMS 원고 20변수(ST_prereg_panel), 4개 질환 = 문헌 주모형 세트 + 강제통제(HM 패널)
  tag <- unname(V2_DZ_TAG[dz]); if (is.na(tag)) stop(sprintf("[v2] '%s'는 v2 주모형이 정의되지 않은 질환 — 기간 동인 분석은 EHEC·이질·장티푸스·노로·살모넬라만", dz), call. = FALSE)
  if (tag == "EHEC") { J <- v2_jkms_main(); return(list(tag = tag, vars = J$vars, panel = "ST_prereg_panel.csv", pc1 = NULL, lab = J$ENG)) }
  S <- jsonlite::fromJSON(file.path(V2_ROOT, "문헌가설_260930", tag, sprintf("sets_%s.json", tag)), simplifyVector = TRUE)
  M <- utils::read.csv(file.path(V2_CACHE, "HM_panel_META.csv"), fileEncoding = "UTF-8-BOM")
  list(tag = tag, vars = unique(c(unlist(S$sets[[S$rules$main_set]]), unlist(S$forced))), panel = sprintf("%s_panel.csv", tag), pc1 = sprintf("%s_kchs_pc1.csv", tag), lab = c(stats::setNames(M$label, M$key), kchs_pc1 = "Handwashing composite (KCHS PC1)"))
}
v2_build_Z_main <- function(dz, y0 = 2011, y1 = 2024, MC = v2_main_cov(dz)) {   # 옛 build_Z(ST_covariates_panel.rds)의 v2 대체: v2 패널에서 창(y0-y1) 안으로 표준화, 평균 채움 없음(결측이면 멈춤)
  d <- .v2_read(MC$panel); d$region <- gsub("\\s+", "", d$region)
  if (!is.null(MC$pc1)) { p1 <- utils::read.csv(file.path(V2_CACHE, MC$pc1), colClasses = c(sigungu_cd = "character")); d$sigungu_cd <- v2_code(d$region); d <- dplyr::left_join(d, p1[, intersect(c("sigungu_cd", "year", "kchs_pc1"), names(p1))], by = intersect(c("sigungu_cd", "year"), names(p1))) }
  b <- d[d$year >= y0 & d$year <= y1, ]
  Z <- data.frame(region = b$region, year = b$year); for (v in MC$vars) { if (is.null(b[[v]])) stop(sprintf("[v2] %s에 열 없음: %s", MC$panel, v), call. = FALSE); Z[[v]] <- v2_zb(b[[v]], v) }
  v2_relabel(Z, sprintf("%s 주모형 공변량(%d-%d)", MC$tag, y0, y1), keep_code = FALSE)
}
v2_key_check <- function(KEY, T2file, outdir, script, base = file.path(V2_ROOT, "_v2_promote_261021", "17_KEY_기준_MASTER_final_table2_260903.csv")) {
  # 542 계열 KEY(원고 고정 4변수)는 바꾸지 않음. 이번 540 결과(T2file)의 1 배제 판정을 원고 기준 Table 2(9/3, base)와 비교해 판정이 바뀐 변수와 KEY 변수를 기록·출력.
  if (!file.exists(T2file)) { message("[v2b KEY 점검] ", basename(T2file), " 없음 — 540을 먼저 돌릴 것"); return(invisible(NULL)) }
  T2 <- utils::read.csv(T2file); B <- if (file.exists(base)) utils::read.csv(base) else NULL
  cv <- union(KEY, T2$cov); R <- data.frame(script = script, cov = cv, in_KEY = cv %in% KEY, stringsAsFactors = FALSE)
  m2 <- match(cv, T2$cov); R$IRR_now <- T2$IRR[m2]; R$lo_now <- T2$lo[m2]; R$hi_now <- T2$hi[m2]; R$cred_now <- as.logical(T2$cred[m2])
  if (!is.null(B)) { mb <- match(cv, B$cov); R$IRR_manuscript <- B$IRR[mb]; R$cred_manuscript <- as.logical(B$cred[mb]) } else { R$IRR_manuscript <- NA; R$cred_manuscript <- NA }
  R$flag <- ifelse(is.na(R$cred_now), "이번 540 Table 2에 없음",
            ifelse(!is.na(R$cred_manuscript) & R$cred_now != R$cred_manuscript, ifelse(R$cred_now, "판정 바뀜: 원고 1 배제 아님 → 이번 1 배제", "판정 바뀜: 원고 1 배제 → 이번 1 배제 아님"), ""))
  R <- R[R$in_KEY | nzchar(R$flag), ]; R$t2_file <- T2file; R$t2_mtime <- format(file.mtime(T2file))
  f <- file.path(outdir, "MASTER_final_KEY_check.csv"); old <- if (file.exists(f)) utils::read.csv(f) else NULL
  if (!is.null(old)) old <- old[old$script != script, , drop = FALSE]; utils::write.csv(rbind(old, R), f, row.names = FALSE)
  bad <- R[nzchar(R$flag), ]
  if (nrow(bad)) message(sprintf("[v2b KEY 점검 %s] ⚠ 원고 대비 1 배제 판정이 바뀐 변수 %d개 — %s. KEY(%s)는 원고 고정 목록이라 그대로 씀; 보고서·원고에 표시할 것 → %s", script, nrow(bad), paste(sprintf("%s%s(%s)", bad$cov, ifelse(bad$in_KEY, "[KEY]", ""), bad$flag), collapse = "; "), paste(KEY, collapse = ", "), f))
  else message(sprintf("[v2b KEY 점검 %s] 원고 대비 1 배제 판정이 바뀐 변수 없음(KEY %s 그대로)", script, paste(KEY, collapse = ", ")))
  invisible(R)
}
