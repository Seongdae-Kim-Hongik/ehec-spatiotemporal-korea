# panel_gate.R — 패널 결합 게이트 (2026-10-21). panel_gate.py와 같은 동작. 결합 뒤, 모형 적합 전에 반드시 통과.
#  실패(stop): 229에 없는 코드·코드화 실패 / (시군구, 연도) 중복 / 지정 열이 229 × 지정 연도 모든 칸에 있고 유한하지 않음, positive 열(인구·E) <= 0
#  경고(warning, strict=TRUE면 stop): 전국 최솟값·최댓값과 같은 값의 연속·반복(채움 의심), 양수 열 같은 값 연속(rule=2 외삽·상수 채움),
#        연도 간 1.3배 초과 변화(positive 열, 그해 전국 중앙 배율 대비; 통합 연도 전후는 ±10%), (선택) area_col이 key 면적과 ±5% 넘게 다름
#  사용: source(".../00_분석전_데이터게이트/panel_gate.R")
#        panel_gate(DD, cols = c("pop", "E_std", "cases"), positive = c("pop", "E_std"), years = 2011:2024, name_col = "region")
local({
  f <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  d <- if (!is.na(f)) dirname(f) else Sys.getenv("EHEC_GATE", unset=file.path(Sys.getenv("EHEC_ROOT", unset="."),"gate"))
  if (!exists("to_code", mode = "function")) source(file.path(d, "sigungu_key.R"), local = globalenv())
})

.pg_runs <- function(ys, flags) {   # 연속 구간 data.frame(a, b, n)
  out <- list(); cur <- NULL
  for (i in seq_along(ys)) {
    if (isTRUE(flags[i]) && !is.null(cur) && ys[i] == cur[2] + 1) cur[2] <- ys[i]
    else if (isTRUE(flags[i])) { if (!is.null(cur)) out[[length(out) + 1]] <- cur; cur <- c(ys[i], ys[i]) }
    else { if (!is.null(cur)) out[[length(out) + 1]] <- cur; cur <- NULL }
  }
  if (!is.null(cur)) out[[length(out) + 1]] <- cur
  if (!length(out)) return(data.frame(a = integer(0), b = integer(0), n = integer(0)))
  m <- do.call(rbind, out); data.frame(a = m[, 1], b = m[, 2], n = m[, 2] - m[, 1] + 1)
}
.pg_fmt <- function(s) { s <- sort(unique(as.integer(s))); r <- s[1]:s[length(s)]; rr <- .pg_runs(r, r %in% s); paste(ifelse(rr$a == rr$b, rr$a, paste0(rr$a, "-", rr$b)), collapse = ", ") }

panel_gate <- function(df, cols, positive = character(0), years = NULL, code_col = NULL, name_col = NULL, year_col = "year",
                       jump = 1.3, run_min = 2, const_run = 3, area_col = NULL, jump_cols = NULL, merge_jump = 1.1, strict = FALSE, label = "", verbose = TRUE, max_show = 60, stop_on_fail = TRUE) {
  df <- as.data.frame(df)
  for (c in unique(c(cols, positive))) if (!c %in% names(df)) stop(sprintf("[%s] 열 없음: %s", label, c))
  T <- sgk_load(); key <- T$key; name_of <- T$name_of; codes229 <- sort(key$sigungu_cd)
  d <- df; fails <- list(); warns <- list()
  if (is.null(code_col)) {
    if (is.null(name_col)) stop("code_col 또는 name_col 하나는 필요")
    d <- tryCatch(attach_code(d, name_col), error = function(e) stop(sprintf("[%s] 실패: 이름을 코드로 못 바꿈\n%s", label, conditionMessage(e)), call. = FALSE))
    code_col <- "sigungu_cd"
  }
  d[[code_col]] <- formatC(sub("\\.0$", "", trimws(as.character(d[[code_col]]))), width = 5, flag = "0")
  d[[year_col]] <- suppressWarnings(as.integer(d[[year_col]]))
  if (anyNA(d[[year_col]])) fails[["연도 결측"]] <- head(d[is.na(d[[year_col]]), , drop = FALSE], max_show)
  yrs <- sort(unique(if (is.null(years)) d[[year_col]][!is.na(d[[year_col]])] else as.integer(years)))
  d <- d[d[[year_col]] %in% yrs, , drop = FALSE]
  bad <- setdiff(unique(d[[code_col]]), codes229)
  if (length(bad)) fails[["229에 없는 코드"]] <- data.frame(sigungu_cd = bad)
  kk <- paste(d[[code_col]], d[[year_col]])
  if (anyDuplicated(kk)) {
    dk <- kk %in% kk[duplicated(kk)]; t <- aggregate(list(행수 = rep(1, sum(dk))), list(sigungu_cd = d[[code_col]][dk], year = d[[year_col]][dk]), sum)
    t$이름 <- name_of[t$sigungu_cd]
    if (!is.null(name_col)) t$원래이름들 <- mapply(function(c, y) paste(sort(unique(as.character(d[[name_col]][d[[code_col]] == c & d[[year_col]] == y]))), collapse = ", "), t$sigungu_cd, t$year)
    fails[["같은 (시군구, 연도) 중복"]] <- t
  }
  g <- expand.grid(code = codes229, year = yrs, stringsAsFactors = FALSE); names(g) <- c(code_col, year_col)
  dd <- d[!duplicated(kk), unique(c(code_col, year_col, cols, positive, area_col)), drop = FALSE]; dd$.row <- 1
  g <- merge(g, dd, by = c(code_col, year_col), all.x = TRUE, sort = TRUE)
  miss <- NULL
  for (c in unique(c(cols, positive))) {
    x <- suppressWarnings(as.numeric(g[[c]]))
    why <- ifelse(is.na(g$.row), "행 없음", ifelse(is.na(x), "결측", ifelse(!is.finite(x), "비유한", ifelse(c %in% positive & x <= 0, "0 이하", ""))))
    if (any(why != "")) miss <- rbind(miss, data.frame(sigungu_cd = g[[code_col]][why != ""], year = g[[year_col]][why != ""], 열 = c, 이유 = why[why != ""], stringsAsFactors = FALSE))
  }
  if (!is.null(miss)) {
    s <- aggregate(year ~ sigungu_cd + 열 + 이유, miss, .pg_fmt); s$이름 <- name_of[s$sigungu_cd]
    fails[[sprintf("빠진 칸 %d개(229 × %d-%d)", nrow(miss), min(yrs), max(yrs))]] <- s[order(s$sigungu_cd), c("sigungu_cd", "이름", "열", "이유", "year")]
  }
  fl <- NULL
  for (c in unique(c(cols, positive))) {
    x <- suppressWarnings(as.numeric(g[[c]])); ok <- is.finite(x); pos <- c %in% positive
    if (sum(ok) < 3) next
    xs <- x[ok & (if (pos) x > 0 else TRUE)]; if (!length(xs)) next
    for (tg in c("전국 최솟값", "전국 최댓값")) {
      v <- if (tg == "전국 최솟값") min(xs) else max(xs)
      if (v == 0 && !pos) next
      hit <- ok & abs(x - v) <= 1e-9 * abs(v)
      if (sum(hit) < 2) next
      regs <- unique(g[[code_col]][hit])
      for (r in regs) {
        sub <- g[g[[code_col]] == r, ]; sub <- sub[order(sub[[year_col]]), ]; xv <- suppressWarnings(as.numeric(sub[[c]]))
        rr <- .pg_runs(sub[[year_col]], !is.na(xv) & abs(xv - v) <= 1e-9 * abs(v)); rr <- rr[rr$n >= run_min, ]
        if (nrow(rr)) fl <- rbind(fl, data.frame(sigungu_cd = r, 이름 = name_of[[r]], 열 = c, 패턴 = sprintf("%s(%.6g)와 같은 값", tg, v), 연도 = paste0(rr$a, "-", rr$b), 칸수 = rr$n))
      }
      if (length(regs) >= 2) fl <- rbind(fl, data.frame(sigungu_cd = "여러", 이름 = sprintf("%d개 지역", length(regs)), 열 = c, 패턴 = sprintf("%s(%.6g)가 %d칸에 반복", tg, v, sum(hit)), 연도 = "", 칸수 = sum(hit)))
    }
    if (pos) for (r in codes229) {
      sub <- g[g[[code_col]] == r, ]; sub <- sub[order(sub[[year_col]]), ]; v <- suppressWarnings(as.numeric(sub[[c]])); ys <- sub[[year_col]]
      same <- c(FALSE, vapply(seq_along(v)[-1], function(i) is.finite(v[i]) && is.finite(v[i - 1]) && abs(v[i] - v[i - 1]) <= 1e-12 * abs(v[i]), TRUE))
      rr <- .pg_runs(ys, same); rr <- rr[rr$n + 1 >= const_run, ]
      if (nrow(rr)) fl <- rbind(fl, data.frame(sigungu_cd = r, 이름 = name_of[[r]], 열 = c, 패턴 = sprintf("같은 값 %.6g이 연속(상수 채움·rule=2 외삽 의심)", v[match(rr$b, ys)]), 연도 = paste0(rr$a - 1, "-", rr$b), 칸수 = rr$n + 1))
    }
  }
  if (!is.null(fl)) warns[["채움 의심"]] <- fl
  jl <- NULL   # 지역 배율 / 그해 전국 중앙 배율 (전국 발생률을 따라 움직이는 E도 지역 고유 급변만 걸림)
  for (c in (if (is.null(jump_cols)) positive else jump_cols)) {
    w <- tapply(suppressWarnings(as.numeric(g[[c]])), list(g[[code_col]], g[[year_col]]), function(v) v[1])
    for (k in seq_along(yrs)[-1]) {
      y0 <- yrs[k - 1]; y1 <- yrs[k]; if (y1 != y0 + 1) next
      a0 <- w[, as.character(y0)]; a1 <- w[, as.character(y1)]; ok <- is.finite(a0) & is.finite(a1) & a0 > 0 & a1 > 0
      if (sum(ok) < 10) next
      rt <- a1[ok] / a0[ok]; med <- stats::median(rt); rel <- rt / med; hit <- names(rt)[rel > jump | rel < 1 / jump]
      if (length(hit)) jl <- rbind(jl, data.frame(sigungu_cd = hit, 이름 = unname(name_of[hit]), 열 = c, 연도 = paste0(y0, "→", y1), 앞 = unname(a0[hit]), 뒤 = unname(a1[hit]), 배율 = round(unname(rt[hit]), 3), 전국중앙배율 = round(med, 3)))
    }
  }
  # 통합 연도 단절: 통합 전 이름의 마지막 해 직전→그해 상대 배율이 ±10% 넘으면 경고(통합 전 단위 누락·이중 합산 의심)
  mg <- T$al[T$al$kind == "통합전" & T$al$valid_to != "", ]; tr <- unique(data.frame(code = mg$sigungu_cd, y1 = as.integer(substr(mg$valid_to, 1, 4))))
  tl <- NULL
  for (c in (if (is.null(jump_cols)) positive else jump_cols)) {
    w <- tapply(suppressWarnings(as.numeric(g[[c]])), list(g[[code_col]], g[[year_col]]), function(v) v[1])
    for (k in seq_len(nrow(tr))) {
      y1 <- tr$y1[k]; y0 <- y1 - 1; cd <- tr$code[k]
      if (!(as.character(y0) %in% colnames(w)) || !(as.character(y1) %in% colnames(w)) || !(cd %in% rownames(w))) next
      a0 <- w[, as.character(y0)]; a1 <- w[, as.character(y1)]; ok <- is.finite(a0) & is.finite(a1) & a0 > 0 & a1 > 0
      if (sum(ok) < 10 || !isTRUE(ok[cd])) next
      med <- stats::median(a1[ok] / a0[ok]); rel <- (a1[cd] / a0[cd]) / med
      if (rel > merge_jump || rel < 1 / merge_jump) tl <- rbind(tl, data.frame(sigungu_cd = cd, 이름 = name_of[[cd]], 열 = c, 연도 = paste0(y0, "→", y1), 앞 = unname(a0[cd]), 뒤 = unname(a1[cd]), 배율 = round(unname(a1[cd] / a0[cd]), 3), 전국중앙배율 = round(med, 3)))
    }
  }
  if (!is.null(tl)) warns[[sprintf("통합 연도 전후 단절(±%d%% 넘음; 통합 전 단위 누락·이중 합산 의심)", round((merge_jump - 1) * 100))]] <- tl
  if (!is.null(jl)) warns[[sprintf("연도 간 %s배 초과 변화(전국 중앙 배율 대비; 실제 변화인지 확인)", jump)]] <- jl
  if (!is.null(area_col)) {
    a <- tapply(suppressWarnings(as.numeric(g[[area_col]])), g[[code_col]], median, na.rm = TRUE)
    ref <- as.numeric(key$area_km2[match(names(a), key$sigungu_cd)]); rat <- a / ref; off <- which(abs(rat - 1) > 0.05)
    if (length(off)) warns[[sprintf("%s이 shp 기하 면적과 ±5%% 넘게 다름(%d곳)", area_col, length(off))]] <- data.frame(sigungu_cd = names(a)[off], 이름 = name_of[names(a)[off]], area = a[off], key_area_km2 = ref[off], 배율 = round(rat[off], 3))
  }
  rep <- list(label = label, n_rows = nrow(df), years = range(yrs), cols = cols, positive = positive, fails = fails, warns = warns,
              passed = !length(fails) && !(strict && length(warns)))
  if (verbose) {
    cat(sprintf("== panel_gate [%s] %d행, 연도 %d-%d, 열 %s (양수 %s)\n결과: %s%s\n", label, nrow(df), min(yrs), max(yrs), paste(cols, collapse = ","), paste(positive, collapse = ","),
                if (length(fails)) "실패" else "통과", if (length(warns)) sprintf(" · 경고 %d종", length(warns)) else " · 경고 없음"))
    for (tg in c("실패", "경고")) { L <- if (tg == "실패") fails else warns
      for (k in names(L)) { cat(sprintf("-- %s: %s (%d줄)\n", tg, k, nrow(L[[k]]))); print(utils::head(L[[k]], max_show), row.names = FALSE) } }
  }
  if (length(warns) && !strict) warning(sprintf("[%s] panel_gate 경고 %d종: %s", label, length(warns), paste(names(warns), collapse = "; ")), call. = FALSE)
  if (stop_on_fail && (length(fails) || (strict && length(warns)))) stop(sprintf("[%s] panel_gate 실패: %s", label, paste(c(names(fails), if (strict) names(warns)), collapse = "; ")), call. = FALSE)
  invisible(rep)
}
