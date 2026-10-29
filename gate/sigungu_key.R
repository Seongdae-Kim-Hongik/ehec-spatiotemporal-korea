# sigungu_key.R — 시군구 코드 결합 공용 도구 (2026-10-21). sigungu_key.py와 같은 동작.
#  이름으로 결합하지 말고 반드시 sigungu_cd(5자리, 229 체계)로 바꾼 뒤 결합한다.
#   to_code(names, year=NULL)               이름(표기 아무거나) -> sigungu_cd. 못 찾으면 stop()(조용한 NA 없음)
#   attach_code(df, name_col, year_col)     데이터프레임에 sigungu_cd 열을 붙임. 실패 시 멈춤
#   read_shp_checked(path)                  final.shp를 읽어 코드 기준으로 canonical_region을 다시 붙이고 코드·이름·면적(±5%) 단언
#   name_kind(names)                        표기 종류(시군구/일반구/통합전/시도/전국/모호/미확인/검토필요)
#   explicit_fill(df, col, how, reason)     결측 채움이 꼭 필요할 때만(채운 칸 기록)
#  사용: source(".../00_분석전_데이터게이트/sigungu_key.R")  — 표는 같은 폴더의 sigungu_key.csv, sigungu_alias.csv
local({
  f <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  assign(".SGK_DIR", if (!is.na(f)) dirname(f) else
    Sys.getenv("EHEC_GATE", unset=file.path(Sys.getenv("EHEC_ROOT", unset="."),"gate")), envir = globalenv())
})
.SGK_SHP <- file.path(Sys.getenv("EHEC_DATA", unset="data"),"FBD_DATA_ZIP","final.shp")
.SGK_SIDO <- list(
  "서울" = c("서울특별시", "서울시", "서울"), "부산" = c("부산광역시", "부산시", "부산"),
  "대구" = c("대구광역시", "대구시", "대구"), "인천" = c("인천광역시", "인천시", "인천"),
  "광주" = c("광주광역시", "광주시", "광주"), "대전" = c("대전광역시", "대전시", "대전"),
  "울산" = c("울산광역시", "울산시", "울산"), "세종" = c("세종특별자치시", "세종시", "세종"),
  "경기" = c("경기도", "경기"), "강원" = c("강원특별자치도", "강원도", "강원"),
  "충북" = c("충청북도", "충북"), "충남" = c("충청남도", "충남"),
  "전북" = c("전북특별자치도", "전라북도", "전북"), "전남" = c("전라남도", "전남"),
  "경북" = c("경상북도", "경북"), "경남" = c("경상남도", "경남"),
  "제주" = c("제주특별자치도", "제주도", "제주"))
.SGK_FORMS <- local({ f <- unlist(.SGK_SIDO); s <- rep(names(.SGK_SIDO), lengths(.SGK_SIDO)); o <- order(-nchar(f)); data.frame(f = unname(f[o]), s = s[o], stringsAsFactors = FALSE) })
.SGK_NATIONAL <- c("전국", "합계", "계", "총계", "전체", "소계", "national", "total")
.SGK_ENV <- new.env()

sgk_norm <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  x <- stringi::stri_trans_nfc(x)
  x <- gsub("\\(\\s*[0-9]+\\s*\\)", "", x, perl = TRUE)
  gsub("[\\s　_-]", "", x, perl = TRUE)
}

sgk_split_sido <- function(n) {   # -> data.frame(s, r): 시도약칭·나머지 후보(겹쳐 적힌 시도는 벗긴 후보도)
  S <- character(0); R <- character(0)
  for (i in seq_len(nrow(.SGK_FORMS))) {
    f <- .SGK_FORMS$f[i]; s <- .SGK_FORMS$s[i]
    if (startsWith(n, f) && nchar(n) > nchar(f)) {
      rest <- substring(n, nchar(f) + 1); S <- c(S, s); R <- c(R, rest)
      repeat {
        ch <- FALSE
        for (f2 in .SGK_SIDO[[s]]) if (startsWith(rest, f2) && nchar(rest) > nchar(f2)) { rest <- substring(rest, nchar(f2) + 1); S <- c(S, s); R <- c(R, rest); ch <- TRUE; break }
        if (!ch) break
      }
    }
  }
  data.frame(s = S, r = R, stringsAsFactors = FALSE)
}
.sgk_strip <- function(r) sub("(시|군|구)$", "", r)
.sgk_yr <- function(d) { d <- as.character(d); ifelse(grepl("^[0-9]{4}", d), as.integer(substr(d, 1, 4)), NA_integer_) }

sgk_load <- function(key_csv = file.path(.SGK_DIR, "sigungu_key.csv"), alias_csv = file.path(.SGK_DIR, "sigungu_alias.csv")) {
  id <- paste(key_csv, alias_csv)
  if (!is.null(.SGK_ENV$T) && identical(.SGK_ENV$id, id)) return(.SGK_ENV$T)
  key <- utils::read.csv(key_csv, colClasses = "character", fileEncoding = "UTF-8-BOM", check.names = FALSE)
  al <- utils::read.csv(alias_csv, colClasses = "character", fileEncoding = "UTF-8-BOM", check.names = FALSE); al[is.na(al)] <- ""
  stopifnot(nrow(key) == 229, !anyDuplicated(key$sigungu_cd))
  ok <- al$needs_review != "1"
  bad <- setdiff(al$sigungu_cd[ok & al$sigungu_cd != ""], key$sigungu_cd); if (length(bad)) stop("alias 표에 229에 없는 코드: ", paste(bad, collapse = ","))
  rv <- al[!ok, ]; al <- al[ok, ]
  al$i <- seq_len(nrow(al)); al$n <- sgk_norm(al$alias); al$g <- sgk_norm(al$sigungu_part)
  al$vf <- .sgk_yr(al$valid_from); al$vt <- .sgk_yr(al$valid_to)
  grp <- function(k, v) split(v, k)
  hasg <- al$g != ""; hass <- hasg & al$sido_part != ""
  T <- list(key = key, al = al,
            exact = grp(al$n, al$i),
            pair = grp(paste(al$sido_part, al$g, sep = "|")[hass], al$i[hass]),
            short = grp(paste(al$sido_part, .sgk_strip(al$g), sep = "|")[hass], al$i[hass]),
            bare = grp(al$g[hasg], al$i[hasg]),
            review = setNames(rv$source_note, sgk_norm(rv$alias)),
            name_of = setNames(key$canonical_region, key$sigungu_cd))
  .SGK_ENV$T <- T; .SGK_ENV$id <- id; T
}

.sgk_resolve <- function(name, T, sido_hint = NULL) {   # -> list(rows=alias 행 번호, how, msg)
  n <- sgk_norm(name)
  if (n == "") return(list(rows = integer(0), how = "빈값", msg = "빈 이름"))
  if (!is.null(sido_hint)) {   # 시도를 따로 아는 경우(shp sigungu_nm 점검 등): (시도, 이름) 조합만 봄
    cd <- sgk_split_sido(n); cd <- rbind(data.frame(s = sido_hint, r = n, stringsAsFactors = FALSE), cd[cd$s == sido_hint, , drop = FALSE])
    rows <- unlist(lapply(seq_len(nrow(cd)), function(j) T$pair[[paste(cd$s[j], cd$r[j], sep = "|")]]))
    return(if (length(rows)) list(rows = rows, how = "시도힌트", msg = "") else list(rows = integer(0), how = "없음", msg = paste0(sido_hint, "|", n, " 조합이 별칭 표에 없음")))
  }
  if (!is.na(T$review[n])) return(list(rows = integer(0), how = "검토필요", msg = paste0("검토 필요 별칭: ", T$review[[n]])))
  if (!is.null(T$exact[[n]])) return(list(rows = T$exact[[n]], how = "정확", msg = ""))
  if (n %in% .SGK_FORMS$f) return(list(rows = integer(0), how = "시도", msg = paste0("시도 이름(시도 합계 행으로 보임)", if (n == "광주시") ". 경기 광주시라면 '경기도 광주시'처럼 시도를 붙일 것" else "")))
  if (tolower(n) %in% .SGK_NATIONAL) return(list(rows = integer(0), how = "전국", msg = "전국/합계 행"))
  cd <- sgk_split_sido(n)
  if (nrow(cd) && any(mapply(function(s, r) r %in% .SGK_SIDO[[s]], cd$s, cd$r))) return(list(rows = integer(0), how = "시도", msg = "시도 이름이 겹쳐 적힌 시도 합계 행"))
  rows <- integer(0)
  for (j in seq_len(nrow(cd))) rows <- c(rows, T$pair[[paste(cd$s[j], cd$r[j], sep = "|")]])
  if (length(rows)) return(list(rows = rows, how = "시도+시군구", msg = ""))
  for (j in seq_len(nrow(cd))) {
    rows <- c(rows, T$short[[paste(cd$s[j], cd$r[j], sep = "|")]])
    if (!length(rows) && .sgk_strip(cd$r[j]) != "") rows <- c(rows, T$short[[paste(cd$s[j], .sgk_strip(cd$r[j]), sep = "|")]])
  }
  if (length(rows)) return(list(rows = rows, how = "시도+시군구(시·군·구 생략)", msg = ""))
  if (!is.null(T$bare[[n]])) return(list(rows = T$bare[[n]], how = "시군구만", msg = ""))
  if (nrow(cd)) return(list(rows = integer(0), how = "없음", msg = paste0("시도-시군구 조합이 별칭 표에 없음: ", paste(paste(cd$s, cd$r, sep = "|"), collapse = ", "))))
  list(rows = integer(0), how = "없음", msg = "별칭 표에 없음")
}

.sgk_pick <- function(rows, year, subdistrict, T) {
  a <- T$al[rows, ]; codes <- sort(unique(a$sigungu_cd))
  if (length(codes) > 1) return(list(code = NA_character_, msg = paste0("모호: ", paste(codes, collapse = ", "))))
  if (all(a$kind == "일반구") && subdistrict != "city")
    return(list(code = NA_character_, msg = paste0("일반구·하위구역 표기(시 ", codes, "의 일부). 시 합계 행과 중복 합산될 수 있어 막음 — 의도했다면 subdistrict='city'")))
  msg <- ""
  if (!is.null(year) && !is.na(year)) {
    y <- as.integer(year); ok <- any((is.na(a$vf) | y >= a$vf) & (is.na(a$vt) | y <= a$vt))
    if (!ok) msg <- paste0("시기불일치: ", paste0(a$alias, "(", ifelse(is.na(a$vf), "", a$vf), "~", ifelse(is.na(a$vt), "", a$vt), ")", collapse = "; "))
  }
  list(code = codes, msg = msg)
}

to_code <- function(names, year = NULL, subdistrict = c("error", "city"), strict_period = FALSE, detail = FALSE, ...) {
  subdistrict <- match.arg(subdistrict); T <- sgk_load(...)
  names <- as.character(names)
  yrs <- if (is.null(year)) rep(NA_integer_, length(names)) else if (length(year) == 1) rep(as.integer(year), length(names)) else as.integer(year)
  if (length(yrs) != length(names)) stop("year 길이가 names와 다름")
  k <- paste(names, yrs, sep = "\r"); u <- !duplicated(k); res <- vector("list", sum(u)); ku <- k[u]; nu <- names[u]; yu <- yrs[u]
  for (i in seq_along(ku)) {
    r <- .sgk_resolve(nu[i], T)
    res[[i]] <- if (length(r$rows)) { p <- .sgk_pick(r$rows, yu[i], subdistrict, T); list(code = p$code, how = r$how, msg = p$msg) } else list(code = NA_character_, how = r$how, msg = r$msg)
  }
  code_u <- vapply(res, `[[`, "", "code"); msg_u <- vapply(res, `[[`, "", "msg"); how_u <- vapply(res, `[[`, "", "how")
  if (detail) return(data.frame(name = nu, year = yu, sigungu_cd = code_u, how = how_u, msg = msg_u, stringsAsFactors = FALSE))
  f <- is.na(code_u)
  if (any(f)) {
    tb <- unique(data.frame(n = nu[f], m = msg_u[f], stringsAsFactors = FALSE))
    ys <- tapply(yu[f], paste(nu[f], msg_u[f]), function(v) paste(head(sort(unique(v[!is.na(v)])), 6), collapse = ","))
    lines <- sprintf("  - '%s' : %s%s", tb$n, tb$m, ifelse(ys[paste(tb$n, tb$m)] != "", paste0(" (연도 ", ys[paste(tb$n, tb$m)], ")"), ""))
    stop(sprintf("시군구 코드화 실패 %d개 이름 — 결합을 멈춤(조용한 NA 금지).\n%s\n별칭이 맞으면 sigungu_alias.csv에 한 줄 추가(출처 기록) 후 다시 실행.",
                 nrow(tb), paste(head(lines, 200), collapse = "\n")), call. = FALSE)
  }
  pm <- startsWith(msg_u, "시기불일치")
  if (any(pm)) {
    m <- paste0("별칭 유효기간 밖 연도(같은 코드로 붙였음; 원자료의 소급 표기인지 확인):\n", paste(sprintf("  - '%s' 연도 %s : %s", nu[pm], yu[pm], msg_u[pm]), collapse = "\n"))
    if (strict_period) stop(m, call. = FALSE) else warning(m, call. = FALSE)
  }
  unname(code_u[match(k, ku)])
}

name_kind <- function(names, ...) {
  T <- sgk_load(...); nm <- as.character(names); u <- unique(nm)
  kk <- vapply(u, function(x) {
    r <- .sgk_resolve(x, T)
    if (!length(r$rows)) return(switch(r$how, "시도" = "시도", "전국" = "전국", "검토필요" = "검토필요", "미확인"))
    a <- T$al[r$rows, ]
    if (length(unique(a$sigungu_cd)) > 1) "모호" else if (all(a$kind == "일반구")) "일반구" else if (all(a$kind == "통합전")) "통합전" else "시군구"
  }, "")
  unname(kk[match(nm, u)])
}

attach_code <- function(df, name_col, year_col = NULL, out_col = "sigungu_cd", subdistrict = c("error", "city"), drop_kinds = NULL, strict_period = FALSE, ...) {
  subdistrict <- match.arg(subdistrict); d <- df
  if (length(drop_kinds)) {
    k <- name_kind(d[[name_col]], ...); dr <- k %in% drop_kinds
    if (any(dr)) message(sprintf("[attach_code] %s: %d행 제외(종류 %s), 예: %s", name_col, sum(dr), paste(sort(unique(k[dr])), collapse = ","), paste(head(unique(d[[name_col]][dr]), 5), collapse = " / ")))
    d <- d[!dr, , drop = FALSE]
  }
  d[[out_col]] <- to_code(d[[name_col]], year = if (!is.null(year_col)) d[[year_col]] else NULL, subdistrict = subdistrict, strict_period = strict_period, ...)
  T <- sgk_load(...); u <- unique(d[[name_col]])
  mg <- u[vapply(u, function(x) { r <- .sgk_resolve(x, T); length(r$rows) > 0 && any(T$al$kind[r$rows] == "통합전") }, TRUE)]
  if (length(mg)) message(sprintf("[attach_code] 통합 전 단위 표기 %d개가 통합 코드로 붙음 — (코드, 연도)별 합산 필요: %s", length(mg), paste(head(mg, 8), collapse = " / ")))
  d
}

explicit_fill <- function(df, col, how = c("zero", "interp"), reason, code_col = "sigungu_cd", year_col = "year", log_path = NULL) {
  how <- match.arg(how); if (missing(reason) || !nzchar(reason)) stop("reason(채우는 이유)을 반드시 적을 것")
  d <- df[order(df[[code_col]], df[[year_col]]), , drop = FALSE]; miss <- is.na(d[[col]])
  if (how == "zero") d[[col]][miss] <- 0 else {
    d[[col]] <- ave(d[[col]], d[[code_col]], FUN = function(v) { ok <- !is.na(v); if (sum(ok) == 0) v else if (sum(ok) == 1) rep(v[ok], length(v)) else stats::approx(which(ok), v[ok], xout = seq_along(v), rule = 2)$y })
  }
  filled <- data.frame(d[miss, c(code_col, year_col)], column = col, value = d[[col]][miss], how = how, reason = reason)
  message(sprintf("[explicit_fill] %s: %d칸 채움(%s) — %s", col, nrow(filled), how, reason)); if (nrow(filled)) print(utils::head(filled, 30))
  if (!is.null(log_path) && nrow(filled)) utils::write.table(filled, log_path, sep = ",", append = file.exists(log_path), col.names = !file.exists(log_path), row.names = FALSE, fileEncoding = "UTF-8")
  if (anyNA(d[[col]])) stop(sprintf("%s: 채운 뒤에도 %d칸 결측(지역 전체가 빈 경우). 원자료에서 확인할 것", col, sum(is.na(d[[col]]))))
  d
}

read_shp_checked <- function(path = .SGK_SHP, crs = 5179, area_tol = 0.05, ...) {
  T <- sgk_load(...); key <- T$key
  g <- sf::st_read(path, quiet = TRUE)
  if (is.na(sf::st_crs(g))) g <- sf::st_set_crs(g, crs) else if (!identical(sf::st_crs(g)$epsg, as.integer(crs))) g <- sf::st_transform(g, crs)
  if (!"sigungu_cd" %in% names(g)) stop("shp에 sigungu_cd 열 없음")
  g$sigungu_cd <- formatC(trimws(as.character(g$sigungu_cd)), width = 5, flag = "0")
  errs <- character(0)
  if (nrow(g) != 229) errs <- c(errs, sprintf("행 수 %d != 229", nrow(g)))
  if (anyDuplicated(g$sigungu_cd)) errs <- c(errs, paste("중복 코드", paste(g$sigungu_cd[duplicated(g$sigungu_cd)], collapse = ",")))
  ex <- setdiff(g$sigungu_cd, key$sigungu_cd); ms <- setdiff(key$sigungu_cd, g$sigungu_cd)
  if (length(ex)) errs <- c(errs, paste("229에 없는 코드", paste(ex, collapse = ",")))
  if (length(ms)) errs <- c(errs, paste("빠진 코드", paste(ms, collapse = ",")))
  ki <- match(g$sigungu_cd, key$sigungu_cd)
  if ("sigungu_nm" %in% names(g)) for (i in which(!is.na(ki))) {
    c0 <- g$sigungu_cd[i]; nm <- g$sigungu_nm[i]
    r1 <- .sgk_resolve(nm, T, sido_hint = key$sido[ki[i]]); r2 <- .sgk_resolve(nm, T)
    if (!(c0 %in% T$al$sigungu_cd[r1$rows]) && !(c0 %in% T$al$sigungu_cd[r2$rows]))
      errs <- c(errs, sprintf("이름 불일치 %s: shp '%s' vs key '%s'", c0, nm, key$canonical_region[ki[i]]))
  }
  a <- as.numeric(sf::st_area(g)) / 1e6; ref <- as.numeric(key$area_km2[ki])
  off <- which(!is.na(ref) & abs(a / ref - 1) > area_tol)
  if (length(off)) errs <- c(errs, sprintf("면적 불일치 %s %s: shp %.1f vs key %.1f km²", g$sigungu_cd[off], key$canonical_region[ki[off]], a[off], ref[off]))
  if (length(errs)) stop("read_shp_checked 실패:\n  ", paste(head(errs, 60), collapse = "\n  "), call. = FALSE)
  g$canonical_region <- key$canonical_region[ki]
  if ("region" %in% names(g)) {
    bad <- which(sgk_norm(g$region) != g$canonical_region)
    if (length(bad)) warning("shp region 라벨이 코드와 다른 행(코드 기준으로 교정함):\n", paste(sprintf("  %s: '%s' -> '%s'", g$sigungu_cd[bad], g$region[bad], g$canonical_region[bad]), collapse = "\n"), call. = FALSE)
    names(g)[names(g) == "region"] <- "region_shp"
  }
  g$region <- g$canonical_region; g$area_km2 <- a
  g
}
