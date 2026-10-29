#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""시군구 코드 결합 공용 도구 (2026-10-21).

이름으로 결합하지 말고 반드시 sigungu_cd(5자리, 229 체계)로 바꾼 뒤 결합한다.
  to_code(names, year=None)              이름(표기 아무거나) -> sigungu_cd. 못 찾으면 예외로 멈춤(조용한 NA 없음).
  attach_code(df, name_col, year_col)    데이터프레임에 sigungu_cd 열을 붙임. 실패 시 멈춤.
  read_shp_checked(path)                 final.shp를 읽어 코드 기준으로 canonical_region을 다시 붙이고 코드·이름·면적(±5%) 단언.
  name_kind(names)                       표기 종류(시군구/일반구/통합전/시도/전국/미확인) 확인용.
  explicit_fill(df, ...)                 결측 채움이 꼭 필요할 때만 쓰는 명시적 채움(채운 칸을 기록으로 남김).

표: 같은 폴더의 sigungu_key.csv(229 코드), sigungu_alias.csv(별칭). 규칙은 README_시군구결합_규칙.md.
명령줄: python3 sigungu_key.py "인천광역시 미추홀구" "경상북도 군위군" [--year 2020]
        python3 sigungu_key.py --file panel.csv --col region [--year-col year]
"""
from __future__ import annotations
import os, re, sys, unicodedata, warnings
from collections import defaultdict
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
KEY_CSV = os.path.join(HERE, "sigungu_key.csv")
ALIAS_CSV = os.path.join(HERE, "sigungu_alias.csv")
DEFAULT_SHP = os.path.join(os.environ.get("EHEC_DATA","data"),"FBD_DATA_ZIP","final.shp")

# 시도 표기 -> 약칭. 긴 표기부터 맞춘다. R 판(sigungu_key.R)과 같은 표.
SIDO_FORMS = {
    "서울": ["서울특별시", "서울시", "서울"], "부산": ["부산광역시", "부산시", "부산"],
    "대구": ["대구광역시", "대구시", "대구"], "인천": ["인천광역시", "인천시", "인천"],
    "광주": ["광주광역시", "광주시", "광주"], "대전": ["대전광역시", "대전시", "대전"],
    "울산": ["울산광역시", "울산시", "울산"], "세종": ["세종특별자치시", "세종시", "세종"],
    "경기": ["경기도", "경기"], "강원": ["강원특별자치도", "강원도", "강원"],
    "충북": ["충청북도", "충북"], "충남": ["충청남도", "충남"],
    "전북": ["전북특별자치도", "전라북도", "전북"], "전남": ["전라남도", "전남"],
    "경북": ["경상북도", "경북"], "경남": ["경상남도", "경남"],
    "제주": ["제주특별자치도", "제주도", "제주"],
}
_FORMS = sorted([(f, s) for s, fs in SIDO_FORMS.items() for f in fs], key=lambda x: -len(x[0]))
NATIONAL = {"전국", "합계", "계", "총계", "전체", "소계", "national", "total"}


class SigunguKeyError(ValueError):
    """시군구 이름을 코드로 바꾸지 못했거나 표·지도 단언이 깨졌을 때."""


def norm(x) -> str:
    """NFC, 괄호 속 숫자 코드 제거, 공백·'-'·'_' 제거."""
    if x is None or (isinstance(x, float) and x != x):
        return ""
    s = unicodedata.normalize("NFC", str(x))
    s = re.sub(r"\(\s*\d+\s*\)", "", s)
    s = re.sub(r"[\s　\-_]", "", s)
    return s


def split_sido(n: str):
    """정규화된 문자열 -> [(시도약칭, 나머지)] 후보. 같은 시도가 두 번 적힌 경우('경기도경기여주군')는 벗긴 후보도 함께 냄.
    벗긴 후보와 안 벗긴 후보를 모두 내므로 '제주도제주시'는 (제주,'제주시')로도 맞음."""
    out = []
    for f, s in _FORMS:
        if n.startswith(f) and len(n) > len(f):
            rest = n[len(f):]
            out.append((s, rest))
            changed = True
            while changed:
                changed = False
                for f2 in SIDO_FORMS[s]:
                    if rest.startswith(f2) and len(rest) > len(f2):
                        rest = rest[len(f2):]; out.append((s, rest)); changed = True; break
    return out


def _strip_suffix(r: str) -> str:
    return re.sub(r"(시|군|구)$", "", r)


_CACHE = {}


def _num(v):
    v = pd.to_numeric(v, errors="coerce")
    return None if pd.isna(v) else int(v)


def _yr(d):
    d = "" if d is None or (isinstance(d, float) and d != d) else str(d).strip()
    return int(d[:4]) if re.match(r"^\d{4}", d) else None


def load_tables(key_csv: str = KEY_CSV, alias_csv: str = ALIAS_CSV):
    k = (key_csv, alias_csv)
    if k in _CACHE:
        return _CACHE[k]
    key = pd.read_csv(key_csv, dtype=str, encoding="utf-8-sig")
    al = pd.read_csv(alias_csv, dtype=str, encoding="utf-8-sig").fillna("")
    assert len(key) == 229 and key.sigungu_cd.is_unique, "sigungu_key.csv는 229개 고유 코드여야 함"
    codes = set(key.sigungu_cd)
    bad = set(al.loc[(al.needs_review != "1") & (al.sigungu_cd != ""), "sigungu_cd"]) - codes
    assert not bad, f"alias 표에 229에 없는 코드: {bad}"
    exact = defaultdict(list); pair = defaultdict(list); bare = defaultdict(list); short = defaultdict(list); review = {}
    for i, r in al.iterrows():
        n = norm(r.alias)
        if r.needs_review == "1":
            review[n] = r.source_note; continue
        row = (r.sigungu_cd, r.kind, _yr(r.valid_from), _yr(r.valid_to), r.alias)
        exact[n].append(row)
        sp, gp = r.sido_part, norm(r.sigungu_part)
        if gp:
            if sp:
                pair[(sp, gp)].append(row)
                short[(sp, _strip_suffix(gp))].append(row)
            bare[gp].append(row)
    T = dict(key=key, alias=al, codes=codes, exact=exact, pair=pair, bare=bare, short=short, review=review,
             name_of=dict(zip(key.sigungu_cd, key.canonical_region)))
    _CACHE[k] = T
    return T


def _resolve(name, T, sido_hint=None):
    """한 이름 -> (rows, 방법, 메시지). rows가 비면 실패."""
    n = norm(name)
    if not n:
        return [], "빈값", "빈 이름"
    if sido_hint is not None:               # 시도를 따로 아는 경우(shp sigungu_nm 점검 등): (시도, 이름) 조합만 봄
        cands = [(sido_hint, n)] + [c for c in split_sido(n) if c[0] == sido_hint]
        rows = [x for s_, r in cands for x in T["pair"].get((s_, r), [])]
        return (rows, "시도힌트", "") if rows else ([], "없음", f"{sido_hint}|{n} 조합이 별칭 표에 없음")
    if n in T["review"]:
        return [], "검토필요", "검토 필요 별칭: " + T["review"][n]
    if n in T["exact"]:
        return T["exact"][n], "정확", ""
    if n in {f for f, _ in _FORMS}:
        return [], "시도", "시도 이름(시도 합계 행으로 보임)" + (". 경기 광주시라면 '경기도 광주시'처럼 시도를 붙일 것" if n == "광주시" else "")
    if n.lower() in NATIONAL:
        return [], "전국", "전국/합계 행"
    cands = split_sido(n)
    if any(r in SIDO_FORMS[sd] for sd, r in cands):
        return [], "시도", "시도 이름이 겹쳐 적힌 시도 합계 행"
    rows = []
    for s, r in cands:                       # 1) 시도+시군구 전체 이름
        rows += T["pair"].get((s, r), [])
    if rows:
        return rows, "시도+시군구", ""
    for s, r in cands:                       # 2) 시도+시군구에서 끝의 시·군·구 생략(식약처형 '경남 마산', '강원 양구')
        rows += T["short"].get((s, r), [])
        if not rows and _strip_suffix(r):
            rows += T["short"].get((s, _strip_suffix(r)), [])
    if rows:
        return rows, "시도+시군구(시·군·구 생략)", ""
    if n in T["bare"]:  # 3) 시도 없는 시군구 이름(같은 이름이 여러 시도에 있으면 모호로 실패)
        return T["bare"][n], "시군구만", ""
    if cands:
        return [], "없음", f"시도-시군구 조합이 별칭 표에 없음: {[f'{s}|{r}' for s, r in cands]}"
    return [], "없음", "별칭 표에 없음"


def _pick(rows, year, subdistrict):
    codes = sorted({r[0] for r in rows})
    if len(codes) > 1:
        return None, "모호: " + ", ".join(codes)
    kinds = {r[1] for r in rows}
    if kinds == {"일반구"} and subdistrict != "city":
        return None, f"일반구·하위구역 표기(시 {codes[0]}의 일부). 시 합계 행과 중복 합산될 수 있어 막음 — 의도했다면 subdistrict='city'"
    msg = ""
    if year is not None and not (isinstance(year, float) and year != year):
        y = int(year)
        ok = any((r[2] is None or y >= r[2]) and (r[3] is None or y <= r[3]) for r in rows)
        if not ok:
            msg = "시기불일치: " + "; ".join(f"{r[4]}({r[2] or ''}~{r[3] or ''})" for r in rows)
    return codes[0], msg


def to_code(names, year=None, subdistrict: str = "error", strict_period: bool = False,
            key_csv: str = KEY_CSV, alias_csv: str = ALIAS_CSV, _detail: bool = False):
    """이름(들) -> sigungu_cd(5자리 문자열) 목록. 하나라도 못 찾으면 SigunguKeyError.

    year: 스칼라 또는 names와 같은 길이. 별칭의 유효기간 밖이면 경고(strict_period=True면 예외).
    subdistrict: 'error'(기본, 일반구·출장소·보건소구역 표기를 막음) | 'city'(소속 시 코드로 붙임; 합산은 직접).
    """
    T = load_tables(key_csv, alias_csv)
    scalar = isinstance(names, (str, bytes)) or not hasattr(names, "__iter__")
    names = [names] if scalar else list(names)
    years = [year] * len(names) if (year is None or not hasattr(year, "__iter__") or isinstance(year, str)) else list(year)
    if len(years) != len(names):
        raise ValueError("year 길이가 names와 다름")
    memo = {}; out = []; fails = defaultdict(set); period = defaultdict(set); detail = []
    for nm, y in zip(names, years):
        k = (nm, y)
        if k not in memo:
            rows, how, m = _resolve(nm, T)
            if rows:
                c, m2 = _pick(rows, y, subdistrict)
                memo[k] = (c, how, m2 if c else m2)
            else:
                memo[k] = (None, how, m)
        c, how, m = memo[k]
        out.append(c); detail.append((nm, y, c, how, m))
        if c is None:
            fails[(str(nm), m)].add(y)
        elif m.startswith("시기불일치"):
            period[(str(nm), m)].add(y)
    if fails:
        lines = [f"  - '{n}' : {m}" + (f" (연도 {sorted(v for v in ys if v is not None)[:6]})" if any(v is not None for v in ys) else "")
                 for (n, m), ys in sorted(fails.items())]
        raise SigunguKeyError(f"시군구 코드화 실패 {len(fails)}개 이름 — 결합을 멈춤(조용한 NA 금지).\n" + "\n".join(lines[:200]) +
                              "\n별칭이 맞으면 sigungu_alias.csv에 한 줄 추가(출처 기록) 후 다시 실행.")
    if period:
        msg = "별칭 유효기간 밖 연도(같은 코드로 붙였음; 원자료의 소급 표기인지 확인):\n" + "\n".join(
            f"  - '{n}' 연도 {sorted(ys)[:8]} : {m}" for (n, m), ys in sorted(period.items()))
        if strict_period:
            raise SigunguKeyError(msg)
        warnings.warn(msg, stacklevel=2)
    if _detail:
        return detail
    return out[0] if scalar else out


def name_kind(names, key_csv: str = KEY_CSV, alias_csv: str = ALIAS_CSV):
    """각 이름의 종류: 시군구 / 일반구(하위구역) / 통합전 / 시도 / 전국 / 모호 / 미확인 / 검토필요."""
    T = load_tables(key_csv, alias_csv); out = []
    for nm in (names if hasattr(names, "__iter__") and not isinstance(names, str) else [names]):
        rows, how, m = _resolve(nm, T)
        if not rows:
            out.append({"시도": "시도", "전국": "전국", "검토필요": "검토필요"}.get(how, "미확인")); continue
        codes = {r[0] for r in rows}; kinds = {r[1] for r in rows}
        out.append("모호" if len(codes) > 1 else ("일반구" if kinds == {"일반구"} else ("통합전" if kinds <= {"통합전"} else "시군구")))
    return out


def attach_code(df: pd.DataFrame, name_col: str, year_col: str | None = None, out_col: str = "sigungu_cd",
                subdistrict: str = "error", drop_kinds=None, strict_period: bool = False, **kw) -> pd.DataFrame:
    """df에 sigungu_cd 열을 붙여 새 df를 돌려줌. 이름 하나라도 실패하면 멈춤.

    drop_kinds: 미리 빼고 싶은 행 종류 목록(예: ['시도','전국']). 뺀 행 수를 출력함(조용히 빼지 않음).
    통합 전 단위(마산·진해·청원·연기 등)가 같은 코드로 붙으면 (코드, 연도) 중복이 생길 수 있으니 합산은 직접 하고 panel_gate로 확인.
    """
    d = df.copy()
    if drop_kinds:
        k = pd.Series(name_kind(d[name_col].tolist(), **kw), index=d.index)
        drop = k.isin(drop_kinds)
        if drop.any():
            print(f"[attach_code] {name_col}: {drop.sum()}행 제외(종류 {sorted(k[drop].unique())}), 예: {d.loc[drop, name_col].drop_duplicates().head(5).tolist()}")
        d = d[~drop].copy()
    yrs = d[year_col].tolist() if year_col else None
    d[out_col] = to_code(d[name_col].tolist(), year=yrs, subdistrict=subdistrict, strict_period=strict_period, **kw)
    det = to_code(d[name_col].drop_duplicates().tolist(), subdistrict=subdistrict, _detail=True, **kw)
    T = load_tables(**kw)
    merged = [(n, c) for n, y, c, h, m in det if any(r[1] == "통합전" for r in _resolve(n, T)[0])]
    if merged:
        print(f"[attach_code] 통합 전 단위 표기 {len(merged)}개가 통합 코드로 붙음 — (코드, 연도)별 합산 필요: {merged[:8]}")
    return d


def explicit_fill(df: pd.DataFrame, col: str, how: str, reason: str, code_col: str = "sigungu_cd", year_col: str = "year",
                  log_path: str | None = None) -> pd.DataFrame:
    """결측 채움이 꼭 필요할 때만. how='zero'(건수형 구조적 0) | 'interp'(같은 지역 안 선형 보간, 양끝은 가까운 값).
    전국 최솟값·평균 같은 '남의 값' 채움은 제공하지 않음. 채운 칸 목록을 출력하고 log_path가 있으면 CSV로 덧붙임."""
    if not reason:
        raise ValueError("reason(채우는 이유)을 반드시 적을 것")
    d = df.sort_values([code_col, year_col]).copy()
    miss = d[col].isna()
    if how == "zero":
        d.loc[miss, col] = 0
    elif how == "interp":
        d[col] = d.groupby(code_col)[col].transform(lambda s: s.interpolate(limit_direction="both"))
    else:
        raise ValueError("how는 'zero' 또는 'interp'")
    filled = d.loc[miss, [code_col, year_col, col]].assign(column=col, how=how, reason=reason)
    print(f"[explicit_fill] {col}: {len(filled)}칸 채움({how}) — {reason}")
    if len(filled):
        print(filled.head(30).to_string(index=False))
    if log_path and len(filled):
        filled.to_csv(log_path, mode="a", header=not os.path.exists(log_path), index=False, encoding="utf-8-sig")
    still = d[col].isna().sum()
    if still:
        raise SigunguKeyError(f"{col}: 채운 뒤에도 {still}칸 결측(지역 전체가 빈 경우). 원자료에서 확인할 것")
    return d


def read_shp_checked(path: str = DEFAULT_SHP, crs: int = 5179, area_tol: float = 0.05, key_csv: str = KEY_CSV,
                     alias_csv: str = ALIAS_CSV):
    """final.shp를 읽어 sigungu_cd 기준으로 canonical_region을 다시 붙이고 단언.
    단언: 행 229, 코드 고유·229와 일치, sigungu_nm이 코드의 이름과 일치, 기하 면적이 key 면적 ±5%.
    shp의 region 라벨은 믿지 않음(다르면 경고로 알리고 코드 기준 값으로 바꿈; 원래 값은 region_shp 열)."""
    import geopandas as gpd
    T = load_tables(key_csv, alias_csv); key = T["key"]
    g = gpd.read_file(path)
    if g.crs is None:
        g = g.set_crs(crs)
    elif g.crs.to_epsg() != crs:
        g = g.to_crs(crs)
    errs = []
    if "sigungu_cd" not in g.columns:
        raise SigunguKeyError("shp에 sigungu_cd 열 없음")
    g["sigungu_cd"] = g["sigungu_cd"].astype(str).str.strip().str.zfill(5)
    if len(g) != 229: errs.append(f"행 수 {len(g)} != 229")
    dup = g.sigungu_cd[g.sigungu_cd.duplicated()].tolist()
    if dup: errs.append(f"중복 코드 {dup}")
    extra = sorted(set(g.sigungu_cd) - T["codes"]); missing = sorted(T["codes"] - set(g.sigungu_cd))
    if extra: errs.append(f"229에 없는 코드 {extra}")
    if missing: errs.append(f"빠진 코드 {missing}")
    k = key.set_index("sigungu_cd")
    if "sigungu_nm" in g.columns:
        for c, nm in zip(g.sigungu_cd, g.sigungu_nm):
            if c not in k.index: continue
            rows, _, _ = _resolve(nm, T, sido_hint=k.at[c, "sido"])
            if c not in {r[0] for r in rows}:
                rows2, _, _ = _resolve(nm, T)
                if c not in {r[0] for r in rows2}:
                    errs.append(f"이름 불일치 {c}: shp '{nm}' vs key '{k.at[c, 'canonical_region']}'")
    a = g.geometry.area / 1e6
    for c, ai in zip(g.sigungu_cd, a):
        if c in k.index:
            ref = float(k.at[c, "area_km2"])
            if abs(ai / ref - 1) > area_tol:
                errs.append(f"면적 불일치 {c} {k.at[c, 'canonical_region']}: shp {ai:.1f} vs key {ref:.1f} km²")
    if errs:
        raise SigunguKeyError("read_shp_checked 실패:\n  " + "\n  ".join(errs[:60]))
    g["canonical_region"] = g.sigungu_cd.map(k["canonical_region"])
    if "region" in g.columns:
        lab = g["region"].astype(str).map(norm)
        bad = g.loc[lab != g.canonical_region, ["sigungu_cd", "region", "canonical_region"]]
        if len(bad):
            warnings.warn("shp region 라벨이 코드와 다른 행(코드 기준으로 교정함):\n" + bad.to_string(index=False), stacklevel=2)
        g = g.rename(columns={"region": "region_shp"})
    g["region"] = g["canonical_region"]
    g["area_km2"] = a.values
    return g


def _main(argv):
    import argparse
    ap = argparse.ArgumentParser(description="시군구 이름 -> sigungu_cd")
    ap.add_argument("names", nargs="*"); ap.add_argument("--year", type=int)
    ap.add_argument("--file"); ap.add_argument("--col", default="region"); ap.add_argument("--year-col")
    ap.add_argument("--subdistrict", default="error")
    a = ap.parse_args(argv)
    try:
        if a.file:
            d = pd.read_csv(a.file, dtype=str, encoding="utf-8-sig")
            d = attach_code(d, a.col, a.year_col, subdistrict=a.subdistrict)
            print(f"통과: {a.file} {a.col} {d[a.col].nunique()}개 이름 -> {d.sigungu_cd.nunique()}개 코드")
        else:
            for n, y, c, h, m in to_code(a.names, year=a.year, subdistrict=a.subdistrict, _detail=True):
                print(f"{n}\t{c}\t{load_tables()['name_of'][c]}\t{h}\t{m}")
    except SigunguKeyError as e:
        print(str(e)); return 1
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
