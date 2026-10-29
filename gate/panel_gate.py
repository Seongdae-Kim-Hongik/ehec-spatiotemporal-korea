#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""패널 결합 게이트 (2026-10-21): 결합을 끝낸 뒤, 모형 적합 전에 반드시 통과시킨다.

실패(멈춤)
  1. 229에 없는 코드, 이름을 코드로 못 바꾼 행
  2. 같은 (시군구, 연도)가 두 번 이상
  3. 지정한 열(인구, 기대건수 E, 결과)이 229 × 지정 연도 모든 칸에 있고 유한한지, positive 열(인구·E)은 > 0인지
경고(멈추지 않음, --strict면 멈춤)
  4. 채움 의심: 값이 그 열의 전국 최솟값·최댓값과 같은 칸이 한 지역에서 여러 해 연속 / 여러 지역에서 반복
     양수 열이 같은 값으로 여러 해 연속(보간 rule=2·상수 채움 흔적)
  5. 연도 간 1.3배 초과 변화(positive 열, 그해 전국 중앙 배율 대비 상대 배율) — 실제 인구 변화(세종·신도시 등)일 수 있으므로 경고만
     통합 연도(청원→청주 2014, 연기→세종 2012, 마산·진해→창원 2010, 북·남제주 2006) 전후는 ±10%로 더 엄하게
  6. (선택) area_col이 sigungu_key.csv 면적과 ±5% 넘게 다름

사용(명령줄): python3 panel_gate.py panel.csv --name-col region --cols pop,E_std,cases --positive pop,E_std --years 2011-2024
사용(코드):   from panel_gate import panel_gate; panel_gate(df, cols=[...], positive=[...], years=range(2011,2025), name_col="region")
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sigungu_key as sk


class PanelGateError(RuntimeError):
    pass


def _runs(years, flags):
    """연속 연도 구간 [(시작, 끝, 길이)] (flags가 True인 연도들)."""
    out = []; cur = None
    for y, f in zip(years, flags):
        if f and cur and y == cur[1] + 1:
            cur[1] = y
        elif f:
            if cur: out.append(cur)
            cur = [y, y]
        else:
            if cur: out.append(cur); cur = None
    if cur: out.append(cur)
    return [(a, b, b - a + 1) for a, b in out]


def _fmt_years(s):
    s = sorted(set(int(y) for y in s)); rng = list(range(s[0], s[-1] + 1))
    return ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b, _ in _runs(rng, [y in s for y in rng]))


def panel_gate(df: pd.DataFrame, cols, positive=None, years=None, code_col: str | None = None, name_col: str | None = None,
               year_col: str = "year", jump: float = 1.3, run_min: int = 2, const_run: int = 3, area_col: str | None = None, jump_cols=None, merge_jump: float = 1.1,
               strict: bool = False, label: str = "", verbose: bool = True, max_show: int = 60, raise_on_fail: bool = True):
    positive = list(positive or []); cols = list(cols)
    for c in cols + positive:
        if c not in df.columns:
            raise PanelGateError(f"[{label}] 열 없음: {c}")
    T = sk.load_tables(); key = T["key"]; name_of = T["name_of"]; codes229 = sorted(T["codes"])
    d = df.copy()
    fails, warns = [], []
    # 1) 코드
    if code_col is None:
        if name_col is None:
            raise PanelGateError("code_col 또는 name_col 하나는 필요")
        try:
            d = sk.attach_code(d, name_col, None)
        except sk.SigunguKeyError as e:
            raise PanelGateError(f"[{label}] 실패: 이름을 코드로 못 바꿈\n{e}")
        code_col = "sigungu_cd"
    d[code_col] = d[code_col].astype(str).str.strip().str.replace(r"\.0$", "", regex=True).str.zfill(5)
    d[year_col] = pd.to_numeric(d[year_col], errors="coerce")
    if d[year_col].isna().any():
        fails.append(("연도 결측", d[d[year_col].isna()].head(max_show)))
    yrs = sorted(int(y) for y in (years if years is not None else d[year_col].dropna().unique()))
    d = d[d[year_col].isin(yrs)].copy(); d[year_col] = d[year_col].astype(int)
    bad = sorted(set(d[code_col]) - set(codes229))
    if bad:
        fails.append(("229에 없는 코드", pd.DataFrame({"sigungu_cd": bad})))
    # 2) 중복
    dup = d[d.duplicated([code_col, year_col], keep=False)]
    if len(dup):
        t = dup.groupby([code_col, year_col]).size().reset_index(name="행수")
        t["이름"] = t[code_col].map(name_of)
        if name_col: t["원래 이름들"] = [", ".join(sorted(set(dup.loc[(dup[code_col] == c) & (dup[year_col] == y), name_col].astype(str)))) for c, y in zip(t[code_col], t[year_col])]
        fails.append(("같은 (시군구, 연도) 중복", t))
    # 3) 완전성
    grid = pd.MultiIndex.from_product([codes229, yrs], names=[code_col, year_col]).to_frame(index=False)
    dd = d.drop_duplicates([code_col, year_col])
    g = grid.merge(dd[[code_col, year_col] + list(dict.fromkeys(cols + positive + ([area_col] if area_col else [])))].assign(_row=1),
                   on=[code_col, year_col], how="left")
    miss = []
    for c in dict.fromkeys(cols + positive):
        x = pd.to_numeric(g[c], errors="coerce")
        why = np.where(g._row.isna(), "행 없음", np.where(x.isna(), "결측", np.where(~np.isfinite(x.fillna(0)), "비유한",
                       np.where((c in positive) & (x <= 0), "0 이하", ""))))
        m = g.loc[why != "", [code_col, year_col]].assign(열=c, 이유=why[why != ""], 값=x[why != ""])
        miss.append(m)
    miss = pd.concat(miss, ignore_index=True)
    if len(miss):
        miss["이름"] = miss[code_col].map(name_of)
        # 연속 연도로 묶어 보기 좋게
        summ = (miss.sort_values([code_col, "열", year_col]).groupby([code_col, "이름", "열", "이유"])[year_col]
                .apply(_fmt_years)
                .reset_index().rename(columns={year_col: "연도"}))
        fails.append((f"빠진 칸 {len(miss)}개(229 × {yrs[0]}-{yrs[-1]})", summ))
    # 4) 채움 의심
    fill = []
    for c in dict.fromkeys(cols + positive):
        x = pd.to_numeric(g[c], errors="coerce")
        ok = np.isfinite(x)
        if ok.sum() < 3: continue
        pos_only = c in positive
        xs = x[ok & ((x > 0) if pos_only else True)]
        if not len(xs): continue
        for tag, v in (("전국 최솟값", xs.min()), ("전국 최댓값", xs.max())):
            if v == 0 and not pos_only: continue
            hit = ok & np.isclose(x, v, rtol=1e-9, atol=0)
            if hit.sum() < 2: continue
            regs = g.loc[hit, code_col].unique()
            for r in regs:
                sub = g[g[code_col] == r].sort_values(year_col)
                rr = _runs(sub[year_col].tolist(), np.isclose(pd.to_numeric(sub[c], errors="coerce"), v, rtol=1e-9, atol=0).tolist())
                for a, b, n in rr:
                    if n >= run_min:
                        fill.append((r, name_of[r], c, f"{tag}({v:.6g})와 같은 값", f"{a}-{b}", n))
            if len(regs) >= 2:
                fill.append(("여러", f"{len(regs)}개 지역", c, f"{tag}({v:.6g})가 {int(hit.sum())}칸에 반복", "", int(hit.sum())))
        if pos_only:
            for r, sub in g.sort_values(year_col).groupby(code_col):
                v = pd.to_numeric(sub[c], errors="coerce").to_numpy(); ys = sub[year_col].tolist()
                same = [False] + [bool(np.isfinite(v[i]) and np.isfinite(v[i - 1]) and np.isclose(v[i], v[i - 1], rtol=1e-12, atol=0)) for i in range(1, len(v))]
                for a, b, n in _runs(ys, same):
                    if n + 1 >= const_run:
                        fill.append((r, name_of[r], c, f"같은 값 {v[ys.index(b)]:.6g}이 연속(상수 채움·rule=2 외삽 의심)", f"{a - 1}-{b}", n + 1))
    if fill:
        warns.append(("채움 의심", pd.DataFrame(fill, columns=["sigungu_cd", "이름", "열", "패턴", "연도", "칸수"])))
    # 5) 급변: 지역 배율을 그해 전국 중앙 배율로 나눈 상대 배율로 판단(전국 발생률을 따라 매년 움직이는 E도 지역 고유 급변만 걸림)
    jl = []
    for c in (positive if jump_cols is None else jump_cols):
        w = g.pivot_table(index=code_col, columns=year_col, values=c, aggfunc="first").reindex(columns=yrs)
        w = w.apply(pd.to_numeric, errors="coerce")
        for y0, y1 in zip(yrs[:-1], yrs[1:]):
            if y1 != y0 + 1: continue
            a0, a1 = w[y0], w[y1]; ok = (a0 > 0) & (a1 > 0) & np.isfinite(a0) & np.isfinite(a1)
            if ok.sum() < 10: continue
            ratio = a1[ok] / a0[ok]; med = float(np.median(ratio)); rel = ratio / med
            for r in rel.index[(rel > jump) | (rel < 1 / jump)]:
                jl.append((r, name_of[r], c, f"{y0}→{y1}", a0[r], a1[r], round(float(ratio[r]), 3), round(med, 3)))
    # 5b) 통합 연도 단절: 통합 전 이름(청원·연기·마산·진해·북/남제주)의 마지막 해 직전→그해 상대 배율이 ±10% 넘으면 경고
    #     (통합 전 단위를 빠뜨려 앞 해가 작거나, 옛·새 이름을 함께 더해 부푼 경우). 1.3배 기준으로는 안 잡히는 청주(청원 누락 -19%) 같은 경우용
    al = T["alias"]; mg = al[(al.kind == "통합전") & (al.needs_review.astype(str) != "1") & (al.valid_to != "")]
    trans = sorted({(c, int(v[:4])) for c, v in zip(mg.sigungu_cd, mg.valid_to)})
    tl = []
    for c in (positive if jump_cols is None else jump_cols):
        w = g.pivot_table(index=code_col, columns=year_col, values=c, aggfunc="first").reindex(columns=yrs).apply(pd.to_numeric, errors="coerce")
        for code, y1 in trans:
            y0 = y1 - 1
            if y0 not in w.columns or y1 not in w.columns or code not in w.index: continue
            a0, a1 = w[y0], w[y1]; ok = (a0 > 0) & (a1 > 0)
            if ok.sum() < 10 or not ok.get(code, False): continue
            med = float(np.median(a1[ok] / a0[ok])); rel = (a1[code] / a0[code]) / med
            if rel > merge_jump or rel < 1 / merge_jump:
                tl.append((code, name_of[code], c, f"{y0}→{y1}", a0[code], a1[code], round(float(a1[code] / a0[code]), 3), round(med, 3)))
    if tl:
        warns.append((f"통합 연도 전후 단절(±{round((merge_jump - 1) * 100)}% 넘음; 통합 전 단위 누락·이중 합산 의심)", pd.DataFrame(tl, columns=["sigungu_cd", "이름", "열", "연도", "앞", "뒤", "배율", "전국중앙배율"])))
    if jl:
        warns.append((f"연도 간 {jump}배 초과 변화(전국 중앙 배율 대비; 실제 변화인지 확인)", pd.DataFrame(jl, columns=["sigungu_cd", "이름", "열", "연도", "앞", "뒤", "배율", "전국중앙배율"])))
    # 6) 면적
    if area_col:
        a = g.groupby(code_col)[area_col].median()
        ref = key.set_index("sigungu_cd").area_km2.astype(float)
        rat = (pd.to_numeric(a, errors="coerce") / ref.reindex(a.index))
        off = rat[(rat - 1).abs() > 0.05]
        if len(off):
            warns.append((f"{area_col}이 shp 기하 면적과 ±5% 넘게 다름({len(off)}곳)", pd.DataFrame({"sigungu_cd": off.index, "이름": off.index.map(name_of),
                          area_col: a.reindex(off.index).values, "key_area_km2": ref.reindex(off.index).values, "배율": off.round(3).values})))
    rep = dict(label=label, n_rows=len(df), years=(yrs[0], yrs[-1]), cols=cols, positive=positive, fails=fails, warns=warns,
               passed=not fails and not (strict and warns))
    if verbose:
        print(report_text(rep, max_show))
    if raise_on_fail and (fails or (strict and warns)):
        raise PanelGateError(f"[{label}] panel_gate 실패: " + "; ".join(k for k, _ in fails + (warns if strict else [])))
    return rep


def report_text(rep, max_show=60, md=False):
    L = [f"== panel_gate [{rep['label']}] {rep['n_rows']}행, 연도 {rep['years'][0]}-{rep['years'][1]}, 열 {rep['cols']} (양수 {rep['positive']})"]
    L.append("결과: " + ("통과" if not rep["fails"] else "실패") + (f" · 경고 {len(rep['warns'])}종" if rep["warns"] else " · 경고 없음"))
    for tag, items in (("실패", rep["fails"]), ("경고", rep["warns"])):
        for k, t in items:
            L.append(f"-- {tag}: {k} ({len(t)}줄)")
            L.append(t.head(max_show).to_markdown(index=False) if md else t.head(max_show).to_string(index=False))
            if len(t) > max_show: L.append(f"   ... {len(t) - max_show}줄 더")
    return "\n".join(L)


def _main(argv):
    ap = argparse.ArgumentParser(description="패널 결합 게이트")
    ap.add_argument("panel"); ap.add_argument("--name-col"); ap.add_argument("--code-col"); ap.add_argument("--year-col", default="year")
    ap.add_argument("--cols", required=True); ap.add_argument("--positive", default=""); ap.add_argument("--years")
    ap.add_argument("--jump", type=float, default=1.3); ap.add_argument("--area-col"); ap.add_argument("--strict", action="store_true")
    a = ap.parse_args(argv)
    d = pd.read_csv(a.panel, encoding="utf-8-sig", low_memory=False)
    yrs = None
    if a.years:
        y0, y1 = a.years.split("-"); yrs = range(int(y0), int(y1) + 1)
    try:
        panel_gate(d, cols=a.cols.split(","), positive=[x for x in a.positive.split(",") if x], years=yrs, code_col=a.code_col,
                   name_col=a.name_col, year_col=a.year_col, jump=a.jump, area_col=a.area_col, strict=a.strict, label=os.path.basename(a.panel))
    except PanelGateError as e:
        print(str(e)); return 1
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
