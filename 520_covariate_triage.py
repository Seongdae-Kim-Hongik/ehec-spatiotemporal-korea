# -*- coding: utf-8 -*-
# 520 — 전수 공변량 1차 분류(triage). 규칙으로 걸러낼 수 있는 것만 걸러내고, 나머지는 수동 판정으로 넘긴다.
#  원칙: 제외는 사유를 남긴다. 결과를 보고 제외하지 않는다(이 단계에서 어떤 IRR도 계산하지 않음).
import pandas as pd, re, unicodedata as ud, os
C="data_cache"; R=pd.read_csv(os.path.join(C,"FULL_COVARIATE_CENSUS.csv"))
Q=R[(R["first"]==2011)&(R.n_years_obs>=10)].copy().reset_index(drop=True)
Q["base"]=Q.col.map(lambda c: re.sub(r"_\d$","",re.sub(r"_(표준화율|조율)(_\d)?$","",c)))

def rule(row):
    c=row.col; f=row.file
    # ── A. 결과변수·불가능 변수 (사유 고정) ──
    if "사망" in c or "사망률" in c: return "EXCLUDE","outcome: mortality (consequence, not exposure)"
    if re.search(r"발생|환자수|이환|유병|감염병_|신고",c) and "식중독" not in c: return "EXCLUDE","outcome: disease incidence (leakage risk)"
    if "식중독" in f: return "EXCLUDE","outcome: foodborne-illness counts (direct leakage)"
    if "청소년" in c or "_1" in c[-2:]: return "EXCLUDE","sub-population/duplicate series"
    if re.search(r"_남자$|_여자$",c): return "EXCLUDE","sex-stratified duplicate of total"
    if "_조율" in c and row.base+"_표준화율" in set(Q.col): return "EXCLUDE","crude duplicate; age-standardized version retained"
    # ── B. 양식·어업 (EHEC 경로 아님; HAV/노로용) ──
    if "양식" in f or "어패류" in f or "해양" in f or "수온" in f: return "EXCLUDE","marine/aquaculture: not an EHEC pathway (retain for HAV/norovirus)"
    # ── C. 비병인 행정지표 ──
    if re.search(r"축제|관광|문화|체육시설|도서관|공원",c): return "EXCLUDE","administrative/leisure: no plausible EHEC pathway"
    if re.search(r"화재|교통사고|범죄|자살|익사|추락",c): return "EXCLUDE","injury/safety: no plausible EHEC pathway"
    # ── D. 병인 경로 후보 (수동 판정 대상) ──
    D=[("D1 reservoir",r"소|한육우|젖소|가축|축산|분뇨|목장|방목|사육"),
       ("D2 agri-land",r"경지|밭|논|임야|농|시설재배|비료|농약|관개|국토|용도"),
       ("D3 water-sanit",r"상수도|하수|지하수|정화조|급수|음용|수질|관로|찌꺼기|슬러지|처리구역"),
       ("D4 hygiene-food",r"손씻|비누|위생|급식|음식점|조리|외식|배달|패스트|주점|식품|영양표시|식생활|아침"),
       ("D5 SES-demog",r"재정|기초생활|수급|소득|교육|학교|대졸|실업|고령|독거|인구|성비|외국인|가구|주택|경제|세|예산"),
       ("D6 behaviour",r"접종|검진|흡연|음주|걷기|운동|신체활동|비만|체중|건강수준|스트레스|수면|영양"),
       ("D7 detection",r"의사|병원|의료|진료|병상|응급|보건소|미충족|의원|약국|인력")]
    for dom,pat in D:
        if re.search(pat,c) or re.search(pat,f): return "CANDIDATE",dom
    return "REVIEW","unclassified: manual review"

Q[["status","reason"]]=Q.apply(lambda r: pd.Series(rule(r)),axis=1)
Q.to_csv(os.path.join(C,"COVARIATE_TRIAGE_v1.csv"),index=False,encoding="utf-8-sig")
print(f"후보 풀 {len(Q):,}개 → 1차 분류")
print(Q.status.value_counts().to_string()); print()
print("── EXCLUDE 사유 ──"); print(Q[Q.status=="EXCLUDE"].reason.value_counts().to_string()); print()
print("── CANDIDATE 도메인 ──"); print(Q[Q.status=="CANDIDATE"].reason.value_counts().to_string()); print()
print(f"── REVIEW(미분류) {int((Q.status=='REVIEW').sum())}개 — 예시 ──")
for c in Q[Q.status=="REVIEW"].col.head(40): print("   ",c)
