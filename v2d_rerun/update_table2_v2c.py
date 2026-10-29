# -*- coding: utf-8 -*-
# update_table2_v2c.py — 새 버전 폴더(JKMS_FBD_v2c)의 ehec_tables.json t2와 content_ehec.py 사본의 Table 2 각주 DIC·φ를 results_ST 최신 540 산출로 바꿈(원본 JKMS_FBD는 건드리지 않음)
import os, json, re, pandas as pd
ROOT = os.environ.get("EHEC_ROOT",".")
N = os.path.join(ROOT, "_v2_rerun_261022/JKMS/JKMS_FBD_v2c"); O = os.path.join(ROOT, "JKMS_FBD")
T = pd.read_csv(os.path.join(ROOT, "results_ST/MASTER_final_table2.csv")).set_index("cov")
OL = pd.read_csv(os.path.join(ROOT, "_v2_promote_261021/17_KEY_기준_MASTER_final_table2_260903.csv")); L2C = dict(zip(OL.label, OL["cov"]))
J = json.load(open(os.path.join(O, "ehec_tables.json"), encoding="utf-8")); new = []
for lab, irr, cri, cred, hdr, unit in J["t2"]:
    if hdr: new.append([lab, irr, cri, cred, hdr, unit]); continue
    r = T.loc[L2C[lab.replace("‡", "").strip().replace("1,000", "1000")]]
    new.append([lab, f"{r.IRR:.3f}", f"{r.lo:.3f}-{r.hi:.3f}", "credible" if str(r.cred).upper() == "TRUE" else "", hdr, unit])
J["t2"] = new; json.dump(J, open(os.path.join(N, "ehec_tables.json"), "w", encoding="utf-8"), ensure_ascii=False)
fit = dict(l.strip().split("=") for l in open(os.path.join(ROOT, "results_ST/final_fit.txt")))
p = os.path.join(N, "content_ehec.py"); s = open(p, encoding="utf-8").read()
s2 = re.sub(r"deviance information criterion \(DIC\) [\d.]+; structured spatial fraction [\d.]+\.", f"deviance information criterion (DIC) {fit['DIC']}; structured spatial fraction {fit['phi']}.", s)
assert s2 != s or f"(DIC) {fit['DIC']};" in s
open(p, "w", encoding="utf-8").write(s2); print("t2 rows", sum(1 for x in new if not x[4]), "| DIC", fit["DIC"], "phi", fit["phi"])
