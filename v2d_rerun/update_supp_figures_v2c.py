# -*- coding: utf-8 -*-
# update_supp_figures_v2c.py — (2026-10-22 v2c) 보충그림 docx 사본에서 재실행으로 다시 그린 그림만 바꿈(그림 폭 유지, 높이는 새 그림 비율로).
#  Fig S14 초과확률(415) · S15 COVID 전후 RR(499) · S16 창별 결정요인(544) · S17 사양 곡선(543, 다중우주 끝난 뒤에만)
#  원본 = _pre_v2_rerun_261022 정본 사본, 출력 = JKMS_FBD_v2c/JKMS_EHEC/Supplementary/Supplementary_Figures_EHEC.docx
import os, re, zipfile, hashlib, glob
from PIL import Image
ROOT = os.environ.get("EHEC_ROOT",".")
RR = os.path.join(ROOT, "_v2_rerun_261022", "JKMS")
SRC = os.path.join(ROOT, "_pre_v2_rerun_261022/JKMS/JKMS_FBD/JKMS_EHEC/Supplementary/Supplementary_Figures_EHEC.docx")
DST = os.path.join(RR, "JKMS_FBD_v2c/JKMS_EHEC/Supplementary/Supplementary_Figures_EHEC.docx")
OLD = os.path.join(ROOT, "_pre_v2_rerun_261022/JKMS/EHEC_top")
md5 = lambda b: hashlib.md5(b).hexdigest()
RS = os.path.join(ROOT, "results_ST")
mv_done = len(glob.glob(os.path.join(RS, "PREREG_MV_results_seed1[01][0-9].csv"))) == 12 and os.path.exists(os.path.join(ROOT, "EHEC/EHEC_Fig_specification_curve.png")) and \
          os.path.getmtime(os.path.join(ROOT, "EHEC/EHEC_Fig_specification_curve.png")) > os.path.getmtime(os.path.join(RS, "PREREG_MV_summary.csv")) - 1 and \
          os.path.getmtime(os.path.join(RS, "PREREG_MV_summary.csv")) > max(os.path.getmtime(f) for f in glob.glob(os.path.join(RS, "PREREG_MV_results_seed1[01][0-9].csv")))
SWAP = ["EHEC_Fig_agestd_exceedance.png", "EHEC_Fig_covid_split_RR.png", "EHEC_Fig_det_windows.png"] + (["EHEC_Fig_specification_curve.png"] if mv_done else [])
zin = zipfile.ZipFile(SRC); files = {n: zin.read(n) for n in zin.namelist()}
oldhash = {md5(open(os.path.join(OLD, f), "rb").read()): f for f in SWAP}
media = {n: oldhash[md5(b)] for n, b in files.items() if n.startswith("word/media/") and md5(b) in oldhash}
rels = files["word/_rels/document.xml.rels"].decode("utf-8"); doc = files["word/document.xml"].decode("utf-8")
log = []
for mname, fig in media.items():
    newb = open(os.path.join(ROOT, "EHEC", fig), "rb").read(); files[mname] = newb
    rid = re.search(r'Id="(rId\d+)"[^>]*Target="media/%s"' % re.escape(os.path.basename(mname)), rels) or re.search(r'Target="media/%s"[^>]*Id="(rId\d+)"' % re.escape(os.path.basename(mname)), rels)
    rid = rid.group(1); w, h = Image.open(os.path.join(ROOT, "EHEC", fig)).size
    # 이 rId를 쓰는 drawing 블록의 extent(cx 유지, cy = cx × h/w)
    for m in re.finditer(r"<w:drawing>.*?</w:drawing>", doc, flags=re.S):
        blk = m.group(0)
        if f'r:embed="{rid}"' not in blk: continue
        cx = int(re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"', blk).group(1)); cy = int(round(cx * h / w))
        nb = re.sub(r'(<wp:extent cx="\d+" cy=")\d+(")', rf"\g<1>{cy}\2", blk); nb = re.sub(r'(<a:ext cx="\d+" cy=")\d+(")', rf"\g<1>{cy}\2", nb)
        doc = doc.replace(blk, nb); break
    log.append(f"{fig} → {mname} ({rid}), {w}x{h}")
files["word/document.xml"] = doc.encode("utf-8")
os.makedirs(os.path.dirname(DST), exist_ok=True)
with zipfile.ZipFile(DST, "w", zipfile.ZIP_DEFLATED) as z:
    for n in zin.namelist(): z.writestr(n, files[n])
print("\n".join(log)); print("다중우주 그림 교체:", mv_done)
