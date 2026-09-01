import openpyxl
import json
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# ============================================================
# 1. Excelから座標データを抽出し、JSONとして整形
# ============================================================
wb = openpyxl.load_workbook('Cadps2_8subgroups_classified_dirmatch.xlsx', data_only=True)
sheets = ['Subgroup1_Aging-specific#1', 'Subgroup2_Aging-specific#2', 'Subgroup3_Aging_and_early_pancr',
          'Subgroup4_Aging_and_longstandin', 'Subgroup5_Pancreatitis-specific', 'Subgroup6_Pancreatitis-specific',
          'Subgroup7_Longstanding_chronic_', 'Subgroup8_Longstanding_chronic_']

# 表現型カテゴリーごとのマーカー遺伝子リスト
markers = {
    "acinar": ["Amy1", "Amy2a5", "Amy2b", "Cela1", "Cela2a", "Cela3b", "Cpa1", "Cpa2", "Cpb1", "Prss1", "Prss2",
               "Prss3", "Ctrb1", "Ctrc", "Pnlip", "Clps", "Manf", "Tmed6", "Tmed11", "Cuzd1", "Serpini2",
               "Serpinb1a", "Retreg1", "Ptf1a", "Nr5a2", "Bhlha15", "Rbpjl"],
    "secretion": ["Cadps", "Cadps2", "Rab27a", "Rab27b", "Syt7", "Syt1", "Snap25", "Stx1a", "Vamp2", "Unc13a",
                  "Unc13b", "Rims1", "Rims2"],
    "fibrosis_immune": ["Col1a1", "Col1a2", "Col3a1", "Col4a1", "Fn1", "Acta2", "Postn", "Pdgfrb", "Tgfb1",
                         "Tgfb2", "Timp1", "Sparc", "Lox", "Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Cd68",
                         "Ptprc", "Itgam", "Ccl2", "Cxcl1", "Il6", "Tnf", "Il1b", "Ccr2"],
    "adm": ["Sox9", "Krt19", "Krt7", "Krt8", "Krt18", "Hnf1b", "Onecut1", "Muc1"],
}


def classify(sym):
    for cat, genes in markers.items():
        if sym in genes:
            return cat
    return None


# 全サブグループシートから Symbol / KO_aging_log2FC / KO_aging_padj を収集
all_genes = []
for sn in sheets:
    ws = wb[sn]
    header = [c.value for c in ws[1]]
    idx = {h: i for i, h in enumerate(header)}
    for row in ws.iter_rows(min_row=2, values_only=True):
        sym = row[idx['Symbol']]
        if sym is None:
            continue
        fc = row[idx['KO_aging_log2FC']]
        padj = row[idx['KO_aging_padj']]
        cat = classify(sym)
        all_genes.append({"sym": sym, "fc": fc, "padj": padj, "cat": cat, "sg": sn})

# 背景用(カテゴリーに該当しない遺伝子)の座標リスト
bg = [
    (round(g['fc'], 2), round(min(-math.log10(max(g['padj'], 1e-300)), 60), 2))
    for g in all_genes if g['cat'] is None
]

# カテゴリー別(表現型マーカーに該当する遺伝子)の座標リスト
cats = {}
for g in all_genes:
    if g['cat']:
        cats.setdefault(g['cat'], []).append({
            "x": round(g['fc'], 2),
            "y": round(min(-math.log10(max(g['padj'], 1e-300)), 60), 2),
            "g": g['sym']
        })

json.dump(bg, open('ko_volcano_bg.json', 'w'))
json.dump(cats, open('ko_volcano_cats.json', 'w'), ensure_ascii=False)

# ============================================================
# 2. ボルケーノプロット作成(Cadps2は除外)
# ============================================================
CAT_COLOR2 = {
    "acinar": "#2a78d6",
    "secretion": "#eb6834",
    "fibrosis_immune": "#e34948",
    "adm": "#1baf7a",
}
CAT_LABEL2 = {
    "acinar": "Acinar / secretory granule",
    "secretion": "Exocytosis machinery",
    "fibrosis_immune": "Fibrosis / immune infiltration",
    "adm": "ADM / ductal",
}
CAT_ORDER2 = ["acinar", "secretion", "fibrosis_immune", "adm"]

fig2, ax2 = plt.subplots(figsize=(7.2, 6.0))

# 背景:KO加齢DEG全体(灰色の点)
bx = [min(p[0], 12) for p in bg]
by = [min(p[1], 60) for p in bg]
ax2.scatter(bx, by, s=8, color="#a3a299", alpha=0.35, linewidths=0, zorder=2)

# ラベルが重ならないようにする微調整オフセット
offsets = {
    "Snap25": (-18, 5), "Rims2": (16, 5), "Serpinb1a": (-30, -2),
    "Cuzd1": (22, 2), "Itgam": (-14, 6), "Ccr2": (-12, -7),
}

# カテゴリーごとに色分け(Cadps2は除外)
for cat in CAT_ORDER2:
    pts = [p for p in cats.get(cat, []) if p["g"] != "Cadps2"]
    if not pts:
        continue
    xs2 = [min(p["x"], 12) for p in pts]
    ys2 = [min(p["y"], 60) for p in pts]
    ax2.scatter(xs2, ys2, s=55, color=CAT_COLOR2[cat], edgecolors="white",
                linewidths=0.7, zorder=4, label=CAT_LABEL2[cat])
    for p in pts:
        x = min(p["x"], 12)
        y = min(p["y"], 60)
        dx, dy = offsets.get(p["g"], (0, 5))
        ax2.annotate(p["g"], (x, y), textcoords="offset points", xytext=(dx, dy),
                     fontsize=7, ha="center", style="italic", color="#3a3936", zorder=5)

ax2.set_xlabel("log2FC (KO 14m vs 6m)", fontsize=9)
ax2.set_ylabel("-log10(padj)", fontsize=9)
ax2.set_title("KO-aging DEGs by phenotype category (Cadps2 excluded)", fontsize=10.5, loc="left", pad=8)
ax2.axvline(0, color="#c3c2b7", linewidth=0.8, zorder=1)
ax2.grid(color="#e1e0d9", linewidth=0.6, zorder=0)
for spine in ["top", "right"]:
    ax2.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax2.spines[spine].set_color("#c3c2b7")
ax2.tick_params(length=0, labelsize=8)

handles2 = [mpatches.Patch(color="#a3a299", alpha=0.5, label="Other KO-aging DEGs")]
handles2 += [mpatches.Patch(color=CAT_COLOR2[c], label=CAT_LABEL2[c]) for c in CAT_ORDER2]
ax2.legend(handles=handles2, loc="upper left", frameon=False, fontsize=7.3)

fig2.savefig("ko_volcano_no_cadps2.pdf", dpi=300, bbox_inches="tight")
fig2.savefig("ko_volcano_no_cadps2.png", dpi=300, bbox_inches="tight")
