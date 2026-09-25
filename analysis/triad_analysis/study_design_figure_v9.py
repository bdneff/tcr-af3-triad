"""
Study design figure - v9.

Refined card-based layout:
- IEDB pipeline laid out as a horizontal filmstrip (5 panels left-to-right)
- Two held-out cards (PDB, CRESTA) sit ABOVE the IEDB card
- Filter 2 connector: dashed arrows from held-out cards DOWN into IEDB's Filter 2 panel
- AF3 classifier + AUC results sit BELOW everything
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
from matplotlib.lines import Line2D

# ---------------- canvas ----------------
fig = plt.figure(figsize=(17, 10.5), dpi=200, facecolor="white")
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_xlim(0, 170)
ax.set_ylim(0, 108)
ax.axis("off")

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
inv_data = ax.transData.inverted()

def text_size_data(s, **kwargs):
    sentinel = "x"
    t = ax.text(0, 0, sentinel + s + sentinel, **kwargs)
    bb = t.get_window_extent(renderer=renderer)
    (x0a, y0), (x1a, y1) = inv_data.transform(((bb.x0, bb.y0), (bb.x1, bb.y1)))
    h = y1 - y0
    t.remove()
    t2 = ax.text(0, 0, sentinel + sentinel, **kwargs)
    bb2 = t2.get_window_extent(renderer=renderer)
    (x0b, _), (x1b, _) = inv_data.transform(((bb2.x0, bb2.y0), (bb2.x1, bb2.y1)))
    w_sent = x1b - x0b
    t2.remove()
    return ((x1a - x0a) - w_sent, h)

# ---------------- palette ----------------
INK     = "#1F2933"
DIM     = "#52606D"
SOFT    = "#8B95A1"
HAIR    = "#D9DEE3"

IEDB_DK = "#1E5F8E";  IEDB_MD = "#3D8DBE"; IEDB_LT = "#D6E7F2"
PDB_DK  = "#9A3F0C";  PDB_MD  = "#C36C2E"; PDB_LT  = "#F2DCC7"
CRES_DK = "#5E2F75";  CRES_MD = "#8A57A6"; CRES_LT = "#E6D6EF"
NON_DK  = "#7C2A2A";  NON_MD  = "#B04545"; NON_LT  = "#FCEFEF"
MODEL_BG = "#1F2933"

# ---------------- helpers ----------------
def text(x, y, s, size=11, color=INK, weight="normal", ha="center", va="center",
         style="normal", zorder=5):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, fontstyle=style, zorder=zorder)

def arrow(x0, y0, x1, y1, color=INK, lw=1.6, head=12, zorder=4, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=head, color=color, linewidth=lw,
                                 zorder=zorder, shrinkA=0, shrinkB=0, linestyle=ls))

def line(x0, y0, x1, y1, color=HAIR, lw=0.8, ls="-", zorder=1):
    ax.add_line(Line2D([x0, x1], [y0, y1], color=color, lw=lw,
                       linestyle=ls, zorder=zorder))

def rounded(x, y, w, h, fc="white", ec=HAIR, lw=1.0, r=0.6, zorder=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0.02,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder))

def chip(cx, cy, label, *, size=10, color=INK, weight="bold",
         fc="white", ec=HAIR, lw=1.0, pad_x=1.3, pad_y=0.6, zorder=3):
    w_text, h_text = text_size_data(label, fontsize=size, fontweight=weight)
    w = w_text + 2 * pad_x
    h = h_text + 2 * pad_y
    rounded(cx - w/2, cy - h/2, w, h, fc=fc, ec=ec, lw=lw, r=h/2, zorder=zorder)
    text(cx, cy, label, size=size, color=color, weight=weight, zorder=zorder + 2)
    return w, h

# =================================================================
# TITLE
# =================================================================
text(85, 104, "Study design & data inclusion", size=22, color=INK, weight="bold")
text(85, 100.5, "Three datasets, three roles  ·  AF3-based TCR : antigen specificity classifier",
     size=11.5, color=DIM, style="italic")

# =================================================================
# TOP ROW: Two HELD-OUT cards side by side
# =================================================================
# Banner header for held-out region
HOLD_BANNER_Y = 95.5
text(85, HOLD_BANNER_Y, "HELD-OUT VALIDATION SETS", size=10, color="#A07414", weight="bold")
text(85, HOLD_BANNER_Y - 1.8,
     "sealed away during classifier development  —  used only for the final test in Fig 4",
     size=9, color=DIM, style="italic")
line(8, 91.5, 162, 91.5, color="#D9C68F", lw=0.7)

# PDB card (left)
PDB_X, PDB_Y, PDB_W, PDB_H = 8, 65, 75, 25   # taller card to fit split + bracket
def heldout_card(x, y, w, h, db, db_dk, db_lt, source_line, class_lbl, n_label,
                 filter_text, noncog_count, fig_tag_text, db_aside=None):
    rounded(x, y, w, h, fc="white", ec=db_dk, lw=1.6, r=1.0, zorder=2)
    # color header
    rounded(x, y + h - 6, w, 6, fc=db_dk, ec=db_dk, lw=0, r=1.0, zorder=3)
    ax.add_patch(Rectangle((x, y + h - 6), w, 3,
                           facecolor=db_dk, edgecolor="none", zorder=3))
    # db name
    text(x + 4, y + h - 2.5, db, size=15, color="white", weight="bold", ha="left", zorder=5)
    # optional aside next to db name (smaller, lighter, italic)
    if db_aside:
        db_w, _ = text_size_data(db, fontsize=15, fontweight="bold")
        text(x + 4 + db_w + 1.8, y + h - 2.6, db_aside, size=10, color="#FFE9D2",
             ha="left", style="italic", zorder=5)
    # source line
    text(x + 4, y + h - 4.8, source_line, size=8.5, color="white", ha="left",
         style="italic", zorder=5)
    # class label on right
    text(x + w - 4, y + h - 3.5, class_lbl, size=9, color="white", weight="bold",
         ha="right", zorder=5)
    # Body: three columns
    # Col 1: starting cognates chip
    body_y_mid = y + (h - 6) / 2
    chip(x + 9, body_y_mid + 2, n_label, size=10, color=db_dk, fc=db_lt, ec=db_dk,
         lw=0.9, pad_x=1.5)
    text(x + 9, body_y_mid - 2, "starting", size=8.5, color=DIM, style="italic")
    text(x + 9, body_y_mid - 3.8, "cognate triads", size=8.5, color=DIM, style="italic")
    # arrow col1 -> col2
    arrow(x + 18, body_y_mid - 0.5, x + 22, body_y_mid - 0.5, color=db_dk, lw=1.4, head=10)
    # Col 2: filter description
    text(x + 24, body_y_mid + 2.5, "Filter 1", size=10, color=INK, weight="bold", ha="left")
    text(x + 24, body_y_mid + 0.3, filter_text, size=8.5, color=DIM, ha="left", style="italic")
    text(x + 24, body_y_mid - 1.8, "no triads removed (none had AF3-training homology)",
         size=8, color=DIM, ha="left", style="italic")
    # arrow col2 -> col3
    col3_x = x + w - 18
    arrow(col3_x - 4, body_y_mid - 0.5, col3_x - 0.5, body_y_mid - 0.5,
          color=NON_MD, lw=1.4, head=10)
    # Col 3: non-cognate addition
    chip(col3_x + 5, body_y_mid + 2, noncog_count, size=10, color=NON_DK, fc=NON_LT,
         ec=NON_DK, lw=0.9, pad_x=1.5)
    text(col3_x + 5, body_y_mid - 2, "non-cognate", size=8.5, color=DIM, style="italic")
    text(col3_x + 5, body_y_mid - 3.8, "controls added", size=8.5, color=DIM, style="italic")
    # fig tag - top-right corner overlay
    if fig_tag_text:
        chip(x + w - 8, y + h - 7.5, fig_tag_text, size=8.5, color=db_dk, fc="#FFFCEB",
             ec=db_dk, lw=1.0, pad_x=1.2, pad_y=0.5, zorder=6)

# ----- Custom PDB card: shows pre/post-cutoff split + Fig 1 bracket -----
def pdb_card(x, y, w, h):
    # outer card
    rounded(x, y, w, h, fc="white", ec=PDB_DK, lw=1.6, r=1.0, zorder=2)
    # color header
    rounded(x, y + h - 6, w, 6, fc=PDB_DK, ec=PDB_DK, lw=0, r=1.0, zorder=3)
    ax.add_patch(Rectangle((x, y + h - 6), w, 3,
                           facecolor=PDB_DK, edgecolor="none", zorder=3))
    text(x + 4, y + h - 2.5, "PDB", size=15, color="white", weight="bold", ha="left", zorder=5)
    db_w, _ = text_size_data("PDB", fontsize=15, fontweight="bold")
    text(x + 4 + db_w + 1.8, y + h - 2.6,
         "(structures confirmed via STCRDab)",
         size=10, color="#FFE9D2", ha="left", style="italic", zorder=5)
    text(x + 4, y + h - 4.8,
         "Solved MHC:peptide:TCR crystal structures  ·  partitioned by AF3 training cutoff",
         size=8.5, color="white", ha="left", style="italic", zorder=5)
    text(x + w - 4, y + h - 3.5, "CLASS I + II", size=9, color="white", weight="bold",
         ha="right", zorder=5)

    # ---------- Body: two-stage layout ----------
    # Body region: from y to y + h - 7 (below header)
    # Top of body: bracket label chip
    # Then bracket arm
    # Then sub-cards
    # Then arrow to held-out chip + non-cognate annotation
    body_top = y + h - 7      # 76 for our default card size
    body_bot = y + 1.5        # 59.5

    # Sub-card geometry: placed lower in the body so bracket above them is INSIDE the body
    sub_w = 18
    sub_h = 8
    sub_y = body_bot + 4       # ~63.5, leaving room below for held-out path? Actually we put held-out alongside, not below
    # Actually let's put the sub-cards center-vertically in the body and bracket above them
    sub_y = body_top - 12     # 64 - sub-cards span y=64 to y=72
    pre_x  = x + 4
    post_x = x + 4 + sub_w + 3

    # bracket region (above sub-cards): label chip at top, then horizontal bar, then arms down to sub-cards
    br_label_y = body_top - 1.7
    br_arm_y   = body_top - 3.2  # horizontal bar y-position
    br_x0 = pre_x
    br_x1 = post_x + sub_w
    # horizontal bar
    line(br_x0, br_arm_y, br_x1, br_arm_y, color=PDB_DK, lw=1.6)
    # vertical arms dropping from the bar down toward the top of the sub-cards
    line(br_x0, br_arm_y, br_x0, sub_y + sub_h + 0.5, color=PDB_DK, lw=1.6)
    line(br_x1, br_arm_y, br_x1, sub_y + sub_h + 0.5, color=PDB_DK, lw=1.6)
    # small terminal ticks at the bottom of each arm (visual finish)
    line(br_x0 - 0.5, sub_y + sub_h + 0.5, br_x0 + 0.5, sub_y + sub_h + 0.5,
         color=PDB_DK, lw=1.6)
    line(br_x1 - 0.5, sub_y + sub_h + 0.5, br_x1 + 0.5, sub_y + sub_h + 0.5,
         color=PDB_DK, lw=1.6)
    # bracket label chip (centered above the bar)
    chip((br_x0 + br_x1) / 2, br_label_y, "Fig 1 — structural benchmarking",
         size=9, color=PDB_DK, fc="#FFFCEB", ec=PDB_DK, lw=1.0, pad_x=1.5, pad_y=0.5, zorder=6)

    # pre-cutoff sub-card
    rounded(pre_x, sub_y, sub_w, sub_h, fc="#FAF0E5", ec=PDB_MD, lw=0.9, r=0.5, zorder=3)
    text(pre_x + sub_w/2, sub_y + sub_h - 2.0, "130 triads", size=11, color=PDB_DK, weight="bold")
    text(pre_x + sub_w/2, sub_y + sub_h - 4.0, "pre-AF3 cutoff", size=8.5, color=DIM, style="italic")
    text(pre_x + sub_w/2, sub_y + sub_h - 5.7, "89 class I  +  41 class II", size=8, color=DIM)
    text(pre_x + sub_w/2, sub_y + 0.7, "(in AF3 training data)", size=7.5, color=DIM, style="italic")

    # post-cutoff sub-card
    rounded(post_x, sub_y, sub_w, sub_h, fc=PDB_LT, ec=PDB_DK, lw=1.0, r=0.5, zorder=3)
    text(post_x + sub_w/2, sub_y + sub_h - 2.0, "21 triads", size=11, color=PDB_DK, weight="bold")
    text(post_x + sub_w/2, sub_y + sub_h - 4.0, "post-AF3 cutoff", size=8.5, color=DIM, style="italic")
    text(post_x + sub_w/2, sub_y + sub_h - 5.7, "16 class I  +  5 class II", size=8, color=DIM)
    text(post_x + sub_w/2, sub_y + 0.7, "(unseen by AF3 training)", size=7.5, color=DIM, style="italic")

    # Filter 1 -> 16 class I held-out for Fig 4a
    # Arrow from post-cutoff sub-card to a "16 class I" result on the right
    held_x = post_x + sub_w + 4
    arrow(post_x + sub_w + 0.3, sub_y + sub_h/2,
          held_x - 0.3, sub_y + sub_h/2,
          color=PDB_DK, lw=1.4, head=10)
    # The held-out subset: 16 class I triads + 181 non-cognates
    h_chip_y = sub_y + sub_h/2 + 1.5
    chip(held_x + 8, h_chip_y, "16 class I", size=10, color=PDB_DK,
         fc=PDB_LT, ec=PDB_DK, lw=1.0, pad_x=1.5)
    text(held_x + 8, h_chip_y - 2.5, "held-out for Fig 4a", size=8, color=DIM, style="italic")
    text(held_x + 8, h_chip_y - 4.2, "+ 181 non-cognate controls", size=8, color=NON_DK,
         style="italic", weight="bold")

# Place the PDB card
pdb_card(PDB_X, PDB_Y, PDB_W, PDB_H)

CRES_X, CRES_Y, CRES_W, CRES_H = 87, 72, 75, 18
heldout_card(CRES_X, CRES_Y, CRES_W, CRES_H,
             "CRESTA", CRES_DK, CRES_LT,
             "Single-cell barcoded-probe assay (Mead et al. 2024)  ·  M. tuberculosis epitopes",
             "CLASS II",
             "205 triads",
             "None of the 8 antigens has prior PDB structures",
             "+ 1,228",
             None,
             db_aside="(antigens never structurally characterized)")

# =================================================================
# MIDDLE: IEDB card with horizontal pipeline
# =================================================================
# (Section banner removed - IEDB card already says "TRAINING SOURCE")

# Big IEDB card
IEDB_X, IEDB_Y, IEDB_W, IEDB_H = 8, 30, 154, 25
rounded(IEDB_X, IEDB_Y, IEDB_W, IEDB_H, fc="white", ec=IEDB_DK, lw=1.6, r=1.0, zorder=2)
# colored header strip across top
rounded(IEDB_X, IEDB_Y + IEDB_H - 5, IEDB_W, 5, fc=IEDB_DK, ec=IEDB_DK, lw=0, r=1.0, zorder=3)
ax.add_patch(Rectangle((IEDB_X, IEDB_Y + IEDB_H - 5), IEDB_W, 2.5,
                       facecolor=IEDB_DK, edgecolor="none", zorder=3))
text(IEDB_X + 4, IEDB_Y + IEDB_H - 2.5, "IEDB",
     size=15, color="white", weight="bold", ha="left", zorder=5)
text(IEDB_X + 18, IEDB_Y + IEDB_H - 2.5,
     "Immune Epitope Database  ·  literature-curated cognate triads, classes I + II",
     size=9.5, color="#CFE2F0", ha="left", zorder=5)
text(IEDB_X + IEDB_W - 4, IEDB_Y + IEDB_H - 2.5, "TRAINING SOURCE",
     size=9, color="white", weight="bold", ha="right", zorder=5)

# Inner pipeline area: 5 horizontal panels
PIPE_TOP = IEDB_Y + IEDB_H - 6.5
PIPE_BOT = IEDB_Y + 1.5
PIPE_H   = PIPE_TOP - PIPE_BOT - 2  # leave headroom

# Five panel centers, evenly distributed across the IEDB inner width
inner_x0 = IEDB_X + 3
inner_x1 = IEDB_X + IEDB_W - 3
n_panels = 5
panel_w = 26
gap = (inner_x1 - inner_x0 - n_panels * panel_w) / (n_panels - 1)
panel_xs = [inner_x0 + i * (panel_w + gap) for i in range(n_panels)]

steps_iedb = [
    {
        "n": 1, "title": "Start",
        "subtitle": "All eligible human triads",
        "count": "9,792 cognate triads",
        "count_sub": "8,104 class I  ·  1,688 class II",
        "fc": IEDB_LT, "ec": IEDB_DK, "circle_c": IEDB_DK,
    },
    {
        "n": 2, "title": "Filter 1",
        "subtitle": "Drop look-alikes of AF3-trained PDB structures",
        "count": "9,336",
        "count_sub": "7,711 class I  ·  1,625 class II",
        "fc": "white", "ec": HAIR, "circle_c": IEDB_DK,
    },
    {
        "n": 3, "title": "Filter 2",
        "subtitle": "Drop look-alikes of held-out PDB / CRESTA",
        "count": "9,327",
        "count_sub": "7,702 class I  ·  1,625 class II",
        "fc": "white", "ec": HAIR, "circle_c": IEDB_DK,
    },
    {
        "n": 4, "title": "+ Non-cognate",
        "subtitle": "Permuted TCRs as negative controls",
        "count": "190,653",
        "count_sub": "157,025 class I  ·  33,628 class II",
        "fc": NON_LT, "ec": NON_MD, "circle_c": NON_DK,
    },
    {
        "n": 5, "title": "Feature selection",
        "subtitle": "Score 59 AF3 features per triad",
        "count": "PTI-PAE",
        "count_sub": "selected as the classifier",
        "fc": IEDB_LT, "ec": IEDB_DK, "circle_c": IEDB_DK,
    },
]

panel_h = PIPE_H
panel_y = PIPE_BOT + 1
for i, step in enumerate(steps_iedb):
    px = panel_xs[i]
    rounded(px, panel_y, panel_w, panel_h, fc=step["fc"], ec=step["ec"], lw=1.1,
            r=0.6, zorder=3)
    # step number circle, top-left
    ax.add_patch(Circle((px + 3, panel_y + panel_h - 3),
                        radius=1.8, facecolor=step["circle_c"],
                        edgecolor="none", zorder=4))
    text(px + 3, panel_y + panel_h - 3, str(step["n"]), size=10, color="white",
         weight="bold", zorder=5)
    # title
    text(px + 6.5, panel_y + panel_h - 3, step["title"],
         size=10.5, color=INK, weight="bold", ha="left", zorder=4)
    # subtitle (italic small)
    text(px + panel_w / 2, panel_y + panel_h - 6.5, step["subtitle"],
         size=8.3, color=DIM, style="italic", zorder=4)
    # count (big, bold, colored)
    count_color = NON_DK if step["circle_c"] == NON_DK else IEDB_DK
    text(px + panel_w / 2, panel_y + panel_h - 11, step["count"],
         size=12, color=count_color, weight="bold", zorder=4)
    # count subtitle
    text(px + panel_w / 2, panel_y + panel_h - 13.5, step["count_sub"],
         size=8, color=DIM, style="italic", zorder=4)

# arrows between panels
for i in range(n_panels - 1):
    x0 = panel_xs[i] + panel_w
    x1 = panel_xs[i + 1]
    y_mid = panel_y + panel_h / 2
    arrow(x0 + 0.3, y_mid, x1 - 0.3, y_mid, color=IEDB_DK, lw=1.6, head=12)

# =================================================================
# Filter 2 connector: dashed lines from held-out cards DOWN into IEDB's Filter 2 panel
# =================================================================
# Filter 2 is panel index 2; its top edge is at panel_y + panel_h
f2_cx = panel_xs[2] + panel_w / 2
f2_top = panel_y + panel_h

# from PDB card's bottom edge down to Filter 2 panel top
pdb_bottom_x = PDB_X + PDB_W * 0.7
pdb_bottom_y = PDB_Y
arrow(pdb_bottom_x, pdb_bottom_y - 0.3,
      f2_cx - 1.5, f2_top + 0.3,
      color=PDB_DK, lw=1.0, head=10, ls=(0, (3, 2)))

cres_bottom_x = CRES_X + CRES_W * 0.3
cres_bottom_y = CRES_Y
arrow(cres_bottom_x, cres_bottom_y - 0.3,
      f2_cx + 1.5, f2_top + 0.3,
      color=CRES_DK, lw=1.0, head=10, ls=(0, (3, 2)))

# small label on the connector midway
mid_y = 58  # between PDB bottom (65) and IEDB top (55)
text(f2_cx, mid_y + 1, "held-out sets", size=8, color=DIM, style="italic", weight="bold")
text(f2_cx, mid_y - 1, "consulted to filter training", size=8, color=DIM, style="italic")

# =================================================================
# BOTTOM ROW: AF3 classifier + AUC results
# =================================================================
# Section banner — left-aligned
text(8, 25, "VALIDATION", size=10, color="#A07414", weight="bold", ha="left")
text(8, 23, "trained classifier scored on each held-out triad  →  per-antigen AUCs (Fig 4)",
     size=9, color=DIM, style="italic", ha="left")

# Classifier block
MODEL_X, MODEL_Y, MODEL_W, MODEL_H = 8, 6, 38, 16
rounded(MODEL_X, MODEL_Y, MODEL_W, MODEL_H, fc=MODEL_BG, ec=MODEL_BG, lw=0, r=1.0, zorder=2)
# left half of classifier: name
text(MODEL_X + 9, MODEL_Y + MODEL_H - 4.5, "AF3",
     size=18, color="white", weight="bold")
text(MODEL_X + 9, MODEL_Y + MODEL_H - 8.5, "PTI-PAE",
     size=12, color="white", weight="bold")
text(MODEL_X + 9, MODEL_Y + MODEL_H - 12, "classifier",
     size=10, color="#9AA5B1")
# divider
line(MODEL_X + 18, MODEL_Y + 2, MODEL_X + 18, MODEL_Y + MODEL_H - 2,
     color="#3E4C59", lw=0.6)
# right half: what it does
text(MODEL_X + 28, MODEL_Y + MODEL_H - 4.5,
     "trained on IEDB feature", size=9, color="#9AA5B1", style="italic")
text(MODEL_X + 28, MODEL_Y + MODEL_H - 7.0,
     "selection (panel 5)", size=9, color="#9AA5B1", style="italic")
text(MODEL_X + 28, MODEL_Y + MODEL_H - 10.5,
     "scored on each held-out", size=9, color="#9AA5B1", style="italic")
text(MODEL_X + 28, MODEL_Y + MODEL_H - 13,
     "triad → outputs AUCs", size=9, color="#9AA5B1", style="italic")

# arrow IEDB feature selection -> classifier
fs_x = panel_xs[4] + panel_w / 2
fs_y = panel_y  # bottom of panel
arrow(fs_x, fs_y - 0.3, MODEL_X + MODEL_W / 2, MODEL_Y + MODEL_H + 0.3,
      color=IEDB_DK, lw=1.8, head=14)

# Held-out data -> AUC pills (the classifier is a step, but each held-out dataset
# ultimately produces a result reported in its own AUC pill).
# Route AROUND the IEDB card (down the outer margins) to avoid crossing it.

# PDB route: from PDB card bottom-left, sweep LEFT to outside the IEDB card,
# DOWN through the channel above the bottom row, then RIGHT under the IEDB card
# to the PDB AUC pill.
PDB_PILL_X = MODEL_X + MODEL_W + 5   # left edge of PDB AUC pill (defined below)
pdb_out_x = 4.5                       # x outside (left of) IEDB card
pdb_channel_y = 24                    # y in the channel between IEDB bottom and bottom row
pdb_start_x = PDB_X + 4               # where the arrow leaves the PDB card
pdb_start_y = PDB_Y                   # bottom edge of PDB card

# Segment 1: down a bit from PDB card
line(pdb_start_x, pdb_start_y, pdb_start_x, pdb_start_y - 2.5,
     color=PDB_MD, lw=1.2, ls=(0, (4, 3)), zorder=3)
# Segment 2: leftward from PDB card to outside IEDB
line(pdb_start_x, pdb_start_y - 2.5, pdb_out_x, pdb_start_y - 2.5,
     color=PDB_MD, lw=1.2, ls=(0, (4, 3)), zorder=3)
# Segment 3: down along outer margin
line(pdb_out_x, pdb_start_y - 2.5, pdb_out_x, pdb_channel_y,
     color=PDB_MD, lw=1.2, ls=(0, (4, 3)), zorder=3)
# Segment 4: right along channel under IEDB to the PDB pill
arrow(pdb_out_x, pdb_channel_y, PDB_PILL_X - 0.3, pdb_channel_y,
      color=PDB_MD, lw=1.2, head=10, ls=(0, (4, 3)))

# CRESTA route: mirror — from CRESTA card bottom-right, sweep RIGHT to outside IEDB card,
# DOWN through channel, then LEFT under IEDB card to CRESTA AUC pill.
CRES_PILL_X = MODEL_X + MODEL_W + 45  # left edge of CRESTA AUC pill
cres_out_x = 165.5                    # x outside (right of) IEDB card
cres_channel_y = 24                   # same channel
cres_start_x = CRES_X + CRES_W - 4
cres_start_y = CRES_Y

line(cres_start_x, cres_start_y, cres_start_x, cres_start_y - 2.5,
     color=CRES_MD, lw=1.2, ls=(0, (4, 3)), zorder=3)
line(cres_start_x, cres_start_y - 2.5, cres_out_x, cres_start_y - 2.5,
     color=CRES_MD, lw=1.2, ls=(0, (4, 3)), zorder=3)
line(cres_out_x, cres_start_y - 2.5, cres_out_x, cres_channel_y,
     color=CRES_MD, lw=1.2, ls=(0, (4, 3)), zorder=3)
arrow(cres_out_x, cres_channel_y, CRES_PILL_X + 36 + 0.3, cres_channel_y,
      color=CRES_MD, lw=1.2, head=10, ls=(0, (4, 3)))

# AUC pills
def auc_pill(x, y, w, h, db, db_dk, auc, sub):
    rounded(x, y, w, h, fc="white", ec=db_dk, lw=1.5, r=0.7, zorder=3)
    ax.add_patch(Rectangle((x, y), 1.5, h, facecolor=db_dk, edgecolor="none", zorder=4))
    # left section: name + AUC label
    text(x + 4, y + h - 3, db, size=12, color=db_dk, weight="bold", ha="left")
    text(x + 4, y + h - 6, "validation AUC", size=8.5, color=DIM, style="italic", ha="left")
    # right section: big AUC + subtitle
    text(x + w - 4, y + h - 5, auc, size=22, color=INK, weight="bold", ha="right")
    text(x + w - 4, y + h - 11, sub, size=8, color=DIM, style="italic", ha="right")

# Two AUC pills to the right of the classifier
auc_pill(MODEL_X + MODEL_W + 5, MOD_Y_TOP := MODEL_Y, 36, MODEL_H,
         "PDB",    PDB_DK,  "0.92", "median, 14 antigens")
auc_pill(MODEL_X + MODEL_W + 45, MODEL_Y, 36, MODEL_H,
         "CRESTA", CRES_DK, "0.78", "median, 8 antigens")

# arrows classifier -> AUC pills
arrow(MODEL_X + MODEL_W + 0.3, MODEL_Y + MODEL_H/2,
      MODEL_X + MODEL_W + 5 - 0.3, MODEL_Y + MODEL_H/2,
      color=PDB_DK, lw=1.6, head=12)
arrow(MODEL_X + MODEL_W + 41 + 0.3, MODEL_Y + MODEL_H/2,
      MODEL_X + MODEL_W + 45 - 0.3, MODEL_Y + MODEL_H/2,
      color=CRES_DK, lw=1.6, head=12)

# Fig 4 tag - top-right corner of BOTH AUC pills (both pills are Fig 4 in the manuscript)
chip(PDB_PILL_X + 36 - 6, MODEL_Y + MODEL_H - 1.5, "Fig 4a",
     size=8.5, color="#333", fc="#FFFCEB", ec="#444", lw=1.0, pad_x=1.2, pad_y=0.5, zorder=7)
chip(CRES_PILL_X + 36 - 6, MODEL_Y + MODEL_H - 1.5, "Fig 4b",
     size=8.5, color="#333", fc="#FFFCEB", ec="#444", lw=1.0, pad_x=1.2, pad_y=0.5, zorder=7)

# =================================================================
# Figure callouts inside the cards
# =================================================================
# (Fig 1 bracket is drawn inside the PDB card itself)

# Fig 2-3 tag — sits ABOVE panel 5, inside the IEDB header strip area
# (positioned left of "TRAINING SOURCE" right-aligned label)
chip(panel_xs[4] + 5, IEDB_Y + IEDB_H - 2.5, "Fig 2–3",
     size=8.5, color=IEDB_DK, fc="#FFFCEB", ec=IEDB_DK, lw=1.0, pad_x=1.2, pad_y=0.5, zorder=6)

# footer
line(2, 3.5, 168, 3.5, color=HAIR, lw=0.5)
text(2, 1.5, "Woods et al — in revision",
     size=9, color=DIM, style="italic", ha="left")

# ---------------- save ----------------
from pathlib import Path
_here = Path(__file__).resolve().parent   # write next to this script
out_png = _here / "study_design_flowchart.png"
out_pdf = _here / "study_design_flowchart.pdf"
plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.25)
plt.savefig(out_pdf, bbox_inches="tight", facecolor="white", pad_inches=0.25)
print("saved:", out_png)
print("saved:", out_pdf)