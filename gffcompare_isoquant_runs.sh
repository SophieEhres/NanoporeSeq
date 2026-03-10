#!/usr/bin/env bash
# compare_gtfs_parallel.sh
# Compares multiple GTF files against a reference using gffcompare,
# outputs a presence/absence matrix, summary table, and publication-quality plots.
#
# Usage: bash compare_gtfs_parallel.sh -r reference.gff3 [-o output.tsv] [-t threads] run1.gtf run2.gtf ... runN.gtf

set -euo pipefail

usage() {
    echo ""
    echo "Usage: $0 -r <reference.gff3> [-o output.tsv] [-t threads] run1.gtf run2.gtf ... runN.gtf"
    echo ""
    echo "  -r   Reference annotation file (GFF3 or GTF)"
    echo "  -o   Output prefix (default: transcript_comparison)"
    echo "  -t   Number of parallel threads (default: 4)"
    echo ""
    echo "Example:"
    echo "  $0 -r gencode.gff3 -o results -t 8 run1.gtf run2.gtf run3.gtf"
    echo ""
    exit 1
}

REFERENCE=""
OUTPUT_PREFIX="transcript_comparison"
THREADS=4

while getopts "r:o:t:h" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

GTF_FILES=("$@")
N=${#GTF_FILES[@]}

if [[ -z "$REFERENCE" || $N -eq 0 ]]; then
    echo "ERROR: Reference and at least one GTF file are required."
    usage
fi

echo "================================================"
echo " GTF Comparison Pipeline"
echo "================================================"
echo "Reference:  $REFERENCE"
echo "GTF files:  $N"
echo "Threads:    $THREADS"
echo "Output:     ${OUTPUT_PREFIX}.*"
echo "================================================"
echo ""

TMPDIR="${OUTPUT_PREFIX}_tmp"
mkdir -p "$TMPDIR"
echo "      Temp dir: $TMPDIR"
echo ""

# --- Step 1: Run gffcompare in parallel ---
echo "[1/4] Running gffcompare in parallel..."

run_gffcompare() {
    i=$1
    gtf=$2
    ref=$3
    tmpdir=$4
    prefix="${tmpdir}/gffcmp_${i}"
    gffcompare -r "$ref" -o "$prefix" "$gtf" 2>/dev/null
    echo "      Done: $(basename $gtf)"
}
export -f run_gffcompare

for i in "${!GTF_FILES[@]}"; do
    echo "$i ${GTF_FILES[$i]}"
done | xargs -P "$THREADS" -n 2 bash -c 'run_gffcompare "$1" "$2" "'"$REFERENCE"'" "'"$TMPDIR"'"' _

echo ""

# --- Step 2: Collect tmap files and labels ---
echo "[2/4] Collecting results..."

declare -a LABELS
declare -a TMAP_FILES
declare -a STATS_FILES

for i in "${!GTF_FILES[@]}"; do
    gtf="${GTF_FILES[$i]}"
    label=$(basename "$gtf" .gtf)
    label=$(basename "$label" .transcript_models)
    LABELS[$i]="$label"

    tmap=$(find "$TMPDIR" -name "gffcmp_${i}.*.tmap" | head -1)
    stats="${TMPDIR}/gffcmp_${i}.stats"
    TMAP_FILES[$i]="$tmap"
    STATS_FILES[$i]="$stats"

    count=$(awk 'NR>1 && !seen[$5]++ {c++} END {print c+0}' "$tmap")
    echo "      ${label}: $count unique transcripts"
done

echo ""

# --- Step 3: Extract transcript sets in parallel ---
echo "[3/4] Extracting transcript sets..."

extract_set() {
    i=$1
    tmap=$2
    tmpdir=$3
    # Single pass: deduplicate transcript IDs (col 5) without sort, extract class codes (col 4)
    awk 'NR > 1 {
        print $4 >> "'"${TMPDIR}"'/codes_" i ".txt"
        if (!seen[$5]++) print $5
    }' i="$i" "$tmap" > "${tmpdir}/set_${i}.txt"
}
export -f extract_set

for i in "${!TMAP_FILES[@]}"; do
    echo "$i ${TMAP_FILES[$i]}"
done | xargs -P "$THREADS" -n 2 bash -c 'extract_set "$1" "$2" "'"$TMPDIR"'"' _

echo ""

# --- Step 4: Build matrix + plots using Python ---
echo "[4/4] Building matrix and generating plots..."

printf '%s\n' "${LABELS[@]}" > "$TMPDIR/labels.txt"
printf '%s\n' "${STATS_FILES[@]}" > "$TMPDIR/stats_files.txt"

python3 - "$TMPDIR" "$N" "$OUTPUT_PREFIX" << 'PYEOF'
import sys
import os
import re
import warnings
warnings.filterwarnings("ignore")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("      WARNING: matplotlib not found. Skipping plots. Install with: pip install matplotlib")

tmpdir  = sys.argv[1]
n       = int(sys.argv[2])
prefix  = sys.argv[3]

with open(os.path.join(tmpdir, "labels.txt")) as f:
    labels = [l.strip() for l in f]

with open(os.path.join(tmpdir, "stats_files.txt")) as f:
    stats_files = [l.strip() for l in f]

# ── Load transcript sets ──────────────────────────────────────────────────────
sets = []
for i in range(n):
    with open(os.path.join(tmpdir, f"set_{i}.txt")) as f:
        sets.append(set(l.strip() for l in f if l.strip()))

union  = sorted(set().union(*sets))
total  = len(union)
print(f"      Total unique transcripts: {total}")

# ── Bitmask counts (used for summary + upset) ─────────────────────────────────
from collections import Counter
bitmask_counts = Counter()
for tx in union:
    mask = sum((1 << i) for i, s in enumerate(sets) if tx in s)
    bitmask_counts[mask] += 1

# ── Write presence/absence matrix ─────────────────────────────────────────────
matrix_path = f"{prefix}.tsv"
with open(matrix_path, "w") as out:
    out.write("transcript_id\t" + "\t".join(labels) + "\tnum_runs_found\tfound_in\n")
    for tx in union:
        presence  = [1 if tx in s else 0 for s in sets]
        count     = sum(presence)
        found_in  = ",".join(labels[i] for i, p in enumerate(presence) if p)
        out.write(tx + "\t" + "\t".join(map(str, presence)) + f"\t{count}\t{found_in}\n")

# ── Write summary table ───────────────────────────────────────────────────────
summary_path = f"{prefix}_summary.tsv"
summary_counts = []
with open(summary_path, "w") as out:
    out.write("num_runs_found\ttranscript_count\tpercent_of_total\n")
    for k in range(1, n + 1):
        # sum counts of all masks that have exactly k bits set
        cnt = sum(v for mask, v in bitmask_counts.items() if bin(mask).count("1") == k)
        pct = f"{cnt/total*100:.1f}" if total > 0 else "0"
        summary_counts.append(cnt)
        out.write(f"{k}\t{cnt}\t{pct}\n")

print(f"      Matrix  → {matrix_path}")
print(f"      Summary → {summary_path}")
print()
print(f"  {'Runs found':<15} {'Count':>10} {'Percent':>10}")
print(f"  {'-'*37}")
for k, cnt in enumerate(summary_counts, 1):
    pct = cnt/total*100 if total > 0 else 0
    label_k = f"{'all' if k==n else str(k)} / {n}"
    print(f"  {label_k:<15} {cnt:>10,} {pct:>9.1f}%")

if not HAS_MPL:
    sys.exit(0)

# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

# ── Colour palette ────────────────────────────────────────────────────────────
BG      = "#0d1117"
PANEL   = "#161b22"
BORDER  = "#30363d"
TEXT    = "#e6edf3"
MUTED   = "#8b949e"
ACCENT  = "#58a6ff"
PALETTE = ["#58a6ff","#3fb950","#f78166","#d2a8ff","#ffa657",
           "#79c0ff","#56d364","#ff7b72","#bc8cff","#ffb86c"]

plt.rcParams.update({
    "figure.facecolor":  BG,
    "axes.facecolor":    PANEL,
    "axes.edgecolor":    BORDER,
    "axes.labelcolor":   TEXT,
    "axes.titlecolor":   TEXT,
    "xtick.color":       MUTED,
    "ytick.color":       MUTED,
    "text.color":        TEXT,
    "grid.color":        BORDER,
    "grid.linewidth":    0.6,
    "font.family":       "DejaVu Sans",
    "font.size":         10,
})

def savefig(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=BG)
    print(f"      Plot    → {path}")

# ── 1. Summary bar chart: transcripts found in k/N runs ───────────────────────
fig, ax = plt.subplots(figsize=(max(6, n * 1.2), 5))
fig.patch.set_facecolor(BG)

x      = np.arange(1, n + 1)
colors = [PALETTE[i % len(PALETTE)] for i in range(n)]
bars   = ax.bar(x, summary_counts, color=colors, width=0.6,
                edgecolor=BG, linewidth=1.5, zorder=3)

for bar, cnt in zip(bars, summary_counts):
    pct = cnt / total * 100 if total > 0 else 0
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + max(summary_counts)*0.01,
            f"{cnt:,}\n({pct:.1f}%)",
            ha="center", va="bottom", fontsize=8, color=TEXT)

ax.set_xticks(x)
ax.set_xticklabels([f"{k}/{n}" for k in x], fontsize=9)
ax.set_xlabel("Number of runs transcript was found in", labelpad=8)
ax.set_ylabel("Number of transcripts", labelpad=8)
ax.set_title("Transcript Support Across Runs", fontsize=13, fontweight="bold", pad=12)
ax.yaxis.grid(True, zorder=0)
ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)
savefig(fig, f"{prefix}_support_barplot.png")
plt.close(fig)

# ── 2. Per-run transcript count bar chart ─────────────────────────────────────
run_counts = [len(s) for s in sets]

fig, ax = plt.subplots(figsize=(max(6, n * 1.4), 5))
fig.patch.set_facecolor(BG)

bars = ax.bar(range(n), run_counts,
              color=[PALETTE[i % len(PALETTE)] for i in range(n)],
              width=0.6, edgecolor=BG, linewidth=1.5, zorder=3)

for bar, cnt in zip(bars, run_counts):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + max(run_counts)*0.01,
            f"{cnt:,}", ha="center", va="bottom", fontsize=9, color=TEXT)

ax.set_xticks(range(n))
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Number of transcripts", labelpad=8)
ax.set_title("Transcripts per Run", fontsize=13, fontweight="bold", pad=12)
ax.yaxis.grid(True, zorder=0)
ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)
savefig(fig, f"{prefix}_per_run_counts.png")
plt.close(fig)

# ── 3. UpSet-style plot ───────────────────────────────────────────────────────
from itertools import combinations

# Reuse bitmask_counts computed above
combo_counts = {}
for mask, cnt in bitmask_counts.items():
    combo = tuple(i for i in range(n) if mask & (1 << i))
    if combo:
        combo_counts[combo] = cnt

# Sort by count descending, take top 20
top_combos = sorted(combo_counts.items(), key=lambda x: -x[1])[:20]
combo_labels_list = [c for c, _ in top_combos]
combo_vals        = [v for _, v in top_combos]

fig_h = max(5, n * 0.5 + 3)
fig   = plt.figure(figsize=(max(10, len(top_combos) * 0.7), fig_h + 2))
fig.patch.set_facecolor(BG)

gs  = gridspec.GridSpec(2, 1, height_ratios=[2, 1], hspace=0.05)
ax_bar = fig.add_subplot(gs[0])
ax_dot = fig.add_subplot(gs[1], sharex=ax_bar)

# Bar chart (top)
bar_colors = [PALETTE[0] if len(c)==1 else PALETTE[min(len(c)-1, len(PALETTE)-1)]
              for c in combo_labels_list]
bars = ax_bar.bar(range(len(top_combos)), combo_vals,
                  color=bar_colors, width=0.6, edgecolor=BG, zorder=3)
for bar, cnt in zip(bars, combo_vals):
    ax_bar.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + max(combo_vals)*0.01,
                f"{cnt:,}", ha="center", va="bottom", fontsize=7, color=TEXT)

ax_bar.set_ylabel("Transcripts", labelpad=6)
ax_bar.set_title("Intersection Sizes (Top 20)", fontsize=13, fontweight="bold", pad=12)
ax_bar.yaxis.grid(True, zorder=0)
ax_bar.set_axisbelow(True)
ax_bar.spines[["top","right"]].set_visible(False)
plt.setp(ax_bar.get_xticklabels(), visible=False)

# Dot matrix (bottom)
ax_dot.set_facecolor(PANEL)
ax_dot.set_ylim(-0.5, n - 0.5)
ax_dot.set_yticks(range(n))
ax_dot.set_yticklabels(labels, fontsize=8)
ax_dot.spines[["top","right","bottom"]].set_visible(False)
ax_dot.xaxis.set_visible(False)

for col_idx, (combo, cnt) in enumerate(top_combos):
    combo_set = set(combo)
    filled = [i for i in range(n) if i in combo_set]
    empty  = [i for i in range(n) if i not in combo_set]
    # empty dots
    for row in empty:
        ax_dot.scatter(col_idx, row, color=BORDER, s=60, zorder=3)
    # filled dots + connecting line
    if filled:
        ax_dot.plot([col_idx, col_idx], [min(filled), max(filled)],
                    color=ACCENT, linewidth=2, zorder=2)
        for row in filled:
            ax_dot.scatter(col_idx, row, color=ACCENT, s=80, zorder=4)

savefig(fig, f"{prefix}_upset.png")
plt.close(fig)

# ── 4. Class code distribution (stacked bar) ──────────────────────────────────
CLASS_LABELS = {
    "=": "Exact match",
    "c": "Contained",
    "j": "Novel isoform",
    "e": "Single exon",
    "i": "Intronic",
    "o": "Generic overlap",
    "p": "Polymerase run",
    "r": "Repeat",
    "u": "Intergenic",
    "x": "Antisense",
    "s": "Antisense intron",
    "k": "Containment",
}
all_codes = list(CLASS_LABELS.keys()) + ["other"]
code_palette = {
    "=": "#3fb950", "c": "#58a6ff", "j": "#d2a8ff", "e": "#ffa657",
    "i": "#ff7b72", "o": "#79c0ff", "p": "#e3b341", "r": "#8b949e",
    "u": "#f78166", "x": "#bc8cff", "s": "#56d364", "k": "#ffb86c",
    "other": "#6e7681"
}

code_data = {}
for i in range(n):
    code_file = os.path.join(tmpdir, f"codes_{i}.txt")
    counts = {}
    with open(code_file) as f:
        for line in f:
            code = line.strip()
            counts[code] = counts.get(code, 0) + 1
    code_data[labels[i]] = counts

present_codes = sorted(
    set(c for d in code_data.values() for c in d),
    key=lambda c: -sum(d.get(c, 0) for d in code_data.values())
)

fig, ax = plt.subplots(figsize=(max(7, n * 1.4), 5))
fig.patch.set_facecolor(BG)

bottoms = np.zeros(n)
x       = np.arange(n)

for code in present_codes:
    vals   = np.array([code_data[labels[i]].get(code, 0) for i in range(n)], dtype=float)
    color  = code_palette.get(code, "#6e7681")
    ax.bar(x, vals, bottom=bottoms, color=color, width=0.6,
           edgecolor=BG, linewidth=0.8, label=CLASS_LABELS.get(code, code), zorder=3)
    bottoms += vals

ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Number of transcripts", labelpad=8)
ax.set_title("Transcript Class Code Distribution per Run", fontsize=13, fontweight="bold", pad=12)
ax.legend(loc="upper right", fontsize=7, framealpha=0.2,
          facecolor=PANEL, edgecolor=BORDER, ncol=2)
ax.yaxis.grid(True, zorder=0)
ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)
savefig(fig, f"{prefix}_class_codes.png")
plt.close(fig)

# ── 5. Parse .stats files and plot Sensitivity / Precision ────────────────────
def parse_stats(filepath):
    metrics = {}
    if not os.path.exists(filepath):
        return metrics
    with open(filepath) as f:
        for line in f:
            m = re.search(r'([\w\s]+level):\s+Sn:\s+([\d.]+)\s+\|\s+Pr:\s+([\d.]+)', line)
            if m:
                level = m.group(1).strip()
                metrics[level] = {
                    "Sensitivity": float(m.group(2)),
                    "Precision":   float(m.group(3))
                }
    return metrics

all_stats = [parse_stats(sf) for sf in stats_files]
levels    = list(all_stats[0].keys()) if all_stats[0] else []

if levels:
    fig, axes = plt.subplots(1, len(levels),
                             figsize=(4.5 * len(levels), 5),
                             sharey=True)
    if len(levels) == 1:
        axes = [axes]
    fig.patch.set_facecolor(BG)

    metric_names = ["Sensitivity", "Precision"]
    bar_w        = 0.75 / n

    for ax, level in zip(axes, levels):
        ax.set_facecolor(PANEL)
        for spine in ax.spines.values():
            spine.set_edgecolor(BORDER)

        for i, (stats, lbl) in enumerate(zip(all_stats, labels)):
            vals   = [stats.get(level, {}).get(m, 0) for m in metric_names]
            xpos   = np.arange(len(metric_names))
            offset = (i - n/2 + 0.5) * bar_w
            bars   = ax.bar(xpos + offset, vals, bar_w,
                            label=lbl,
                            color=PALETTE[i % len(PALETTE)],
                            edgecolor=BG, linewidth=1, zorder=3)
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 0.5,
                        f"{val:.1f}", ha="center", va="bottom", fontsize=7, color=TEXT)

        ax.set_title(level, fontweight="bold", pad=8)
        ax.set_xticks(np.arange(len(metric_names)))
        ax.set_xticklabels(metric_names, fontsize=9)
        ax.set_ylim(0, 108)
        ax.set_ylabel("%" if ax == axes[0] else "", labelpad=6)
        ax.legend(fontsize=7, framealpha=0.2,
                  facecolor=PANEL, edgecolor=BORDER)
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)

    fig.suptitle("Sensitivity & Precision vs Reference",
                 fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    savefig(fig, f"{prefix}_sn_pr.png")
    plt.close(fig)

print()
print("  All done!")
PYEOF

echo ""
echo "================================================"
echo " Output files"
echo "================================================"
echo "  ${OUTPUT_PREFIX}.tsv                  — presence/absence matrix"
echo "  ${OUTPUT_PREFIX}_summary.tsv           — transcripts found in k/N runs"
echo "  ${OUTPUT_PREFIX}_support_barplot.png   — bar chart: support across runs"
echo "  ${OUTPUT_PREFIX}_per_run_counts.png    — bar chart: transcripts per run"
echo "  ${OUTPUT_PREFIX}_upset.png             — UpSet intersection plot"
echo "  ${OUTPUT_PREFIX}_class_codes.png       — class code distribution"
echo "  ${OUTPUT_PREFIX}_sn_pr.png             — sensitivity & precision"
echo "================================================"
