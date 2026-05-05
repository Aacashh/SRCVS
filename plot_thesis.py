"""
plot_thesis.py
==============
Generates publication-quality thesis plots from Phase 7-9 MD pipeline outputs.

CSV structure from MDAnalysis RMSD (5 columns, 4-name header → read by position):
  col0: frame_idx  col1: time_ps  col2: backbone_rmsd_A  col3: ca_rmsd_A  col4: ligand_rmsd_A

RMSF CSV: residue, rmsf (Angstroms)
MMGBSA CSV: lead, seed, frame, dG_bind (kcal/mol)
"""

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy.ndimage import uniform_filter1d

# ─── Aesthetics ───────────────────────────────────────────────────────────────
sns.set_theme(style="ticks", context="talk", font_scale=1.15)
plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,          # preview; bumped to 600 on final save
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.35,
    "grid.linestyle": "--",
})

COLORS  = {"44259_01": "#0d7377", "72307_01": "#e05c5c"}
LABELS  = {"44259_01": "Compound 44259_01", "72307_01": "Compound 72307_01"}
LEADS   = ["44259_01", "72307_01"]
ANA_DIR = Path("out/08_analysis")
MM_DIR  = Path("out/09_mmgbsa_post")
OUT_DIR = Path("out/thesis_plots")

# ─── Loaders ─────────────────────────────────────────────────────────────────
def load_rmsd(lead: str) -> pd.DataFrame:
    """
    MDAnalysis writes 5 columns but only 4 header names.
    We read with header=None and name the columns explicitly.
    Columns: frame_idx | time_ps | backbone_rmsd | ca_rmsd | ligand_rmsd
    """
    p = ANA_DIR / lead / "seed1" / "rmsd.csv"
    if not p.exists():
        return pd.DataFrame()
    # skip the malformed header row
    df = pd.read_csv(p, header=0, names=["frame_idx","time_ps","backbone_rmsd","ca_rmsd","ligand_rmsd"])
    df["time_ns"] = df["time_ps"] / 1000.0
    df["lead"]    = lead
    return df


def load_rmsf(lead: str) -> pd.DataFrame:
    p = ANA_DIR / lead / "seed1" / "rmsf.csv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p)
    df["lead"] = lead
    return df


def load_mmgbsa() -> pd.DataFrame:
    p = MM_DIR / "mmgbsa_post.csv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p)
    # map frame → time_ns  (frames_start_pct=50, stride=10, trajectory_interval_ps=10 ps)
    # frame number from the mmgbsa script is the trajectory frame index
    df["time_ns"] = df["frame"] * 10.0 / 1000.0
    return df


# ─── Smoothing helper ─────────────────────────────────────────────────────────
def smooth(y, w=9):
    return uniform_filter1d(y, size=w, mode="nearest")


# ─── Plot functions ───────────────────────────────────────────────────────────

def plot_ligand_rmsd_time(rmsd_dfs):
    fig, ax = plt.subplots(figsize=(11, 5))
    for lead, df in rmsd_dfs.items():
        t   = df["time_ns"].values
        raw = df["ligand_rmsd"].values
        smt = smooth(raw, w=15)
        ax.plot(t, raw, color=COLORS[lead], alpha=0.25, lw=0.8)
        ax.plot(t, smt, color=COLORS[lead], lw=2.2, label=LABELS[lead])
        ax.fill_between(t, smt, alpha=0.08, color=COLORS[lead])
    ax.set_xlabel("Simulation Time (ns)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("Ligand Positional Stability over 2 ns MD Simulation")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "01_ligand_rmsd_time.png")
    plt.close(fig)
    print("  ✓ 01_ligand_rmsd_time.png")


def plot_protein_rmsd_time(rmsd_dfs):
    fig, ax = plt.subplots(figsize=(11, 5))
    for lead, df in rmsd_dfs.items():
        t   = df["time_ns"].values
        raw = df["backbone_rmsd"].values
        smt = smooth(raw, w=15)
        ax.plot(t, raw, color=COLORS[lead], alpha=0.25, lw=0.8)
        ax.plot(t, smt, color=COLORS[lead], lw=2.2, label=LABELS[lead])
        ax.fill_between(t, smt, alpha=0.08, color=COLORS[lead])
    ax.set_xlabel("Simulation Time (ns)")
    ax.set_ylabel("Protein Backbone RMSD (Å)")
    ax.set_title("Receptor Backbone Structural Drift over 2 ns MD Simulation")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "02_protein_rmsd_time.png")
    plt.close(fig)
    print("  ✓ 02_protein_rmsd_time.png")


def plot_ligand_rmsd_dist(rmsd_dfs):
    fig, ax = plt.subplots(figsize=(8, 5))
    for lead, df in rmsd_dfs.items():
        # Use only the production half
        half = len(df) // 2
        vals = df["ligand_rmsd"].values[half:]
        sns.kdeplot(vals, ax=ax, color=COLORS[lead], fill=True, alpha=0.35,
                    linewidth=2.2, label=LABELS[lead])
    ax.set_xlabel("Ligand RMSD (Å)")
    ax.set_ylabel("Probability Density")
    ax.set_title("Ligand RMSD Distribution (Production Phase)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "03_ligand_rmsd_dist.png")
    plt.close(fig)
    print("  ✓ 03_ligand_rmsd_dist.png")


def plot_protein_rmsf(rmsf_dfs):
    fig, ax = plt.subplots(figsize=(13, 5))
    for lead, df in rmsf_dfs.items():
        ax.plot(df["residue"], df["rmsf"],
                color=COLORS[lead], lw=1.4, alpha=0.85, label=LABELS[lead])
    ax.set_xlabel("Residue Index")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Per-Residue Protein Flexibility (RMSF)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "04_protein_rmsf.png")
    plt.close(fig)
    print("  ✓ 04_protein_rmsf.png")


def plot_rmsf_difference(rmsf_dfs):
    if len(rmsf_dfs) < 2:
        return
    df44 = rmsf_dfs["44259_01"].set_index("residue")["rmsf"]
    df72 = rmsf_dfs["72307_01"].set_index("residue")["rmsf"]
    common = df44.index.intersection(df72.index)
    diff = df44[common] - df72[common]

    fig, ax = plt.subplots(figsize=(13, 5))
    pos = diff.clip(lower=0)
    neg = diff.clip(upper=0)
    ax.fill_between(diff.index, pos, 0, color=COLORS["44259_01"], alpha=0.6,
                    label="44259_01 more flexible")
    ax.fill_between(diff.index, neg, 0, color=COLORS["72307_01"], alpha=0.6,
                    label="72307_01 more flexible")
    ax.axhline(0, color="black", lw=0.8, linestyle="--")
    ax.set_xlabel("Residue Index")
    ax.set_ylabel("ΔRMSF (Å)  [44259_01 − 72307_01]")
    ax.set_title("Differential Per-Residue Flexibility")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "05_rmsf_difference.png")
    plt.close(fig)
    print("  ✓ 05_rmsf_difference.png")


def plot_mmgbsa_time(mm_df):
    fig, ax = plt.subplots(figsize=(10, 5))
    for lead in LEADS:
        d = mm_df[mm_df["lead"] == lead].sort_values("time_ns")
        ax.scatter(d["time_ns"], d["dG_bind"],
                   color=COLORS[lead], s=40, alpha=0.5, zorder=3)
        if len(d) >= 3:
            smt = smooth(d["dG_bind"].values, w=3)
            ax.plot(d["time_ns"], smt,
                    color=COLORS[lead], lw=2.5, label=LABELS[lead], zorder=4)
    ax.set_xlabel("Simulation Time (ns)")
    ax.set_ylabel("Binding Free Energy ΔG (kcal/mol)")
    ax.set_title("MM-GBSA Binding Free Energy over Production Trajectory")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "06_mmgbsa_energy_time.png")
    plt.close(fig)
    print("  ✓ 06_mmgbsa_energy_time.png")


def plot_mmgbsa_violin(mm_df):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.violinplot(data=mm_df, x="lead", y="dG_bind", hue="lead",
                   palette=COLORS, inner="quartile", linewidth=1.6,
                   legend=False, ax=ax)
    sns.swarmplot(data=mm_df, x="lead", y="dG_bind",
                  color="white", size=5, edgecolor="black", linewidth=0.8, ax=ax)

    # Annotate mean
    for lead in LEADS:
        mean_val = mm_df[mm_df["lead"] == lead]["dG_bind"].mean()
        x_pos = LEADS.index(lead)
        ax.text(x_pos, mean_val + 1.5, f"μ = {mean_val:.1f}", ha="center",
                fontsize=12, fontweight="bold", color=COLORS[lead])

    ax.set_xlabel("Compound")
    ax.set_ylabel("Binding Free Energy ΔG (kcal/mol)")
    ax.set_title("MM-GBSA Binding Free Energy Distribution")
    ax.set_xticklabels([LABELS[l] for l in LEADS])
    fig.tight_layout()
    fig.savefig(OUT_DIR / "07_mmgbsa_violin.png")
    plt.close(fig)
    print("  ✓ 07_mmgbsa_violin.png")


def plot_mmgbsa_boxplot(mm_df):
    fig, ax = plt.subplots(figsize=(7, 6))
    bp = ax.boxplot(
        [mm_df[mm_df["lead"] == lead]["dG_bind"].values for lead in LEADS],
        labels=[LABELS[l] for l in LEADS],
        patch_artist=True,
        medianprops=dict(color="white", linewidth=2.5),
        widths=0.5,
    )
    for patch, lead in zip(bp["boxes"], LEADS):
        patch.set_facecolor(COLORS[lead])
        patch.set_alpha(0.75)
    # jitter
    for i, lead in enumerate(LEADS):
        vals = mm_df[mm_df["lead"] == lead]["dG_bind"].values
        jitter = np.random.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   color=COLORS[lead], s=35, zorder=5, edgecolors="black", linewidths=0.5)
    ax.set_ylabel("Binding Free Energy ΔG (kcal/mol)")
    ax.set_title("MM-GBSA Energy — Box Plot Comparison")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "08_mmgbsa_boxplot.png")
    plt.close(fig)
    print("  ✓ 08_mmgbsa_boxplot.png")


def plot_combined_rmsd_panel(rmsd_dfs):
    """4-panel figure combining ligand + protein RMSD for both compounds."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharey="row")
    fig.suptitle("Molecular Dynamics RMSD Analysis", fontweight="bold", fontsize=16)

    row_labels = ["Ligand RMSD (Å)", "Protein Backbone RMSD (Å)"]
    rmsd_cols  = ["ligand_rmsd", "backbone_rmsd"]

    for col_i, lead in enumerate(LEADS):
        df = rmsd_dfs[lead]
        t  = df["time_ns"].values
        for row_i, (col, ylabel) in enumerate(zip(rmsd_cols, row_labels)):
            ax  = axes[row_i][col_i]
            raw = df[col].values
            smt = smooth(raw, w=15)
            ax.plot(t, raw, color=COLORS[lead], alpha=0.2, lw=0.8)
            ax.plot(t, smt, color=COLORS[lead], lw=2.2)
            ax.fill_between(t, smt, alpha=0.1, color=COLORS[lead])
            ax.set_xlabel("Time (ns)")
            if col_i == 0:
                ax.set_ylabel(ylabel)
            ax.set_title(LABELS[lead])

    fig.tight_layout()
    fig.savefig(OUT_DIR / "09_combined_rmsd_panel.png")
    plt.close(fig)
    print("  ✓ 09_combined_rmsd_panel.png")


def plot_energy_vs_ligand_rmsd(rmsd_dfs, mm_df):
    """Scatter: ΔG vs ligand RMSD per frame (if frames align)."""
    frames_merged = []
    for lead in LEADS:
        df_r = rmsd_dfs[lead][["time_ns", "ligand_rmsd"]].copy()
        df_r["time_ns"] = df_r["time_ns"].round(3)
        df_m = mm_df[mm_df["lead"] == lead][["time_ns", "dG_bind"]].copy()
        df_m["time_ns"] = df_m["time_ns"].round(3)
        merged = pd.merge(df_r, df_m, on="time_ns")
        merged["lead"] = lead
        frames_merged.append(merged)

    all_merged = pd.concat(frames_merged, ignore_index=True)
    if all_merged.empty:
        print("  ✗ 10_energy_vs_rmsd skipped — no overlapping frames")
        return

    fig, ax = plt.subplots(figsize=(9, 6))
    for lead in LEADS:
        d = all_merged[all_merged["lead"] == lead]
        ax.scatter(d["ligand_rmsd"], d["dG_bind"],
                   color=COLORS[lead], s=55, alpha=0.75, edgecolors="black",
                   linewidths=0.5, label=LABELS[lead])
        # Fit a linear trend
        if len(d) >= 3:
            z = np.polyfit(d["ligand_rmsd"], d["dG_bind"], 1)
            x_range = np.linspace(d["ligand_rmsd"].min(), d["ligand_rmsd"].max(), 80)
            ax.plot(x_range, np.polyval(z, x_range),
                    color=COLORS[lead], lw=1.8, linestyle="--", alpha=0.7)

    ax.set_xlabel("Ligand RMSD (Å)")
    ax.set_ylabel("Binding Free Energy ΔG (kcal/mol)")
    ax.set_title("Binding Affinity vs. Ligand Positional Stability")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "10_energy_vs_rmsd.png")
    plt.close(fig)
    print("  ✓ 10_energy_vs_rmsd.png")


# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.random.seed(42)

    rmsd_dfs = {l: load_rmsd(l) for l in LEADS}
    rmsf_dfs = {l: load_rmsf(l) for l in LEADS}
    mm_df    = load_mmgbsa()

    # Drop empty
    rmsd_dfs = {k: v for k, v in rmsd_dfs.items() if not v.empty}
    rmsf_dfs = {k: v for k, v in rmsf_dfs.items() if not v.empty}

    print(f"\nGenerating plots → {OUT_DIR}/\n")

    if rmsd_dfs:
        plot_ligand_rmsd_time(rmsd_dfs)
        plot_protein_rmsd_time(rmsd_dfs)
        plot_ligand_rmsd_dist(rmsd_dfs)
        plot_combined_rmsd_panel(rmsd_dfs)

    if rmsf_dfs:
        plot_protein_rmsf(rmsf_dfs)
        plot_rmsf_difference(rmsf_dfs)

    if not mm_df.empty:
        plot_mmgbsa_time(mm_df)
        plot_mmgbsa_violin(mm_df)
        plot_mmgbsa_boxplot(mm_df)

    if rmsd_dfs and not mm_df.empty:
        plot_energy_vs_ligand_rmsd(rmsd_dfs, mm_df)

    print(f"\nDone. All plots saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
