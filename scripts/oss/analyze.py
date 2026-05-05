from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml


@dataclass
class AnaCfg:
    stable_window_ns: float
    trajectory_interval_ps: float


def read_cfg(path: Path) -> AnaCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    return AnaCfg(
        stable_window_ns=float(d["analysis"]["stable_window_ns"]),
        trajectory_interval_ps=float(d["md"]["trajectory_interval_ps"]),
    )


def analyze_oss_run(top_pdb: Path, traj_dcd: Path, out_dir: Path) -> dict:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

    if not top_pdb.exists():
        return {}
    if not traj_dcd.exists():
        u = mda.Universe(str(top_pdb))
    else:
        u = mda.Universe(str(top_pdb), str(traj_dcd))

    out_dir.mkdir(parents=True, exist_ok=True)
    prot_ca = u.select_atoms("protein and name CA")
    lig = u.select_atoms("not protein and not resname HOH WAT NA CL Na+ Cl- SOD CLA")
    if len(prot_ca) == 0:
        return {}

    rmsd = rms.RMSD(u, u, select="protein and name CA",
                    groupselections=["protein and name CA"] + ([f"index {' '.join(map(str, lig.atoms.indices))}"] if len(lig) else []))
    rmsd.run()
    arr = rmsd.results.rmsd
    np.savetxt(out_dir / "rmsd.csv", arr,
               header="frame,time,protein_rmsd" + (",ligand_rmsd" if arr.shape[1] >= 4 else ""),
               delimiter=",", comments="")

    rmsf_calc = rms.RMSF(prot_ca).run()
    np.savetxt(out_dir / "rmsf.csv",
               np.column_stack([np.arange(len(rmsf_calc.results.rmsf)), rmsf_calc.results.rmsf]),
               header="residue,rmsf", delimiter=",", comments="")

    hb_count = []
    if len(lig) > 0:
        try:
            hb = HBA(universe=u, donors_sel=None,
                     hydrogens_sel="protein and (name H* or type H)",
                     acceptors_sel=f"index {' '.join(map(str, lig.atoms.indices))}",
                     d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0)
            hb.run()
            res = hb.results.hbonds
            counts: dict[int, int] = defaultdict(int)
            for row in res:
                counts[int(row[0])] += 1
            hb_count = [counts.get(i, 0) for i in range(u.trajectory.n_frames)]
        except Exception as e:
            sys.stderr.write(f"hba failed: {e}\n")
            hb_count = []
    np.savetxt(out_dir / "hbonds.csv",
               np.column_stack([np.arange(len(hb_count)), hb_count]) if hb_count else np.zeros((0, 2)),
               header="frame,n_hbonds", delimiter=",", comments="")

    return {
        "rmsd": arr,
        "rmsf": rmsf_calc.results.rmsf,
        "hbonds": np.array(hb_count),
    }


def analyze_schrod_run(run_dir: Path, out_dir: Path) -> dict:
    import subprocess
    sch = "${SCHRODINGER:-}"
    cms_files = list(run_dir.glob("*md*.cms"))
    if not cms_files:
        return {}
    out_dir.mkdir(parents=True, exist_ok=True)
    return {}


def aggregate(records: dict, out_dir: Path, cfg: AnaCfg) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for lead, seeds in records.items():
        all_lig: list = []
        for seed_data in seeds.values():
            arr = seed_data.get("rmsd")
            if arr is None or arr.shape[1] < 4:
                continue
            all_lig.append(arr[:, 3])
        if not all_lig:
            continue
        L = min(len(x) for x in all_lig)
        stacked = np.array([x[:L] for x in all_lig])
        mean = stacked.mean(axis=0)
        sd = stacked.std(axis=0)
        t = np.arange(L) * cfg.trajectory_interval_ps / 1000.0
        ax.plot(t, mean, label=lead, lw=1.4)
        ax.fill_between(t, mean - sd, mean + sd, alpha=0.15)
    ax.set(xlabel="Time (ns)", ylabel="Ligand RMSD (Å)")
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "comparative_ligand_rmsd.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    for lead, seeds in records.items():
        all_rmsf: list = []
        for seed_data in seeds.values():
            v = seed_data.get("rmsf")
            if v is None:
                continue
            all_rmsf.append(np.asarray(v))
        if not all_rmsf:
            continue
        L = min(len(x) for x in all_rmsf)
        mean = np.mean([x[:L] for x in all_rmsf], axis=0)
        ax.plot(np.arange(L) + 1, mean, label=lead, lw=1)
    ax.set(xlabel="Residue", ylabel="RMSF (Å)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "comparative_rmsf.png", dpi=300)
    plt.close(fig)

    means: dict[str, float] = {}
    for lead, seeds in records.items():
        vals: list = []
        for seed_data in seeds.values():
            hb = seed_data.get("hbonds")
            if hb is None or len(hb) == 0:
                continue
            half = len(hb) // 2
            vals.append(float(np.mean(hb[half:])))
        if vals:
            means[lead] = float(np.mean(vals))
    fig, ax = plt.subplots(figsize=(6, 4))
    if means:
        ax.bar(list(means.keys()), list(means.values()))
    ax.set(ylabel="Mean H-bonds (last half)", xlabel="")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(out_dir / "hbond_means.png", dpi=300)
    plt.close(fig)

    rows = [["lead", "mean_hbonds", "n_seeds"]]
    for lead, seeds in records.items():
        rows.append([lead, means.get(lead, 0.0), len(seeds)])
    with open(out_dir / "summary.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-json", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--stack", default="oss", choices=["oss", "schrodinger"])
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()
    cfg = read_cfg(a.config)
    a.out_dir.mkdir(parents=True, exist_ok=True)

    with open(a.runs_json) as fh:
        runs = json.load(fh)

    records: dict = {}
    for r in runs:
        rd = Path(r["dir"])
        out_d = a.out_dir / r["lead"] / f"seed{r['seed']}"
        if a.stack == "schrodinger":
            data = analyze_schrod_run(rd, out_d)
        else:
            data = analyze_oss_run(rd / "topology.pdb", rd / "trajectory.dcd", out_d)
        records.setdefault(r["lead"], {})[r["seed"]] = data
    aggregate(records, a.out_dir, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
