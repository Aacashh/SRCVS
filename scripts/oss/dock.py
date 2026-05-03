from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class DockCfg:
    htvs_keep_pct: float
    sp_keep_pct: float
    xp_top_n: int
    exhaustiveness: int
    num_modes: int
    cpu: int


def read_cfg(path: Path) -> DockCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    dc = d["docking"]
    return DockCfg(
        htvs_keep_pct=float(dc["htvs_keep_pct"]),
        sp_keep_pct=float(dc["sp_keep_pct"]),
        xp_top_n=int(dc["xp_top_n"]),
        exhaustiveness=int(dc["exhaustiveness"]),
        num_modes=int(dc["num_modes"]),
        cpu=int(dc["cpu"]),
    )


def find_docking_binary() -> str:
    for name in ("vina", "smina", "qvina2", "qvina-w"):
        p = shutil.which(name)
        if p:
            return p
    sys.stderr.write("no docking binary found in PATH (install via: sudo apt install -y autodock-vina)\n")
    sys.exit(2)


def receptor_to_pdbqt(receptor: Path, out_pdbqt: Path) -> None:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which("obabel") is None:
        sys.stderr.write("obabel not found in PATH\n")
        sys.exit(2)
    subprocess.run(
        ["obabel", str(receptor), "-O", str(out_pdbqt), "-xr", "--partialcharge", "gasteiger"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def parse_score(pdbqt_path: Path) -> float | None:
    if not pdbqt_path.exists():
        return None
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        return float(parts[3])
                    except ValueError:
                        return None
            if line.startswith("REMARK minimizedAffinity"):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        return float(parts[2])
                    except ValueError:
                        return None
    return None


def dock_one(binary: str, receptor_pdbqt: Path, ligand_pdbqt: Path, out_pdbqt: Path,
             center: list, size: list, cfg: DockCfg, exh: int) -> float | None:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        binary,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", str(exh),
        "--num_modes", str(cfg.num_modes),
        "--out", str(out_pdbqt),
        "--cpu", "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return None
    return parse_score(out_pdbqt)


def write_csv(path: Path, header: list, rows: list) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def pose_pdbqt_to_sdf(pdbqt: Path, out_sdf: Path) -> None:
    if shutil.which("obabel") is None or not pdbqt.exists():
        return
    subprocess.run(
        ["obabel", str(pdbqt), "-O", str(out_sdf), "-l", "1"],
        check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, type=Path)
    ap.add_argument("--grid", required=True, type=Path)
    ap.add_argument("--ligands", required=True, type=Path)
    ap.add_argument("--pdbqt-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()

    cfg = read_cfg(a.config)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    poses_dir = a.out_dir / "poses"
    poses_dir.mkdir(exist_ok=True)

    with open(a.grid) as fh:
        grid = json.load(fh)
    center = grid["center"]
    size = grid["outer"]

    binary = find_docking_binary()

    receptor_pdbqt = a.out_dir / "receptor.pdbqt"
    if not receptor_pdbqt.exists():
        receptor_to_pdbqt(a.receptor, receptor_pdbqt)

    ligand_files = sorted(a.pdbqt_dir.glob("*.pdbqt"))
    if not ligand_files:
        sys.stderr.write(f"no ligand pdbqt files in {a.pdbqt_dir}\n")
        return 3

    htvs_rows: list = []
    htvs_dir = poses_dir / "_htvs"
    for lp in ligand_files:
        out_p = htvs_dir / f"{lp.stem}.pdbqt"
        score = dock_one(binary, receptor_pdbqt, lp, out_p, center, size, cfg,
                         exh=max(4, cfg.exhaustiveness // 2))
        if score is not None:
            htvs_rows.append([lp.stem, score, str(out_p)])
    htvs_rows.sort(key=lambda r: r[1])
    write_csv(a.out_dir / "htvs_scores.csv", ["title", "score", "pose"], htvs_rows)

    n_sp = max(1, int(len(htvs_rows) * cfg.htvs_keep_pct / 100.0))
    sp_inputs = htvs_rows[:n_sp]
    sp_rows: list = []
    sp_dir = poses_dir / "_sp"
    for title, _, _ in sp_inputs:
        lp = a.pdbqt_dir / f"{title}.pdbqt"
        out_p = sp_dir / f"{title}.pdbqt"
        score = dock_one(binary, receptor_pdbqt, lp, out_p, center, size, cfg,
                         exh=cfg.exhaustiveness)
        if score is not None:
            sp_rows.append([title, score, str(out_p)])
    sp_rows.sort(key=lambda r: r[1])
    write_csv(a.out_dir / "sp_scores.csv", ["title", "score", "pose"], sp_rows)

    n_xp_pct = max(1, int(len(sp_rows) * cfg.sp_keep_pct / 100.0))
    n_xp = max(n_xp_pct, cfg.xp_top_n) if cfg.xp_top_n else n_xp_pct
    xp_inputs = sp_rows[:n_xp]
    xp_rows: list = []
    for title, _, _ in xp_inputs:
        lp = a.pdbqt_dir / f"{title}.pdbqt"
        out_p = poses_dir / f"{title}.pdbqt"
        score = dock_one(binary, receptor_pdbqt, lp, out_p, center, size, cfg,
                         exh=cfg.exhaustiveness * 2)
        if score is not None:
            sdf_p = poses_dir / f"{title}.sdf"
            pose_pdbqt_to_sdf(out_p, sdf_p)
            xp_rows.append([title, score, str(out_p), str(sdf_p)])
    xp_rows.sort(key=lambda r: r[1])
    write_csv(a.out_dir / "xp_scores.csv", ["title", "score", "pose", "sdf"], xp_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
