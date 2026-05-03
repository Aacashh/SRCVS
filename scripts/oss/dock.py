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
from rdkit import Chem


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


def receptor_to_pdbqt(receptor: Path, out_pdbqt: Path) -> None:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which("obabel") is None:
        sys.stderr.write("obabel not found in PATH\n")
        sys.exit(2)
    subprocess.run(
        ["obabel", str(receptor), "-O", str(out_pdbqt), "-xr", "--partialcharge", "gasteiger"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def dock_one(vina, ligand_pdbqt: Path, out_pdbqt: Path, cfg: DockCfg, exh: int) -> float | None:
    vina.set_ligand_from_file(str(ligand_pdbqt))
    vina.dock(exhaustiveness=exh, n_poses=cfg.num_modes)
    energies = vina.energies(n_poses=cfg.num_modes)
    if len(energies) == 0:
        return None
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    vina.write_poses(str(out_pdbqt), n_poses=cfg.num_modes, overwrite=True)
    return float(energies[0][0])


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def pose_pdbqt_to_sdf(pdbqt: Path, out_sdf: Path) -> None:
    if shutil.which("obabel") is None:
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

    receptor_pdbqt = a.out_dir / "receptor.pdbqt"
    if not receptor_pdbqt.exists():
        receptor_to_pdbqt(a.receptor, receptor_pdbqt)

    from vina import Vina
    v = Vina(sf_name="vina", cpu=cfg.cpu, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.compute_vina_maps(center=grid["center"], box_size=grid["outer"])

    ligand_files = sorted(a.pdbqt_dir.glob("*.pdbqt"))
    if not ligand_files:
        sys.stderr.write(f"no ligand pdbqt files in {a.pdbqt_dir}\n")
        return 3

    htvs_rows: list[list] = []
    for lp in ligand_files:
        out_p = poses_dir / "_htvs" / f"{lp.stem}.pdbqt"
        score = None
        try:
            score = dock_one(v, lp, out_p, cfg, exh=max(4, cfg.exhaustiveness // 2))
        except Exception as e:
            sys.stderr.write(f"htvs dock failed for {lp.name}: {e}\n")
        if score is not None:
            htvs_rows.append([lp.stem, score, str(out_p)])

    htvs_rows.sort(key=lambda r: r[1])
    write_csv(a.out_dir / "htvs_scores.csv", ["title", "score", "pose"], htvs_rows)

    n_sp = max(1, int(len(htvs_rows) * cfg.htvs_keep_pct / 100.0))
    sp_inputs = htvs_rows[:n_sp]
    sp_rows: list[list] = []
    for title, _, _ in sp_inputs:
        lp = a.pdbqt_dir / f"{title}.pdbqt"
        out_p = poses_dir / "_sp" / f"{title}.pdbqt"
        try:
            score = dock_one(v, lp, out_p, cfg, exh=cfg.exhaustiveness)
        except Exception as e:
            sys.stderr.write(f"sp dock failed for {title}: {e}\n")
            continue
        if score is not None:
            sp_rows.append([title, score, str(out_p)])
    sp_rows.sort(key=lambda r: r[1])
    write_csv(a.out_dir / "sp_scores.csv", ["title", "score", "pose"], sp_rows)

    n_xp = max(1, int(len(sp_rows) * cfg.sp_keep_pct / 100.0))
    xp_inputs = sp_rows[:max(n_xp, cfg.xp_top_n)] if cfg.xp_top_n else sp_rows[:n_xp]
    xp_rows: list[list] = []
    for title, _, _ in xp_inputs:
        lp = a.pdbqt_dir / f"{title}.pdbqt"
        out_p = poses_dir / f"{title}.pdbqt"
        try:
            score = dock_one(v, lp, out_p, cfg, exh=cfg.exhaustiveness * 2)
        except Exception as e:
            sys.stderr.write(f"xp dock failed for {title}: {e}\n")
            continue
        if score is not None:
            sdf_p = poses_dir / f"{title}.sdf"
            pose_pdbqt_to_sdf(out_p, sdf_p)
            xp_rows.append([title, score, str(out_p), str(sdf_p)])
    xp_rows.sort(key=lambda r: r[1])
    write_csv(a.out_dir / "xp_scores.csv", ["title", "score", "pose", "sdf"], xp_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
