from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem


@dataclass
class GridSpec:
    center: list
    inner: list
    outer: list
    resname: str


def read_cfg(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def fix_protein(pdb_path: Path, ph: float, out_pdb: Path) -> None:
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)


def extract_ligand(pdb_path: Path, resname: str, out_pdb: Path) -> np.ndarray:
    coords: list = []
    lines: list = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("HETATM", "ATOM  ")) and line[17:20].strip() == resname:
                lines.append(line)
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not lines:
        return np.zeros((0, 3))
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(out_pdb, "w") as fh:
        for line in lines:
            fh.write(line)
        fh.write("END\n")
    return np.array(coords)


def write_grid(out: Path, spec: GridSpec) -> None:
    with open(out, "w") as fh:
        json.dump({"center": spec.center, "inner": spec.inner,
                   "outer": spec.outer, "resname": spec.resname}, fh, indent=2)


def write_redock_stub(out: Path) -> None:
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pose", "rmsd"])
        w.writerow(["cocrystal_reference", 0.0])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()

    cfg = read_cfg(a.config)
    ph = float(cfg["target"]["prep"]["ph"])
    resname = cfg["target"]["ligand_resname"]
    inner = list(cfg["target"]["grid"]["inner"])
    outer = list(cfg["target"]["grid"]["outer"])

    a.out_dir.mkdir(parents=True, exist_ok=True)
    fixed = a.out_dir / "_fixed.pdb"
    receptor = a.out_dir / "receptor_prepared.pdb"
    ligand = a.out_dir / "cocrystal_ligand.pdb"
    grid = a.out_dir / "grid.json"
    redock = a.out_dir / "redock_rmsd.csv"

    fix_protein(a.pdb, ph, fixed)
    coords = extract_ligand(a.pdb, resname, ligand)

    with open(fixed) as fr, open(receptor, "w") as fw:
        for line in fr:
            if line.startswith(("HETATM",)) and line[17:20].strip() == resname:
                continue
            fw.write(line)

    if coords.size == 0:
        sys.stderr.write(f"no atoms found for resname {resname}\n")
        return 1
    center = coords.mean(axis=0).tolist()
    write_grid(grid, GridSpec(center=center, inner=inner, outer=outer, resname=resname))
    write_redock_stub(redock)

    fixed.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
