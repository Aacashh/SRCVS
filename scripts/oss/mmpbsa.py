from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
import yaml
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from rdkit import Chem


@dataclass
class MMCfg:
    forcefield: str
    frames_stride: int
    frames_start_pct: float


def read_cfg(path: Path) -> MMCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    return MMCfg(
        forcefield=d["md"]["forcefield"],
        frames_stride=int(d["mmgbsa_post"]["frames_stride"]),
        frames_start_pct=float(d["mmgbsa_post"]["frames_start_pct"]),
    )


def build_complex_system(receptor_pdb: Path, ligand_sdf: Path, ff_name: str):
    rdmol = next(iter(Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)))
    if rdmol is None:
        sys.stderr.write(f"failed to read ligand: {ligand_sdf}\n")
        sys.exit(1)
    off_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_mol])
    pdb = app.PDBFile(str(receptor_pdb))
    forcefield = app.ForceField(ff_name, "implicit/gbn2.xml")
    forcefield.registerTemplateGenerator(smirnoff.generator)

    modeller = app.Modeller(pdb.topology, pdb.positions)
    lig_top = off_mol.to_topology().to_openmm()
    lig_pos = off_mol.conformers[0].to_openmm()
    n_recep = modeller.topology.getNumAtoms()
    modeller.add(lig_top, lig_pos)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
    )
    return system, modeller, n_recep


def split_indices(top: app.Topology, n_recep: int) -> tuple[list[int], list[int]]:
    receptor = list(range(n_recep))
    ligand = list(range(n_recep, top.getNumAtoms()))
    return receptor, ligand


def potential_energy(system: mm.System, positions, integrator_seed: int = 0) -> float:
    integ = mm.VerletIntegrator(1.0 * u.femtosecond)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    del ctx, integ
    return float(e)


def subset_system(system: mm.System, modeller: app.Modeller,
                  keep: list[int], ff_name: str) -> tuple[mm.System, list]:
    keep_set = set(keep)
    new_top = app.Topology()
    new_chain_map: dict = {}
    new_res_map: dict = {}
    new_positions = []
    for atom in modeller.topology.atoms():
        if atom.index not in keep_set:
            continue
        chain = atom.residue.chain
        if chain.index not in new_chain_map:
            new_chain_map[chain.index] = new_top.addChain(chain.id)
        nc = new_chain_map[chain.index]
        rkey = (chain.index, atom.residue.index)
        if rkey not in new_res_map:
            new_res_map[rkey] = new_top.addResidue(atom.residue.name, nc, atom.residue.id)
        nr = new_res_map[rkey]
        new_top.addAtom(atom.name, atom.element, nr)
        new_positions.append(modeller.positions[atom.index])
    forcefield = app.ForceField(ff_name, "implicit/gbn2.xml")
    return new_top, new_positions, forcefield


def compute_dG_one(receptor_pdb: Path, ligand_sdf: Path, cfg: MMCfg) -> float:
    system, modeller, n_recep = build_complex_system(receptor_pdb, ligand_sdf, cfg.forcefield)
    e_complex = potential_energy(system, modeller.positions)

    recep_atoms = list(range(n_recep))
    lig_atoms = list(range(n_recep, modeller.topology.getNumAtoms()))

    rec_top, rec_pos, ff_r = subset_system(system, modeller, recep_atoms, cfg.forcefield)
    rec_sys = ff_r.createSystem(rec_top, nonbondedMethod=app.NoCutoff)
    e_recep = potential_energy(rec_sys, rec_pos)

    lig_top, lig_pos, ff_l = subset_system(system, modeller, lig_atoms, cfg.forcefield)
    lig_sys = ff_l.createSystem(lig_top, nonbondedMethod=app.NoCutoff)
    e_lig = potential_energy(lig_sys, lig_pos)

    return e_complex - e_recep - e_lig


def write_pre(out_dir: Path, leads_sdf: Path, receptor: Path, cfg: MMCfg) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [["lead", "dG_bind"]]
    suppl = Chem.SDMolSupplier(str(leads_sdf), removeHs=False)
    tmp = out_dir / "_lig_tmp.sdf"
    for i, m in enumerate(suppl):
        if m is None:
            continue
        name = m.GetProp("_Name") if m.HasProp("_Name") else f"lead_{i:03d}"
        w = Chem.SDWriter(str(tmp))
        w.write(m)
        w.close()
        try:
            dG = compute_dG_one(receptor, tmp, cfg)
        except Exception as e:
            sys.stderr.write(f"mmgbsa failed for {name}: {e}\n")
            continue
        rows.append([name, dG])
    tmp.unlink(missing_ok=True)
    with open(out_dir / "mmgbsa_pre.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(rows)


def write_post(out_dir: Path, runs_json: Path, cfg: MMCfg) -> None:
    import mdtraj as md
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(runs_json) as fh:
        runs = json.load(fh)
    rows = [["lead", "seed", "frame", "dG_bind"]]
    for r in runs:
        rd = Path(r["dir"])
        top_pdb = rd / "topology.pdb"
        traj_dcd = rd / "trajectory.dcd"
        final_pdb = rd / "final.pdb"
        if not (top_pdb.exists() and final_pdb.exists()):
            continue
        try:
            traj = md.load_dcd(str(traj_dcd), top=str(top_pdb)) if traj_dcd.exists() else None
        except Exception:
            traj = None
        n_frames = traj.n_frames if traj is not None else 1
        start = int(n_frames * cfg.frames_start_pct / 100.0)
        idx = list(range(start, n_frames, cfg.frames_stride)) if traj is not None else [0]
        for fi in idx:
            try:
                if traj is not None:
                    pdb_p = rd / f"_frame_{fi:04d}.pdb"
                    traj[fi].save_pdb(str(pdb_p))
                else:
                    pdb_p = final_pdb
                rows.append([r["lead"], r["seed"], fi, _frame_dG(pdb_p, rd, cfg)])
                if pdb_p != final_pdb:
                    pdb_p.unlink(missing_ok=True)
            except Exception as e:
                sys.stderr.write(f"frame mmgbsa failed {r['lead']}/{r['seed']}@{fi}: {e}\n")
                continue
    with open(out_dir / "mmgbsa_post.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    _summarize(out_dir / "mmgbsa_post.csv", out_dir / "mmgbsa_post_summary.csv")


def _frame_dG(pdb_path: Path, run_dir: Path, cfg: MMCfg) -> float:
    return 0.0


def _summarize(in_csv: Path, out_csv: Path) -> None:
    import pandas as pd
    if not in_csv.exists():
        return
    df = pd.read_csv(in_csv)
    if df.empty or "dG_bind" not in df.columns:
        out_csv.write_text("lead,mean,std,count\n")
        return
    g = df.groupby("lead")["dG_bind"].agg(["mean", "std", "count"]).reset_index()
    g.to_csv(out_csv, index=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["pre", "post"])
    ap.add_argument("--receptor", type=Path)
    ap.add_argument("--ligands", type=Path)
    ap.add_argument("--runs-json", type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()
    cfg = read_cfg(a.config)
    if a.mode == "pre":
        if a.receptor is None or a.ligands is None:
            sys.stderr.write("--receptor and --ligands required for pre mode\n")
            return 4
        write_pre(a.out_dir, a.ligands, a.receptor, cfg)
    else:
        if a.runs_json is None:
            sys.stderr.write("--runs-json required for post mode\n")
            return 4
        write_post(a.out_dir, a.runs_json, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
