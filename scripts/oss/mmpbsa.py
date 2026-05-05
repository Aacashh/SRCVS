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

# MONKEY PATCH FOR FAST CHARGES
_orig_assign = Molecule.assign_partial_charges
def _fast_assign(self, partial_charge_method=None, **kwargs):
    _orig_assign(self, partial_charge_method="gasteiger", **kwargs)
Molecule.assign_partial_charges = _fast_assign


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
    return system, modeller, n_recep, smirnoff


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
                  keep: list[int], ff_name: str) -> tuple[app.Topology, list, app.ForceField]:
    keep_set = set(keep)
    new_top = app.Topology()
    new_chain_map: dict = {}
    new_res_map: dict = {}
    new_positions = []
    old_to_new_atom = {}
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
        new_atom = new_top.addAtom(atom.name, atom.element, nr)
        old_to_new_atom[atom] = new_atom
        new_positions.append(modeller.positions[atom.index])
    
    for bond in modeller.topology.bonds():
        if bond[0] in old_to_new_atom and bond[1] in old_to_new_atom:
            new_top.addBond(old_to_new_atom[bond[0]], old_to_new_atom[bond[1]], type=bond.type, order=bond.order)
            
    forcefield = app.ForceField(ff_name, "implicit/gbn2.xml")
    return new_top, new_positions, forcefield


def compute_dG_one(receptor_pdb: Path, ligand_sdf: Path, cfg: MMCfg) -> float:
    system, modeller, n_recep, smirnoff = build_complex_system(receptor_pdb, ligand_sdf, cfg.forcefield)
    e_complex = potential_energy(system, modeller.positions)

    recep_atoms = list(range(n_recep))
    lig_atoms = list(range(n_recep, modeller.topology.getNumAtoms()))

    rec_top, rec_pos, ff_r = subset_system(system, modeller, recep_atoms, cfg.forcefield)
    rec_sys = ff_r.createSystem(rec_top, nonbondedMethod=app.NoCutoff)
    e_recep = potential_energy(rec_sys, rec_pos)

    lig_top, lig_pos, ff_l = subset_system(system, modeller, lig_atoms, cfg.forcefield)
    ff_l.registerTemplateGenerator(smirnoff.generator)
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


_SOLVENT_RESNAMES = {"HOH", "WAT", "TIP", "TIP3", "NA", "CL", "Na+", "Cl-", "SOD", "CLA"}


def _find_lead_sdf(run_dir: Path, lead: str) -> Path | None:
    candidates = [
        run_dir.parent.parent / "_leads" / f"{lead}.sdf",
        run_dir.parent / f"{lead}.sdf",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def _build_run_systems(receptor_pdb: Path, ligand_sdf: Path, ff_name: str):
    sys_full, modeller, n_recep, smirnoff = build_complex_system(receptor_pdb, ligand_sdf, ff_name)
    n_total = modeller.topology.getNumAtoms()
    rec_atoms = list(range(n_recep))
    lig_atoms = list(range(n_recep, n_total))

    rec_top, _, ff_r = subset_system(sys_full, modeller, rec_atoms, ff_name)
    rec_sys = ff_r.createSystem(rec_top, nonbondedMethod=app.NoCutoff)

    lig_top, _, ff_l = subset_system(sys_full, modeller, lig_atoms, ff_name)
    ff_l.registerTemplateGenerator(smirnoff.generator)
    lig_sys = ff_l.createSystem(lig_top, nonbondedMethod=app.NoCutoff)

    plat = mm.Platform.getPlatformByName("CPU")
    ctx_full = mm.Context(sys_full, mm.VerletIntegrator(1.0 * u.femtosecond), plat)
    ctx_rec = mm.Context(rec_sys, mm.VerletIntegrator(1.0 * u.femtosecond), plat)
    ctx_lig = mm.Context(lig_sys, mm.VerletIntegrator(1.0 * u.femtosecond), plat)
    return ctx_full, ctx_rec, ctx_lig, n_recep, n_total


def _ctx_energy(ctx: mm.Context) -> float:
    return float(
        ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    )


def _frame_solute_positions(traj_frame, n_solute: int):
    xyz_nm = traj_frame.xyz[0, :n_solute, :]
    return [mm.Vec3(float(x), float(y), float(z)) for x, y, z in xyz_nm] * u.nanometer


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
        lead_sdf = _find_lead_sdf(rd, r["lead"])
        if lead_sdf is None:
            sys.stderr.write(f"missing lead sdf for {r['lead']}; skipping\n")
            continue

        receptor_only_pdb = rd / "_receptor_only.pdb"
        if not receptor_only_pdb.exists():
            _write_receptor_only(top_pdb, receptor_only_pdb)

        try:
            ctx_full, ctx_rec, ctx_lig, n_rec, n_tot = _build_run_systems(
                receptor_only_pdb, lead_sdf, cfg.forcefield
            )
        except Exception as e:
            sys.stderr.write(f"system build failed for {r['lead']}/{r['seed']}: {e}\n")
            continue

        try:
            top = md.load_topology(str(top_pdb))
            solute_idx = top.select(
                "not water and not resname NA and not resname CL "
                "and not resname 'Na+' and not resname 'Cl-' "
                "and not resname SOD and not resname CLA"
            )
            traj = md.load_dcd(str(traj_dcd), top=str(top_pdb), atom_indices=solute_idx) \
                if traj_dcd.exists() else None
        except Exception as e:
            sys.stderr.write(f"trajectory load failed for {r['lead']}/{r['seed']}: {e}\n")
            traj = None

        n_frames = traj.n_frames if traj is not None else 1
        start = int(n_frames * cfg.frames_start_pct / 100.0)
        idx = list(range(start, n_frames, cfg.frames_stride)) if traj is not None else [0]
        if not idx:
            idx = [n_frames - 1] if n_frames else [0]

        for fi in idx:
            try:
                if traj is not None:
                    full_pos = _frame_solute_positions(traj[fi], n_tot)
                else:
                    md_final = md.load_pdb(str(final_pdb))
                    sub = md_final.atom_slice(md_final.topology.select(
                        "not water and not resname NA and not resname CL"
                    ))
                    full_pos = _frame_solute_positions(sub, n_tot)
                rec_pos = full_pos[:n_rec]
                lig_pos = full_pos[n_rec:n_tot]
                ctx_full.setPositions(full_pos)
                ctx_rec.setPositions(rec_pos)
                ctx_lig.setPositions(lig_pos)
                dG = _ctx_energy(ctx_full) - _ctx_energy(ctx_rec) - _ctx_energy(ctx_lig)
                rows.append([r["lead"], r["seed"], fi, dG])
            except Exception as e:
                sys.stderr.write(f"frame mmgbsa failed {r['lead']}/{r['seed']}@{fi}: {e}\n")
                continue

        del ctx_full, ctx_rec, ctx_lig

    with open(out_dir / "mmgbsa_post.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    _summarize(out_dir / "mmgbsa_post.csv", out_dir / "mmgbsa_post_summary.csv")


def _write_receptor_only(topology_pdb: Path, out_pdb: Path) -> None:
    with open(topology_pdb) as fh, open(out_pdb, "w") as fw:
        for line in fh:
            rec = line[:6]
            if rec == "ATOM  ":
                resn = line[17:20].strip()
                if resn in _SOLVENT_RESNAMES:
                    continue
                fw.write(line)
            elif rec in ("TER   ", "END   ", "ENDMDL"):
                fw.write(line)
            elif rec.startswith(("CRYST", "MODEL", "HEADER", "TITLE", "REMARK")):
                fw.write(line)


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
