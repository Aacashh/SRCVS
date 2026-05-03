from __future__ import annotations

import argparse
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
class MDCfg:
    box_buffer_nm: float
    salt_M: float
    temperature_K: float
    pressure_bar: float
    production_ns: float
    equilibration_ns: float
    timestep_fs: float
    trajectory_interval_ps: float
    forcefield: str
    water_model: str


def read_cfg(path: Path) -> MDCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    m = d["md"]
    return MDCfg(
        box_buffer_nm=float(m["box_buffer_nm"]),
        salt_M=float(m["ion_concentration_M"]),
        temperature_K=float(m["temperature_K"]),
        pressure_bar=float(m["pressure_bar"]),
        production_ns=float(m["production_ns"]),
        equilibration_ns=float(m["equilibration_ns"]),
        timestep_fs=float(m["timestep_fs"]),
        trajectory_interval_ps=float(m["trajectory_interval_ps"]),
        forcefield=m["forcefield"],
        water_model=m["water_model"],
    )


def build_system(receptor_pdb: Path, ligand_sdf: Path, cfg: MDCfg, work: Path):
    rdmol = next(iter(Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)))
    if rdmol is None:
        sys.stderr.write(f"failed to read ligand sdf: {ligand_sdf}\n")
        sys.exit(1)
    if rdmol.GetNumConformers() == 0:
        sys.stderr.write(f"ligand has no conformer: {ligand_sdf}\n")
        sys.exit(1)
    off_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_mol])

    pdb = app.PDBFile(str(receptor_pdb))
    forcefield = app.ForceField(cfg.forcefield, cfg.water_model)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    modeller = app.Modeller(pdb.topology, pdb.positions)

    lig_top = off_mol.to_topology().to_openmm()
    lig_pos = off_mol.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)

    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=cfg.box_buffer_nm * u.nanometer,
        ionicStrength=cfg.salt_M * u.molar,
        positiveIon="Na+",
        negativeIon="Cl-",
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * u.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )
    work.mkdir(parents=True, exist_ok=True)
    with open(work / "topology.pdb", "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    return system, modeller


def restrain_heavy(system: mm.System, modeller: app.Modeller, k: float) -> int:
    force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k * u.kilojoule_per_mole / u.nanometer ** 2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    pos = modeller.positions
    for i, atom in enumerate(modeller.topology.atoms()):
        if atom.element is None or atom.element.symbol == "H":
            continue
        if atom.residue.name in ("HOH", "WAT", "NA", "CL", "Na+", "Cl-"):
            continue
        x, y, z = pos[i].value_in_unit(u.nanometer)
        force.addParticle(i, [x, y, z])
    return system.addForce(force)


def run_simulation(receptor_pdb: Path, ligand_sdf: Path, out_dir: Path,
                   seed: int, cfg: MDCfg) -> None:
    system, modeller = build_system(receptor_pdb, ligand_sdf, cfg, out_dir)
    restraint_idx = restrain_heavy(system, modeller, k=1000.0)

    barostat = mm.MonteCarloBarostat(cfg.pressure_bar * u.bar,
                                     cfg.temperature_K * u.kelvin, 25)
    barostat_idx = system.addForce(barostat)

    integrator = mm.LangevinMiddleIntegrator(cfg.temperature_K * u.kelvin,
                                             1.0 / u.picosecond,
                                             cfg.timestep_fs * u.femtosecond)
    integrator.setRandomNumberSeed(int(seed))

    platform = mm.Platform.getPlatformByName(_pick_platform())
    sim = app.Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)

    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(cfg.temperature_K * u.kelvin, int(seed))

    eq_steps = max(1, int(cfg.equilibration_ns * 1000.0 / cfg.timestep_fs * 1000.0))
    sim.step(eq_steps // 2)

    system.removeForce(restraint_idx)
    sim.context.reinitialize(preserveState=True)
    sim.step(eq_steps // 2)

    interval_steps = max(1, int(cfg.trajectory_interval_ps * 1000.0 / cfg.timestep_fs))
    prod_steps = max(1, int(cfg.production_ns * 1000.0 * 1000.0 / cfg.timestep_fs))

    dcd = app.DCDReporter(str(out_dir / "trajectory.dcd"), interval_steps)
    state_log = app.StateDataReporter(str(out_dir / "state.log"), interval_steps,
                                      step=True, potentialEnergy=True,
                                      temperature=True, volume=True)
    sim.reporters.append(dcd)
    sim.reporters.append(state_log)
    sim.step(prod_steps)

    state = sim.context.getState(getPositions=True, getVelocities=True,
                                 getForces=False, getEnergy=True,
                                 enforcePeriodicBox=True)
    with open(out_dir / "final.pdb", "w") as fh:
        app.PDBFile.writeFile(sim.topology, state.getPositions(), fh, keepIds=True)
    with open(out_dir / "state.xml", "w") as fh:
        fh.write(mm.XmlSerializer.serialize(state))


def _pick_platform() -> str:
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            mm.Platform.getPlatformByName(name)
            return name
        except Exception:
            continue
    return "CPU"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, type=Path)
    ap.add_argument("--ligand", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--seed", required=True, type=int)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()
    cfg = read_cfg(a.config)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    run_simulation(a.receptor, a.ligand, a.out_dir, a.seed, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
