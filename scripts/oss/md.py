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

# MONKEY PATCH FOR FAST CHARGES
_orig_assign = Molecule.assign_partial_charges
def _fast_assign(self, partial_charge_method=None, **kwargs):
    _orig_assign(self, partial_charge_method="gasteiger", **kwargs)
Molecule.assign_partial_charges = _fast_assign


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
    implicit_solvent: str


def read_cfg(path: Path) -> MDCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    m = d["md"]
    return MDCfg(
        box_buffer_nm=float(m.get("box_buffer_nm", 1.0)),
        salt_M=float(m.get("ion_concentration_M", 0.15)),
        temperature_K=float(m["temperature_K"]),
        pressure_bar=float(m["pressure_bar"]),
        production_ns=float(m["production_ns"]),
        equilibration_ns=float(m["equilibration_ns"]),
        timestep_fs=float(m["timestep_fs"]),
        trajectory_interval_ps=float(m["trajectory_interval_ps"]),
        forcefield=m["forcefield"],
        implicit_solvent=m.get("implicit_solvent", "implicit/obc2.xml"),
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
    forcefield = app.ForceField(cfg.forcefield, cfg.implicit_solvent)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    modeller = app.Modeller(pdb.topology, pdb.positions)

    lig_top = off_mol.to_topology().to_openmm()
    lig_pos = off_mol.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        rigidWater=False,
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

    platform_name = _pick_platform()
    print(f"Initializing OpenMM simulation using platform: {platform_name}")
    platform = mm.Platform.getPlatformByName(platform_name)

    # ---------------------------------------------------------------
    # Stage 1: Fast cascaded minimization + NVT bumps
    #   Phase A: restrain all heavy atoms, relax water (2000 iter + 2000 NVT steps)
    #   Phase B: restrain backbone, relax sidechains/ligand (2000 iter + 2000 NVT steps)
    #   Phase C: fully free final minimisation (2000 iter)
    # ---------------------------------------------------------------
    def _add_restraint(k_val, mask_fn):
        force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", k_val * u.kilojoule_per_mole / u.nanometer**2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        pos = modeller.positions
        for i, atom in enumerate(modeller.topology.atoms()):
            if mask_fn(atom):
                x, y, z = pos[i].value_in_unit(u.nanometer)
                force.addParticle(i, [x, y, z])
        return system.addForce(force)

    def _is_heavy_non_water(a):
        return (a.element is not None
                and a.element.symbol != "H"
                and a.residue.name not in ("HOH", "WAT", "NA", "CL", "Na+", "Cl-"))

    def _is_backbone(a):
        return _is_heavy_non_water(a) and a.name in ("CA", "C", "N", "O")

    # Use a dummy integrator for minimization only (no dynamics between phases)
    integ_min = mm.LangevinMiddleIntegrator(
        5.0 * u.kelvin, 1.0 / u.picosecond, 0.001 * u.femtosecond)
    integ_min.setRandomNumberSeed(int(seed))

    # Phase A — restrain ALL heavy atoms strongly, minimize water only
    rid_a = _add_restraint(5000.0, _is_heavy_non_water)
    sim = app.Simulation(modeller.topology, system, integ_min, platform)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=5000)  # pure minimization, no dynamics
    system.removeForce(rid_a)
    sim.context.reinitialize(preserveState=True)

    # Phase B — restrain backbone only, relax sidechains + ligand
    rid_b = _add_restraint(1000.0, _is_backbone)
    sim.context.reinitialize(preserveState=True)
    sim.minimizeEnergy(maxIterations=5000)  # pure minimization, no dynamics
    system.removeForce(rid_b)
    sim.context.reinitialize(preserveState=True)

    # Phase C — fully free final minimisation
    sim.minimizeEnergy(maxIterations=5000)
    print("Minimization done.")

    # ---------------------------------------------------------------
    # Stage 2: Ultra-gentle NVT warm-up 0 K → T_target
    #   Start at 0.01 fs timestep to survive any residual clashes,
    #   ramp temperature and timestep together over ~30 ps total.
    # ---------------------------------------------------------------
    t_target = cfg.temperature_K
    # Phase 1: 0→50 K at 0.01 fs (2000 steps = 0.02 ps)
    integ_warmup = mm.LangevinMiddleIntegrator(
        1.0 * u.kelvin, 1.0 / u.picosecond, 0.01 * u.femtosecond)
    integ_warmup.setRandomNumberSeed(int(seed))
    state_min = sim.context.getState(getPositions=True)
    sim_w = app.Simulation(modeller.topology, system, integ_warmup, platform)
    sim_w.context.setPositions(state_min.getPositions())
    sim_w.context.setVelocitiesToTemperature(1 * u.kelvin, int(seed))
    for t in [10, 30, 50]:
        sim_w.integrator.setTemperature(t * u.kelvin)
        sim_w.step(2000)
    # Phase 2: 50→T_target at 0.1 fs (1000 steps per ramp = ~0.1 ps each)
    sim_w.integrator.setStepSize(0.1 * u.femtosecond)
    for t in [100, 150, 200, 250, int(t_target)]:
        sim_w.integrator.setTemperature(t * u.kelvin)
        sim_w.step(1000)
    # Phase 3: T_target at 0.5 fs (10 ramps × 4000 steps = 20 ps)
    sim_w.integrator.setStepSize(0.5 * u.femtosecond)
    for i in range(1, 11):
        sim_w.integrator.setTemperature(t_target * u.kelvin)
        sim_w.step(4000)
    print("NVT warm-up done.")

    # --- Stage 3: NVT equilibration (unrestrained, full timestep) ---
    integrator = mm.LangevinMiddleIntegrator(cfg.temperature_K * u.kelvin,
                                             1.0 / u.picosecond,
                                             cfg.timestep_fs * u.femtosecond)
    integrator.setRandomNumberSeed(int(seed))
    state_nvt = sim_w.context.getState(getPositions=True, getVelocities=True)
    sim2 = app.Simulation(modeller.topology, system, integrator, platform)
    sim2.context.setPositions(state_nvt.getPositions())
    sim2.context.setVelocities(state_nvt.getVelocities())

    eq_steps = max(1, int(cfg.equilibration_ns * 1e6 / cfg.timestep_fs))
    sim2.step(eq_steps)
    print("Equilibration done.")

    # --- Stage 4: Production ---
    interval_steps = max(1, int(cfg.trajectory_interval_ps * 1000.0 / cfg.timestep_fs))
    prod_steps = max(1, int(cfg.production_ns * 1e6 / cfg.timestep_fs))

    dcd = app.DCDReporter(str(out_dir / "trajectory.dcd"), interval_steps)
    state_log = app.StateDataReporter(str(out_dir / "state.log"), interval_steps,
                                      step=True, potentialEnergy=True,
                                      temperature=True, volume=True)
    sim2.reporters.append(dcd)
    sim2.reporters.append(state_log)
    sim2.step(prod_steps)

    state = sim2.context.getState(getPositions=True, getVelocities=True,
                                  getForces=False, getEnergy=True)
    with open(out_dir / "final.pdb", "w") as fh:
        app.PDBFile.writeFile(sim2.topology, state.getPositions(), fh, keepIds=True)
    with open(out_dir / "state.xml", "w") as fh:
        fh.write(mm.XmlSerializer.serialize(state))


def _pick_platform() -> str:
    """Return the first platform that actually works (not just registered)."""
    import openmm as _mm
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            plat = _mm.Platform.getPlatformByName(name)
            # Probe with a trivial 1-particle system
            sys_test = _mm.System()
            sys_test.addParticle(1.0)
            integ_test = _mm.VerletIntegrator(0.001)
            ctx = _mm.Context(sys_test, integ_test, plat)
            del ctx
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
