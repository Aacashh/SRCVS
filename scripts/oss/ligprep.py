from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass
class LigCfg:
    ph: float
    pht: float
    max_stereo: int
    max_tautomers: int
    max_input: int


def read_cfg(path: Path) -> LigCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    lp = d["library"]["ligprep"]
    return LigCfg(
        ph=float(lp["ph"]),
        pht=float(lp["pht"]),
        max_stereo=int(lp["max_stereo"]),
        max_tautomers=int(lp["max_tautomers"]),
        max_input=int(d["library"].get("max_input_compounds", 0)),
    )


def enumerate_states(mol: Chem.Mol, cfg: LigCfg) -> list:
    out = []
    opts = StereoEnumerationOptions(maxIsomers=cfg.max_stereo, onlyUnassigned=True, unique=True)
    iso_list = list(EnumerateStereoisomers(mol, options=opts))
    if not iso_list:
        iso_list = [mol]
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(cfg.max_tautomers)
    seen = set()
    for iso in iso_list:
        try:
            tauts = list(enumerator.Enumerate(iso))
        except Exception:
            tauts = [iso]
        for t in tauts:
            smi = Chem.MolToSmiles(t, canonical=True)
            if smi in seen:
                continue
            seen.add(smi)
            out.append(t)
    return out

def embed3d(mol: Chem.Mol, seed: int = 42) -> Chem.Mol | None:
    try:
        # validate molecule
        Chem.SanitizeMol(mol)

        m = Chem.AddHs(mol)

        # embed 3D coordinates
        if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
            return None

        # optimize geometry
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:
            AllChem.UFFOptimizeMolecule(m, maxIters=500)

        return m

    except Exception:
        # skip invalid molecules
        return None



def to_pdbqt(mol: Chem.Mol, out_path: Path) -> bool:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)
    # Handle both list and generator returns across Meeko versions
    if isinstance(mol_setups, list):
        if not mol_setups:
            return False
        setup = mol_setups[0]
    else:
        try:
            setup = next(iter(mol_setups))
        except (StopIteration, TypeError):
            return False
    pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(setup)
    if not is_ok:
        return False
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(pdbqt_string)
    return True


def safe_name(name: str, idx: int) -> str:
    base = name.strip() or f"lig_{idx:06d}"
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in base)[:80]


def iter_input(paths: list[Path], limit: int):
    n = 0
    for p in paths:
        # Use sanitize=False to avoid silently dropping molecules;
        # we sanitize manually so failures are logged.
        suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=False)
        for idx, m in enumerate(suppl):
            if m is None:
                sys.stderr.write(f"warning: molecule #{idx} in {p.name} could not be read, skipping\n")
                continue
            try:
                Chem.SanitizeMol(m)
            except Exception as e:
                sys.stderr.write(f"warning: molecule #{idx} in {p.name} failed sanitization: {e}\n")
                continue
            yield p.stem, m
            n += 1
            if limit and n >= limit:
                return


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--inputs", nargs="+", required=True, type=Path)
    a = ap.parse_args()

    cfg = read_cfg(a.config)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_dir = a.out_dir / "pdbqt"
    pdbqt_dir.mkdir(exist_ok=True)
    out_sdf = a.out_dir / "library_prepared.sdf"
    summary_csv = a.out_dir / "ligprep_summary.csv"

    sd_writer = Chem.SDWriter(str(out_sdf))
    rows = [["source", "input_name", "out_name", "smiles", "pdbqt"]]
    counter: dict[str, int] = {}
    mol_count = 0
    skip_count = 0
    t0 = time.time()
    for src, mol in iter_input(a.inputs, cfg.max_input):
        mol_count += 1
        title = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if mol_count % 25 == 0:
            elapsed = time.time() - t0
            sys.stderr.write(f"[ligprep] processed {mol_count} molecules ({skip_count} skipped) in {elapsed:.1f}s\n")
        for state in enumerate_states(mol, cfg):
            embedded = embed3d(state, seed=42)
            if embedded is None:
                skip_count += 1
                continue
            counter[title] = counter.get(title, 0) + 1
            sname = safe_name(title, counter[title])
            uname = f"{sname}_{counter[title]:02d}"
            embedded.SetProp("_Name", uname)
            sd_writer.write(embedded)
            pdbqt_path = pdbqt_dir / f"{uname}.pdbqt"
            ok = False
            try:
                ok = to_pdbqt(embedded, pdbqt_path)
            except Exception as e:
                sys.stderr.write(f"pdbqt failed for {uname}: {e}\n")
            smi = Chem.MolToSmiles(Chem.RemoveHs(embedded))
            rows.append([src, title, uname, smi, str(pdbqt_path) if ok else ""])
    sd_writer.close()
    elapsed = time.time() - t0
    sys.stderr.write(f"[ligprep] done: {mol_count} input molecules, {len(rows)-1} outputs, {skip_count} skipped in {elapsed:.1f}s\n")
    with open(summary_csv, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
