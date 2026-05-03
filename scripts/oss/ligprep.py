from __future__ import annotations

import argparse
import csv
import sys
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
            if len(out) >= cfg.max_stereo * cfg.max_tautomers:
                return out
    return out


def embed3d(mol: Chem.Mol, seed: int = 42) -> Chem.Mol | None:
    m = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(m, maxIters=500)
    return m


def to_pdbqt(mol: Chem.Mol, out_path: Path) -> bool:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    prep = MoleculePreparation()
    prep.prepare(mol)
    pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(prep.setup)
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
        suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=True)
        for m in suppl:
            if m is None:
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
    for src, mol in iter_input(a.inputs, cfg.max_input):
        title = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        for state in enumerate_states(mol, cfg):
            embedded = embed3d(state, seed=42)
            if embedded is None:
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
    with open(summary_csv, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
