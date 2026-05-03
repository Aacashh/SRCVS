from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, FilterCatalog, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalogParams


@dataclass
class AdmetCfg:
    lipinski_max_violations: int
    pains_filter: bool
    bbb_logbb_min: float
    bbb_logbb_max: float
    herg_min: float
    caco2_min: float
    mdck_min: float


def read_cfg(path: Path) -> AdmetCfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    a = d["admet"]
    return AdmetCfg(
        lipinski_max_violations=int(a["lipinski_max_violations"]),
        pains_filter=bool(a["pains_filter"]),
        bbb_logbb_min=float(a["bbb_logbb_min"]),
        bbb_logbb_max=float(a["bbb_logbb_max"]),
        herg_min=float(a["herg_logic50_min"]),
        caco2_min=float(a["caco2_min"]),
        mdck_min=float(a["mdck_min"]),
    )


def lipinski_violations(m: Chem.Mol) -> int:
    v = 0
    if Descriptors.MolWt(m) > 500:
        v += 1
    if Crippen.MolLogP(m) > 5:
        v += 1
    if Lipinski.NumHDonors(m) > 5:
        v += 1
    if Lipinski.NumHAcceptors(m) > 10:
        v += 1
    return v


def veber_pass(m: Chem.Mol) -> bool:
    return (rdMolDescriptors.CalcTPSA(m) <= 140) and (Lipinski.NumRotatableBonds(m) <= 10)


def boiled_egg_logbb(m: Chem.Mol) -> float:
    logp = Crippen.MolLogP(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    return 0.152 * logp - 0.0148 * tpsa + 0.139


def caco2_proxy(m: Chem.Mol) -> float:
    tpsa = rdMolDescriptors.CalcTPSA(m)
    logp = Crippen.MolLogP(m)
    return max(0.0, 100.0 - tpsa) * (1.0 + 0.1 * logp)


def mdck_proxy(m: Chem.Mol) -> float:
    return caco2_proxy(m) * 0.9


def herg_proxy(m: Chem.Mol) -> float:
    logp = Crippen.MolLogP(m)
    base = -3.0 - 0.5 * max(0.0, logp - 3.0)
    return base


def pains_match(catalog, m: Chem.Mol) -> bool:
    return catalog.HasMatch(m)


def smiles_index(sdf_path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=True)
    for m in suppl:
        if m is None:
            continue
        title = m.GetProp("_Name") if m.HasProp("_Name") else ""
        if not title:
            continue
        out[title] = Chem.MolToSmiles(m)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xp-csv", required=True, type=Path)
    ap.add_argument("--ligand-sdf", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()

    cfg = read_cfg(a.config)
    a.out_dir.mkdir(parents=True, exist_ok=True)

    if not a.xp_csv.exists() or a.xp_csv.stat().st_size == 0:
        sys.stderr.write(f"missing or empty xp csv: {a.xp_csv}\n")
        return 3
    xp = pd.read_csv(a.xp_csv)
    smi_map = smiles_index(a.ligand_sdf)

    params = FilterCatalogParams()
    for pf in [FilterCatalogParams.FilterCatalogs.PAINS_A,
               FilterCatalogParams.FilterCatalogs.PAINS_B,
               FilterCatalogParams.FilterCatalogs.PAINS_C]:
        params.AddCatalog(pf)
    pains = FilterCatalog.FilterCatalog(params)

    rows = []
    for _, r in xp.iterrows():
        title = str(r["title"])
        smi = smi_map.get(title)
        if smi is None:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        hbd = Lipinski.NumHDonors(m)
        hba = Lipinski.NumHAcceptors(m)
        tpsa = rdMolDescriptors.CalcTPSA(m)
        rot = Lipinski.NumRotatableBonds(m)
        viol = lipinski_violations(m)
        veb = veber_pass(m)
        logbb = boiled_egg_logbb(m)
        caco2 = caco2_proxy(m)
        mdck = mdck_proxy(m)
        herg = herg_proxy(m)
        pa = not pains_match(pains, m) if cfg.pains_filter else True
        rows.append({
            "name": title,
            "smiles": smi,
            "score": float(r.get("score", 0.0)),
            "mw": mw, "logp": logp, "hbd": hbd, "hba": hba,
            "tpsa": tpsa, "rotors": rot, "lipinski_violations": viol,
            "veber_pass": veb, "logbb": logbb, "caco2": caco2,
            "mdck": mdck, "herg_logic50": herg, "pains_pass": pa,
        })

    df = pd.DataFrame(rows)
    df.to_csv(a.out_dir / "admet_descriptors.csv", index=False)

    mask = (
        (df["lipinski_violations"] <= cfg.lipinski_max_violations)
        & df["veber_pass"]
        & df["logbb"].between(cfg.bbb_logbb_min, cfg.bbb_logbb_max)
        & (df["caco2"] >= cfg.caco2_min)
        & (df["mdck"] >= cfg.mdck_min)
        & (df["herg_logic50"] >= cfg.herg_min)
        & df["pains_pass"]
    )
    passed = df[mask].sort_values("score").reset_index(drop=True)
    passed.to_csv(a.out_dir / "admet_passed.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
