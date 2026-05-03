from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina


@dataclass
class Cfg:
    root: Path
    out_root: Path
    log_dir: Path
    n_leads: int
    cutoff: float
    fp_radius: int
    fp_bits: int


def load_cfg(path: Path) -> Cfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    root = Path(d["project"]["root"]).resolve()
    sel = d["selection"]
    return Cfg(
        root=root,
        out_root=root / d["project"]["out_dir"],
        log_dir=root / d["project"]["log_dir"],
        n_leads=int(sel["n_leads"]),
        cutoff=float(sel["cluster_cutoff"]),
        fp_radius=int(sel["fp_radius"]),
        fp_bits=int(sel["fp_bits"]),
    )


def _setup_log(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=str(p), level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s", force=True,
    )


def _read_manifest(p: Path) -> dict:
    if not p.exists():
        sys.stderr.write(f"missing dependency manifest: {p}\n")
        sys.exit(3)
    with open(p) as fh:
        return json.load(fh)


def write_manifest(out_dir: Path, payload: dict) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    mp = out_dir / "manifest.json"
    with open(mp, "w") as fh:
        json.dump(payload, fh, indent=2)
    return mp


def cluster_select(df: pd.DataFrame, cfg: Cfg) -> pd.DataFrame:
    mols = []
    keep_idx = []
    for i, smi in enumerate(df["smiles"].tolist()):
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if m is None:
            continue
        mols.append(m)
        keep_idx.append(i)
    if not mols:
        return df.head(0)
    df = df.iloc[keep_idx].reset_index(drop=True)
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, cfg.fp_radius, cfg.fp_bits) for m in mols]
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1.0 - x for x in sims])
    if n == 1:
        return df.head(1)
    clusters = Butina.ClusterData(dists, n, cfg.cutoff, isDistData=True)
    score_col = "score" if "score" in df.columns else df.select_dtypes("number").columns[0]
    reps = []
    for cluster in sorted(clusters, key=len, reverse=True)[: cfg.n_leads]:
        sub = df.iloc[list(cluster)].sort_values(score_col)
        reps.append(sub.iloc[0])
    return pd.DataFrame(reps).reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stack", default="oss", choices=["oss", "schrodinger"])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    cfg = load_cfg(a.config)
    out_dir = cfg.out_root / "05_select"
    log_path = cfg.log_dir / "phase5.log"
    mp = out_dir / "manifest.json"
    if mp.exists() and not a.force:
        sys.stdout.write(str(mp) + "\n")
        return 0

    _setup_log(log_path)
    logging.info("phase5 start")

    admet_mf = _read_manifest(cfg.out_root / "04_admet" / "manifest.json")
    passed_csv = Path(admet_mf["outputs"]["passed"])
    if not passed_csv.exists():
        sys.stderr.write(f"missing admet passed csv: {passed_csv}\n")
        return 3

    df = pd.read_csv(passed_csv)
    if "smiles" not in df.columns:
        sys.stderr.write("admet csv lacks 'smiles' column\n")
        return 4

    leads = cluster_select(df, cfg)
    out_dir.mkdir(parents=True, exist_ok=True)
    leads_csv = out_dir / "final_leads.csv"
    leads_sdf = out_dir / "final_leads.sdf"
    leads.to_csv(leads_csv, index=False)

    writer = Chem.SDWriter(str(leads_sdf))
    for _, row in leads.iterrows():
        m = Chem.MolFromSmiles(row["smiles"])
        if m is None:
            continue
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(m)
        m.SetProp("_Name", str(row.get("name", row.get("title", ""))))
        for col in row.index:
            v = row[col]
            if pd.notna(v):
                m.SetProp(col, str(v))
        writer.write(m)
    writer.close()

    payload = {
        "phase": 5,
        "name": "lead_selection",
        "stack": a.stack,
        "outputs": {
            "leads_csv": str(leads_csv),
            "leads_sdf": str(leads_sdf),
        },
    }
    p = write_manifest(out_dir, payload)
    sys.stdout.write(str(p) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
