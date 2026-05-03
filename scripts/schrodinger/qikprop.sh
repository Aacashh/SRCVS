#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SCHRODINGER:-}" ]]; then
    echo "SCHRODINGER environment variable not set" >&2
    exit 2
fi

xp_csv="$1"
ligands="$2"
out_dir="$3"
config="$4"

mkdir -p "$out_dir"
cd "$out_dir"

cp "$ligands" xp_hits.maegz
"$SCHRODINGER/qikprop" \
    -fast \
    -outname src_qikprop \
    xp_hits.maegz \
    -WAIT

"$SCHRODINGER/run" python3 - "$out_dir" "$xp_csv" "$config" <<'PYEOF'
import sys, csv, os
import pandas as pd, yaml
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

out_dir, xp_csv, cfg_path = sys.argv[1], sys.argv[2], sys.argv[3]
with open(cfg_path) as fh: cfg = yaml.safe_load(fh)
qcfg = cfg["admet"]

qp = pd.read_csv(os.path.join(out_dir, "src_qikprop.csv"))
glide = pd.read_csv(xp_csv)
df = qp.merge(glide, left_on="molecule", right_on="title", how="inner")

df["lipinski_pass"] = df["#stars"] <= qcfg["lipinski_max_violations"]
df["bbb_ok"] = df["QPlogBB"].between(qcfg["bbb_logbb_min"], qcfg["bbb_logbb_max"])
df["herg_ok"] = df["QPlogHERG"] > qcfg["herg_logic50_min"]
df["caco2_ok"] = df["QPPCaco"] > qcfg["caco2_min"]
df["mdck_ok"] = df["QPPMDCK"] > qcfg["mdck_min"]

params = FilterCatalogParams()
for pf in [FilterCatalogParams.FilterCatalogs.PAINS_A,
           FilterCatalogParams.FilterCatalogs.PAINS_B,
           FilterCatalogParams.FilterCatalogs.PAINS_C]:
    params.AddCatalog(pf)
pains = FilterCatalog.FilterCatalog(params)

def pains_clean(smi):
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    return False if m is None else not pains.HasMatch(m)
df["pains_pass"] = df.get("SMILES", df.get("smiles", pd.Series([""]*len(df)))).apply(pains_clean)

df.to_csv(os.path.join(out_dir, "admet_descriptors.csv"), index=False)
passed = df[df["lipinski_pass"] & df["bbb_ok"] & df["herg_ok"] & df["caco2_ok"] & df["mdck_ok"] & df["pains_pass"]]
passed = passed.rename(columns={"SMILES":"smiles","molecule":"name","glide_score":"score"})
passed.to_csv(os.path.join(out_dir, "admet_passed.csv"), index=False)
PYEOF
