#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SCHRODINGER:-}" ]]; then
    echo "SCHRODINGER environment variable not set" >&2
    exit 2
fi

pdb="$1"
out_dir="$2"
config="$3"

mkdir -p "$out_dir"

ph=$(awk '/^  prep:/{f=1;next} f&&/^    ph:/{print $2;exit}' "$config")
rmsd=$(awk '/^  prep:/{f=1;next} f&&/^    rmsd:/{print $2;exit}' "$config")
ff=$(awk '/^  prep:/{f=1;next} f&&/^    forcefield:/{print $2;exit}' "$config")
lig_resname=$(awk '/^target:/{f=1;next} f&&/^  ligand_resname:/{print $2;exit}' "$config")

prepared_mae="$out_dir/receptor_prepared.mae"
prepared_pdb="$out_dir/receptor_prepared.pdb"

"$SCHRODINGER/utilities/prepwizard" \
    -fillsidechains \
    -fillloops \
    -disulfides \
    -mse \
    -propka_pH "$ph" \
    -minimize_adj_h \
    -f "$ff" \
    -rmsd "$rmsd" \
    -WAIT \
    "$pdb" "$prepared_mae"

"$SCHRODINGER/utilities/structconvert" \
    -imae "$prepared_mae" \
    -opdb "$prepared_pdb"

"$SCHRODINGER/run" split_structure.py \
    -m ligand \
    -groupby molecule \
    "$prepared_mae" \
    "$out_dir/split.mae"

ligand_mae=$(ls "$out_dir"/split-lig-*.mae 2>/dev/null | head -1 || true)
if [[ -z "$ligand_mae" ]]; then
    echo "no co-crystal ligand extracted" >&2
    exit 1
fi
"$SCHRODINGER/utilities/structconvert" \
    -imae "$ligand_mae" \
    -opdb "$out_dir/cocrystal_ligand.pdb"

cat > "$out_dir/grid.in" <<EOF
GRID_CENTER_USING_LIG  True
LIGAND_INDEX           1
LIGAND_FILE            $ligand_mae
RECEP_FILE             $prepared_mae
INNERBOX               12, 12, 12
OUTERBOX               24, 24, 24
GRIDFILE               $out_dir/grid.zip
FORCEFIELD             OPLS_2005
EOF

"$SCHRODINGER/glide" "$out_dir/grid.in" -OVERWRITE -WAIT -HOST localhost

"$SCHRODINGER/run" python3 - "$ligand_mae" "$out_dir" <<'PYEOF'
import sys, json
from schrodinger import structure
lig = next(structure.StructureReader(sys.argv[1]))
xs = [a.x for a in lig.atom]; ys = [a.y for a in lig.atom]; zs = [a.z for a in lig.atom]
center = [sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)]
data = {"center": center, "inner": [12,12,12], "outer": [24,24,24], "resname": lig.atom[1].pdbres.strip()}
with open(f"{sys.argv[2]}/grid.json","w") as fh: json.dump(data, fh, indent=2)
PYEOF

cat > "$out_dir/redock.in" <<EOF
GRIDFILE   $out_dir/grid.zip
LIGANDFILE $ligand_mae
PRECISION  XP
EOF
"$SCHRODINGER/glide" "$out_dir/redock.in" -OVERWRITE -WAIT -HOST localhost || true

"$SCHRODINGER/run" python3 - "$out_dir" "$lig_resname" <<'PYEOF'
import sys, csv, glob
from schrodinger import structure
from schrodinger.structutils import rmsd as srmsd
out_dir = sys.argv[1]
ref_path = glob.glob(f"{out_dir}/split-lig-*.mae")[0]
ref = next(structure.StructureReader(ref_path))
pv = glob.glob(f"{out_dir}/redock_pv.maegz")
rows = [["pose","rmsd"]]
if pv:
    reader = structure.StructureReader(pv[0])
    next(reader)
    for i, st in enumerate(reader):
        try:
            r = srmsd.calculate_in_place_rmsd(ref, list(range(1, ref.atom_total+1)), st, list(range(1, st.atom_total+1)))
            rows.append([i, r])
        except Exception:
            continue
with open(f"{out_dir}/redock_rmsd.csv","w",newline="") as fh:
    csv.writer(fh).writerows(rows)
PYEOF
