#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SCHRODINGER:-}" ]]; then
    echo "SCHRODINGER environment variable not set" >&2
    exit 2
fi

receptor="$1"
ligands="$2"
out_dir="$3"
config="$4"

mkdir -p "$out_dir/poses"
target_dir="$(dirname "$receptor")"
grid_zip="$target_dir/grid.zip"
if [[ ! -s "$grid_zip" ]]; then
    echo "missing grid: $grid_zip" >&2
    exit 3
fi

vsw_inp="$out_dir/vsw.inp"
cat > "$vsw_inp" <<EOF
[SET:ORIGINAL_LIGANDS]
    VARCLASS   Structures
    FILES      $ligands,
[STAGE:HTVS]
    STAGECLASS         glide.DockingStage
    INPUTS             ORIGINAL_LIGANDS,
    OUTPUTS            HTVS_OUT,
    RECOMBINE          NO
    PRECISION          HTVS
    GRIDFILE           $grid_zip
    POSE_OUTTYPE       ligandlib_sd
    DOCKING_METHOD     confgen
    POSES_PER_LIG      1
[STAGE:SCORE_FILTER_HTVS]
    STAGECLASS         glide.ScoreFilterStage
    INPUTS             HTVS_OUT,
    OUTPUTS            SCORE_FILTER_HTVS_OUT,
    FILTER             "r_i_docking_score < 0.0"
    NUMBER_OUT         10%
[STAGE:SP]
    STAGECLASS         glide.DockingStage
    INPUTS             SCORE_FILTER_HTVS_OUT,
    OUTPUTS            SP_OUT,
    PRECISION          SP
    GRIDFILE           $grid_zip
    POSES_PER_LIG      3
[STAGE:SCORE_FILTER_SP]
    STAGECLASS         glide.ScoreFilterStage
    INPUTS             SP_OUT,
    OUTPUTS            SCORE_FILTER_SP_OUT,
    NUMBER_OUT         10%
[STAGE:XP]
    STAGECLASS         glide.DockingStage
    INPUTS             SCORE_FILTER_SP_OUT,
    OUTPUTS            XP_OUT,
    PRECISION          XP
    GRIDFILE           $grid_zip
    POSES_PER_LIG      5
    WRITE_XP_DESC      YES
[USEROUTS]
    USEROUTS  XP_OUT,
    STRUCTOUT XP_OUT
EOF

cd "$out_dir"
"$SCHRODINGER/vsw" "$vsw_inp" \
    -host_glide "localhost:8" \
    -adjust \
    -OVERWRITE \
    -WAIT

"$SCHRODINGER/run" python3 - "$out_dir" <<'PYEOF'
import sys, csv, glob, os
from schrodinger import structure
out_dir = sys.argv[1]
def emit(stage, out_csv):
    rows = [["title","glide_score","glide_emodel","glide_evdw","glide_ecoul","ligand_efficiency"]]
    for f in sorted(glob.glob(os.path.join(out_dir, f"*{stage}*pv*.maegz"))):
        reader = structure.StructureReader(f)
        try:
            next(reader)
        except StopIteration:
            continue
        for st in reader:
            rows.append([st.title,
                         st.property.get("r_i_docking_score"),
                         st.property.get("r_i_glide_emodel"),
                         st.property.get("r_i_glide_evdw"),
                         st.property.get("r_i_glide_ecoul"),
                         st.property.get("r_i_glide_ligand_efficiency")])
    rows[1:] = sorted(rows[1:], key=lambda r: (r[1] if r[1] is not None else 0))
    with open(os.path.join(out_dir, out_csv), "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
emit("HTVS_OUT", "htvs_scores.csv")
emit("SP_OUT",   "sp_scores.csv")
emit("XP_OUT",   "xp_scores.csv")
PYEOF
