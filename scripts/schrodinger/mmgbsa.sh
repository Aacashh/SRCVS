#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SCHRODINGER:-}" ]]; then
    echo "SCHRODINGER environment variable not set" >&2
    exit 2
fi

mode=""
receptor=""
ligands=""
runs_json=""
out_dir=""
config=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) mode="$2"; shift 2 ;;
        --receptor) receptor="$2"; shift 2 ;;
        --ligands) ligands="$2"; shift 2 ;;
        --runs-json) runs_json="$2"; shift 2 ;;
        --out-dir) out_dir="$2"; shift 2 ;;
        --config) config="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 4 ;;
    esac
done

if [[ -z "$mode" || -z "$out_dir" || -z "$config" ]]; then
    echo "missing required arg" >&2
    exit 4
fi

mkdir -p "$out_dir"
cd "$out_dir"

if [[ "$mode" == "pre" ]]; then
    "$SCHRODINGER/run" python3 - "$receptor" "$ligands" "$out_dir" <<'PYEOF'
import sys
from schrodinger import structure
recep = next(structure.StructureReader(sys.argv[1]))
out = sys.argv[3] + "/leads_pv.maegz"
w = structure.StructureWriter(out)
w.append(recep)
for st in structure.StructureReader(sys.argv[2]):
    w.append(st)
w.close()
PYEOF
    "$SCHRODINGER/prime_mmgbsa" \
        -ligand "(res.ptype UNK) OR (chain L)" \
        -out_type COMPLEX \
        -job_type REAL_MIN \
        -csv_output yes \
        -jobname mmgbsa_pre \
        leads_pv.maegz \
        -HOST localhost \
        -WAIT
    cp mmgbsa_pre-out.csv mmgbsa_pre.csv
elif [[ "$mode" == "post" ]]; then
    "$SCHRODINGER/run" python3 - "$runs_json" "$out_dir" "$config" <<'PYEOF'
import sys, json, os, csv, glob, subprocess
import yaml, pandas as pd
runs = json.load(open(sys.argv[1]))
out_dir = sys.argv[2]
cfg = yaml.safe_load(open(sys.argv[3]))
stride = int(cfg["mmgbsa_post"]["frames_stride"])
start_pct = float(cfg["mmgbsa_post"]["frames_start_pct"])
sch = os.environ["SCHRODINGER"]
all_rows = []
for r in runs:
    rd = r["dir"]; lead = r["lead"]; seed = r["seed"]
    cms = glob.glob(os.path.join(rd, "*md*.cms"))
    if not cms: continue
    trj = glob.glob(os.path.join(rd, "*trj"))
    if not trj: continue
    sub = os.path.join(rd, "frames.maegz")
    n_frames = 1000
    start = int(n_frames * start_pct / 100)
    subprocess.run([f"{sch}/run","trj_extract_subsystem.py", cms[0],
                    "-t", trj[0], "-ASL", "all",
                    "-frames", f"{start}:{n_frames}:{stride}",
                    "-out", sub], check=False)
    subprocess.run([f"{sch}/prime_mmgbsa","-ligand","res.ptype UNK",
                    "-job_type","ENERGY","-out_type","COMPLEX","-csv_output","yes",
                    "-jobname", f"{lead}_s{seed}_post", sub,
                    "-HOST","localhost","-WAIT"], check=False, cwd=rd)
    csv_p = os.path.join(rd, f"{lead}_s{seed}_post-out.csv")
    if os.path.exists(csv_p):
        df = pd.read_csv(csv_p)
        df["lead"] = lead; df["seed"] = seed
        all_rows.append(df)
if all_rows:
    df = pd.concat(all_rows, ignore_index=True)
    df.to_csv(os.path.join(out_dir, "mmgbsa_post.csv"), index=False)
    summary = df.groupby("lead")["r_psp_MMGBSA_dG_Bind"].agg(["mean","std","count"]).reset_index()
    summary.to_csv(os.path.join(out_dir, "mmgbsa_post_summary.csv"), index=False)
else:
    open(os.path.join(out_dir,"mmgbsa_post.csv"),"w").write("lead,seed,frame,dG_bind\n")
    open(os.path.join(out_dir,"mmgbsa_post_summary.csv"),"w").write("lead,mean,std,count\n")
PYEOF
else
    echo "invalid --mode: $mode" >&2
    exit 4
fi
