#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cfg="$root/tests/smoke_config.yaml"

mkdir -p "$root/tests"
python - "$root/config.yaml" "$cfg" <<'PYEOF'
import sys, yaml
src, dst = sys.argv[1], sys.argv[2]
with open(src) as fh: d = yaml.safe_load(fh)
d["library"]["max_input_compounds"] = 50
d["library"]["include_controls"] = True
d["docking"]["xp_top_n"] = 10
d["docking"]["exhaustiveness"] = 4
d["selection"]["n_leads"] = 3
d["md"]["production_ns"] = 0.05
d["md"]["equilibration_ns"] = 0.01
d["md"]["replicates"] = 1
d["md"]["trajectory_interval_ps"] = 5
d["mmgbsa_post"]["frames_stride"] = 5
d["project"]["out_dir"] = "out_smoke"
d["project"]["log_dir"] = "out_smoke/logs"
with open(dst, "w") as fh: yaml.safe_dump(d, fh)
PYEOF

bash "$root/run.sh" --config "$cfg" --force >/dev/null

for p in 01_target 02_library 03_dock 04_admet 05_select 06_mmgbsa_pre 07_md 08_analysis 09_mmgbsa_post; do
    mf="$root/out_smoke/$p/manifest.json"
    if [[ ! -s "$mf" ]]; then
        echo "missing manifest: $mf" >&2
        exit 1
    fi
done

printf '%s\n' "$root/out_smoke/run_manifest.json"
