#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SCHRODINGER:-}" ]]; then
    echo "SCHRODINGER environment variable not set" >&2
    exit 2
fi

out_dir="$1"
config="$2"
shift 2
inputs=("$@")

mkdir -p "$out_dir"

ph=$(awk '/^  ligprep:/{f=1;next} f&&/^    ph:/{print $2;exit}' "$config")
pht=$(awk '/^  ligprep:/{f=1;next} f&&/^    pht:/{print $2;exit}' "$config")
stereo=$(awk '/^  ligprep:/{f=1;next} f&&/^    max_stereo:/{print $2;exit}' "$config")
taut=$(awk '/^  ligprep:/{f=1;next} f&&/^    max_tautomers:/{print $2;exit}' "$config")
ffid=$(awk '/^  ligprep:/{f=1;next} f&&/^    forcefield_id:/{print $2;exit}' "$config")

per_outs=()
for sdf in "${inputs[@]}"; do
    base=$(basename "$sdf" .sdf)
    out_mae="$out_dir/${base}_prep.maegz"
    "$SCHRODINGER/ligprep" \
        -isd "$sdf" \
        -omae "$out_mae" \
        -epik \
        -ph "$ph" \
        -pht "$pht" \
        -bff "$ffid" \
        -s "$stereo" \
        -t "$taut" \
        -i 0 \
        -nt \
        -HOST localhost \
        -WAIT
    per_outs+=("$out_mae")
done

"$SCHRODINGER/utilities/structcat" \
    -imae "${per_outs[@]}" \
    -omae "$out_dir/library_prepared.maegz"

mkdir -p "$out_dir/pdbqt"
echo "title,n_states" > "$out_dir/ligprep_summary.csv"
"$SCHRODINGER/run" python3 - "$out_dir/library_prepared.maegz" "$out_dir/ligprep_summary.csv" <<'PYEOF'
import sys, csv, collections
from schrodinger import structure
counts = collections.Counter()
for st in structure.StructureReader(sys.argv[1]):
    counts[st.title] += 1
with open(sys.argv[2],"a",newline="") as fh:
    w = csv.writer(fh)
    for t, n in counts.items():
        w.writerow([t, n])
PYEOF
