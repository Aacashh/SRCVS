#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="$root/data/raw"
mkdir -p "$raw_dir"
log="$raw_dir/fetch.log"
checksums_log="$raw_dir/checksums.computed"
ua="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) SRCVS/1.0"

: > "$log"

fetch_with_retry() {
    local url="$1"
    local out="$2"
    local i=0
    while [[ $i -lt 3 ]]; do
        if curl -sSL --fail \
                --user-agent "$ua" \
                --retry 3 --retry-delay 2 \
                --connect-timeout 30 --max-time 600 \
                -o "$out" "$url" 2>>"$log"; then
            return 0
        fi
        i=$((i + 1))
        sleep $((i * 2))
    done
    return 1
}

verify_sha256() {
    local file="$1"
    local expected="${2:-}"
    local actual
    actual=$(sha256sum "$file" | awk '{print $1}')
    printf '%s  %s\n' "$actual" "$file" >> "$checksums_log"
    if [[ -n "$expected" && "$actual" != "$expected" ]]; then
        echo "sha256 mismatch for $file" >&2
        return 1
    fi
}

pdb_id="3G5D"
pdb_out="$raw_dir/${pdb_id}.pdb"
if [[ ! -s "$pdb_out" ]]; then
    if ! fetch_with_retry "https://files.rcsb.org/download/${pdb_id}.pdb" "$pdb_out"; then
        echo "failed to fetch PDB $pdb_id (see $log)" >&2
        exit 1
    fi
fi
verify_sha256 "$pdb_out" "" || true

declare -A controls=(
    [dasatinib]=3062316
    [bosutinib]=5328940
    [saracatinib]=10302451
    [pp1]=4878
    [pp2]=4879
)
for name in "${!controls[@]}"; do
    cid="${controls[$name]}"
    out="$raw_dir/control_${name}_${cid}.sdf"
    if [[ ! -s "$out" ]]; then
        url_3d="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d"
        url_2d="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF"
        if ! fetch_with_retry "$url_3d" "$out"; then
            sleep 2
            if ! fetch_with_retry "$url_2d" "$out"; then
                echo "failed to fetch PubChem CID $cid (see $log)" >&2
                exit 1
            fi
        fi
    fi
    verify_sha256 "$out" "" || true
done

dude_actives="$raw_dir/dude_src_actives.ism"
dude_decoys="$raw_dir/dude_src_decoys.ism"
if [[ ! -s "$dude_actives" ]]; then
    fetch_with_retry "https://dude.docking.org/targets/src/actives_final.ism" "$dude_actives" || true
fi
if [[ ! -s "$dude_decoys" ]]; then
    fetch_with_retry "https://dude.docking.org/targets/src/decoys_final.ism" "$dude_decoys" || true
fi
[[ -s "$dude_actives" ]] && { verify_sha256 "$dude_actives" "" || true; }
[[ -s "$dude_decoys" ]]  && { verify_sha256 "$dude_decoys"  "" || true; }

nubbe_path="$raw_dir/nubbe_full.sdf"
nubbe_urls=(
    "https://nubbe.iq.unesp.br/portal/files/nubbedb.sdf"
    "https://nubbe.iq.unesp.br/portal/nubbedb.sdf"
    "https://nubbe.iq.unesp.br/portal/files/nubbe-search-results.sdf"
)
if [[ ! -s "$nubbe_path" ]]; then
    for u in "${nubbe_urls[@]}"; do
        if fetch_with_retry "$u" "$nubbe_path"; then
            if [[ -s "$nubbe_path" ]] && head -c 200 "$nubbe_path" | grep -qi "<\|html\|<!DOCTYPE"; then
                rm -f "$nubbe_path"
                continue
            fi
            break
        fi
    done
fi
if [[ ! -s "$nubbe_path" ]]; then
    echo "auto-fetch failed for NuBBE; download nubbe_full.sdf from https://nubbe.iq.unesp.br/portal/nubbedb.html and place it at $nubbe_path" >&2
    exit 3
fi
verify_sha256 "$nubbe_path" "" || true

enamine_path="$raw_dir/enamine_np.sdf"
if [[ ! -s "$enamine_path" ]]; then
    echo "warning: $enamine_path not present (optional - request library from https://enamine.net)" >&2
else
    verify_sha256 "$enamine_path" "" || true
fi
