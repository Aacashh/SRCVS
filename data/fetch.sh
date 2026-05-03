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
                --connect-timeout 30 --max-time 900 \
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

is_html() {
    [[ -s "$1" ]] && head -c 256 "$1" | grep -qi -e "<html" -e "<!DOCTYPE" -e "<head"
}

is_valid_sdf() {
    [[ -s "$1" ]] || return 1
    if is_html "$1"; then return 1; fi
    grep -q "^M  END\|^\\\$\\\$\\\$\\\$" "$1"
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
[[ -s "$dude_actives" ]] || fetch_with_retry "https://dude.docking.org/targets/src/actives_final.ism" "$dude_actives" || true
[[ -s "$dude_decoys"  ]] || fetch_with_retry "https://dude.docking.org/targets/src/decoys_final.ism"  "$dude_decoys"  || true

nubbe_path="$raw_dir/nubbe_full.sdf"

if [[ ! -s "$nubbe_path" ]]; then
    coconut_zip="$raw_dir/coconut.sdf.zip"
    coconut_urls=(
        "https://zenodo.org/records/13692394/files/coconut_complete-10-2024.sdf.zip"
        "https://zenodo.org/record/13692394/files/coconut_complete-10-2024.sdf.zip"
    )
    for u in "${coconut_urls[@]}"; do
        if fetch_with_retry "$u" "$coconut_zip"; then
            if unzip -o -q "$coconut_zip" -d "$raw_dir" 2>>"$log"; then
                cand=$(ls "$raw_dir"/coconut*.sdf "$raw_dir"/COCONUT*.sdf 2>/dev/null | head -1 || true)
                if [[ -n "$cand" ]] && is_valid_sdf "$cand"; then
                    mv "$cand" "$nubbe_path"
                    rm -f "$coconut_zip"
                    break
                fi
            fi
            rm -f "$coconut_zip"
        fi
    done
fi

if [[ ! -s "$nubbe_path" ]]; then
    cids_file="$root/data/np_seed_cids.txt"
    if [[ -s "$cids_file" ]]; then
        tmp_dir="$raw_dir/_pubchem_np"
        mkdir -p "$tmp_dir"
        > "$nubbe_path"
        batch_size=50
        cids_clean=$(grep -E '^[0-9]+' "$cids_file" | sort -u)
        n_total=$(printf '%s\n' "$cids_clean" | wc -l)
        printf '%s\n' "$cids_clean" | awk -v b=$batch_size 'NR%b==1{c++} {print > ("'"$tmp_dir"'/batch_"c)}'
        ok_n=0
        for bf in "$tmp_dir"/batch_*; do
            [[ -s "$bf" ]] || continue
            ids=$(tr '\n' ',' < "$bf" | sed 's/,$//')
            chunk_out="$bf.sdf"
            url_3d="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${ids}/SDF?record_type=3d"
            url_2d="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${ids}/SDF"
            if fetch_with_retry "$url_3d" "$chunk_out" || fetch_with_retry "$url_2d" "$chunk_out"; then
                if is_valid_sdf "$chunk_out"; then
                    cat "$chunk_out" >> "$nubbe_path"
                    ok_n=$((ok_n + 1))
                fi
            fi
            sleep 1
        done
        rm -rf "$tmp_dir"
        if [[ $ok_n -eq 0 || ! -s "$nubbe_path" ]]; then
            rm -f "$nubbe_path"
        fi
    fi
fi

if [[ ! -s "$nubbe_path" ]]; then
    echo "natural-products library auto-fetch failed; pipeline will run on the 5 PubChem controls only" >&2
else
    verify_sha256 "$nubbe_path" "" || true
fi

enamine_path="$raw_dir/enamine_np.sdf"
if [[ ! -s "$enamine_path" ]]; then
    echo "warning: enamine_np.sdf not present (optional - request library from https://enamine.net)" >&2
else
    verify_sha256 "$enamine_path" "" || true
fi
