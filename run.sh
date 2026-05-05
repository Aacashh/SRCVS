#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
usage: bash run.sh [--config PATH] [--phase N] [--stack schrodinger|oss] [--force] [--dry]
EOF
    exit 4
}

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config="$root/config.yaml"
phase=""
stack="oss"
force=0
dry=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) config="${2:-}"; shift 2 ;;
        --phase) phase="${2:-}"; shift 2 ;;
        --stack) stack="${2:-}"; shift 2 ;;
        --force) force=1; shift ;;
        --dry) dry=1; shift ;;
        -h|--help) usage ;;
        *) echo "unknown argument: $1" >&2; exit 4 ;;
    esac
done

if [[ "$stack" != "schrodinger" && "$stack" != "oss" ]]; then
    echo "invalid --stack: $stack" >&2
    exit 4
fi
if [[ ! -f "$config" ]]; then
    echo "config not found: $config" >&2
    exit 4
fi

ts() { date +"%H:%M:%S"; }
log()  { printf '[%s] %s\n' "$(ts)" "$*"; }
hr()   { printf '\n========================================================\n'; }

export PROJECT_ROOT="$root"

log "srcvs pipeline starting (stack=$stack force=$force dry=$dry)"
log "config: $config"
log "root:   $root"

log "ensuring conda env (env/setup.sh)..."
bash "$root/env/setup.sh" >/dev/null
log "env/setup.sh done"

env_bin=""
if [[ -d "$HOME/miniforge3/envs/srcvs/bin" ]]; then
    env_bin="$HOME/miniforge3/envs/srcvs/bin"
elif [[ -d "$HOME/miniconda/envs/srcvs/bin" ]]; then
    env_bin="$HOME/miniconda/envs/srcvs/bin"
elif command -v conda >/dev/null 2>&1; then
    env_bin="$(conda info --base 2>/dev/null)/envs/srcvs/bin"
fi
if [[ -n "$env_bin" && -x "$env_bin/python" ]]; then
    export PATH="$env_bin:$PATH"
    log "python: $env_bin/python"
fi

if [[ "$dry" == "0" ]]; then
    log "fetching reference data (data/fetch.sh)..."
    t0=$(date +%s)
    bash "$root/data/fetch.sh"
    log "data fetch done in $(( $(date +%s) - t0 ))s"
fi

phases=(1 2 3 4 5 6 7 8 9)
if [[ -n "$phase" ]]; then
    if ! [[ "$phase" =~ ^[1-9]$ ]]; then
        echo "invalid --phase: $phase" >&2
        exit 4
    fi
    phases=("$phase")
fi

declare -A names=(
    [1]=phase1_target
    [2]=phase2_library
    [3]=phase3_dock
    [4]=phase4_admet
    [5]=phase5_select
    [6]=phase6_mmgbsa_pre
    [7]=phase7_md
    [8]=phase8_analyze
    [9]=phase9_mmgbsa_post
)

declare -A descriptions=(
    [1]="target preparation (PDBFixer, grid)"
    [2]="ligand library prep (LigPrep / RDKit + Meeko)"
    [3]="virtual screening (Glide / Vina HTVS+SP+XP)"
    [4]="ADMET filtering (QikProp / RDKit)"
    [5]="lead selection (Butina clustering)"
    [6]="pre-MD MM-GBSA"
    [7]="molecular dynamics"
    [8]="trajectory analysis (RMSD/RMSF/H-bonds)"
    [9]="post-MD MM-GBSA"
)

declare -a outputs=()
for p in "${phases[@]}"; do
    script="$root/scripts/${names[$p]}.py"
    if [[ ! -f "$script" ]]; then
        echo "missing phase script: $script" >&2
        exit 1
    fi
    desc="${descriptions[$p]}"
    if [[ "$dry" == "1" ]]; then
        log "[phase $p / dry] would run ${names[$p]} — $desc"
        outputs+=("$script --config $config --stack $stack")
        continue
    fi
    hr
    log "[phase $p / 9] starting: $desc"
    log "[phase $p] script:  ${names[$p]}.py"
    log "[phase $p] log:     out/logs/phase${p}.log"
    # Truncate the log file so old errors don't persist
    : > "$root/out/logs/phase${p}.log"
    extra=()
    [[ "$force" == "1" ]] && extra+=("--force")
    t0=$(date +%s)
    rc=0
    out_path=$(python "$script" --config "$config" --stack "$stack" "${extra[@]}" 2>&1) || rc=$?
    if [[ "$rc" -eq 0 ]]; then
        dt=$(( $(date +%s) - t0 ))
        log "[phase $p] done in ${dt}s — manifest: $out_path"
        outputs+=("$out_path")
    else
        log "[phase $p] FAILED (exit $rc) — see out/logs/phase${p}.log"
        echo "$out_path" >&2
        exit $rc
    fi
done

hr
mkdir -p "$root/out"
run_manifest="$root/out/run_manifest.json"
{
    echo "{"
    echo "  \"config\": \"$config\","
    echo "  \"stack\": \"$stack\","
    echo "  \"force\": $force,"
    echo "  \"dry\": $dry,"
    echo "  \"phases\": ["
    n=${#outputs[@]}
    i=0
    for o in "${outputs[@]}"; do
        i=$((i + 1))
        sep=","
        [[ $i -eq $n ]] && sep=""
        echo "    \"$o\"$sep"
    done
    echo "  ]"
    echo "}"
} > "$run_manifest"

log "all phases complete"
log "run manifest: $run_manifest"
printf '%s\n' "$run_manifest"
