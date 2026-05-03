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

export PROJECT_ROOT="$root"

bash "$root/env/setup.sh" >/dev/null

env_bin=""
if [[ -d "$HOME/miniforge3/envs/srcvs/bin" ]]; then
    env_bin="$HOME/miniforge3/envs/srcvs/bin"
elif command -v conda >/dev/null 2>&1; then
    env_bin="$(conda info --base 2>/dev/null)/envs/srcvs/bin"
fi
if [[ -n "$env_bin" && -x "$env_bin/python" ]]; then
    export PATH="$env_bin:$PATH"
fi

if [[ "$dry" == "0" ]]; then
    bash "$root/data/fetch.sh"
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

declare -a outputs=()
for p in "${phases[@]}"; do
    script="$root/scripts/${names[$p]}.py"
    if [[ ! -f "$script" ]]; then
        echo "missing phase script: $script" >&2
        exit 1
    fi
    if [[ "$dry" == "1" ]]; then
        outputs+=("$script --config $config --stack $stack")
        continue
    fi
    extra=()
    [[ "$force" == "1" ]] && extra+=("--force")
    out_path=$(python "$script" --config "$config" --stack "$stack" "${extra[@]}")
    outputs+=("$out_path")
done

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

printf '%s\n' "$run_manifest"
