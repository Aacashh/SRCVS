#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
prefix="${MINIFORGE_HOME:-$HOME/miniforge3}"
env_name="srcvs"

ensure_miniforge() {
    if command -v conda >/dev/null 2>&1; then
        prefix="$(conda info --base)"
        return 0
    fi
    if [[ -x "$prefix/bin/conda" ]]; then
        return 0
    fi
    local installer="/tmp/miniforge_installer.sh"
    local url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    if ! curl -sSL --fail --retry 3 -o "$installer" "$url"; then
        echo "failed to download miniforge" >&2
        exit 2
    fi
    bash "$installer" -b -p "$prefix" >/dev/null
    rm -f "$installer"
}

ensure_env() {
    local conda="$prefix/bin/conda"
    local mamba="$prefix/bin/mamba"
    if [[ ! -x "$conda" ]]; then
        echo "conda not executable: $conda" >&2
        exit 2
    fi
    if "$conda" env list | awk 'NR>2 {print $1}' | grep -qx "$env_name"; then
        return 0
    fi
    if [[ -x "$mamba" ]]; then
        "$mamba" env create -f "$root/environment.yml" -n "$env_name" >/dev/null
    else
        "$conda" env create -f "$root/environment.yml" -n "$env_name" >/dev/null
    fi
}

ensure_miniforge
ensure_env
