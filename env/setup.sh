#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
prefix="${MINIFORGE_HOME:-$HOME/miniforge3}"
env_name="srcvs"

ts() { date +"%H:%M:%S"; }
log() { printf '[%s] env/setup: %s\n' "$(ts)" "$*"; }

detect_installer_url() {
    local os arch
    os="$(uname -s)"
    arch="$(uname -m)"
    case "$os" in
        Linux)
            case "$arch" in
                x86_64)  echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" ;;
                aarch64) echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh" ;;
                *) echo ""; return 1 ;;
            esac ;;
        Darwin)
            case "$arch" in
                arm64)  echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh" ;;
                x86_64) echo "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh" ;;
                *) echo ""; return 1 ;;
            esac ;;
        MINGW*|MSYS*|CYGWIN*)
            echo "windows-manual"; return 0 ;;
        *) echo ""; return 1 ;;
    esac
}

ensure_miniforge() {
    if command -v conda >/dev/null 2>&1; then
        prefix="$(conda info --base)"
        log "found existing conda at $prefix"
        return 0
    fi
    if [[ -x "$prefix/bin/conda" ]]; then
        log "using local miniforge at $prefix"
        return 0
    fi
    local url
    if ! url="$(detect_installer_url)" || [[ -z "$url" ]]; then
        echo "unsupported platform: $(uname -s) $(uname -m); install Miniforge manually" >&2
        exit 2
    fi
    if [[ "$url" == "windows-manual" ]]; then
        echo "windows detected: install Miniforge from https://conda-forge.org/download/ then re-run" >&2
        exit 2
    fi
    log "no conda found, downloading miniforge..."
    local installer="/tmp/miniforge_installer.sh"
    if ! curl -sSL --fail --retry 3 -o "$installer" "$url"; then
        echo "failed to download miniforge from $url" >&2
        exit 2
    fi
    log "installing miniforge to $prefix..."
    bash "$installer" -b -p "$prefix"
    rm -f "$installer"
}

env_exists() {
    local conda="$1"
    "$conda" env list 2>/dev/null | awk 'NR>2 {print $1}' | grep -qx "$env_name"
}

try_create() {
    local conda="$1"
    local mode="$2"
    local cmd=()
    case "$mode" in
        mamba)    cmd=("$prefix/bin/mamba" env create -f "$root/environment.yml" -n "$env_name") ;;
        libmamba) cmd=("$conda" env create -f "$root/environment.yml" -n "$env_name" --solver=libmamba) ;;
        classic)  cmd=("$conda" env create -f "$root/environment.yml" -n "$env_name") ;;
    esac
    log "attempting env create with: $mode"
    if "${cmd[@]}"; then
        log "env create succeeded with: $mode"
        return 0
    fi
    log "env create FAILED with: $mode (continuing to next strategy)"
    return 1
}

clean_caches() {
    local conda="$1"
    log "cleaning corrupt conda package caches..."
    "$conda" clean --packages --tarballs --index-cache --yes >/dev/null 2>&1 || true
    "$conda" env remove -n "$env_name" --yes >/dev/null 2>&1 || true
}

ensure_env() {
    local conda="$prefix/bin/conda"
    local mamba="$prefix/bin/mamba"
    if [[ ! -x "$conda" ]]; then
        echo "conda not executable: $conda" >&2
        exit 2
    fi
    if env_exists "$conda"; then
        log "env '$env_name' already exists, skipping create"
        return 0
    fi
    log "creating env '$env_name' (this may take 5-15 minutes)..."
    free -h 2>/dev/null | head -2 | sed 's/^/[setup mem] /' || true

    if [[ -x "$mamba" ]] && try_create "$conda" mamba; then return 0; fi
    if try_create "$conda" libmamba; then return 0; fi
    if try_create "$conda" classic; then return 0; fi

    log "first pass failed; cleaning caches and retrying once..."
    clean_caches "$conda"
    if [[ -x "$mamba" ]] && try_create "$conda" mamba; then return 0; fi
    if try_create "$conda" libmamba; then return 0; fi
    if try_create "$conda" classic; then return 0; fi

    echo "all env-create strategies failed; check memory (try: free -h) and disk space" >&2
    exit 1
}

ensure_miniforge
ensure_env
