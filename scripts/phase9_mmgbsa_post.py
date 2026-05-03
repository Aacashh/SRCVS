from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class Cfg:
    root: Path
    out_root: Path
    log_dir: Path


def load_cfg(path: Path) -> Cfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    root = Path(d["project"]["root"]).resolve()
    return Cfg(
        root=root,
        out_root=root / d["project"]["out_dir"],
        log_dir=root / d["project"]["log_dir"],
    )


def _setup_log(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=str(p), level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s", force=True,
    )


def _run(cmd: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as fh:
        subprocess.run(cmd, check=True, stdout=fh, stderr=fh)


def _read_manifest(p: Path) -> dict:
    if not p.exists():
        sys.stderr.write(f"missing dependency manifest: {p}\n")
        sys.exit(3)
    with open(p) as fh:
        return json.load(fh)


def write_manifest(out_dir: Path, payload: dict) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    mp = out_dir / "manifest.json"
    with open(mp, "w") as fh:
        json.dump(payload, fh, indent=2)
    return mp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stack", default="oss", choices=["oss", "schrodinger"])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    cfg = load_cfg(a.config)
    out_dir = cfg.out_root / "09_mmgbsa_post"
    log_path = cfg.log_dir / "phase9.log"
    mp = out_dir / "manifest.json"
    if mp.exists() and not a.force:
        sys.stdout.write(str(mp) + "\n")
        return 0

    _setup_log(log_path)
    logging.info("phase9 start stack=%s", a.stack)

    md_mf = _read_manifest(cfg.out_root / "07_md" / "manifest.json")
    runs_json = out_dir / "runs.json"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(runs_json, "w") as fh:
        json.dump(md_mf["outputs"]["runs"], fh)

    if a.stack == "schrodinger":
        script = cfg.root / "scripts" / "schrodinger" / "mmgbsa.sh"
        _run(["bash", str(script),
              "--mode", "post",
              "--runs-json", str(runs_json),
              "--out-dir", str(out_dir),
              "--config", str(a.config)], log_path)
    else:
        script = cfg.root / "scripts" / "oss" / "mmpbsa.py"
        _run([sys.executable, str(script),
              "--mode", "post",
              "--runs-json", str(runs_json),
              "--out-dir", str(out_dir),
              "--config", str(a.config)], log_path)

    payload = {
        "phase": 9,
        "name": "mmgbsa_post_md",
        "stack": a.stack,
        "outputs": {
            "csv": str(out_dir / "mmgbsa_post.csv"),
            "summary": str(out_dir / "mmgbsa_post_summary.csv"),
        },
    }
    p = write_manifest(out_dir, payload)
    sys.stdout.write(str(p) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
