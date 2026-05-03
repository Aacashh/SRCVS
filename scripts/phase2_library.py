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
    raw_dir: Path
    sources: list
    include_controls: bool
    max_input_compounds: int


def load_cfg(path: Path) -> Cfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    root = Path(d["project"]["root"]).resolve()
    return Cfg(
        root=root,
        out_root=root / d["project"]["out_dir"],
        log_dir=root / d["project"]["log_dir"],
        raw_dir=root / d["project"]["raw_dir"],
        sources=d["library"]["sources"],
        include_controls=bool(d["library"].get("include_controls", True)),
        max_input_compounds=int(d["library"].get("max_input_compounds", 0)),
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
    out_dir = cfg.out_root / "02_library"
    log_path = cfg.log_dir / "phase2.log"
    mp = out_dir / "manifest.json"
    if mp.exists() and not a.force:
        sys.stdout.write(str(mp) + "\n")
        return 0

    _setup_log(log_path)
    logging.info("phase2 start stack=%s", a.stack)

    sdf_inputs: list[Path] = []
    for src in cfg.sources:
        p = cfg.raw_dir / src["file"]
        if not p.exists():
            sys.stderr.write(f"missing library: {p}\n")
            return 3
        sdf_inputs.append(p)
    if cfg.include_controls:
        sdf_inputs.extend(sorted(cfg.raw_dir.glob("control_*.sdf")))

    out_dir.mkdir(parents=True, exist_ok=True)
    if a.stack == "schrodinger":
        script = cfg.root / "scripts" / "schrodinger" / "ligprep.sh"
        _run(["bash", str(script), str(out_dir), str(a.config), *map(str, sdf_inputs)], log_path)
        prepared_sdf = out_dir / "library_prepared.maegz"
    else:
        script = cfg.root / "scripts" / "oss" / "ligprep.py"
        _run([sys.executable, str(script),
              "--out-dir", str(out_dir),
              "--config", str(a.config),
              "--inputs", *map(str, sdf_inputs)], log_path)
        prepared_sdf = out_dir / "library_prepared.sdf"

    payload = {
        "phase": 2,
        "name": "library_preparation",
        "stack": a.stack,
        "outputs": {
            "prepared": str(prepared_sdf),
            "pdbqt_dir": str(out_dir / "pdbqt"),
            "summary": str(out_dir / "ligprep_summary.csv"),
        },
    }
    p = write_manifest(out_dir, payload)
    sys.stdout.write(str(p) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
