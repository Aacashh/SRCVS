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
    pdb_id: str
    redock_threshold: float


def load_cfg(path: Path) -> Cfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    root = Path(d["project"]["root"]).resolve()
    return Cfg(
        root=root,
        out_root=root / d["project"]["out_dir"],
        log_dir=root / d["project"]["log_dir"],
        raw_dir=root / d["project"]["raw_dir"],
        pdb_id=d["target"]["pdb_id"],
        redock_threshold=float(d["target"]["redock_rmsd_threshold"]),
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
    out_dir = cfg.out_root / "01_target"
    log_path = cfg.log_dir / "phase1.log"
    mp = out_dir / "manifest.json"
    if mp.exists() and not a.force:
        sys.stdout.write(str(mp) + "\n")
        return 0

    _setup_log(log_path)
    logging.info("phase1 start stack=%s", a.stack)

    pdb_path = cfg.raw_dir / f"{cfg.pdb_id}.pdb"
    if not pdb_path.exists():
        sys.stderr.write(f"missing pdb: {pdb_path}\n")
        return 3

    out_dir.mkdir(parents=True, exist_ok=True)
    if a.stack == "schrodinger":
        script = cfg.root / "scripts" / "schrodinger" / "prep.sh"
        _run(["bash", str(script), str(pdb_path), str(out_dir), str(a.config)], log_path)
        receptor = out_dir / "receptor_prepared.mae"
    else:
        script = cfg.root / "scripts" / "oss" / "prep.py"
        _run([sys.executable, str(script),
              "--pdb", str(pdb_path),
              "--out-dir", str(out_dir),
              "--config", str(a.config)], log_path)
        receptor = out_dir / "receptor_prepared.pdb"

    payload = {
        "phase": 1,
        "name": "target_preparation",
        "stack": a.stack,
        "outputs": {
            "receptor": str(receptor),
            "ligand": str(out_dir / "cocrystal_ligand.pdb"),
            "grid": str(out_dir / "grid.json"),
            "redock": str(out_dir / "redock_rmsd.csv"),
        },
    }
    p = write_manifest(out_dir, payload)
    sys.stdout.write(str(p) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
