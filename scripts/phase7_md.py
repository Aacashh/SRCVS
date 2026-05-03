from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml
from rdkit import Chem


@dataclass
class Cfg:
    root: Path
    out_root: Path
    log_dir: Path
    replicates: int


def load_cfg(path: Path) -> Cfg:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    root = Path(d["project"]["root"]).resolve()
    return Cfg(
        root=root,
        out_root=root / d["project"]["out_dir"],
        log_dir=root / d["project"]["log_dir"],
        replicates=int(d["md"]["replicates"]),
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


def split_leads(sdf_path: Path, out_dir: Path) -> list[tuple[str, Path]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    pairs: list[tuple[str, Path]] = []
    for i, m in enumerate(suppl):
        if m is None:
            continue
        name = (m.GetProp("_Name") if m.HasProp("_Name") else "").strip() or f"lead_{i:03d}"
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
        out_p = out_dir / f"{safe}.sdf"
        w = Chem.SDWriter(str(out_p))
        w.write(m)
        w.close()
        pairs.append((safe, out_p))
    return pairs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stack", default="oss", choices=["oss", "schrodinger"])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    cfg = load_cfg(a.config)
    out_dir = cfg.out_root / "07_md"
    log_path = cfg.log_dir / "phase7.log"
    mp = out_dir / "manifest.json"
    if mp.exists() and not a.force:
        sys.stdout.write(str(mp) + "\n")
        return 0

    _setup_log(log_path)
    logging.info("phase7 start stack=%s", a.stack)

    target_mf = _read_manifest(cfg.out_root / "01_target" / "manifest.json")
    select_mf = _read_manifest(cfg.out_root / "05_select" / "manifest.json")

    out_dir.mkdir(parents=True, exist_ok=True)
    leads = split_leads(Path(select_mf["outputs"]["leads_sdf"]), out_dir / "_leads")
    if not leads:
        sys.stderr.write("no leads to simulate\n")
        return 1

    runs = []
    for name, sdf in leads:
        for seed in range(1, cfg.replicates + 1):
            run_dir = out_dir / name / f"seed{seed}"
            run_dir.mkdir(parents=True, exist_ok=True)
            done_marker = run_dir / "_done"
            if done_marker.exists() and not a.force:
                runs.append({"lead": name, "seed": seed, "dir": str(run_dir)})
                continue
            if a.stack == "schrodinger":
                script = cfg.root / "scripts" / "schrodinger" / "desmond.sh"
                _run(["bash", str(script),
                      "--receptor", target_mf["outputs"]["receptor"],
                      "--ligand", str(sdf),
                      "--out-dir", str(run_dir),
                      "--seed", str(seed),
                      "--config", str(a.config)], log_path)
            else:
                script = cfg.root / "scripts" / "oss" / "md.py"
                _run([sys.executable, str(script),
                      "--receptor", target_mf["outputs"]["receptor"],
                      "--ligand", str(sdf),
                      "--out-dir", str(run_dir),
                      "--seed", str(seed),
                      "--config", str(a.config)], log_path)
            done_marker.touch()
            runs.append({"lead": name, "seed": seed, "dir": str(run_dir)})

    payload = {
        "phase": 7,
        "name": "molecular_dynamics",
        "stack": a.stack,
        "outputs": {"runs": runs},
    }
    p = write_manifest(out_dir, payload)
    sys.stdout.write(str(p) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
