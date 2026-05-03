## Prereqs
- Ubuntu 22.04 or later
- bash, curl, git
- Internet access for first run (PDB, PubChem, DUD-E)
- Optional: NVIDIA GPU + CUDA driver for OSS molecular dynamics
- Optional: Schrödinger Suite 2025-4 with `$SCHRODINGER` set, for `--stack schrodinger`
- Manual downloads required: `data/raw/nubbe_full.sdf`, `data/raw/enamine_np.sdf`

## Usage
```
bash run.sh
bash run.sh --phase 3
bash run.sh --stack schrodinger
bash run.sh --force
bash run.sh --dry
bash tests/smoke.sh
```

## Exit codes
- 0 success
- 1 generic failure
- 2 missing dependency
- 3 missing data file
- 4 bad config
