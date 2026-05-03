# Research Directive (v2): Virtual Screening and Molecular Dynamics Investigation of Natural Compounds Against SRC Kinase for Chronic Pain Treatment

> **Aligned to the methodology of Akhilesh et al., 2022 (SphK2 / NuBBE) and Akhilesh et al., 2024 (Autotaxin / Enamine) — Tiwari Lab, IIT-BHU.**
> This is your master plan from data acquisition through thesis defense.

---

## 0. Why This Version Exists

Your two reference papers come from the same group at IIT-BHU and share an **identical methodological skeleton**:

| Step | Paper 1 (SphK2, 2022) | Paper 2 (Autotaxin, 2024) | Your Thesis (SRC) |
|------|----------------------|---------------------------|-------------------|
| Target | SphK2 | Autotaxin (ATX) | **SRC kinase** |
| Database | NuBBE (~2,221 Brazilian NPs) | Enamine phenotypic screening library | **NuBBE + Enamine** (recommended hybrid) |
| Protein prep | Schrödinger Protein Prep Wizard | Schrödinger Protein Prep Wizard | Same |
| Ligand prep | LigPrep / Epik (pH 7.0±2.0) | LigPrep / Epik (pH 7.0±2.0) | Same |
| Docking | Glide HTVS → SP → XP | Glide HTVS → SP → XP | Same |
| Filters | Lipinski, ADMET (QikProp), PAINS | Lipinski, ADMET (QikProp), PAINS | Same |
| Free energy | Prime MM-GBSA | Prime MM-GBSA | Same |
| MD | Desmond, 100 ns, OPLS, TIP3P | Desmond, 100 ns, OPLS, TIP3P | Same (with **OPLS4**, current) |

Your thesis is essentially a **target-translated extension** of this template. Your novelty comes from:
1. Choosing SRC kinase, with its distinct mechanism in NMDAR-mediated central sensitization
2. Hit chemistry (different scaffolds emerging from a different binding pocket)
3. CNS/BBB-aware lead prioritization (more important for SRC analgesia than SphK2/ATX)
4. Methodological refinement: triplicate MD, longer trajectory window, OPLS4 force field

---

## 1. Comprehensive Data & Resources Manifest

### 1.1 Software & Licenses

| Item | Version | Source | Required? |
|------|---------|--------|-----------|
| **Schrödinger Suite** | **2025-4** (latest) or 2024-x | IIT-BHU institutional license | **Mandatory** (used in both reference papers) |
| Maestro (GUI) | bundled | — | Mandatory |
| LigPrep + Epik | bundled | — | Mandatory |
| Glide (HTVS/SP/XP) | bundled | — | Mandatory |
| QikProp (ADMET) | bundled | — | Mandatory |
| Prime (MM-GBSA) | bundled | — | Mandatory |
| Desmond (MD) | bundled | — | Mandatory |
| **OPLS4** force field | bundled | — | Mandatory (upgrade from OPLS_2005 used in 2022 paper) |
| **PyMOL** | 3.x | open-source / academic | Recommended (figures) |
| **Discovery Studio Visualizer** | 2024 | BIOVIA, free | Recommended (2D interaction plots) |
| **SwissADME** | web | http://www.swissadme.ch | Recommended (BOILED-Egg BBB) |
| **ADMETlab 3.0** | web | https://admetlab3.scbdd.com | Recommended (cross-validation of QikProp) |
| **ProTox 3.0** | web | https://tox.charite.de | Recommended (toxicity classes) |
| **Python 3.10+ with RDKit** | conda-forge | open-source | For supplementary cheminformatics |

> **Action item:** Confirm with your supervisor / IIT-BHU IT that you have a **named user license seat** on the Schrödinger floating license server. Token contention with senior students kills weekends of work.

### 1.2 Computational Hardware

| Resource | Minimum | Recommended | Notes |
|----------|---------|-------------|-------|
| CPU cores | 16 | 32–64 | Glide HTVS scales linearly; more cores = faster |
| GPU | 1× RTX 3060 (12 GB) | RTX 4090 / A100 (40 GB) | **Mandatory for Desmond MD** |
| RAM | 32 GB | 64–128 GB | Large MD systems need ≥64 GB |
| Storage | 500 GB SSD | 2 TB NVMe | Trajectories are 5–20 GB each |
| Network | — | HPC access | IIT-BHU's PARAM Shivay or local lab cluster |

> **Time budget:** ~2 weeks of compute including queue waits, on a single GPU workstation.

### 1.3 Target (Receptor) Data

**Protein:** Proto-oncogene tyrosine-protein kinase Src (Human)
**UniProt:** P12931
**Gene:** SRC (chromosome 20q11.23)
**Length:** 536 residues; kinase domain: ~270–520

**Recommended PDB structures** (you will benchmark and pick one):

| PDB | Resolution | State | Co-crystal ligand | Notes |
|-----|-----------|-------|--------------------|-------|
| **3G5D** | 1.95 Å | Active | Dasatinib (1N1) | **Primary recommendation** — active conformation + clinical inhibitor for redock |
| **1Y57** | 1.91 Å | Active | None | Apo active form |
| **2SRC** | 1.50 Å | Inactive (autoinhibited) | AMP-PNP | Reference for full-length |
| **4MXO** | 2.22 Å | Active | ATP analog | ATP-pocket reference |
| **3F3V** | 2.30 Å | Active | PP2 | Tool compound benchmark |

### 1.4 Ligand Library Data

**Recommendation: Use NuBBE as your primary library (matching Paper 1) with Enamine as supplementary (matching Paper 2).** This dual-library approach gives you strongest narrative continuity with both papers.

| Database | Size | Coverage | URL | License |
|----------|------|----------|-----|---------|
| **NuBBE** | ~2,221 | Brazilian/Latin American NPs | https://nubbe.iq.unesp.br/portal/nubbe-search.html | Free, open |
| **Enamine NP screening library** | ~25,000 | Commercial NP-derived | https://enamine.net (request) | Free academic with request |
| COCONUT 2.0 | ~700,000 | Aggregator of 50+ NP DBs | https://coconut.naturalproducts.net | Free |
| IMPPAT 2.0 | ~17,967 | Indian medicinal plants | https://cb.imsc.res.in/imppat/ | Free, open |
| LOTUS | ~750,000 | Open, taxonomy-tagged | https://lotus.naturalproducts.net | Free |
| ZINC22 NP subset | ~250,000 | Vendor-purchasable | https://zinc22.docking.org | Free |
| NPASS v2.0 | ~30,000 | Bioactivity-annotated | http://bidd.group/NPASS/ | Free |

**Storage:** NuBBE SDF ~30 MB; Enamine NP library ~500 MB; if you also pull COCONUT, plan for ~2–4 GB.

### 1.5 Positive Controls (Reference Compounds)

You must include known SRC inhibitors as benchmarks. Dock them, MD them, MM-GBSA them — your hits should be **comparable or better**.

| Compound | Type | PubChem CID | Use |
|----------|------|-------------|-----|
| **Dasatinib** | FDA-approved SRC/Abl inhibitor | 3062316 | Re-docking validation + MD reference |
| **Bosutinib** | FDA-approved SRC/Abl | 5328940 | MD reference |
| **Saracatinib (AZD0530)** | Phase II clinical | 10302451 | MD reference (CNS-tested) |
| **PP1** | Tool compound | 4878 | Selectivity benchmark |
| **PP2** | Tool compound | 4879 | Selectivity benchmark (3F3V) |

Pull SDF files for all five from PubChem before you begin Phase 4.

### 1.6 Validation Datasets

| Dataset | Use |
|---------|-----|
| **DUD-E SRC subset** (http://dude.docking.org) | Active vs decoy ROC-AUC validation of your docking pipeline |
| **PDBbind general set** | Benchmark MM-GBSA correlation with experimental Ki/Kd |
| **ChEMBL CHEMBL267 (SRC)** | ~10,000 SRC bioassay records — confirm your hits aren't already known with poor activity |

### 1.7 Storage & Reproducibility Plan

```
~/thesis/                                    # 200–500 GB total
├── 00_admin/         (license logs, license_acks, weekly notes)
├── 01_target/        (PDBs, prepared .mae, redock validation)
├── 02_library/       (raw SDFs, prepared .maegz with states)
├── 03_glide/         (grid, dock outputs, _pv.maegz pose viewer files)
├── 04_qikprop/       (CSV outputs)
├── 05_mmgbsa/        (Prime output, decomposition CSVs)
├── 06_desmond/       (MD systems, trajectories ~10 GB each)
├── 07_analysis/      (SID reports, figures)
├── 08_figures/       (publication-ready PNGs/PDFs)
├── 09_thesis/        (LaTeX/Word, bibliography, supplements)
├── scripts/          (Python automation; this directive's code)
└── logs/             (job stdout/stderr)
```

Back up `01_target`, `02_library`, `scripts`, `09_thesis` to **two external drives + cloud**. Trajectories go to one external; you can regenerate them if lost.

End-of-thesis: archive everything (minus trajectories) to **Zenodo** for a citable DOI in your manuscript.

---

## 2. Pipeline at a Glance

```
                                                      ┌──────────────────────────┐
                                                      │   POSITIVE CONTROLS:     │
                                                      │  Dasatinib, Bosutinib,   │
                                                      │  Saracatinib, PP1, PP2   │
                                                      └──────────┬───────────────┘
                                                                 │
[ PDB 3G5D ]──► Protein Prep Wizard ──► Receptor Grid ──┐        │
                                                        ▼        ▼
[ NuBBE +Enamine ]──► LigPrep (Epik) ──► HTVS ──► SP ──► XP (top 0.5–1%)
                                                                 │
                                                                 ▼
                                       QikProp ADMET + Lipinski + PAINS
                                                                 │
                                                                 ▼
                                              Prime MM-GBSA on top 20–50
                                                                 │
                                                                 ▼
                                       Top 5–10 leads + 1 control selected
                                                                 │
                                                                 ▼
                                       Desmond MD (100 ns × 3 replicates)
                                                                 │
                                                                 ▼
                                       SID analysis + RMSD/RMSF/H-bonds + PCA
                                                                 │
                                                                 ▼
                                       Final MM-GBSA on MD ensemble
                                                                 │
                                                                 ▼
                                                Thesis figures + tables
```

---

## 3. Phase 1 — Project Setup & Environment

```bash
mkdir -p ~/thesis/{00_admin,01_target,02_library,03_glide,04_qikprop,05_mmgbsa,06_desmond,07_analysis,08_figures,09_thesis,scripts,logs}
cd ~/thesis

# Confirm Schrödinger is on PATH
echo $SCHRODINGER
$SCHRODINGER/licadmin STAT     # check available tokens
$SCHRODINGER/run python3 -c "import schrodinger; print(schrodinger.__version__)"
```

Create `scripts/env.sh`:

```bash
#!/bin/bash
export SCHRODINGER=/opt/schrodinger/2025-4    # adjust to your install
export PROJECT=$HOME/thesis
export PYTHONPATH=$PROJECT/scripts:$PYTHONPATH
alias srun="$SCHRODINGER/run"
alias spython="$SCHRODINGER/run python3"
```

Source it: `echo "source ~/thesis/scripts/env.sh" >> ~/.bashrc`

---

## 4. Phase 2 — Target Preparation

### 4.1 Strategy
Match the reference papers' Protein Preparation Wizard parameters: assign bond orders, add hydrogens, create disulfide bonds, fill missing side chains/loops via Prime, optimize H-bond network using PROPKA at pH 7.0, and minimize with OPLS4 (RMSD convergence 0.30 Å).

### 4.2 Script: `scripts/01_prepare_target.sh`

```bash
#!/bin/bash
# Download SRC kinase active structure (with dasatinib for redock validation)
cd ~/thesis/01_target

PDB="3G5D"
wget -q "https://files.rcsb.org/download/${PDB}.pdb" -O ${PDB}.pdb

# Run Schrödinger Protein Preparation Wizard
$SCHRODINGER/utilities/prepwizard \
    -fillsidechains \
    -fillloops \
    -disulfides \
    -mse \
    -propka_pH 7.0 \
    -minimize_adj_h \
    -f S-OPLS \
    -rmsd 0.30 \
    -WAIT \
    ${PDB}.pdb \
    ${PDB}_prep.mae

# Extract co-crystal dasatinib (resname 1N1) as redock reference
$SCHRODINGER/run split_structure.py \
    -m ligand \
    -groupby molecule \
    ${PDB}_prep.mae \
    ${PDB}_split.mae

# Build apo receptor (drop ligand, keep waters? — match Paper 1 by removing all waters)
spython << 'EOF'
from schrodinger import structure
st = next(structure.StructureReader("3G5D_prep.mae"))
# Remove HETATM (ligand + waters) but keep ions if catalytically relevant
to_delete = [a.index for a in st.atom
             if a.pdbres.strip() in ("HOH", "1N1", "GOL", "SO4")]
st.deleteAtoms(to_delete)
st.write("3G5D_apo.mae")
print(f"Apo SRC: {st.atom_total} atoms")
EOF

echo "Done. Files in ~/thesis/01_target/"
ls -la
```

### 4.3 Receptor Grid Generation: `scripts/02_grid.sh`

```bash
#!/bin/bash
cd ~/thesis/01_target

# Define box centered on co-crystal dasatinib (autodetected if specified)
cat > grid.in << 'EOF'
GRID_CENTER_USING_LIG  True
LIGAND_INDEX           1
LIGAND_FILE            3G5D_split-lig-1.mae
RECEP_FILE             3G5D_apo.mae
INNERBOX               12, 12, 12
OUTERBOX               24, 24, 24
GRIDFILE               src_grid.zip
FORCEFIELD             OPLS_2005
EOF

$SCHRODINGER/glide grid.in -OVERWRITE -WAIT -HOST localhost
echo "Grid created: src_grid.zip"
```

### 4.4 Re-docking Validation (mandatory before screening)

```bash
cat > redock.in << 'EOF'
GRIDFILE   src_grid.zip
LIGANDFILE 3G5D_split-lig-1.mae
PRECISION  XP
EOF

$SCHRODINGER/glide redock.in -OVERWRITE -WAIT
# Inspect: open redock_pv.maegz in Maestro, overlay docked vs crystal pose
# Required: heavy-atom RMSD ≤ 2.0 Å between top docked pose and crystal
```

**Quality gates before proceeding:**
- [ ] Catalytic residues present and correctly protonated: K298, E310, D386, F405 (DFG)
- [ ] Hinge region intact: E339, M341, T338
- [ ] Activation loop visible (residues 404–432) — if missing, used Prime to fill
- [ ] Re-dock RMSD ≤ 2.0 Å
- [ ] Grid box visually encompasses ATP pocket in Maestro

### 4.5 Thesis hooks (Chapter 4)
- Figure 4.1: Cartoon of prepared SRC with annotated catalytic residues
- Figure 4.2: Re-docking overlay (crystal vs predicted dasatinib pose)
- Table 4.1: PDB stats (resolution, R-factor, missing residues filled)

---

## 5. Phase 3 — Ligand Library Curation

### 5.1 Acquire Libraries

```bash
cd ~/thesis/02_library

# 1. NuBBE — register at https://nubbe.iq.unesp.br to download the SDF
#    Save the master SDF as: nubbe_full.sdf (~2,221 compounds)

# 2. Enamine NP screening library — request via https://enamine.net
#    Save as: enamine_np.sdf

# 3. (Optional) COCONUT for breadth
wget -q "https://coconut.naturalproducts.net/download/sdf" -O coconut_full.sdf
```

### 5.2 LigPrep with Epik (matches both papers' protocol)

`scripts/03_ligprep.sh`:

```bash
#!/bin/bash
cd ~/thesis/02_library

for INPUT in nubbe_full.sdf enamine_np.sdf; do
    BASE=$(basename $INPUT .sdf)
    $SCHRODINGER/ligprep \
        -isd $INPUT \
        -omae ${BASE}_prep.maegz \
        -epik \
        -ph 7.0 \
        -pht 2.0 \
        -bff 16 \
        -s 8 \
        -t 4 \
        -i 0 \
        -nt \
        -HOST localhost:8 \
        -WAIT
done
# -ph 7.0 -pht 2.0  →  ionize at pH 7.0 ± 2.0 (paper protocol)
# -bff 16           →  OPLS4 (use 14 for OPLS_2005 if matching 2022 paper exactly)
# -s 8              →  retain top 8 stereoisomers per input
# -t 4              →  generate up to 4 tautomers
# -i 0              →  do not desalt twice
```

Combine and report:

```bash
$SCHRODINGER/utilities/structcat \
    -imae nubbe_full_prep.maegz enamine_np_prep.maegz \
    -omae library_prepared.maegz

spython - << 'EOF'
from schrodinger import structure
n = sum(1 for _ in structure.StructureReader("library_prepared.maegz"))
print(f"Total prepared ligand states: {n}")
EOF
```

Typical output: ~7,000–9,000 states from ~5,000 unique input compounds.

### 5.3 Quality checks
- [ ] LigPrep log shows no parse errors
- [ ] Quercetin, luteolin, curcumin, resveratrol, EGCG appear in the prepared library (sanity check for known SRC binders)
- [ ] State count is reasonable (1.5–2× input compound count)

---

## 6. Phase 4 — Virtual Screening (Glide HTVS → SP → XP)

### 6.1 The Three-Tier Approach
Both reference papers use this exact funnel: HTVS keeps the top 10% by GlideScore, those are re-docked at SP, top 10% of SP go to XP.

### 6.2 Single-command Virtual Screening Workflow (`vsw`)

`scripts/04_screen.sh`:

```bash
#!/bin/bash
cd ~/thesis/03_glide

# Symbolic-link the inputs
ln -sf ../01_target/src_grid.zip .
ln -sf ../02_library/library_prepared.maegz .

cat > vsw.inp << 'EOF'
[SET:ORIGINAL_LIGANDS]
    VARCLASS   Structures
    FILES      library_prepared.maegz,
[STAGE:LIGPREP]
    STAGECLASS  ligprep.LigPrepStage
    INPUTS      ORIGINAL_LIGANDS,
    OUTPUTS     LIGPREP_OUT,
    RECOMBINE   YES
    RETITLE     YES
    SKIP_BAD_LIGANDS  YES
    USE_EPIK    YES
    METAL_BINDING NO
    PH          7.0
    PHT         2.0
    NEUTRALIZE  NO
    GENERATE_TAUTOMERS  YES
    NUM_STEREOISOMERS   8
[STAGE:HTVS]
    STAGECLASS         glide.DockingStage
    INPUTS             LIGPREP_OUT,
    OUTPUTS            HTVS_OUT,
    RECOMBINE          NO
    PRECISION          HTVS
    GRIDFILE           src_grid.zip
    LIGANDFILE         &LIGPREP_OUT;
    POSE_OUTTYPE       ligandlib_sd
    DOCKING_METHOD     confgen
    POSES_PER_LIG      1
[STAGE:SCORE_FILTER_HTVS]
    STAGECLASS         glide.ScoreFilterStage
    INPUTS             HTVS_OUT,
    OUTPUTS            SCORE_FILTER_HTVS_OUT,
    FILTER             "r_i_docking_score < 0.0"
    NUMBER_OUT         10%
[STAGE:SP]
    STAGECLASS         glide.DockingStage
    INPUTS             SCORE_FILTER_HTVS_OUT,
    OUTPUTS            SP_OUT,
    PRECISION          SP
    GRIDFILE           src_grid.zip
    POSES_PER_LIG      3
[STAGE:SCORE_FILTER_SP]
    STAGECLASS         glide.ScoreFilterStage
    INPUTS             SP_OUT,
    OUTPUTS            SCORE_FILTER_SP_OUT,
    NUMBER_OUT         10%
[STAGE:XP]
    STAGECLASS         glide.DockingStage
    INPUTS             SCORE_FILTER_SP_OUT,
    OUTPUTS            XP_OUT,
    PRECISION          XP
    GRIDFILE           src_grid.zip
    POSES_PER_LIG      5
    WRITE_XP_DESC      YES
[USEROUTS]
    USEROUTS XP_OUT,
    STRUCTOUT XP_OUT
EOF

$SCHRODINGER/vsw vsw.inp \
    -host_glide localhost:16 \
    -adjust \
    -OVERWRITE \
    -WAIT
```

> **Compute estimate:** ~9,000 compounds × HTVS (1 s) + 900 × SP (15 s) + 90 × XP (90 s) on 16 cores ≈ 4–6 hours.

### 6.3 Score Aggregation: `scripts/05_score_table.py`

```python
"""Run with: $SCHRODINGER/run python3 scripts/05_score_table.py"""
from pathlib import Path
import csv
from schrodinger import structure

OUT = Path("03_glide/xp_scores.csv")
inputs = list(Path("03_glide").glob("*XP_OUT*pv*.maegz"))

rows = []
for f in inputs:
    reader = structure.StructureReader(str(f))
    next(reader)            # first entry is the receptor
    for st in reader:
        rows.append({
            "title":           st.title,
            "glide_score":     st.property.get("r_i_docking_score"),
            "glide_emodel":    st.property.get("r_i_glide_emodel"),
            "glide_gscore":    st.property.get("r_i_glide_gscore"),
            "glide_evdw":      st.property.get("r_i_glide_evdw"),
            "glide_ecoul":     st.property.get("r_i_glide_ecoul"),
            "glide_eint":      st.property.get("r_i_glide_einternal"),
            "glide_ligand_efficiency": st.property.get("r_i_glide_ligand_efficiency"),
        })

rows.sort(key=lambda r: r["glide_score"])
with open(OUT, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
    writer.writeheader(); writer.writerows(rows)
print(f"Wrote {len(rows)} XP-ranked compounds → {OUT}")
print("Top 10:")
for r in rows[:10]: print(f"  {r['title']:<30s} {r['glide_score']:.3f}")
```

### 6.4 Quality Checks
- [ ] XP scores span ~−12 to −5 kcal/mol (broader than SP)
- [ ] Top hits visually inspected — H-bonds with hinge (E339 backbone, M341 backbone NH), salt bridge with K298 if charged ligand
- [ ] Dasatinib control re-ranks in top ~5% (validates funnel)
- [ ] DFG-in conformation preserved in receptor

### 6.5 Thesis Hooks (Chapter 6)
- Figure 6.1: Histogram of Glide scores at each tier (HTVS, SP, XP)
- Figure 6.2: Top-50 docked compounds with score, MW, source organism (NuBBE provides taxonomy!)
- Figure 6.3: 2D ligand interaction diagram (Maestro Ligand Interaction Diagram tool) for top 5
- Table 6.1: Top 20 hits — title, source organism, GlideScore, ΔE_vdW, ΔE_coul, ligand efficiency

---

## 7. Phase 5 — ADMET Profiling (QikProp + External)

### 7.1 QikProp on XP Hits

`scripts/06_qikprop.sh`:

```bash
#!/bin/bash
cd ~/thesis/04_qikprop
ln -sf ../03_glide/*XP_OUT*pv*.maegz xp_hits.maegz

$SCHRODINGER/qikprop \
    -fast \
    -outname src_qikprop \
    xp_hits.maegz \
    -WAIT
# Output: src_qikprop.csv with 50+ ADMET descriptors per compound
```

### 7.2 Filter Combining Lipinski, Veber, PAINS, CNS-MPO

`scripts/07_filter_admet.py`:

```python
"""Run with: spython scripts/07_filter_admet.py"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog, AllChem
from rdkit.Chem.FilterCatalog import FilterCatalogParams

qp = pd.read_csv("04_qikprop/src_qikprop.csv")
glide = pd.read_csv("03_glide/xp_scores.csv")
df = qp.merge(glide, left_on="molecule", right_on="title")

# QikProp acceptable ranges (from Schrödinger documentation)
# These match what Akhilesh et al. used as primary ADMET filters
df["lipinski_pass"]   = df["#stars"] <= 2
df["ro5_pass"]        = (df["mol_MW"] <= 500) & (df["QPlogPo/w"] <= 5) & \
                        (df["accptHB"] <= 10) & (df["donorHB"] <= 5)
df["veber_pass"]      = (df["PSA"] <= 140) & (df["#rotor"] <= 10)
df["pK_HSA_ok"]       = df["QPlogKhsa"].between(-1.5, 1.5)    # serum-binding
df["log_BB_ok"]       = df["QPlogBB"].between(-3.0, 1.2)      # CNS — IMPORTANT for analgesic
df["caco2_ok"]        = df["QPPCaco"] > 25                    # GI absorption
df["hERG_ok"]         = df["QPlogHERG"] > -5                  # cardiac safety
df["mdck_ok"]         = df["QPPMDCK"] > 25                    # BBB proxy
df["cns_score"]       = (df["log_BB_ok"].astype(int)
                       + (df["mol_MW"] <= 360).astype(int)
                       + (df["QPlogPo/w"] <= 3).astype(int)
                       + (df["PSA"] <= 90).astype(int)
                       + (df["donorHB"] <= 1).astype(int))

# PAINS filter via RDKit
params = FilterCatalogParams()
for pf in [FilterCatalogParams.FilterCatalogs.PAINS_A,
           FilterCatalogParams.FilterCatalogs.PAINS_B,
           FilterCatalogParams.FilterCatalogs.PAINS_C]:
    params.AddCatalog(pf)
pains = FilterCatalog.FilterCatalog(params)

def pains_clean(smi):
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    return False if m is None else not pains.HasMatch(m)
df["pains_pass"] = df["SMILES"].apply(pains_clean)

passed = df[df["lipinski_pass"] & df["ro5_pass"] & df["veber_pass"]
          & df["pK_HSA_ok"] & df["log_BB_ok"] & df["caco2_ok"]
          & df["hERG_ok"] & df["mdck_ok"] & df["pains_pass"]]
passed = passed.sort_values("glide_score").head(50)
passed.to_csv("04_qikprop/admet_passed.csv", index=False)
print(f"Passed: {len(passed)}/{len(df)}")
print(passed[["molecule", "glide_score", "QPlogBB", "QPlogHERG", "cns_score"]].head(15))
```

### 7.3 Cross-validation with web servers (recommended)

Submit `admet_passed.csv` SMILES to:
- **SwissADME** — for BOILED-Egg plot (BBB + GI panels) → Figure 7.1
- **ADMETlab 3.0** — for AMES, hepatotoxicity, BBB classifier (cross-checks QikProp)
- **ProTox 3.0** — predicted LD50 and toxicity classes (organ-specific)

Add results as additional columns to `admet_passed.csv` manually or via APIs.

### 7.4 Diversity Selection of 10 Final Leads (Butina clustering)

```python
# scripts/08_select_leads.py
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina

df = pd.read_csv("04_qikprop/admet_passed.csv").reset_index(drop=True)
mols = [Chem.MolFromSmiles(s) for s in df["SMILES"]]
fps  = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in mols]

n = len(fps); dists = []
for i in range(1, n):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend([1 - x for x in sims])
clusters = Butina.ClusterData(dists, n, 0.4, isDistData=True)

# Pick the best (lowest GlideScore) representative from each of the top 10 largest clusters
representatives = []
for cluster in sorted(clusters, key=len, reverse=True)[:10]:
    best = df.iloc[list(cluster)].sort_values("glide_score").iloc[0]
    representatives.append(best)
final = pd.DataFrame(representatives).reset_index(drop=True)
final.to_csv("04_qikprop/final_leads.csv", index=False)
print(final[["molecule", "glide_score", "QPlogBB", "cns_score"]])
```

### 7.5 Thesis Hooks (Chapter 7)
- Figure 7.1: SwissADME BOILED-Egg plot of leads
- Figure 7.2: Radar chart of normalized physicochemical properties
- Table 7.1: Full ADMET profile per lead
- Table 7.2: Diversity / Tanimoto similarity matrix among 10 finalists

---

## 8. Phase 6 — Prime MM-GBSA (Pre-MD Free Energy)

### 8.1 Why First Round Pre-MD
Both reference papers run MM-GBSA on docked poses **before** committing to MD — it's a cheap re-rank that catches bad poses Glide ranked highly.

### 8.2 Script: `scripts/09_prime_mmgbsa.sh`

```bash
#!/bin/bash
cd ~/thesis/05_mmgbsa

# Extract pose-viewer files for the final 10 leads
spython - << 'EOF'
import pandas as pd
from schrodinger import structure
leads = pd.read_csv("../04_qikprop/final_leads.csv")["molecule"].tolist()
pv_in  = "../03_glide/XP_OUT_pv.maegz"
pv_out = "leads_pv.maegz"
writer = structure.StructureWriter(pv_out)
for i, st in enumerate(structure.StructureReader(pv_in)):
    if i == 0 or st.title in leads:
        writer.append(st)
writer.close()
EOF

$SCHRODINGER/prime_mmgbsa \
    -ligand "(res.ptype UNK) OR (res.ptype LIG)" \
    -out_type COMPLEX \
    -job_type REAL_MIN \
    -csv_output yes \
    -jobname src_mmgbsa_pre \
    leads_pv.maegz \
    -HOST localhost:8 \
    -WAIT

# Output: src_mmgbsa_pre-out.csv with MMGBSA dG_Bind, vdW, Coulomb, etc.
```

Inspect: `src_mmgbsa_pre-out.csv` reports `r_psp_MMGBSA_dG_Bind`. Akhilesh et al. (2024) reported ΔG values from ~−25 to −75 kcal/mol — this is the **gas-phase + implicit-solvent** value, which is far more negative than experimental ΔG and only useful for rank-ordering.

### 8.3 Re-rank and Confirm Final Leads

If MM-GBSA disagrees significantly with GlideScore, drop the worst MM-GBSA performers (likely strained poses) and bring up the next-ranked candidates from Phase 5 to maintain N=5 leads + dasatinib control going into MD.

---

## 9. Phase 7 — Molecular Dynamics with Desmond

### 9.1 Match the Lab's Exact Protocol

From the reference papers' methods sections:
- **Box:** Orthorhombic, 10 Å buffer (10 × 10 × 10 Å)
- **Solvent:** TIP3P water
- **Force field:** OPLS_2005 in 2022; **upgrade to OPLS4** for your work
- **Neutralization:** Na+ / Cl− to 0.15 M
- **Length:** 100 ns production
- **Output:** Default Desmond .cms + .trj

### 9.2 System Builder + Multisim: `scripts/10_md_setup.py`

Run **once per lead + control = 6 systems**.

```python
"""
spython scripts/10_md_setup.py NP_001
Builds Desmond system from the MM-GBSA output, sets up relaxation + production.
"""
import sys, subprocess
from pathlib import Path
from schrodinger import structure
from schrodinger.application.desmond.packages import system_builder_inp

LEAD = sys.argv[1]
WORK = Path(f"06_desmond/{LEAD}"); WORK.mkdir(parents=True, exist_ok=True)

# 1. Pull complex from MM-GBSA output
src_pv = "05_mmgbsa/src_mmgbsa_pre-out.maegz"
out_complex = WORK / f"{LEAD}_complex.mae"
with structure.StructureWriter(str(out_complex)) as w:
    for st in structure.StructureReader(src_pv):
        if st.title == LEAD:
            w.append(st); break

# 2. System Builder configuration (matches papers; OPLS4 upgrade)
msj_path = WORK / "system_build.msj"
msj_path.write_text(f"""\
task {{
  task = "desmond:auto"
}}
build_geometry {{
  box = {{
    shape = orthorhombic
    size = [10.0 10.0 10.0]
    size_type = buffer
  }}
  solvent = TIP3P
  add_counterion = {{ ion = Na species = Cl }}
  salt = {{ concentration = 0.15 negative_ion = Cl positive_ion = Na }}
  rezero_system = false
}}
assign_forcefield {{
  forcefield = OPLS_2005    # change to "S-OPLS" for OPLS4 if licensed
}}
""")

# 3. Production MD .msj — 100 ns, NPT, 300 K, 1.01325 bar
prod_msj = WORK / "prod.msj"
prod_msj.write_text("""\
task { task = "desmond:auto" }
simulate {
  title = "Brownian Minimization"
  annealing = off
  time = 100
  timestep = [0.001 0.001 0.003]
  temperature = 10.0
  ensemble = { class = "NVT" method = "Brownie" brownie = { delta_max = 0.1 } }
  restraints.new = [{ name = posre_harm atoms = solute_heavy_atom force_constants = 50.0 }]
}
simulate {
  title = "NVT, T = 10K, restraints solute_heavy"
  effect_if = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
  time = 12
  temperature = 10.0
  restraints.new = [{ name = posre_harm atoms = solute_heavy_atom force_constants = 50.0 }]
  ensemble = { class = NVT method = Berendsen thermostat.tau = 0.1 }
}
simulate {
  title = "NPT, T = 10K, restraints solute_heavy"
  effect_if = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = MTK']
  time = 12
  pressure = [1.01325 isotropic]
  temperature = 10.0
  restraints.new = [{ name = posre_harm atoms = solute_heavy_atom force_constants = 50.0 }]
  ensemble = { class = NPT method = Berendsen thermostat.tau = 0.1 barostat.tau = 50.0 }
}
simulate {
  title = "NPT, T = 300K, restraints solute_heavy"
  effect_if = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = MTK']
  time = 12
  temperature = 300.0
  restraints.new = [{ name = posre_harm atoms = solute_heavy_atom force_constants = 50.0 }]
  ensemble = { class = NPT method = Berendsen thermostat.tau = 0.1 barostat.tau = 50.0 }
}
simulate {
  title = "NPT, no restraints, equilibration"
  effect_if = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = MTK']
  time = 24
  temperature = 300.0
  ensemble = { class = NPT method = Berendsen thermostat.tau = 0.1 barostat.tau = 2.0 }
}
simulate {
  title = "Production: 100 ns NPT"
  cpu = 1
  dir = "."
  time = 100000.0     # 100 ns in ps
  elapsed_time = 0
  temperature = 300.0
  pressure = [1.01325 isotropic]
  ensemble = { class = NPT method = MTK thermostat.tau = 1.0 barostat.tau = 2.0 }
  trajectory.center = solute
  trajectory.first = 0.0
  trajectory.format = dtr
  trajectory.frames_per_file = 250
  trajectory.interval = 100.0       # write every 100 ps -> 1000 frames
  trajectory.periodicfix = true
  trajectory.write_velocity = false
}
""")

# 4. Build & launch
import os
os.chdir(WORK)
subprocess.run([f"{os.environ['SCHRODINGER']}/utilities/multisim",
                "-JOBNAME", f"{LEAD}_md",
                "-m", "system_build.msj",
                "-c", "system_build.cfg" if Path("system_build.cfg").exists() else None,
                "-i", str(out_complex.name),
                "-o", f"{LEAD}_system-out.cms",
                "-HOST", "localhost",
                "-maxjob", "1",
                "-WAIT"], check=False)

subprocess.run([f"{os.environ['SCHRODINGER']}/utilities/multisim",
                "-JOBNAME", f"{LEAD}_md_prod",
                "-m", "prod.msj",
                "-i", f"{LEAD}_system-out.cms",
                "-o", f"{LEAD}_md-out.cms",
                "-HOST", "localhost:1:gpgpu=1",
                "-maxjob", "1",
                "-WAIT"], check=False)
print(f"MD complete for {LEAD}")
```

> **Compute:** ~12–24 hr per 100 ns on RTX 3080+, ~6 systems = 3–6 days.
> **Replicates (strongly recommended):** Run each system 3× with different random seeds for thesis-defensible statistics. That's ~2 weeks of GPU.

### 9.3 Trajectory Triage Loop

```bash
# scripts/11_run_all_md.sh
#!/bin/bash
LEADS=( "NUBBE_xxx" "NUBBE_yyy" "ENAMINE_zzz" "ENAMINE_aaa" "ENAMINE_bbb" "DASATINIB" )
for L in "${LEADS[@]}"; do
    for SEED in 1 2 3; do
        export DESMOND_SEED=$SEED
        spython scripts/10_md_setup.py "${L}_seed${SEED}"
    done
done
```

---

## 10. Phase 8 — Trajectory Analysis (Simulation Interactions Diagram + Custom)

### 10.1 SID — Built-in Schrödinger Analysis

`scripts/12_sid.sh`:

```bash
#!/bin/bash
for LEAD_DIR in 06_desmond/*/; do
    L=$(basename $LEAD_DIR)
    cd $LEAD_DIR
    $SCHRODINGER/run event_analysis.py analyze \
        -p ${L}_md-out.cms \
        -t ${L}_md-out_trj \
        -o ${L}_sid.eaf \
        -ASL_protein "protein" \
        -ASL_ligand "res.ptype UNK or (chain L)"
    $SCHRODINGER/run event_analysis.py report \
        ${L}_sid.eaf \
        -pdf ${L}_sid_report.pdf
    cd -
done
```

This produces the same RMSD, RMSF, P-L contacts, ligand torsion plots that appear in Figures 4–7 of both reference papers.

### 10.2 Custom Cross-Lead Comparison: `scripts/13_compare.py`

```python
"""
spython scripts/13_compare.py
Aggregates RMSD, RMSF, H-bonds, Rg across all leads and replicates → comparative figures.
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt, seaborn as sns
from schrodinger.application.desmond.packages import analysis, traj, topo

sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)

LEADS = ["NUBBE_xxx", "NUBBE_yyy", "ENAMINE_zzz", "ENAMINE_aaa", "ENAMINE_bbb", "DASATINIB"]
SEEDS = [1, 2, 3]

records = {}
for L in LEADS:
    rmsd_runs = []; rmsf_runs = []; hb_runs = []
    for s in SEEDS:
        cms = f"06_desmond/{L}_seed{s}/{L}_seed{s}_md-out.cms"
        if not Path(cms).exists(): continue
        msys, cms_model = topo.read_cms(cms)
        tr = traj.read_traj(cms.replace("-out.cms", "-out_trj"))

        # RMSD: ligand
        lig_aids = cms_model.select_atom("res.ptype UNK")
        prot_ca  = cms_model.select_atom("protein and atom.ptype CA")
        rmsd_calc = analysis.RMSD(msys, cms_model, prot_ca, fit_aids=prot_ca,
                                  ref_aids=prot_ca, ref_pos=cms_model.getXYZ())
        rmsd_lig  = analysis.RMSD(msys, cms_model, lig_aids, fit_aids=prot_ca,
                                  ref_aids=lig_aids, ref_pos=cms_model.getXYZ())
        rmsd_p, rmsd_l = analysis.analyze(tr, rmsd_calc, rmsd_lig)
        rmsd_runs.append(np.array(rmsd_l))

        # RMSF
        rmsf_calc = analysis.ProteinRMSF(msys, cms_model, prot_ca, fit_aids=prot_ca,
                                         ref_pos=cms_model.getXYZ())
        rmsf = analysis.analyze(tr, rmsf_calc)[0]
        rmsf_runs.append(np.array(rmsf))

        # H-bonds
        hb = analysis.HydrogenBondFinder(msys, cms_model,
                                         aids1=prot_ca, aids2=lig_aids)
        hb_n = [len(x) for x in analysis.analyze(tr, hb)[0]]
        hb_runs.append(np.array(hb_n))

    records[L] = {"rmsd": np.mean(rmsd_runs, axis=0),
                  "rmsd_sd": np.std(rmsd_runs, axis=0),
                  "rmsf": np.mean(rmsf_runs, axis=0),
                  "hb": np.mean(hb_runs, axis=0)}

# Comparative ligand-RMSD plot
fig, ax = plt.subplots(figsize=(8, 4.5))
for L, d in records.items():
    t = np.arange(len(d["rmsd"])) * 0.1   # 100 ps stride → ns
    ax.plot(t, d["rmsd"], label=L, lw=1.4)
    ax.fill_between(t, d["rmsd"]-d["rmsd_sd"], d["rmsd"]+d["rmsd_sd"], alpha=0.15)
ax.set(xlabel="Time (ns)", ylabel="Ligand RMSD (Å)",
       title="Comparative ligand stability across leads (mean ± SD, n=3)")
ax.legend(loc="upper right", ncol=2, fontsize=8)
fig.tight_layout(); fig.savefig("07_analysis/comparative_ligand_rmsd.png", dpi=300)

# Comparative RMSF
fig, ax = plt.subplots(figsize=(10, 4))
for L, d in records.items():
    ax.plot(np.arange(len(d["rmsf"])) + 1, d["rmsf"], label=L, lw=1)
for r, lab in [(298, "K298"), (310, "E310"), (339, "E339"), (386, "D386(DFG)"), (404, "AL")]:
    ax.axvline(r, ls="--", c="grey", lw=0.5)
    ax.text(r, ax.get_ylim()[1]*0.85, lab, fontsize=7, rotation=90)
ax.set(xlabel="Residue", ylabel="RMSF (Å)", title="Per-residue flexibility")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("07_analysis/comparative_rmsf.png", dpi=300)

# H-bond bar chart (mean over last 50 ns)
hb_means = {L: float(np.mean(d["hb"][len(d["hb"])//2:])) for L, d in records.items()}
fig, ax = plt.subplots(figsize=(6, 4))
sns.barplot(x=list(hb_means.keys()), y=list(hb_means.values()), ax=ax)
ax.set(ylabel="Mean # H-bonds (last 50 ns)", xlabel="")
plt.xticks(rotation=30, ha="right")
fig.tight_layout(); fig.savefig("07_analysis/hbond_means.png", dpi=300)
```

### 10.3 Quality Checks
- [ ] Ligand RMSD plateau ≤ 3 Å in last 50 ns for at least 3 leads
- [ ] Protein Cα RMSD plateau ≤ 3 Å for all systems
- [ ] H-bond pattern includes hinge residues (E339, M341 backbone)
- [ ] Activation loop RMSF < 4 Å (suggests inhibition of aberrant motion)
- [ ] DFG-in conformation maintained throughout

### 10.4 Thesis Hooks (Chapter 8)
- Figure 8.1–8.5: SID-style 4-panel per lead (RMSD, RMSF, P-L contacts, ligand torsions)
- Figure 8.6: Comparative ligand RMSD across all leads + apo
- Figure 8.7: Comparative RMSF
- Figure 8.8: H-bond occupancy heatmap
- Table 8.1: Mean RMSD/RMSF/H-bonds in stable window with SD across replicates

---

## 11. Phase 9 — Post-MD MM-GBSA & Per-Residue Decomposition

### 11.1 Re-run MM-GBSA on MD Frames

`scripts/14_mmgbsa_md.sh`:

```bash
#!/bin/bash
for LEAD_DIR in 06_desmond/*/; do
    L=$(basename $LEAD_DIR)
    OUT=07_analysis/mmgbsa_md/$L
    mkdir -p $OUT
    # Extract frames every 1 ns from last 50 ns (50 frames)
    $SCHRODINGER/run trj_extract_subsystem.py \
        $LEAD_DIR/${L}_md-out.cms \
        -t $LEAD_DIR/${L}_md-out_trj \
        -ASL "all" \
        -frames 500:1000:10 \
        -out $OUT/${L}_frames.maegz
    # Run MM-GBSA on the frame ensemble
    $SCHRODINGER/prime_mmgbsa \
        -ligand "res.ptype UNK" \
        -job_type ENERGY \
        -out_type COMPLEX \
        -csv_output yes \
        -jobname ${L}_mmgbsa_post \
        $OUT/${L}_frames.maegz \
        -HOST localhost:8 -WAIT
done
```

### 11.2 Aggregation: `scripts/15_summary.py`

```python
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from pathlib import Path
sns.set_theme(style="whitegrid", context="paper")

rows = []
for csv in Path("07_analysis/mmgbsa_md").rglob("*post-out.csv"):
    df = pd.read_csv(csv)
    L = csv.stem.replace("_mmgbsa_post-out", "")
    rows.append({
        "lead": L,
        "dG_mean": df["r_psp_MMGBSA_dG_Bind"].mean(),
        "dG_sd":   df["r_psp_MMGBSA_dG_Bind"].std(),
        "dG_vdW":  df["r_psp_MMGBSA_dG_Bind_vdW"].mean(),
        "dG_Coul": df["r_psp_MMGBSA_dG_Bind_Coulomb"].mean(),
        "dG_Solv_GB": df["r_psp_MMGBSA_dG_Bind_Solv_GB"].mean(),
        "dG_HB":   df.get("r_psp_MMGBSA_dG_Bind_Hbond", pd.Series([0])).mean()
    })
summary = pd.DataFrame(rows).sort_values("dG_mean")
summary.to_csv("07_analysis/mmgbsa_summary.csv", index=False)

fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(summary["lead"], summary["dG_mean"],
       yerr=summary["dG_sd"], color="steelblue", edgecolor="black")
ax.set(ylabel=r"$\Delta G_{bind}$ (MM-GBSA, kcal/mol)", xlabel="")
plt.xticks(rotation=30, ha="right")
fig.tight_layout(); fig.savefig("07_analysis/mmgbsa_bar.png", dpi=300)
print(summary)
```

### 11.3 Quality Checks
- [ ] At least 3 of your 5 leads have mean ΔG within ±10 kcal/mol of dasatinib
- [ ] Per-frame SD < 8 kcal/mol (indicates stable binding)
- [ ] Top contributing residues (per-residue decomposition) include hinge (E339, M341), gatekeeper (T338), DFG (D386, F405)

### 11.4 Thesis Hooks (Chapter 9)
- Figure 9.1: MM-GBSA ΔG bar chart with SD error bars
- Figure 9.2: Per-residue decomposition heatmap
- Table 9.1: Energy term breakdown (vdW, Coul, GB-solvation, H-bond, lipophilic, packing) per lead
- Table 9.2: Comparison to literature SRC inhibitor experimental Ki/Kd

---

## 12. Phase 10 — Mechanism, Figures, Thesis Compilation

### 12.1 Mechanistic Discussion Framework

For each top lead, address in your Discussion (Chapter 10):

1. **Pocket occupancy:** Hinge engagement — does it form ≥2 H-bonds with E339/M341 backbone? (Type I behavior)
2. **DFG state:** Maintained DFG-in throughout, or did the lead trigger DFG-out (Type II)?
3. **αC-helix:** Glu310 stable in salt bridge with Lys298?
4. **Activation loop:** Y419 phosphorylation site stabilized? (Tie this to the Salter/Pitcher 2012 mechanism — SRC-NMDAR-GluN2B Y1472 axis.)
5. **Selectivity hypothesis:** Compare hit's binding-pocket interactions to those of dasatinib (broad), saracatinib (more SRC-selective), and PP2 (selective tool). Is your hit hinting at SRC-family selectivity?
6. **Predicted analgesic mechanism:** Bridge from kinase inhibition → reduced GluN2B Y1472 phosphorylation → attenuated central sensitization → predicted analgesic effect.

### 12.2 Final Figure Production

`scripts/16_pymol_renders.pml`:

```python
# pymol -cq scripts/16_pymol_renders.pml
load 06_desmond/NP_001/cluster_centroid.pdb, complex
hide everything
show cartoon, polymer
color slate, polymer
select pocket, byres polymer within 5 of resn UNK
show sticks, pocket
color grey80, pocket
util.cnc pocket
show sticks, resn UNK
color magenta, resn UNK
util.cnc resn UNK
# Annotate key residues
label resi 298 and name CA, "K298"
label resi 339 and name CA, "E339"
label resi 341 and name CA, "M341"
label resi 386 and name CA, "D386"
set ray_opaque_background, off
bg_color white
ray 1600, 1200
png 08_figures/NP_001_complex.png, dpi=300
```

For 2D ligand-interaction diagrams, open the SID `.eaf` in Maestro → "Ligand Interactions Diagram" → export as PNG; or use **Discovery Studio Visualizer** → Receptor-Ligand Interactions → Show 2D Diagram.

### 12.3 Mandatory Thesis Figures Checklist
- [ ] Pipeline schematic (Methods chapter, Figure 3.1 — TOC graphic)
- [ ] SRC structure with annotated pocket residues
- [ ] Library property distributions (pre/post LigPrep)
- [ ] Re-docking validation overlay
- [ ] HTVS→SP→XP attrition funnel
- [ ] Top-hits table + score histogram
- [ ] BOILED-Egg ADMET plot
- [ ] Radar chart of physicochemical properties
- [ ] SID 4-panel per lead (×5)
- [ ] Comparative RMSD/RMSF/Rg across leads
- [ ] H-bond occupancy heatmap
- [ ] MM-GBSA ΔG bar chart
- [ ] Per-residue decomposition heatmap
- [ ] 2D interaction diagram per lead
- [ ] PyMOL rendering of best lead bound

---

## 13. Thesis Structure (Mapping)

| Chapter | Title | Phase | Pages (target) |
|---------|-------|-------|----------------|
| 1 | Introduction: Chronic Pain — Burden & Unmet Need | — | 8–12 |
| 2 | SRC Kinase: Structure, Function, and Role in Nociception | — | 12–15 |
| 3 | Computational Drug Discovery — Rationale and Methods Overview | All | 10–12 |
| 4 | Target Preparation and Validation | 4 | 6–8 |
| 5 | Natural Compound Library Curation | 5 | 6–8 |
| 6 | Structure-Based Virtual Screening | 6 | 10–14 |
| 7 | ADMET Profiling and Lead Selection | 7 | 8–10 |
| 8 | Molecular Dynamics Investigation | 8 | 12–18 |
| 9 | MM-GBSA Free-Energy Analysis | 9, 11 | 8–10 |
| 10 | Discussion: Mechanistic Insights and Therapeutic Implications | 12 | 10–14 |
| 11 | Conclusions and Future Directions | — | 4–6 |
| 12 | References, Supplementary Information | — | — |

**Total target:** 100–130 pages excluding supplementary information.

### 13.1 Citing Your Reference Papers

You will cite Akhilesh et al. 2022 and 2024 throughout — particularly in your Methods chapter (justify each step by citing precedent), Discussion (compare your SRC findings to their SphK2/ATX findings as part of a "growing portfolio of computational analgesic discovery in this lab"), and Introduction (establish lineage).

---

## 14. Timeline (6-Month Plan)

| Month | Weeks | Deliverable |
|-------|-------|-------------|
| 1 | 1 | Literature review; environment setup; license confirmed |
| 1 | 2 | Target prep; redock validation; pilot Glide run |
| 1 | 3–4 | Library acquisition; LigPrep |
| 2 | 5–6 | Full HTVS→SP→XP screen |
| 2 | 7–8 | QikProp + ADMET + lead selection |
| 3 | 9–10 | Pre-MD MM-GBSA; system building for MD |
| 3 | 11–12 | First MD replicate set (5 leads + control × 1 seed) |
| 4 | 13–14 | Replicate seeds 2 and 3; SID analysis |
| 4 | 15–16 | Post-MD MM-GBSA; comparative analysis |
| 5 | 17–18 | All figures finalized; first thesis draft (Chapters 4–9) |
| 5 | 19–20 | Discussion (Ch 10), Introduction (Ch 1–2) drafts |
| 6 | 21–22 | Revisions, supervisor feedback iterations |
| 6 | 23–24 | Defense rehearsal; final binding & submission |

---

## 15. Risk Register & Mitigations

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Schrödinger token contention | Medium | High | Run jobs at night/weekends; coordinate with lab; maintain open-source backup pipeline |
| MD instability | Low–Medium | High | Run 3 replicates; the Desmond default protocol is robust |
| Hits show no SRC selectivity | Medium | Medium | Frame thesis around analgesic-relevant kinome; include selectivity discussion |
| Re-dock RMSD > 2 Å | Low | Critical | Try alternative PDB (1Y57, 4MXO); revisit Protein Prep settings |
| Library too small after filters | Low | Medium | Loosen Veber/CNS thresholds; supplement with COCONUT |
| GPU unavailability | Medium | High | Reduce replicate count; use CPU Desmond (~10× slower) |
| Reviewer demands experimental validation | High | Medium | Pre-emptive "Future Directions" plan: kinase activity assay + DRG neuron pY1472 western |

---

## 16. Open-Source Mirror (Backup Plan)

If Schrödinger access falls through mid-project, you can pivot to an open-source equivalent that still maps cleanly onto the methods sections of your reference papers:

| Schrödinger Module | Open-Source Equivalent |
|---------------------|------------------------|
| Protein Preparation Wizard | PDBFixer + PROPKA3 + PyMOL |
| LigPrep (Epik) | RDKit + Dimorphite-DL |
| Glide (HTVS/SP/XP) | AutoDock Vina + Smina or DiffDock |
| QikProp | RDKit descriptors + ADMETlab 3.0 (web) |
| Prime MM-GBSA | gmx_MMPBSA + AMBER MMPBSA.py |
| Desmond | GROMACS (AMBER ff14SB / OpenFF) or OpenMM |
| Simulation Interactions Diagram | MDAnalysis + ProLIF + matplotlib |

**The directive I produced earlier (`research_directive.md`) covers this open-source path in detail and remains valid as a fallback.**

---

## 17. Reading List

**Primary references (cite extensively):**
- Akhilesh, Baidya, Uniyal, Das, Kumar, Tiwari (2022). *Structure-based virtual screening and molecular dynamics simulation for the identification of sphingosine kinase-2 inhibitors as potential analgesics.* J Biomol Struct Dyn, 40(23):12472–12490.
- Akhilesh, Menon, Agrawal, Chouhan, Gadepalli, Das, Kumar, Singh, Tiwari (2024). *Virtual screening and molecular dynamics investigations using natural compounds against autotaxin for the treatment of chronic pain.* J Biomol Struct Dyn, 43(11):5372–5392.

**SRC kinase + pain biology:**
- Roskoski (2015). *Src protein-tyrosine kinase structure, mechanism, and small molecule inhibitors.* Pharmacological Research, 94:9–25.
- Salter & Pitcher (2012). *Dysregulated Src upregulation of NMDA receptor activity: a common link in chronic pain.* FEBS Journal, 279(1):2–11.
- Liu et al. (2008). *Treatment of inflammatory and neuropathic pain by uncoupling Src from the NMDA receptor complex.* Nature Medicine, 14(12):1325–1332.

**Methods:**
- Friesner et al. (2006). *Extra precision Glide: docking and scoring incorporating a model of hydrophobic enclosure for protein-ligand complexes.* J Med Chem, 49(21):6177–6196.
- Lu et al. (2021). *OPLS4: Improving Force Field Accuracy on Challenging Regimes of Chemical Space.* J Chem Theory Comput, 17(7):4291–4300.
- Bowers et al. (2006). *Scalable Algorithms for Molecular Dynamics Simulations on Commodity Clusters.* SC '06.
- Genheden & Ryde (2015). *The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities.* Expert Opin Drug Discov, 10(5):449–461.

---

## 18. Final Word

This is a **template-faithful, novel-target** thesis. Your job is *not* to invent a new pipeline — your job is to apply this lab's validated computational pipeline to a high-impact, mechanistically distinct target (SRC kinase) and to produce a defensible portfolio of natural compound leads with demonstrated pocket binding stability and credible analgesic mechanism.

Stay rigorous on the things that actually matter to your defense:
1. Re-docking validation (≤ 2 Å RMSD).
2. Triplicate MD with reported SD.
3. CNS/BBB relevance for every lead (this is where your work goes *beyond* SphK2/ATX — those targets are peripheral; SRC's analgesic relevance requires CNS access).
4. Mechanistic discussion that bridges kinase pocket → NMDAR-Y1472 → central sensitization.
5. Honest limitations: docking is correlative, MM-GBSA neglects entropy unless explicitly computed, and *in silico* hits require *in vitro* validation.

You have a clear precedent, a defined target, and a known methodology. Execute.
