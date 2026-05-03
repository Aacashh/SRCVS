#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SCHRODINGER:-}" ]]; then
    echo "SCHRODINGER environment variable not set" >&2
    exit 2
fi

receptor=""
ligand=""
out_dir=""
seed=""
config=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --receptor) receptor="$2"; shift 2 ;;
        --ligand) ligand="$2"; shift 2 ;;
        --out-dir) out_dir="$2"; shift 2 ;;
        --seed) seed="$2"; shift 2 ;;
        --config) config="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 4 ;;
    esac
done

mkdir -p "$out_dir"
cd "$out_dir"

prod_ns=$(awk '/^md:/{f=1;next} f&&/^  production_ns:/{print $2;exit}' "$config")
temp_k=$(awk '/^md:/{f=1;next} f&&/^  temperature_K:/{print $2;exit}' "$config")
press=$(awk '/^md:/{f=1;next} f&&/^  pressure_bar:/{print $2;exit}' "$config")
buffer=$(awk '/^md:/{f=1;next} f&&/^  box_buffer_nm:/{print $2;exit}' "$config")
salt=$(awk '/^md:/{f=1;next} f&&/^  ion_concentration_M:/{print $2;exit}' "$config")
trj_int=$(awk '/^md:/{f=1;next} f&&/^  trajectory_interval_ps:/{print $2;exit}' "$config")
buffer_a=$(python3 -c "print(${buffer:-1.0}*10)")
prod_ps=$(python3 -c "print(${prod_ns:-100.0}*1000)")

"$SCHRODINGER/run" python3 - "$receptor" "$ligand" "$out_dir/complex.mae" <<'PYEOF'
import sys
from schrodinger import structure
recep = next(structure.StructureReader(sys.argv[1]))
lig = next(structure.StructureReader(sys.argv[2]))
out = structure.StructureWriter(sys.argv[3])
recep.extend(lig)
out.append(recep)
out.close()
PYEOF

cat > system_build.msj <<EOF
task { task = "desmond:auto" }
build_geometry {
  box = { shape = orthorhombic size = [$buffer_a $buffer_a $buffer_a] size_type = buffer }
  solvent = TIP3P
  add_counterion = { ion = Na species = Cl }
  salt = { concentration = $salt negative_ion = Cl positive_ion = Na }
  rezero_system = false
}
assign_forcefield { forcefield = OPLS_2005 }
EOF

"$SCHRODINGER/utilities/multisim" \
    -JOBNAME md_build \
    -m system_build.msj \
    -i complex.mae \
    -o system-out.cms \
    -HOST localhost \
    -maxjob 1 \
    -WAIT

cat > prod.msj <<EOF
task { task = "desmond:auto" }
simulate {
  title = "Brownian Minimization"
  time = 100
  timestep = [0.001 0.001 0.003]
  temperature = 10.0
  ensemble = { class = "NVT" method = "Brownie" brownie = { delta_max = 0.1 } }
  restraints.new = [{ name = posre_harm atoms = solute_heavy_atom force_constants = 50.0 }]
}
simulate {
  title = "NPT equilibration"
  time = 1000.0
  temperature = $temp_k
  pressure = [$press isotropic]
  ensemble = { class = NPT method = MTK thermostat.tau = 1.0 barostat.tau = 2.0 }
  restraints.new = [{ name = posre_harm atoms = solute_heavy_atom force_constants = 5.0 }]
}
simulate {
  title = "Production"
  time = $prod_ps
  temperature = $temp_k
  pressure = [$press isotropic]
  ensemble = { class = NPT method = MTK thermostat.tau = 1.0 barostat.tau = 2.0 }
  randomize_velocity.seed = $seed
  trajectory.center = solute
  trajectory.format = dtr
  trajectory.interval = $trj_int
  trajectory.periodicfix = true
}
EOF

"$SCHRODINGER/utilities/multisim" \
    -JOBNAME md_prod \
    -m prod.msj \
    -i system-out.cms \
    -o md-out.cms \
    -HOST "localhost:1:gpgpu=1" \
    -maxjob 1 \
    -WAIT
