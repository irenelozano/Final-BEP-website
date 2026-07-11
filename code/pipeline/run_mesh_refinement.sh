#!/bin/bash
# =========================================================================
# run_mesh_refinement.sh — Sphere Surface Mesh Convergence Study
# =========================================================================

set -e

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS & WORST-CASE PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
K_B=1.380649e-23
M_H2=3.346e-27
T_GAS=293.15

D_WORST=10e-6
V_WORST=700
PRESSURE=100
TEMPERATURE=293

# ─────────────────────────────────────────────────────────────────────────────
# LOCKED NUMERICAL PARAMETERS (Derived from verified convergence data)
# ─────────────────────────────────────────────────────────────────────────────
GRID_LOCKED=40
TOTAL_STEPS=5000
TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))

# Domain: L = 50 × d
L_DOMAIN=$(awk -v d="$D_WORST" 'BEGIN {printf "%.8e", 50 * d}')

# CFL timestep for N=40 grid (Using optimal CFL = 0.5 verified from time study)
DX=$(awk -v L="$L_DOMAIN" -v n="$GRID_LOCKED" 'BEGIN {printf "%.8e", L / n}')
V_TH=$(awk -v k="$K_B" -v t="$T_GAS" -v m="$M_H2" 'BEGIN {printf "%.8e", sqrt(3*k*t/m)}')
V_MAX=$(awk -v vt="$V_TH" -v v="$V_WORST" 'BEGIN {printf "%.8e", vt + v}')
DT_FIXED=$(awk -v dx="$DX" -v vm="$V_MAX" 'BEGIN {printf "%.8e", 0.5 * dx / vm}')

# SERVER-SAFE FNUM:
# n_real_total = 3.08e12 molecules
# FNUM = 5,000,000 -> ~616K active particles -> PPC ≈ 9.6
# This perfectly hits your ~12 minute runtime budget!
FNUM_LOCKED=5000000

# ─────────────────────────────────────────────────────────────────────────────
# MESH FILES (relative to pipeline/)
# ─────────────────────────────────────────────────────────────────────────────
declare -A MESH_FILES
MESH_FILES[192]="../global_assets/data.sphere"
MESH_FILES[1200]="../global_assets/data_1200.sphere"
MESH_FILES[4800]="../global_assets/data_2400.sphere3d"

MESH_ORDER=(192 1200 4800)

# ─────────────────────────────────────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────────────────────────────────────
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
mkdir -p ../runtime_output

CSV_FILE="../runtime_output/mesh_convergence_${TIMESTAMP}.csv"
echo "Ntriangles,Mesh_File,Drag_Fx_N,Drag_Fx_pN,PPC_actual" > "$CSV_FILE"

echo "========================================================="
echo " SPHERE SURFACE MESH CONVERGENCE STUDY                   "
echo "========================================================="
echo "Worst-case: d=${D_WORST} m, v=${V_WORST} m/s"
echo "Grid: ${GRID_LOCKED}^3 | dt: ${DT_FIXED} s | fnum: ${FNUM_LOCKED}"
echo "Steps: ${TOTAL_STEPS} (${TRANSIENT_STEPS}+${STEADY_STEPS})"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# MAIN LOOP
# ─────────────────────────────────────────────────────────────────────────────

PREVIOUS_DRAG=""

for N_TRI in "${MESH_ORDER[@]}"; do
    MESH_PATH="${MESH_FILES[$N_TRI]}"

    echo "--- Mesh: ${N_TRI} triangles (${MESH_PATH}) ---"

    LOG_FILE="../runtime_output/mesh_${N_TRI}tri_${TIMESTAMP}.log"
    DATA_FILE="../runtime_output/force_mesh_${N_TRI}tri_${TIMESTAMP}.txt"

         \
        -var vel_val "$V_WORST" \
        -var diam_scale "$D_WORST" \
        -var fnum_val "$FNUM_LOCKED" \
        -var grid_res "$GRID_LOCKED" \
        -var t_step "$DT_FIXED" \
        -var TRANSIENT_STEPS "$TRANSIENT_STEPS" \
        -var STEADY_STEPS "$STEADY_STEPS" \
        -var mesh_file "$MESH_PATH" \
        -var out_name "$DATA_FILE" \
        < in.drag > "$LOG_FILE" 2>&1

    if [ $? -ne 0 ]; then
        echo "ERROR: SPARTA failed for mesh ${N_TRI}. Check ${LOG_FILE}"
        exit 1
    fi

    # Extract total drag from cumulative dump
    DRAG_FX=$(awk '
        BEGIN {sum = 0}
        /^ITEM:/ {next}
        /^#/ {next}
        NF==4 {sum += $2}
        END {if (sum < 0) sum = -sum; printf "%.8e", sum}
    ' "$DATA_FILE")

    DRAG_PN=$(awk -v drag="$DRAG_FX" 'BEGIN {printf "%.4f", drag * 1e12}')

    # Compute actual PPC from log
    TOTAL_CELLS=$(( GRID_LOCKED * GRID_LOCKED * GRID_LOCKED ))
    NP_FINAL=$(grep -A 1 "Step CPU Np" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    PPC=$(awk -v np="$NP_FINAL" -v nc="$TOTAL_CELLS" 'BEGIN {printf "%.2f", np/nc}')

    # Percent change
    if [ -z "$PREVIOUS_DRAG" ]; then
        PCT_CHANGE="—"
    else
        PCT_CHANGE=$(awk -v prev="$PREVIOUS_DRAG" -v curr="$DRAG_FX" \
            'BEGIN {change = (curr - prev) / prev * 100; printf "%.3f", change}')
        PCT_CHANGE="${PCT_CHANGE}%"
    fi

    echo "${N_TRI},${MESH_PATH},${DRAG_FX},${DRAG_PN},${PPC}" >> "$CSV_FILE"

    echo "  Drag: ${DRAG_FX} N = ${DRAG_PN} pN"
    echo "  Particles: ${NP_FINAL} | PPC: ${PPC}"
    echo "  Change from previous mesh: ${PCT_CHANGE}"
    echo ""

    PREVIOUS_DRAG="$DRAG_FX"
done

echo "========================================================="
echo " MESH CONVERGENCE STUDY COMPLETE"
echo " Results written to: ${CSV_FILE}"
echo "========================================================="