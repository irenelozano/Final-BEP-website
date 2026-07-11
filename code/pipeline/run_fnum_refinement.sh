#!/bin/bash
# =========================================================================
# run_fnum_refinement.sh - Statistical Particle Convergence Study
# =========================================================================
# Varies fnum (particle weighting) while locking all other parameters.
# Uses worst-case (d=10µm, v=700 m/s) with converged N=40 grid.
#
# Key physics for N=40 at 10µm:
#   n_real_total = nrho * V_domain = 2.47e22 * (5e-4)^3 = 3.08e12 molecules
#   n_real_per_cell = n_real_total / N^3 = 3.08e12 / 64000 ≈ 48.2 million
#
# SERVER-SAFE FNUM VALUES:
#   fnum=10000000 (1e7) -> PPC ≈ 4.8  (~300K particles)  [Runs in seconds]
#   fnum=5000000  (5e6) -> PPC ≈ 9.6  (~600K particles)
#   fnum=2000000  (2e6) -> PPC ≈ 24.1 (~1.5M particles)
#   fnum=1000000  (1e6) -> PPC ≈ 48.2 (~3.0M particles)  [Target sweet spot]
#   fnum=500000   (5e5) -> PPC ≈ 96.4 (~6.1M particles)  [Heavy reference]
# =========================================================================

set -e

# ─────────────────────────────────────────────────────────────────────────────
# FIXED PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
GRID_LOCKED=40
V_WORST=700
D_WORST=10e-6

# Physical constants
K_B=1.380649e-23
M_H2=3.346e-27
T_GAS=293.15
L_DOMAIN=$(awk -v d="$D_WORST" 'BEGIN {print 50 * d}')
N_RHO=2.47e22
TOTAL_CELLS=$(( GRID_LOCKED * GRID_LOCKED * GRID_LOCKED ))

# CFL-based timestep for N=40 grid
DX=$(awk -v L="$L_DOMAIN" -v n="$GRID_LOCKED" 'BEGIN {printf "%.8e", L / n}')
V_TH=$(awk -v k="$K_B" -v t="$T_GAS" -v m="$M_H2" 'BEGIN {printf "%.8e", sqrt(3*k*t/m)}')
V_MAX=$(awk -v vt="$V_TH" -v v="$V_WORST" 'BEGIN {printf "%.8e", vt + v}')
DT_FIXED=$(awk -v dx="$DX" -v vm="$V_MAX" 'BEGIN {printf "%.8e", 0.9 * dx / vm}')

# Domain volume and total real molecules
V_DOMAIN=$(awk -v L="$L_DOMAIN" 'BEGIN {printf "%.8e", L^3}')
N_REAL_TOTAL=$(awk -v nr="$N_RHO" -v vd="$V_DOMAIN" 'BEGIN {printf "%.8e", nr * vd}')

TOTAL_STEPS=5000
TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))

# ─────────────────────────────────────────────────────────────────────────────
# FNUM VALUES: Reversed so the fastest simulations run first
# ─────────────────────────────────────────────────────────────────────────────
fnum_steps=(10000000 5000000 2000000 1000000 500000)

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
mkdir -p ../runtime_output
CSV_OUT="../runtime_output/fnum_convergence_${TIMESTAMP}.csv"
echo "Fnum,Total_Simulated_Particles,Drag_Fx" > "$CSV_OUT"

echo "========================================================="
echo " FNUM (PARTICLE) CONVERGENCE STUDY                       "
echo " Grid: ${GRID_LOCKED}^3 = ${TOTAL_CELLS} cells          "
echo " dt: ${DT_FIXED} s | Steps: ${TOTAL_STEPS}               "
echo " n_real_total = ${N_REAL_TOTAL} molecules               "
echo "========================================================="
echo ""

for F in "${fnum_steps[@]}"
do
    N_SIM=$(awk -v nrt="$N_REAL_TOTAL" -v f="$F" 'BEGIN {printf "%.0f", nrt / f}')
    PPC_EXPECTED=$(awk -v ns="$N_SIM" -v nc="$TOTAL_CELLS" 'BEGIN {printf "%.1f", ns / nc}')

    echo "Fnum: ${F} | expected particles: ${N_SIM} | PPC: ${PPC_EXPECTED}"

    LOG_FILE="../runtime_output/fnum_study_F${F}_${TIMESTAMP}.log"
    DATA_FILE="../runtime_output/force_F${F}_${TIMESTAMP}.txt"

    mpirun -np 4 /home/ilozano/sparta/src/spa_mpi \
        -var vel_val "$V_WORST" \
        -var diam_scale "$D_WORST" \
        -var fnum_val "$F" \
        -var grid_res "$GRID_LOCKED" \
        -var t_step "$DT_FIXED" \
        -var TRANSIENT_STEPS "$TRANSIENT_STEPS" \
        -var STEADY_STEPS "$STEADY_STEPS" \
        -var out_name "$DATA_FILE" \
        < in.drag > "$LOG_FILE"

    TOTAL_NP=$(grep -A 1 "Step CPU Np" "$LOG_FILE" | tail -n 1 | awk '{print $3}')

    DRAG_FX=$(awk '
        BEGIN {sum=0}
        /^ITEM:/ {next}
        /^#/ {next}
        NF==4 {sum += $2}
        END {if (sum < 0) sum = -sum; printf "%e", sum}
    ' "$DATA_FILE")

    echo "${F},${TOTAL_NP},${DRAG_FX}" >> "$CSV_OUT"
    echo "  Actual particles: ${TOTAL_NP} | Drag: ${DRAG_FX} N"
    echo "---------------------------------------------------------"
done

echo ""
echo "FNUM CONVERGENCE STUDY COMPLETE. Results saved to: $CSV_OUT"