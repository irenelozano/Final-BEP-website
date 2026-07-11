#!/bin/bash
# =========================================================================
# run_time_refinement.sh - CFL-Bounded Temporal Step Independence Study
# =========================================================================

GRID_LOCKED=40
V_WORST=700
D_WORST=10e-6

# Physical constants for CFL computation
K_B=1.380649e-23
M_H2=3.346e-27
T_GAS=293.15
L_DOMAIN=$(awk -v d="$D_WORST" 'BEGIN {print 50 * d}')

# CFL reference timestep for N=40 grid
DX_REF=$(awk -v L="$L_DOMAIN" -v n="$GRID_LOCKED" 'BEGIN {printf "%.8e", L / n}')
V_TH=$(awk -v k="$K_B" -v t="$T_GAS" -v m="$M_H2" 'BEGIN {printf "%.8e", sqrt(3*k*t/m)}')
V_MAX=$(awk -v vt="$V_TH" -v v="$V_WORST" 'BEGIN {printf "%.8e", vt + v}')
DT_CFL=$(awk -v dx="$DX_REF" -v vm="$V_MAX" 'BEGIN {printf "%.8e", 0.9 * dx / vm}')

# ─────────────────────────────────────────────────────────────────────────────
# LOCKED FNUM: Safe Production Sweet Spot (PPC ≈ 48, Particles ≈ 3 Million)
# ─────────────────────────────────────────────────────────────────────────────
FNUM_LOCKED=1000000

# Timestep values: Reversed so the fastest (fewest steps) run first
time_steps=(
    $(awk -v dt="$DT_CFL" 'BEGIN {printf "%.4e", dt*4}')
    $(awk -v dt="$DT_CFL" 'BEGIN {printf "%.4e", dt*2}')
    $(awk -v dt="$DT_CFL" 'BEGIN {printf "%.4e", dt}')
    $(awk -v dt="$DT_CFL" 'BEGIN {printf "%.4e", dt/2}')
    $(awk -v dt="$DT_CFL" 'BEGIN {printf "%.4e", dt/4}')
)

# Reference: 5000 steps at DT_CFL → constant physical time base
PHYS_TIME_BASE=$(awk -v dt="$DT_CFL" 'BEGIN {printf "%.8e", 5000 * dt}')

echo "Grid: ${GRID_LOCKED}^3, Domain: ${L_DOMAIN} m"
echo "CFL dt: ${DT_CFL} s, Physical time base: ${PHYS_TIME_BASE} s"
echo ""

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
mkdir -p ../runtime_output

CSV_OUT="../runtime_output/time_convergence_${TIMESTAMP}.csv"
echo "Timestep_s,Total_Steps,Drag_Fx" > "$CSV_OUT"

echo "========================================================="
echo " CFL-BOUNDED TIMESTEP CONVERGENCE STUDY                  "
echo " Grid Mesh: ${GRID_LOCKED}^3 | Fnum: ${FNUM_LOCKED}      "
echo " Δt_CFL = ${DT_CFL} s                                    "
echo " T_total = ${PHYS_TIME_BASE} s (constant for all dt)     "
echo "========================================================="

for DT in "${time_steps[@]}"
do
    TOTAL_STEPS=$(awk -v tt="$PHYS_TIME_BASE" -v dt="$DT" 'BEGIN {n=int(tt/dt); if(n%2!=0) n++; n=(n<1000?1000:n); printf "%.0f", n}')
    TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
    STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))

    echo "Timestep: ${DT} s | Steps: ${TOTAL_STEPS} (${TRANSIENT_STEPS}+${STEADY_STEPS}) | CFL ratio: $(awk -v dt="$DT" -v cfl="$DT_CFL" 'BEGIN {printf "%.2f", dt/cfl}')"

    LOG_FILE="../runtime_output/time_study_DT${DT}_${TIMESTAMP}.log"
    DATA_FILE="../runtime_output/force_DT${DT}_${TIMESTAMP}.txt"

    mpirun -np 4 /home/ilozano/sparta/src/spa_mpi \
        -var vel_val "$V_WORST" \
        -var diam_scale "$D_WORST" \
        -var fnum_val "$FNUM_LOCKED" \
        -var grid_res "$GRID_LOCKED" \
        -var t_step "$DT" \
        -var TRANSIENT_STEPS "$TRANSIENT_STEPS" \
        -var STEADY_STEPS "$STEADY_STEPS" \
        -var out_name "$DATA_FILE" \
        < in.drag > "$LOG_FILE"

    DRAG_FX=$(awk '
        BEGIN {sum = 0}
        /^ITEM:/ {next}
        /^#/ {next}
        NF==4 {sum += $2}
        END {if (sum < 0) sum = -sum; printf "%e", sum}
    ' "$DATA_FILE")

    echo "${DT},${TOTAL_STEPS},${DRAG_FX}" >> "$CSV_OUT"
    echo "  Drag: ${DRAG_FX} N"
    echo "---------------------------------------------------------"
done

echo "TIMESTEP CONVERGENCE STUDY COMPLETE. Results saved to: $CSV_OUT"