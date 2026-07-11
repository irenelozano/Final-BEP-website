#!/bin/bash
# =========================================================================
# run_steadystate_proof.sh
#
# PURPOSE:
#   Prove steady state for the TWO worst-case scenarios in your matrix:
#     Case A: d=100nm, v=700 m/s  (highest Kn, hardest to converge)
#     Case B: d=10um,  v=700 m/s  (highest Re, strongest transient)
#
#   This script uses in.steadystate_proof, which records the drag force
#   time series from step 0 — including the full transient phase.
#   The output .txt file is what you hand to the MATLAB post-processing
#   script to generate the steady-state proof figure for your report.
#
# HOW IT IS DIFFERENT FROM run_steadystate_test.sh:
#   Old script:  Only 5000 transient + 5000 steady = 5 data points.
#                Records ONLY during steady window. Transient invisible.
#   This script: Reuses your production optimizer step counts.
#                Records from step 0. Transient IS visible in the plot.
#                Runs both your hardest cases automatically.
#
# OUTPUT FILES (per case):
#   *_timeseries.txt      <- time series of Fx from step 0 [PRIMARY PLOT]
#   *_surface_cumul.txt   <- final averaged surface force [DRAG VALUE]
#   *.log                 <- full SPARTA log
#
# USAGE:
#   bash run_steadystate_proof.sh
# =========================================================================

set -euo pipefail

SPARTA_EXE="/home/ilozano/sparta/src/spa_mpi"
MPI_CORES=4
INPUT_FILE="in.steadystate_proof"
OUTDIR="../runtime_output/steadystate_proof"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

mkdir -p "$OUTDIR"

# =========================================================================
# THE DYNAMIC OPTIMIZER (same as run_30case_matrix.sh — kept in sync)
# If you update the optimizer in the matrix script, update it here too.
# =========================================================================
compute_optimal_params() {
    local d=$1 v=$2

    local n_inf=2.47e22
    local d_ref=2.72e-10
    local T_ref=273.15
    local omega=0.67
    local T_inf=293.15
    local m=3.346e-27
    local k=1.380649e-23

    local lambda_inf=$(awk -v n="$n_inf" -v dr="$d_ref" -v w="$omega" -v tr="$T_ref" -v ti="$T_inf" \
        'BEGIN {print 1 / (sqrt(2) * n * atan2(0,-1) * dr^2 * (ti/tr)^(w-0.5))}')

    local Kn=$(awk -v l="$lambda_inf" -v dd="$d" 'BEGIN {print l / dd}')

    local N_grad
    if   (( $(awk -v kn="$Kn" 'BEGIN{print (kn<0.1)}') )); then N_grad=20
    elif (( $(awk -v kn="$Kn" 'BEGIN{print (kn<10)}')  )); then N_grad=10
    else N_grad=5; fi

    local dx_grad=$(awk -v dd="$d" -v ng="$N_grad" 'BEGIN {print dd / ng}')
    local dx_bird=$(awk -v l="$lambda_inf" 'BEGIN {print l / 3}')

    local dx
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn<0.1)}') )); then
        dx=$dx_grad
    else
        dx=$(awk -v b="$dx_bird" -v g="$dx_grad" 'BEGIN {print (b<g ? b : g)}')
    fi

    local L=$(awk -v dd="$d" 'BEGIN {print 50 * dd}')
    local N
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn>10)}') )); then
        local d_worst=100e-9 N_worst=40
        local N_scaled=$(awk -v dw="$d_worst" -v nw="$N_worst" -v dd="$d" \
            'BEGIN {n=int(nw * (dw/dd)^(1/3)); print (n>nw ? nw : n)}')
        N=$(awk -v ns="$N_scaled" 'BEGIN {print (ns<10 ? 10 : ns)}')
    else
        N=$(awk -v Lx="$L" -v dx="$dx" 'BEGIN {n=int(Lx/dx); print (n<10 ? 10 : n)}')
    fi
    local N_min=$(awk -v Lx="$L" -v dd="$d" 'BEGIN {n=int(Lx / dd); print (n<10 ? 10 : n)}')
    N=$(( N > N_min ? N : N_min ))
    N=$(( (N + 1) / 2 * 2 ))
    local dx_actual=$(awk -v Lx="$L" -v nn="$N" 'BEGIN {print Lx / nn}')

    local v_th=$(awk -v k="$k" -v ti="$T_inf" -v m="$m" 'BEGIN {print sqrt(3*k*ti/m)}')
    local v_max=$(awk -v vt="$v_th" -v vv="$v" 'BEGIN {print vt + vv}')
    local dt=$(awk -v dx="$dx_actual" -v vm="$v_max" 'BEGIN {print 0.9 * dx / vm}')

    local V_cell=$(awk -v dx="$dx_actual" 'BEGIN {print dx^3}')
    local n_real=$(awk -v n="$n_inf" -v vc="$V_cell" 'BEGIN {print n * vc}')
    local target_ppc=25
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn>10)}') )); then target_ppc=50; fi
    local fnum=$(awk -v nr="$n_real" -v tp="$target_ppc" \
        'BEGIN {fn=nr/tp; fn=(fn<1?1:int(fn)); print fn}')

    local N_total=$(( N * N * N ))
    local MAX_TOTAL_PARTICLES=1500000
    local fnum_min=$(awk -v nr="$n_real" -v nt="$N_total" -v mx="$MAX_TOTAL_PARTICLES" \
        'BEGIN {f=nr*nt/mx; f=(f<1?1:int(f)+1); print f}')
    if [ "$fnum" -lt "$fnum_min" ]; then fnum=$fnum_min; fi
    local actual_ppc=$(awk -v nr="$n_real" -v fn="$fnum" 'BEGIN {print int(nr/fn)}')

    local MIN_STEPS=50000
    local t_steady=$(awk -v Lx="$L" -v vv="$v" 'BEGIN {print 3 * Lx / vv}')
    local N_steps=$(awk -v ts="$t_steady" -v dt="$dt" -v mn="$MIN_STEPS" \
        'BEGIN {n=int(ts/dt)+1; print (n<mn ? mn : n)}')

    echo "$N $fnum $dt $N_steps $actual_ppc $Kn $lambda_inf"
}

# =========================================================================
# CASE: d=10um, v=700 m/s
#
# This is the case with the strongest transient and highest Reynolds number
# in your parameter matrix, making it the most demanding case to prove
# steady state. It also has enough surface collisions per averaging window
# to produce a clean non-zero force signal in the timeseries.
#
# Kn for 10um is ~12 (low free-molecular), so the transient startup IS
# visible before the running mean converges — exactly what you need for
# the proof figure in your report.
# =========================================================================
D_CASE="10e-6"
V_CASE="700"
LABEL_CASE="d10um_v700"

# Recording interval: how many steps between force samples.
# Smaller = more data points = better resolution of transient.
# Keep <= 500 so you get at least 100 points over a 50k-step run.
RECORD_EVERY=200

# =========================================================================
# SUMMARY LOG FILE
# =========================================================================
SUMMARY="$OUTDIR/steadystate_proof_summary_${TIMESTAMP}.txt"
echo "================================================================" | tee "$SUMMARY"
echo "  STEADY-STATE PROOF RUN — ${TIMESTAMP}"                          | tee -a "$SUMMARY"
echo "  SPARTA: $SPARTA_EXE"                                            | tee -a "$SUMMARY"
echo "  MPI cores: $MPI_CORES"                                          | tee -a "$SUMMARY"
echo "================================================================" | tee -a "$SUMMARY"

# =========================================================================
# HELPER: run one case
# =========================================================================
run_case() {
    local LABEL=$1
    local D_VAL=$2
    local V_VAL=$3

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  CASE: ${LABEL}  |  d = ${D_VAL} m  |  v = ${V_VAL} m/s"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # --- Get optimizer parameters ---
    PARAMS=($(compute_optimal_params "$D_VAL" "$V_VAL"))
    GRID_RES=${PARAMS[0]}
    FNUM_VAL=${PARAMS[1]}
    DT=${PARAMS[2]}
    TOTAL_STEPS_HALF=${PARAMS[3]}   # This is the production steady window length
    ACTUAL_PPC=${PARAMS[4]}
    KN=${PARAMS[5]}

    # For the proof run: use 2x the production steady-state steps.
    #   First half  = transient phase (we watch the force rise and settle)
    #   Second half = steady-state sampling (matches production run exactly)
    # This ensures the transient IS captured in the time series.
    TRANSIENT_STEPS=$TOTAL_STEPS_HALF
    STEADY_STEPS=$TOTAL_STEPS_HALF
    TOTAL_STEPS=$(( TRANSIENT_STEPS + STEADY_STEPS ))

    # --- Output files ---
    OUT_BASE="${OUTDIR}/ss_proof_${LABEL}_${TIMESTAMP}"
    LOG_FILE="${OUT_BASE}.log"

    echo "  Kn          = ${KN}"
    echo "  Grid        = ${GRID_RES}^3"
    echo "  fnum        = ${FNUM_VAL}"
    echo "  dt          = ${DT} s"
    echo "  Transient   = ${TRANSIENT_STEPS} steps"
    echo "  Steady      = ${STEADY_STEPS} steps"
    echo "  Total       = ${TOTAL_STEPS} steps"
    echo "  PPC (target)= ${ACTUAL_PPC}"
    echo "  Record every= ${RECORD_EVERY} steps"
    echo "  Output base = ${OUT_BASE}"
    echo ""

    # --- Log to summary ---
    {
        echo ""
        echo "CASE: ${LABEL}"
        echo "  d=${D_VAL} m | v=${V_VAL} m/s | Kn=${KN}"
        echo "  Grid=${GRID_RES}^3 | fnum=${FNUM_VAL} | dt=${DT} s"
        echo "  Transient=${TRANSIENT_STEPS} | Steady=${STEADY_STEPS} | Total=${TOTAL_STEPS}"
        echo "  Output: ${OUT_BASE}_timeseries.txt"
    } >> "$SUMMARY"

    # --- Run SPARTA ---
    echo "  [$(date +"%H:%M:%S")] Starting SPARTA..."

    mpirun -np "$MPI_CORES" "$SPARTA_EXE" \
        -var diam_scale    "$D_VAL"           \
        -var vel_val       "$V_VAL"           \
        -var grid_res      "$GRID_RES"        \
        -var fnum_val      "$FNUM_VAL"        \
        -var t_step        "$DT"              \
        -var TRANSIENT_STEPS "$TRANSIENT_STEPS" \
        -var STEADY_STEPS  "$STEADY_STEPS"   \
        -var TOTAL_STEPS   "$TOTAL_STEPS"    \
        -var RECORD_EVERY  "$RECORD_EVERY"   \
        -var out_name      "$OUT_BASE"        \
        < "$INPUT_FILE" > "$LOG_FILE" 2>&1

    echo "  [$(date +"%H:%M:%S")] SPARTA finished."

    # --- Quick sanity check on output ---
    TIMESERIES_FILE="${OUT_BASE}_timeseries.txt"
    if [ ! -f "$TIMESERIES_FILE" ]; then
        echo "  ⚠️  ERROR: timeseries file not found! Check log: ${LOG_FILE}" | tee -a "$SUMMARY"
        return 1
    fi

    N_DATAPOINTS=$(grep -c "^[0-9]" "$TIMESERIES_FILE" 2>/dev/null || echo 0)
    echo "  ✓ Timeseries file written: ${N_DATAPOINTS} data points"

    # --- Extract mean drag from last 50% of timeseries (steady window only) ---
    # The timeseries file is a surf dump: one ITEM: TIMESTEP block per RECORD_EVERY steps.
    # We sum Fx (col 2) across all triangles per block, then average blocks in steady window.
    MEAN_FX=$(awk -v half="$TRANSIENT_STEPS" '
        /^ITEM: TIMESTEP/ { in_header=1; current_step=-1; block_fx=0; count=0; next }
        in_header==1 { current_step=$1; in_header=2; next }
        in_header==2 { in_header=3; next }   # NUMBER OF ATOMS line
        in_header==3 { in_header=4; next }   # atom count line
        in_header==4 { in_header=0; next }   # ITEM: ATOMS header line
        NF==4 && current_step > half {
            fx = ($2 < 0 ? -$2 : $2)
            block_fx += fx
            count++
        }
        /^ITEM: TIMESTEP/ || END {
            if (current_step > half && count > 0) {
                sum += block_fx
                n_blocks++
            }
        }
        END {
            if (n_blocks > 0) printf "%.6e", sum/n_blocks
            else print "N/A"
        }
    ' "$TIMESERIES_FILE")

    echo "  Mean |Fx| in steady window = ${MEAN_FX} N"
    echo "  Timeseries: ${TIMESERIES_FILE}"
    echo ""

    echo "  Mean |Fx| in steady window = ${MEAN_FX} N" >> "$SUMMARY"
    echo "  Timeseries: ${TIMESERIES_FILE}" >> "$SUMMARY"

    # --- Write a small config file for MATLAB ---
    # This tells the MATLAB post-processor where the transient ends,
    # what dt is, and what the case label is — no manual editing needed.
    CONFIG_FILE="${OUT_BASE}_config.txt"
    {
        echo "label=${LABEL}"
        echo "diameter_m=${D_VAL}"
        echo "velocity_ms=${V_VAL}"
        echo "Kn=${KN}"
        echo "dt_s=${DT}"
        echo "transient_steps=${TRANSIENT_STEPS}"
        echo "steady_steps=${STEADY_STEPS}"
        echo "total_steps=${TOTAL_STEPS}"
        echo "record_every=${RECORD_EVERY}"
        echo "timeseries_file=${TIMESERIES_FILE}"
        echo "mean_fx_steady_N=${MEAN_FX}"
    } > "$CONFIG_FILE"
    echo "  Config written: ${CONFIG_FILE}"
}

# =========================================================================
# RUN THE CASE
# =========================================================================
run_case "$LABEL_CASE" "$D_CASE" "$V_CASE"

# =========================================================================
# FINAL SUMMARY
# =========================================================================
echo ""
echo "================================================================" | tee -a "$SUMMARY"
echo "  CASE COMPLETE"                                                  | tee -a "$SUMMARY"
echo "  Summary log: ${SUMMARY}"                                        | tee -a "$SUMMARY"
echo ""
echo "  NEXT STEP: open MATLAB and run:"                                | tee -a "$SUMMARY"
echo "    steadystate_proof_plot('${OUTDIR}')"                          | tee -a "$SUMMARY"
echo "================================================================" | tee -a "$SUMMARY"
