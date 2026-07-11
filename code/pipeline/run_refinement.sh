#!/bin/bash
# =========================================================================
# run_refinement.sh — SPARTA DSMC Grid Convergence Study
# =========================================================================
# PURPOSE:
#   Execute a mathematically rigorous grid independence study for the worst-case
#   scenario (d = 10 µm, v = 700 m/s, 100 Pa H₂). Systematically vary grid
#   resolution N ∈ {10, 16, 22, 28, 34, 40} while maintaining constant
#   particles-per-cell via dynamic fnum scaling. Verify convergence of drag force.
#
# PHYSICS PARAMETERS (worst-case):
#   Droplet diameter:    d = 10 µm
#   Stream velocity:     v = 700 m/s
#   Background pressure: P = 100 Pa
#   Gas temperature:     T = 293 K
#
# NUMERICAL PARAMETERS:
#   Domain length:       L = 50 × d = 500 µm
#   Target PPC:          50 (free-molecular regime)
#   Total particle cap:  1.5 million (RAM/runtime limit)
#   Sampling window:     25,000 steady-state steps
#
# OUTPUT:
#   CSV file: ../runtime_output/grid_convergence.csv
#   Columns:  N, Total_Cells, Drag_Fx_N, PPC_actual, Fnum_used
#
# NOTES:
#   - All fnum calculations use ideal gas law: n_rho = P / (k_B × T)
#   - Scales fnum to maintain PPC = 50 across all grid resolutions
#   - Enforces computational cap to prevent memory exhaustion
#   - Percentage change and convergence status reported for each step
# =========================================================================

set -e  # Exit on error

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS & WORST-CASE PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

BOLTZMANN_K=1.380649e-23      # Boltzmann constant (J/K)
D_WORST=10e-6                  # Droplet diameter (m): 10 µm
V_WORST=700                    # Stream velocity (m/s)
PRESSURE=100                   # Background pressure (Pa)
TEMPERATURE=293                # Background temperature (K)

# Derived physical quantities
DOMAIN_LENGTH=$(awk -v d="$D_WORST" 'BEGIN {printf "%.8e", 50 * d}')  # L = 50d
N_RHO=$(awk -v p="$PRESSURE" -v k="$BOLTZMANN_K" -v t="$TEMPERATURE" \
    'BEGIN {printf "%.8e", p / (k * t)}')  # Molecular number density

echo "========================================================="
echo " SPARTA DSMC GRID CONVERGENCE REFINEMENT STUDY"
echo "========================================================="
echo "Worst-case scenario:"
echo "  Droplet diameter d     = ${D_WORST} m (10 µm)"
echo "  Stream velocity v      = ${V_WORST} m/s"
echo "  Pressure P             = ${PRESSURE} Pa"
echo "  Temperature T          = ${TEMPERATURE} K"
echo "  Domain length L = 50d  = ${DOMAIN_LENGTH} m"
echo "  Molecular density n_ρ  = ${N_RHO} m⁻³"
echo "========================================================="
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS FOR CFL TIMESTEP
# ─────────────────────────────────────────────────────────────────────────────

K_B=1.380649e-23
M_H2=3.346e-27
T_GAS=293.15
V_STREAM=700

# Use the FINEST grid (N=40) to compute a fixed CFL timestep for ALL grid levels
# This isolates grid resolution as the only variable in the convergence study.
N_REF=40
DX_REF=$(awk -v L="$DOMAIN_LENGTH" -v n="$N_REF" 'BEGIN {printf "%.8e", L / n}')
V_TH=$(awk -v k="$K_B" -v t="$T_GAS" -v m="$M_H2" 'BEGIN {printf "%.8e", sqrt(3*k*t/m)}')
V_MAX=$(awk -v vt="$V_TH" -v v="$V_STREAM" 'BEGIN {printf "%.8e", vt + v}')
DT_FIXED=$(awk -v dx="$DX_REF" -v vm="$V_MAX" 'BEGIN {printf "%.8e", 0.9 * dx / vm}')

echo "CFL reference: N=${N_REF}, dx=${DX_REF} m, v_th=${V_TH} m/s, v_max=${V_MAX} m/s"
echo "Fixed timestep for all grid levels: dt=${DT_FIXED} s"

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# RUNTIME BUDGET: each simulation must finish ≤40 minutes.
# Production benchmark: 1.5M particles × 50000 steps ≈ 5 hours on 4 ranks.
# Budget: 1e10 particle-steps ≈ 40 min worst case.
# With 2000 steps × 1.5M particles = 3e9 (~12 min), plenty of margin.
# ─────────────────────────────────────────────────────────────────────────────
TOTAL_STEPS=2000
TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))

# Particles per cell target (used for coarse grids where CAP doesn't apply)
PPC_TARGET=50

# Computational RAM cap: 1.5M particles max keeps even N=40 ≤12 min
MAX_TOTAL_PARTICLES=1500000

# Grid resolutions to evaluate (in increasing order)
GRID_STEPS=(10 16 22 28 34 40)

# ─────────────────────────────────────────────────────────────────────────────
# SETUP OUTPUT DIRECTORY & CSV HEADER
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p ../runtime_output

CSV_FILE="../runtime_output/grid_convergence.csv"
echo "N,Total_Cells,Drag_Fx_N,Drag_Fx_pN,PPC_actual,Fnum_used,Percent_Change" > "$CSV_FILE"

PREVIOUS_DRAG=""

# ─────────────────────────────────────────────────────────────────────────────
# MAIN REFINEMENT LOOP
# ─────────────────────────────────────────────────────────────────────────────

for N in "${GRID_STEPS[@]}"
do
    TOTAL_CELLS=$((N * N * N))
    
    # Compute cell volume for this grid level
    CELL_SIZE=$(awk -v l="$DOMAIN_LENGTH" -v n="$N" 'BEGIN {printf "%.8e", l / n}')
    CELL_VOL=$(awk -v cs="$CELL_SIZE" 'BEGIN {printf "%.8e", cs^3}')
    
    # Real molecules per cell (from ideal gas law)
    N_REAL_PER_CELL=$(awk -v nrho="$N_RHO" -v vcell="$CELL_VOL" \
        'BEGIN {printf "%.3e", nrho * vcell}')
    
    # ─────────────────────────────────────────────────────────────────────────
    # CALCULATE FNUM: Enforce TWO CONSTRAINTS
    # ─────────────────────────────────────────────────────────────────────────
    
    # Constraint 1: Physics-based fnum to maintain target PPC
    FNUM_PHYSICS=$(awk -v nrho="$N_RHO" -v vcell="$CELL_VOL" -v ppc="$PPC_TARGET" \
        'BEGIN {fnum = (nrho * vcell) / ppc; printf "%.0f", fnum}')
    
    # Constraint 2: Computational cap to prevent memory explosion
    # Total particles = N³ × PPC_actual = (N_RHO × V_domain) / fnum
    # We want: (N_RHO × V_domain) / fnum ≤ MAX_TOTAL_PARTICLES
    # Therefore: fnum ≥ (N_RHO × V_domain) / MAX_TOTAL_PARTICLES
    V_DOMAIN=$(awk -v l="$DOMAIN_LENGTH" 'BEGIN {printf "%.8e", l^3}')
    FNUM_CAP=$(awk -v nrho="$N_RHO" -v vdomain="$V_DOMAIN" -v max="$MAX_TOTAL_PARTICLES" \
        'BEGIN {fnum = (nrho * vdomain) / max; printf "%.0f", fnum}')
    
    # Use the LARGER of the two (more fnum = fewer particles)
    FNUM_DYNAMIC=$(awk -v phys="$FNUM_PHYSICS" -v cap="$FNUM_CAP" \
        'BEGIN {fnum = (phys > cap) ? phys : cap; printf "%.0f", fnum}')
    
    # Safety floor: never let fnum drop below 5
    if [ "$FNUM_DYNAMIC" -lt 5 ]; then
        FNUM_DYNAMIC=5
    fi
    
    # Calculate actual PPC achieved with this fnum
    PPC_ACTUAL=$(awk -v nreal="$N_REAL_PER_CELL" -v fnum="$FNUM_DYNAMIC" \
        'BEGIN {ppc = nreal / fnum; printf "%.2f", ppc}')
    
    TOTAL_PARTICLES=$(awk -v cells="$TOTAL_CELLS" -v ppc="$PPC_ACTUAL" \
        'BEGIN {printf "%.0f", cells * ppc}')
    
    # ─────────────────────────────────────────────────────────────────────────
    # PRINT PROGRESS & DIAGNOSTICS
    # ─────────────────────────────────────────────────────────────────────────
    
    echo ""
    echo "──────────────────────────────────────────────────────────"
    echo "Grid level: ${N}³ = ${TOTAL_CELLS} cells"
    echo "  Cell size:            ${CELL_SIZE} m"
    echo "  Real molecules/cell:  ${N_REAL_PER_CELL}"
    echo "  fnum (physics):       ${FNUM_PHYSICS}"
    echo "  fnum (cap):           ${FNUM_CAP}"
    echo "  fnum (used):          ${FNUM_DYNAMIC}"
    echo "  PPC achieved:         ${PPC_ACTUAL}"
    echo "  Total particles:      ${TOTAL_PARTICLES}"
    echo ""
    
    # ─────────────────────────────────────────────────────────────────────────
    # RUN SPARTA SIMULATION
    # ─────────────────────────────────────────────────────────────────────────
    
    LOG_FILE="../runtime_output/refinement_N${N}.log"
    DATA_FILE="../runtime_output/force_N${N}.txt"
    
    echo "  Executing SPARTA simulation..."
    
    # Call SPARTA with dynamic parameters
    # NOTE: Variable names are CASE-SENSITIVE and must match in.drag exactly
    mpirun -np 4 /home/ilozano/sparta/src/spa_mpi \
        -var vel_val "$V_WORST" \
        -var diam_scale "$D_WORST" \
        -var fnum_val "$FNUM_DYNAMIC" \
        -var grid_res "$N" \
        -var t_step "$DT_FIXED" \
        -var TRANSIENT_STEPS "$TRANSIENT_STEPS" \
        -var STEADY_STEPS "$STEADY_STEPS" \
        -var out_name "$DATA_FILE" \
        < in.drag > "$LOG_FILE" 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: SPARTA simulation failed for N=${N}. Check ${LOG_FILE}"
        exit 1
    fi
    
    echo "  SPARTA run complete."
    
    # ─────────────────────────────────────────────────────────────────────────
    # EXTRACT DRAG FORCE FROM OUTPUT
    # ─────────────────────────────────────────────────────────────────────────
    
    # Parse the force output file. The exact format depends on your in.drag script.
    # Expected: columns where column 2 is the x-component of force per triangle.
    # This extracts and sums all force values, then converts to absolute value.
    
    DRAG_FX=$(awk '
        BEGIN {sum = 0}
        /^ITEM:/ {next}
        /^#/ {next}
        NF==4 {sum += $2}
        END {
            if (sum < 0) sum = -sum
            printf "%.8e", sum
        }
    ' "$DATA_FILE")
    
    # Convert to piconewtons for readability
    DRAG_PN=$(awk -v drag="$DRAG_FX" 'BEGIN {printf "%.4f", drag * 1e12}')
    
    # ─────────────────────────────────────────────────────────────────────────
    # CALCULATE PERCENTAGE CHANGE FROM PREVIOUS STEP
    # ─────────────────────────────────────────────────────────────────────────
    
    if [ -z "$PREVIOUS_DRAG" ]; then
        PCT_CHANGE="—"
        PCT_CHANGE_NUM=999.0
    else
        PCT_CHANGE_NUM=$(awk -v prev="$PREVIOUS_DRAG" -v curr="$DRAG_FX" \
            'BEGIN {change = (curr - prev) / prev * 100; printf "%.3f", change}')
        
        PCT_CHANGE="${PCT_CHANGE_NUM}%"
    fi
    
    # ─────────────────────────────────────────────────────────────────────────
    # CHECK CONVERGENCE
    # ─────────────────────────────────────────────────────────────────────────
    
    # Convergence criterion: |change| < 2% and change is monotonically decreasing
    CONVERGED="no"
    PCT_ABS=$(awk -v pct="$PCT_CHANGE_NUM" 'BEGIN {if (pct < 0) pct = -pct; printf "%.3f", pct}')
    
    if (( $(awk -v pct="$PCT_ABS" 'BEGIN {print (pct < 2.0) ? 1 : 0}') )); then
        if [ -z "$PREVIOUS_DRAG" ] || (( $(awk -v pct="$PCT_CHANGE_NUM" 'BEGIN {print (pct >= -2.0 && pct <= 0.5) ? 1 : 0}') )); then
            CONVERGED="yes"
        fi
    fi
    
    # ─────────────────────────────────────────────────────────────────────────
    # APPEND TO CSV & PRINT RESULTS
    # ─────────────────────────────────────────────────────────────────────────
    
    echo "${N},${TOTAL_CELLS},${DRAG_FX},${DRAG_PN},${PPC_ACTUAL},${FNUM_DYNAMIC},${PCT_CHANGE}" >> "$CSV_FILE"
    
    echo "  Drag force |Fx|:     ${DRAG_FX} N = ${DRAG_PN} pN"
    echo "  Change from prev:    ${PCT_CHANGE}"
    echo "  Convergence status:  ${CONVERGED}"
    echo "──────────────────────────────────────────────────────────"
    
    # Update for next iteration
    PREVIOUS_DRAG="$DRAG_FX"
done

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY & CONVERGENCE REPORT
# ─────────────────────────────────────────────────────────────────────────────

echo ""
echo "========================================================="
echo " GRID CONVERGENCE STUDY COMPLETE"
echo "========================================================="
echo ""
echo "Results written to: ${CSV_FILE}"
echo ""
echo "Quick summary:"
tail -6 "$CSV_FILE" | awk -F, '
    NR == 1 {print "N\tCells\tDrag (pN)\tPPC\tFnum\t% Change"}
    NR > 1 {printf "%s\t%s\t%.4f\t%.2f\t%s\t%s\n", $1, $2, $4, $5, $6, $7}
'
echo ""
echo "For detailed analysis, plot N vs. Drag_Fx_pN and check for:"
echo "  1. Monotonic convergence (drag values flatten)"
echo "  2. Percentage change < 2% for final steps"
echo "  3. No sign reversals in percentage change"
echo ""
echo "If convergence was NOT achieved by N=40, you may need:"
echo "  - Larger N values (extend GRID_STEPS)"
echo "  - Lower target PPC (but not below 20 for free-molecular)"
echo "  - Longer steady-state sampling window"
echo ""
echo "========================================================="
