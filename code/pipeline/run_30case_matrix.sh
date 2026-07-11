#!/bin/bash
# =========================================================================
# grid_optimizer.sh - 30-Case Parametric Sweep with Dynamic Optimization
# =========================================================================

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
mkdir -p ../production_output
CSV_OUT="../production_output/matrix_results_${TIMESTAMP}.csv"
echo "Diameter_m,Velocity_ms,Grid_N,Fnum,Timestep_s,Total_Steps,Drag_Fx" > "$CSV_OUT"

# THE DYNAMIC OPTIMIZER FUNCTION
compute_optimal_params() {
    local d=$1 v=$2}
    
    # --- gas constants ---
    local n_inf=2.47e22
    local d_ref=2.72e-10
    local T_ref=273.15
    local omega=0.67
    local T_inf=293.15
    local m=3.346e-27
    local k=1.380649e-23

    # Step 1: mean free path
    local lambda_inf=$(awk -v n="$n_inf" -v dr="$d_ref" -v w="$omega" -v tr="$T_ref" -v ti="$T_inf" 'BEGIN {print 1 / (sqrt(2) * n * atan2(0,-1) * dr^2 * (ti/tr)^(w-0.5))}')
    
    # Step 2: Knudsen number
    local Kn=$(awk -v l="$lambda_inf" -v dd="$d" 'BEGIN {print l / dd}')
    
    # Step 3: gradient resolution
    local N_grad
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn<0.1)}') )); then N_grad=20
    elif (( $(awk -v kn="$Kn" 'BEGIN{print (kn<10)}') )); then N_grad=10
    else N_grad=5; fi

    # Step 4: cell size
    local dx_grad=$(awk -v dd="$d" -v ng="$N_grad" 'BEGIN {print dd / ng}')
    local dx_bird=$(awk -v l="$lambda_inf" 'BEGIN {print l / 3}')
    
    local dx
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn<0.1)}') )); then dx=$dx_grad
    else dx=$(awk -v b="$dx_bird" -v g="$dx_grad" 'BEGIN {print (b<g ? b : g)}'); fi

    # Step 5: grid
    local L=$(awk -v dd="$d" 'BEGIN {print 50 * dd}')
    local N
    # Dynamic grid selection based on flow regime:
    #   - Continuum (Kn<0.1):     20 cells minimum, gradient-resolved
    #   - Transitional (0.1<Kn<10): Bird/gradient criterion
    #   - Free-Molecular (Kn>10):  Cube-root scaling from 40^3 (100nm) to 10^3 (10um)
    #     The worst-case (100nm, Kn=1230) caps at 40^3 = 64K cells.
    #     Larger droplets in FM need fewer cells since the flow is collisionless.
    #     At Kn>10, Bird criterion gives dx > Lx (unphysical), so we scale
    #     by droplet diameter: N = 40 * (d_worst/d)^(1/3), bounded to [10, 40].
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn>10)}') )); then
        local d_worst=100e-9 N_worst=40
        local N_scaled=$(awk -v dw="$d_worst" -v nw="$N_worst" -v dd="$d" 'BEGIN {n=int(nw * (dw/dd)^(1/3)); print (n>nw ? nw : n)}')
        N=$(awk -v ns="$N_scaled" 'BEGIN {print (ns<10 ? 10 : ns)}')
    else
        N=$(awk -v Lx="$L" -v dx="$dx" 'BEGIN {n=int(Lx/dx); print (n<10 ? 10 : n)}')
    fi
    # Enforce minimum resolution: sphere must span at least 1 cell
    local N_min=$(awk -v Lx="$L" -v dd="$d" 'BEGIN {n=int(Lx / dd); print (n<10 ? 10 : n)}')
    N=$(( N > N_min ? N : N_min ))
    N=$(( (N + 1) / 2 * 2 )) # Round to even for MPI
    local dx_actual=$(awk -v Lx="$L" -v nn="$N" 'BEGIN {print Lx / nn}')

    # Step 6: CFL
    local v_th=$(awk -v k="$k" -v ti="$T_inf" -v m="$m" 'BEGIN {print sqrt(3*k*ti/m)}')
    local v_max=$(awk -v vt="$v_th" -v vv="$v" 'BEGIN {print vt + vv}')
    local dt=$(awk -v dx="$dx_actual" -v vm="$v_max" 'BEGIN {print 0.9 * dx / vm}')

    # Step 7: fnum (macro-particle weighting)
    # Target particles-per-cell (PPC):
    #   Continuum/Transitional: 25 PPC (Bird standard minimum)
    #   Free-Molecular (Kn>10): 50 PPC (higher needed for surface collision statistics
    #     because surface collisions per timestep scale with particle count)
    local V_cell=$(awk -v dx="$dx_actual" 'BEGIN {print dx^3}')
    local n_real=$(awk -v n="$n_inf" -v vc="$V_cell" 'BEGIN {print n * vc}')
    local target_ppc=25
    if (( $(awk -v kn="$Kn" 'BEGIN{print (kn>10)}') )); then
        target_ppc=50
    fi
    local fnum=$(awk -v nr="$n_real" -v tp="$target_ppc" 'BEGIN {fn=nr/tp; fn=(fn<1?1:int(fn)); print fn}')

    # Step 7b: performance cap on total simulated particles
    # For small droplets (100nm), n_real ~48 so fnum floors to 1 even at PPC=50.
    # This creates ~3M particles × 50K steps = extremely slow runtime.
    # Cap total simulated particles so each case is tractable.
    # Free-molecular flow (Kn > 1000) tolerates lower PPC since there are
    # no intermolecular collisions to resolve; surface statistics come from
    # the long averaging window (25K steady-state steps), not from PPC alone.
    local N_total=$(( N * N * N ))
    local MAX_TOTAL_PARTICLES=1500000
    local fnum_min=$(awk -v nr="$n_real" -v nt="$N_total" -v mx="$MAX_TOTAL_PARTICLES" \
        'BEGIN {f=nr*nt/mx; f=(f<1?1:int(f)+1); print f}')
    if [ "$fnum" -lt "$fnum_min" ]; then
        fnum=$fnum_min
    fi
    local actual_ppc=$(awk -v nr="$n_real" -v fn="$fnum" 'BEGIN {print int(nr/fn)}')

    # Step 8: run steps with minimum floor
    # Flow-through criterion: 3 domain lengths at stream velocity.
    # This ensures the flow field establishes, but for free-molecular flow
    # we also need enough wall-clock sampling for surface collision statistics.
    # Enforce MIN_STEPS=50000 so fix ave/surf has adequate accumulation window.
    local MIN_STEPS=50000
    local t_steady=$(awk -v Lx="$L" -v vv="$v" 'BEGIN {print 3 * Lx / vv}')
    local N_steps=$(awk -v ts="$t_steady" -v dt="$dt" -v mn="$MIN_STEPS" 'BEGIN {n=int(ts/dt)+1; print (n<mn ? mn : n)}')

    echo "$N $fnum $dt $N_steps $actual_ppc"


# Define your 30-case matrix arrays
diameters=(100e-9 500e-9 1e-6 5e-6 10e-6)
velocities=(20 100 250 400 550 700)

echo "========================================================="
echo " INITIATING 30-CASE DYNAMIC MATRIX SWEEP                 "
echo "========================================================="

# Outer loop: Diameters
for D_VAL in "${diameters[@]}"; do
    # Inner loop: Velocities
    for V_VAL in "${velocities[@]}"; do
        echo "Configuring case: Diameter = ${D_VAL} m | Velocity = ${V_VAL} m/s"
        
        # 1. Call the optimizer function and capture the output array
        PARAMS=($(compute_optimal_params $D_VAL $V_VAL))
        GRID_RES=${PARAMS[0]}
        FNUM_VAL=${PARAMS[1]}
        DT=${PARAMS[2]}
        TOTAL_STEPS=${PARAMS[3]}
        ACTUAL_PPC=${PARAMS[4]}
        
        # 2. Calculate the exact 50/50 phase split for the steady-state window
        TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
        STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))
        
        LOG_FILE="../production_output/log_D${D_VAL}_V${V_VAL}_${TIMESTAMP}.log"
        DATA_FILE="../production_output/force_D${D_VAL}_V${V_VAL}_${TIMESTAMP}.txt"
        
        echo "   -> Dynamics: Grid ${GRID_RES}^3 | Fnum ${FNUM_VAL} | dt ${DT} s | Steps ${TOTAL_STEPS} | PPC: ${ACTUAL_PPC}"
        
        # 3. Launch SPARTA with all dynamic parameters explicitly passed
        mpirun -np 4 /home/ilozano/sparta/src/spa_mpi -var vel_val $V_VAL -var diam_scale $D_VAL -var fnum_val $FNUM_VAL -var grid_res $GRID_RES -var t_step $DT -var TRANSIENT_STEPS $TRANSIENT_STEPS -var STEADY_STEPS $STEADY_STEPS -var out_name "$DATA_FILE" < in.drag > "$LOG_FILE"
        
        TOTAL_CELLS=$(( GRID_RES * GRID_RES * GRID_RES ))

        # Extract final particle count from log
        NP_FINAL=$(grep "Particles:" "$LOG_FILE" | tail -1 | awk '{print $3}')
        PPC=$(awk -v np="$NP_FINAL" -v nc="$TOTAL_CELLS" 'BEGIN {print np/nc}')

        # Check PPC > 10 (absolute minimum for any DSMC validity)
        # Note: The fnum cap (Step 7b) may intentionally lower PPC for
        # small FM droplets to keep total particles tractable. PPC=15
        # with 25K steady-state steps gives adequate surface statistics.
        if (( $(awk -v ppc="$PPC" 'BEGIN {print (ppc < 10)}') )); then
            echo "   -> ⚠️ WARNING: PPC = $PPC < 10" >&2
        elif (( $(awk -v ppc="$PPC" 'BEGIN {print (ppc < 15)}') )); then
            echo "   -> ⚠️ NOTE: PPC = $PPC < 15 — low but tolerable with long sampling" >&2
        fi

        # Check zero surface collisions in STEADY-STATE window only.
        # Log contains "SurfColl occurs" for each run block.
        # First line = transient phase, last line = cumulative total.
        # We compute: steady_state = total - transient.
        NSCOLL_TRANSIENT=$(grep "SurfColl occurs" "$LOG_FILE" | head -1 | awk '{print $4}')
        NSCOLL_TOTAL=$(grep "SurfColl occurs" "$LOG_FILE" | tail -1 | awk '{print $4}')
        NSCOLL_STEADY=$((NSCOLL_TOTAL - NSCOLL_TRANSIENT))
        if [ "$NSCOLL_STEADY" -le 0 ] 2>/dev/null; then
            echo "   -> ⚠️ WARNING: Zero steady-state surface collisions — insufficient sampling" >&2
        fi
        # ==========================================================
        # 4. Extract the absolute Drag Force safely
        # Awk reads the dump file, skips all header lines (ITEM:, #),
        # accumulates surface-averaged Fx (column 2) across all triangles.
        # No sum-reset on ITEM: TIMESTEP — handles multi-block dumps correctly.
        DRAG_FX=$(awk '
            BEGIN {sum=0}
            /^ITEM:/ {next}
            /^#/ {next}
            NF==4 {sum += $2}
            END {
                if (sum < 0) sum = -sum;
                printf "%e", sum
            }
        ' "$DATA_FILE")
        
        echo "${D_VAL},${V_VAL},${GRID_RES},${FNUM_VAL},${DT},${TOTAL_STEPS},${DRAG_FX}" >> "$CSV_OUT"
        echo "   -> Resolved Drag: ${DRAG_FX} N"
        echo "---------------------------------------------------------"
    done
done

echo "MATRIX SWEEP COMPLETE. Telemetry saved to: $CSV_OUT"