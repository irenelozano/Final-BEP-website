#!/bin/bash
# =========================================================================
# run_matrix.sh - Parametric Simulation Sweep Pipeline
# Execute from inside the sparta/simulations/tin_drag/pipeline/ directory
# =========================================================================

# Define the matrices
velocities=(20 60 100 250 500 700)
diameters=(100e-9 500e-9 1e-6 5e-6 10e-6)

# WARNING: This script is deprecated. Use run_30case_matrix.sh instead.
# Issues with this script:
#   1. No -var t_step passed (uses broken default 1e-9, CFL-violating for d<1um)
#   2. Fixed N=40 under-resolves 5/10 um spheres (10 cells in 500 um domain)
#   3. Fixed fnum scaling (10000 * vol_ratio) ignores PPC targets
#   4. No transient/steady phase separation in force output
# Kept only for reference / backward compatibility.

# Fixed grid resolution (see run_30case_matrix.sh for dynamic scaling)
GRID_RES=40

# Two-phase steady-state sampling window
TOTAL_STEPS=50000
TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))
SAMPLE_INTERVAL=1000
SAMPLE_BLOCKS=$((STEADY_STEPS / SAMPLE_INTERVAL))

# Ensure the output directory exists
mkdir -p ../runtime_output

# Main orchestration loop
for v in "${velocities[@]}"
do
    for d in "${diameters[@]}"
    do
        echo "Launching: Velocity = $v m/s | Diameter = $d m"
        
        # Calculate dynamic fnum weighting factor on the fly
        # V_ratio scales from 100nm baseline up to the targeted size
        vol_ratio=$(awk -v d="$d" 'BEGIN {print (d / 100e-9)^3}')
        fnum=$(awk -v vr="$vol_ratio" 'BEGIN {print 10000 * vr}')
        
        # Define clean, explicit output filenames
        log_file="../runtime_output/sim_v${v}_d${d}.log"
        data_file="../runtime_output/data_v${v}_d${d}.txt"
        
        # Execute binary: inject parameters as runtime arguments using -var flag
        mpirun -np 4 ../../../src/spa_mpi \
            -var vel_val "$v" \
            -var diam_scale "$d" \
            -var fnum_val "$fnum" \
            -var grid_res $GRID_RES \
            -var t_step 1e-9 \
            -var TRANSIENT_STEPS $TRANSIENT_STEPS \
            -var STEADY_STEPS $STEADY_STEPS \
            -var out_name "$data_file" \
            < in.drag > "$log_file"
            
        echo "Case complete. Telemetry piped to $log_file"
        echo "-----------------------------------------------------"
    done
done
echo "All parametric matrix sweeps complete!"