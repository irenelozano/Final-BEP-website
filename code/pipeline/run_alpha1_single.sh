#!/bin/bash
# =========================================================================
# run_alpha1_single.sh - Single alpha = 1 CLL control simulation
# Execute from inside:
# sparta/simulations/tin_drag/pipeline/
# =========================================================================

# -------------------------------------------------------------------------
# Selected control case
# -------------------------------------------------------------------------
# d = 1 micrometer
# v = 700 m/s
# N = 40
#
# Reason:
# - Kn is still clearly free-molecular.
# - Surface statistics are better than for 100 nm.
# - It avoids the more transitional 10 um case.
# - High velocity makes finite-speed effects relevant.
# -------------------------------------------------------------------------

VEL_VAL=700
DIAM_SCALE=1e-6
GRID_RES=34

# Match the same fnum scaling used in the original matrix script.
# For d = 1 um:
# vol_ratio = (1e-6 / 100e-9)^3 = 1000
# fnum = 10000 * 1000 = 1e7
FNUM_VAL=2380

# CFL-safe timestep for d = 1 um, N = 40.
# dx = 1.25 um and v_th + v0 ≈ 2606 m/s,
# so dt_CFL ≈ 4.8e-10 s.
T_STEP=4.3e-10

TOTAL_STEPS=50000
TRANSIENT_STEPS=$((TOTAL_STEPS / 2))
STEADY_STEPS=$((TOTAL_STEPS - TRANSIENT_STEPS))

# Output base name
OUT_NAME="../runtime_output/alpha1_d1um_v700"

# Ensure output directory exists
mkdir -p ../runtime_output

echo "Launching alpha = 1 control case"
echo "Velocity      = ${VEL_VAL} m/s"
echo "Diameter      = ${DIAM_SCALE} m"
echo "Grid          = ${GRID_RES}^3"
echo "Fnum          = ${FNUM_VAL}"
echo "Timestep      = ${T_STEP} s"
echo "Output prefix = ${OUT_NAME}"
echo "-----------------------------------------------------"

mpirun -np 4 /home/ilozano/sparta/src/spa_mpi \
    -var vel_val "$VEL_VAL" \
    -var diam_scale "$DIAM_SCALE" \
    -var fnum_val "$FNUM_VAL" \
    -var grid_res "$GRID_RES" \
    -var t_step "$T_STEP" \
    -var TRANSIENT_STEPS "$TRANSIENT_STEPS" \
    -var STEADY_STEPS "$STEADY_STEPS" \
    -var out_name "$OUT_NAME" \
    -var mesh_file ../global_assets/data.sphere \
    < validation_in.drag > "../runtime_output/alpha1_d1um_v700.log"

echo "Alpha = 1 control case complete."
echo "Log file:"
echo "../runtime_output/alpha1_d1um_v700.log"
echo "Global cumulative force:"
echo "../runtime_output/alpha1_d1um_v700_global_cumul.txt"
echo "Global block force history:"
echo "../runtime_output/alpha1_d1um_v700_global_blocks.txt"