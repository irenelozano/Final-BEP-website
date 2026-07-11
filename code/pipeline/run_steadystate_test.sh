#!/bin/bash

##################################################################
# RUN SCRIPT: One-case steady-state verification
##################################################################

# -----------------------------------------------------------------
# Case definition
# -----------------------------------------------------------------

DIAM_SCALE=1e-5
VEL_VAL=700
GRID_RES=40
FNUM_VAL=2e6
T_STEP=4.3e-9

TRANSIENT_STEPS=5000
STEADY_STEPS=5000

OUT_NAME="../runtime_output/steady_state_d10um_v700"

# Make sure output folder exists
mkdir -p ../runtime_output

# -----------------------------------------------------------------
# Run SPARTA with MPI
# -----------------------------------------------------------------

mpirun -np 4 /home/ilozano/sparta/src/spa_mpi \
    -var diam_scale ${DIAM_SCALE} \
    -var vel_val ${VEL_VAL} \
    -var grid_res ${GRID_RES} \
    -var fnum_val ${FNUM_VAL} \
    -var t_step ${T_STEP} \
    -var TRANSIENT_STEPS ${TRANSIENT_STEPS} \
    -var STEADY_STEPS ${STEADY_STEPS} \
    -var out_name ${OUT_NAME} \
    < in.steadystate_test