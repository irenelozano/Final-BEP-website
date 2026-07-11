#!/bin/bash
# =========================================================================
# run_null.sh
#   Runs v=0 for all 5 diameters to measure F₀(d) directly.
#   Also runs 100nm at v=700 to check if F₀ changes with velocity.
# =========================================================================

SPARTA=~/sparta/build/sparta
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
mkdir -p ../production_output
OUT_CSV="../production_output/null_results_${TIMESTAMP}.csv"
echo "Diameter_m,Velocity_ms,Grid_N,Fnum,Timestep_s,Total_Steps,Drag_Fx" > "$OUT_CSV"

# (diameter, grid_N, fnum, t_step)
cases=(
    "100e-9 40 3 5.84431e-11"
    "500e-9 24 558 4.87027e-10"
    "1e-6   18 10588 1.29874e-09"
    "5e-6   10 7718750 1.16886e-08"
    "10e-6  10 61750000 2.33772e-08"
)

for case in "${cases[@]}"; do
    read -r d N fnum t_step <<< "$case"
    diam_scale=$d
    out_name="../production_output/null_d${d}_v0"

    echo "--- d=$d @ v=0 (N=$N, fnum=$fnum) ---"
    mpirun -np 1 "$SPARTA" -var out_name "$out_name" \
        -var diam_scale "$diam_scale" \
        -var grid_res "$N" \
        -var fnum_val "$fnum" \
        -var vel_val 0 \
        -var t_step "$t_step" \
        -in in.drag 2>&1 | tail -3

    drag_val=$(awk 'BEGIN {sum=0; n=0}
        /^ITEM:/ {next} /^#/ {next}
        NF==4 {sum += $2; n++}
        END {if (n>0 && sum!=0) print sum; else print "0"}' \
        "$out_name" 2>/dev/null)

    echo "$d,0,$N,$fnum,$t_step,50000,$drag_val" >> "$OUT_CSV"
    echo "  F₀($d) = $drag_val N"
done

# ---- Check velocity dependence: 100nm at v=700 ----
echo ""
echo "--- 100nm @ v=700 (same N=40, fnum=3) ---"
out_name="../production_output/null_d100e-9_v700"
mpirun -np 1 "$SPARTA" -var out_name "$out_name" \
    -var diam_scale 100e-9 \
    -var grid_res 40 \
    -var fnum_val 3 \
    -var vel_val 700 \
    -var t_step 5.84431e-11 \
    -in in.drag 2>&1 | tail -3

drag_val_700=$(awk 'BEGIN {sum=0; n=0}
    /^ITEM:/ {next} /^#/ {next}
    NF==4 {sum += $2; n++}
    END {if (n>0) print sum; else print "0"}' \
    "../production_output/null_d100e-9_v700" 2>/dev/null)

echo "  100nm v=700: F_raw = $drag_val_700 N"
F0_100nm=$(tail -1 "$OUT_CSV" | cut -d',' -f7 | head -1)
echo "  100nm v=0:   F₀     = $(awk -F, 'NR==2{print $7}' $OUT_CSV) N"
echo "  Difference (≈F_drag) = $(awk -F, -v v700="$drag_val_700" 'NR==2{printf "%.4e", v700-$7}' $OUT_CSV) N"
echo ""
echo "Results in $OUT_CSV"
