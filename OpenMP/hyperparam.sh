#!/bin/bash
REFERENCE_FILE=true.gal
INPUT_FILE=../input_data/ellipse_N_03000.gal


total_simulation_time=$(echo "20000 * 1e-8" | bc -l)

for theta_max in 0.02 0.05 0.1 0.2 0.3 0.4 0.5
do
  for delta_t in 1e-3 1e-4 1e-5
  do
    steps=$(echo "$total_simulation_time / $delta_t" | bc)

    ./galsim 3000 $INPUT_FILE $steps $delta_t $theta_max 0 6
    # results error
    pos_maxdiff=$(../compare_gal_files/compare_gal_files 3000 result.gal $REFERENCE_FILE | grep "pos_maxdiff" | awk '{print $3}')

    # save results
    echo "$theta_max $delta_t $pos_maxdiff" >> results.txt
  done
done
