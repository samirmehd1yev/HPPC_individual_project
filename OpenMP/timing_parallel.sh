#!/bin/bash
OUTPUT_FILE=../timing_parallel.txt
echo -n "" > $OUTPUT_FILE


N_VALUES=(10 20 30 40 50 60 70 80 90 100 150 200 300 400 500 600 700 800 900 1000 1500 2000 3000 4000 5000 6000 7000 8000 9000 10000)


THREAD_COUNTS=(2 4 6 8)


for N in "${N_VALUES[@]}"
do

  N_PADDED=$(printf "%05d" $N)


  for THREADS in "${THREAD_COUNTS[@]}"
  do
    output=$(./galsim $N ../input_data/ellipse_N_${N_PADDED}.gal 2000 1e-5 0.5 0 $THREADS)
    total_time=$(echo "$output" | grep "Total time" | awk '{print $3}')
    echo "$N $THREADS $total_time" >> $OUTPUT_FILE
  done
done
