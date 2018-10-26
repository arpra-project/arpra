#!/bin/bash

echo Begin: ${0}
make experiments/experiment2

for i in {0..99}; do
    echo i: ${i}

    rm -f *.dat
    ./experiments/experiment2
    outdir=experiment_2_out/i_${i}

    echo saving to ${outdir}
    mkdir -p ${outdir}
    mv *.dat ${outdir}
done

python3 ./tools/arpra_mpfr_2d.py
