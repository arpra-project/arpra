#!/bin/bash

echo Begin: ${0}
cd $(dirname ${0})/..
make experiments/experiment_2
cd experiments
outdir=experiment_2_out
rm -rf ${outdir}

for i in {0..99}; do
    echo i: ${i}

    outsubdir=${outdir}/i_${i}
    echo Changing to ${outsubdir}
    mkdir -p ${outsubdir}
    cd ${outsubdir}
    ../../experiment_2
    cd ../..
done

python3 ../tools/arpra_mpfr_2d.py
