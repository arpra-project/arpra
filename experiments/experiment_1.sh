#!/bin/bash

echo Begin: ${0}
cd $(dirname ${0})/..
make experiments/experiment_1
cd experiments
outdir=experiment_1_out
rm -rf ${outdir}

for inputs in {0..50..2}; do
    echo inputs: ${inputs}

    for frequency in {0..50..2}; do
        echo frequency: ${frequency}

        outsubdir=${outdir}/in_${inputs}_freq_${frequency}
        echo Changing to ${outsubdir}
        mkdir -p ${outsubdir}
        cd ${outsubdir}
        ../../experiment_1 ${inputs} ${frequency} 2>/dev/null
        cd ../..
    done
done

python3 ./experiment_1.py
