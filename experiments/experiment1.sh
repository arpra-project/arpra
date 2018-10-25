#!/bin/bash

echo Begin: ${0}
make experiments/experiment1

for inputs in {0..50..2}; do
    echo inputs: ${inputs}

    for frequency in {0..50..2}; do
        echo frequency: ${frequency}

        rm -f *.dat
        ./experiments/experiment1 ${inputs} ${frequency}
        outdir=experiment_1_out/in_${inputs}_freq_${frequency}

        echo saving to ${outdir}
        mkdir -p ${outdir}
        mv *.dat ${outdir}
    done
done
