#!/bin/bash

echo Begin: ${0}
cd $(dirname ${0})/..
make experiments/experiment_2
make extra/morris_lecar
cd experiments
outdir=experiment_2_out
rm -rf ${outdir}

for i in {0..99}; do
    echo i: ${i}

    outsubdir=${outdir}/i_${i}
    echo Changing to ${outsubdir}
    mkdir -p ${outsubdir}
    cd ${outsubdir}
    ../../experiment_2 0 2>/dev/null
    cd ../..
done

echo ascending
outsubdir=${outdir}/ascending
echo Changing to ${outsubdir}
mkdir -p ${outsubdir}
cd ${outsubdir}
../../experiment_2 1 2>/dev/null
cd ../..

echo descending
outsubdir=${outdir}/descending
echo Changing to ${outsubdir}
mkdir -p ${outsubdir}
cd ${outsubdir}
../../experiment_2 2 2>/dev/null
cd ../..

echo Arpra sim
outsubdir=${outdir}/arpra
echo Changing to ${outsubdir}
mkdir -p ${outsubdir}
cd ${outsubdir}
../../../extra/morris_lecar
cd ../..

python3 ./experiment_2.py
