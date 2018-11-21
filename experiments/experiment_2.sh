#!/bin/bash

echo Begin: ${0}
cd $(dirname ${0})/..
make experiments/experiment_1
make experiments/experiment_2
cd experiments
outdir=experiment_2_out
rm -rf ${outdir}

# for i in {0..99}; do
#     echo MPFR ${i}
#     outsubdir=${outdir}/i_${i}
#     echo Changing to ${outsubdir}
#     mkdir -p ${outsubdir}
#     cd ${outsubdir}
#     cmd="../../experiment_2 24 0"
#     echo ${cmd}
#     ${cmd} 2>/dev/null
#     cd ../..
# done

# echo MPFR ascending
# outsubdir=${outdir}/ascending
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../experiment_2 24 1"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

# echo MPFR descending
# outsubdir=${outdir}/descending
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../experiment_2 24 2"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

# echo MPFR high precision
# outsubdir=${outdir}/high_prec
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../experiment_2 2048 0"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

echo Arpra bounds
outsubdir=${outdir}/arpra
echo Changing to ${outsubdir}
mkdir -p ${outsubdir}
cd ${outsubdir}
cmd="../../experiment_1 24 2048 500 50"
echo ${cmd}
${cmd} 2>/dev/null
cd ../..

python3 ./experiment_2.py
