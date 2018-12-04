#!/bin/bash

fig="fig5"

echo Begin: ${0}
cd $(dirname ${0})/../..
make experiments/${fig}/ml-arpra
make experiments/${fig}/ml-mpfr
make experiments/${fig}/ml-mpfi
cd experiments/${fig}
outdir=${fig}_out
rm -rf ${outdir}

# for i in {0..99}; do
#     echo MPFR ${i}
#     outsubdir=${outdir}/i_${i}
#     echo Changing to ${outsubdir}
#     mkdir -p ${outsubdir}
#     cd ${outsubdir}
#     cmd="../../ml-mpfr 53 0"
#     echo ${cmd}
#     ${cmd} 2>/dev/null
#     cd ../..
# done

# echo MPFR ascending
# outsubdir=${outdir}/ascending
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../ml-mpfr 53 1"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

# echo MPFR descending
# outsubdir=${outdir}/descending
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../ml-mpfr 53 2"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

# echo MPFR high precision
# outsubdir=${outdir}/high_prec
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../ml-mpfr 2048 0"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

# echo MPFI interval bounds
# outsubdir=${outdir}/mpfi
# echo Changing to ${outsubdir}
# mkdir -p ${outsubdir}
# cd ${outsubdir}
# cmd="../../ml-mpfi"
# echo ${cmd}
# ${cmd} 2>/dev/null
# cd ../..

echo Arpra bounds
outsubdir=${outdir}/arpra
echo Changing to ${outsubdir}
mkdir -p ${outsubdir}
cd ${outsubdir}
cmd="../../ml-arpra 53 2048 500 50"
echo ${cmd}
${cmd} 2>/dev/null
cd ../..

python3 ./${fig}.py
