#! /bin/sh
pardir=$(pwd) # just in case $PWD doesnt exist on your system

# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings
for i in 0.1 0.3 0.5 0.7 0.9; do
    args=("seed=12345" \
        "simName=RGSCintensity${i}" \
        "selectRGSC=${i}"
    )
    # echo ${args[1]} ${args[2]}
    Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" &
done



for i in 0.01 0.05 0.1 0.2 0.3 0.5 0.7 0.9; do
    args=("seed=12345" \
        "simName=RGSCintensity${i}" \
        "selectRGSC=${i}"
        "nYr=30"
    )
    # echo ${args[1]} ${args[2]}
    Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" &
done

for i in 0.01 0.05 0.1 0.2 0.3 0.5 0.7 0.9; do
    args=("seed=12345" \
        "simName=RGSCintensity${i}_nNucl1000" \
        "selectRGSC=${i}"
        "nYr=30" \
        "nNuclear=1000"
    )
    # echo ${args[1]} ${args[2]}
    Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" &
done
