#!/bin/bash

projName=truncation

# rgscint=(0.05 0.1 0.2 0.3 0.4 0.5)
nFam=(4 10 20 40)
famSize=(25 50 100 200)

reps=1
nFounderPops=100
nThreads=35
nYr=20
QTL=1000QTL

# Removed the low heritability trial
selectTrials=("c(0.5, 0.5, 0.4, 0.5)" \
              "c(0.5, 0.4, 0.3, 1/3)" \
              "c(0.25, 0.2, 0.4, 0.25)" \
              "c(0.25, 0.25, 0.2, 0.2)" )
# truncation selection
for j in 0 1 2; do
    for i in {0..3}; do 
        tmpSimName="${projName}${nYr}yr${QTL}_trunc_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
        args=("nThreads=${nThreads}" \
            "nYr=${nYr}" \
            "simName=$tmpSimName" \
            "projName=${projName}" \
            "reps=${reps}" \
            "nFounderPops=${nFounderPops}" \
            "useTruth=${j}" \
            "nFam=${nFam[$i]}" \
            "famSize=${famSize[$i]}"\
            "selectTrials=${selectTrials[$i]}")
        Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
     done
done
