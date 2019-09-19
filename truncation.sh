#!/bin/bash

projName=truncation

nFam=(10 20 40)
famSize=(50 75 100)

reps=1
nFounderPops=100
nThreads=25
nYr=30
QTL=1000QTL
lgen=4

selectTrials=("c(0.5, 0.2, 0.5, 0.4)" \
              "c(0.5, 0.10, 0.4, 1/3)" \
              "c(0.5, 0.05, 0.4, 0.25)" )

# truncation selection
for j in 0; do
    for i in 1 0 2; do 
        for r in 0.2 0.3 0.1; do 
            tmpSimName="${projName}${nYr}yr${QTL}_trunc_truth${j}_rgsc${r}_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=${nThreads}" \
                "lgen=$lgen" \
                "nYr=${nYr}" \
                "simName=$tmpSimName" \
                "projName=${projName}" \
                "reps=${reps}" \
                "nFounderPops=${nFounderPops}" \
                "selectRGSC=${r}" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
         done
    done
done

for j in 0; do
    for i in 2 0 1; do 
        for r in 0.5 0.4; do 
            tmpSimName="${projName}${nYr}yr${QTL}_trunc_truth${j}_rgsc${r}_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=${nThreads}" \
                "lgen=$lgen" \
                "nYr=${nYr}" \
                "simName=$tmpSimName" \
                "projName=${projName}" \
                "reps=${reps}" \
                "nFounderPops=${nFounderPops}" \
                "selectRGSC=${r}" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
         done
    done
done
