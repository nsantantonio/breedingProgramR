#!/bin/bash

projName=traditional

# rgscint=(0.05 0.1 0.2 0.3 0.4 0.5)
nFam=(4 10 20 40)
famSize=(25 50 100 200)

reps=1
nFounderPops=100
nThreads=25
nYr=20
QTL=1000QTL

# Removed the low heritability trial
selectTrials=("c(0.5, 0.5, 0.4, 0.5)" \
              "c(0.5, 0.4, 0.3, 1/3)" \
              "c(0.25, 0.2, 0.4, 0.25)" \
              "c(0.25, 0.25, 0.2, 0.2)" )
# simple truncation, no effort to not mating within family
for i in {0..3}; do 
    for j in 2 1 0; do
        for t in 1 2; do
            tmpSimName="truncVSquadprog${nYr}yr${QTL}_trad${t}_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=${nThreads}" \
                "nYr=${nYr}" \
                "simName=$tmpSimName" \
                "projName=${projName}" \
                "reps=${reps}" \
                "nFounderPops=${nFounderPops}" \
                "selectRGSC=0.2" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                "traditional=${t}"\
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
        done
    done
done

# add 300 additinal lgenotyping lines to select from: as comparison to separate RGSC (300 individuals per year)
for i in {0..3}; do 
    for j in 2 1 0; do
        for t in 1 2; do
            tmpSimName="truncVSquadprog${nYr}yr${QTL}_trad${t}_add300_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=${nThreads}" \
                "nYr=${nYr}" \
                "simName=$tmpSimName" \
                "projName=${projName}" \
                "reps=${reps}" \
                "nFounderPops=${nFounderPops}" \
                "selectRGSC=0.2" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                "traditional=${t}"\
                "withinFamInt=300"\
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
        done
    done
done

# within family , 
for i in {0..2}; do 
    for j in 2 1 0; do
        for t in 1 2; do
            for k in 0.2 1; do
                for l in 0.5 1; do
            tmpSimName="truncVSquadprog${nYr}yr${QTL}_trad${t}_intWithin${k}_intAcross${l}_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=${nThreads}" \
                "nYr=${nYr}" \
                "simName=$tmpSimName" \
                "projName=${projName}" \
                "reps=${reps}" \
                "nFounderPops=${nFounderPops}" \
                "selectRGSC=0.2" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                "traditional=${t}"\
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
        done
    done
done


for i in {0..2}; do 
    for j in 2 1 0; do
        for t in 1 2; do
            tmpSimName="truncVSquadprog${nYr}yr${QTL}_trad${t}_add300_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=${nThreads}" \
                "nYr=${nYr}" \
                "simName=$tmpSimName" \
                "projName=${projName}" \
                "reps=${reps}" \
                "nFounderPops=${nFounderPops}" \
                "selectRGSC=0.2" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                "traditional=${t}"\
                "withinFamInt=300"\
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
        done
    done
done

