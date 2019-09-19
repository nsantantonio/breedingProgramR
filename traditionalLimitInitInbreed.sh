#!/bin/bash

projName=traditionalLimitI

nFam=(10 20 40)
famSize=(50 75 100)

reps=1
nFounderPops=100
nThreads=100
nYr=30
QTL=1000QTL
lgen=4

selectTrials=("c(0.5, 0.2, 0.5, 0.4)" \
              "c(0.5, 0.10, 0.4, 1/3)" \
              "c(0.5, 0.05, 0.4, 0.25)" )

# test with phenotypes only!
# for i in 1; do 
#     for j in 0; do
#         for t in 2; do
#             for k in 0.2; do
#                 for l in 0.5; do
#                     nElite=$((${nFam[$i]} * 2))
#             tmpSimName="${projName}${nYr}yr${QTL}_trad${t}_intWithin${k}_intAcross${l}_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
#             args=("nThreads=${nThreads}" \
#                 "nYr=${nYr}" \
#                 "lgen=$lgen" \
#                 "simName=$tmpSimName" \
#                 "projName=${projName}" \
#                 "reps=${reps}" \
#                 "nFounderPops=${nFounderPops}" \
#                 "selectRGSC=0.2" \
#                 "useTruth=${j}" \
#                 "useIn=pheno" \
#                 "useOut=pheno" \
#                 "nFam=${nFam[$i]}" \
#                 "famSize=${famSize[$i]}"\
#                 "traditional=${t}"\
#                 "selectTrials=${selectTrials[$i]}" \
#                 "intWithin=${k}" \
#                 "intAcross=${l}"\
#                 "nElite=${nElite}")
#             Rscript $(pwd)/singleTraitProgram.R "${args[@]}" 
#                 done
#             done
#         done
#     done
# done


# test
for i in 1; do 
    for j in 0; do
        for t in 2; do
            for k in 0.2; do
                for l in 0.5; do
                    nElite=$((${nFam[$i]} * 2))
                    tmpSimName="${projName}${nYr}yr${QTL}_trad${t}_intWithin${k}_intAcross${l}_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
                    args=("nThreads=${nThreads}" \
                        "nYr=${nYr}" \
                        "lgen=$lgen" \
                        "simName=$tmpSimName" \
                        "projName=${projName}" \
                        "reps=${reps}" \
                        "nFounderPops=${nFounderPops}" \
                        "selectRGSC=0.2" \
                        "useTruth=${j}" \
                        "nFam=${nFam[$i]}" \
                        "famSize=${famSize[$i]}"\
                        "traditional=${t}"\
                        "selectTrials=${selectTrials[$i]}" \
                        "intWithin=${k}" \
                        "intAcross=${l}"\
                        "nElite=${nElite}")
                    echo "nElite: $nElite"
                    Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                done
            done
        done
    done
done


# traditional selection,
for i in {0..2}; do 
    for j in 0; do
        for t in 2 3; do
            for k in 0.2 0.5 1; do
                for l in 0.5 1; do
                    nElite=$((${nFam[$i]} * 2))
                    tmpSimName="${projName}${nYr}yr${QTL}_trad${t}_intWithin${k}_intAcross${l}_truth${j}_vdp${nFam[$i]}x${famSize[$i]}"
                    args=("nThreads=${nThreads}" \
                        "nYr=${nYr}" \
                        "lgen=$lgen" \
                        "simName=$tmpSimName" \
                        "projName=${projName}" \
                        "reps=${reps}" \
                        "nFounderPops=${nFounderPops}" \
                        "selectRGSC=0.2" \
                        "useTruth=${j}" \
                        "nFam=${nFam[$i]}" \
                        "famSize=${famSize[$i]}"\
                        "traditional=${t}"\
                        "selectTrials=${selectTrials[$i]}" \
                        "intWithin=${k}" \
                        "intAcross=${l}"\
                        "nElite=${nElite}")
                    Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                done
            done
        done
    done
done

