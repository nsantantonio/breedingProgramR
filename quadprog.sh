#!/bin/bash

projName=quadprog

# rgscint=(0.05 0.1 0.2 0.3 0.4 0.5)
nFam=(4 10 20 40)
famSize=(25 50 100 200)

reps=1
nFounderPops=100
nThreads=50
nYr=30
QTL=1000QTL

# Removed the low heritability trial
selectTrials=("c(0.5, 0.5, 0.4, 0.5)" \
              "c(0.5, 0.4, 0.3, 1/3)" \
              "c(0.25, 0.2, 0.4, 0.25)" \
              "c(0.25, 0.25, 0.2, 0.2)" )

# use solqp, vary fthreshout, 
for i in {0..2}; do 
    for k in 0.005 0.01 0.05; do
        for j in 2 1 0; do
            for l in 0.05 0.1 0.2; do
                for n in 1; do
                    for p in 0; do
                        tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
                            "limitN=${n}" \
                            "pullCycle=${p}" \
                            "selFuncOut=(solqpOut)" \
                            "selFuncIn=(solqp)" \
                            "fthresh=${k}" \
                            "fthreshOut=${l}" \
                            "selectTrials=${selectTrials[$i]}")
                        Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                    done
                done
            done
        done
    done
done


# use solqp, vary limitN
for i in 1; do 
    for k in 0.01; do
        for j in 2 1 0; do
            for l in 0.1 0.2; do
                for n in 0 1 2; do
                    for p in 0; do
                        tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
                            "limitN=${n}" \
                            "pullCycle=${p}" \
                            "selFuncOut=(solqpOut)" \
                            "selFuncIn=(solqp)" \
                            "fthresh=${k}" \
                            "fthreshOut=${l}" \
                            "selectTrials=${selectTrials[$i]}")
                        Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                    done
                done
            done
        done
    done
done



# use solqp, pull generation
for i in {0..2}; do 
    for k in 0.05 0.005; do
        for j in 2 1 0; do
            for l in 0.1 0.2; do
                for n in 1; do
                    for p in {0..2}; do
                        tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
                            "limitN=${n}" \
                            "pullCycle=${p}" \
                            "selFuncOut=(solqpOut)" \
                            "selFuncIn=(solqp)" \
                            "fthresh=${k}" \
                            "fthreshOut=${l}" \
                            "selectTrials=${selectTrials[$i]}")
                        Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                    done
                done
            done
        done
    done
done


# use solqp, dont branch....
for i in {0..2}; do 
    for k in 0.05 0.01 0.005; do
        for j in 2 0 1; do
            tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_pull3_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
                "selFuncIn=(solqp)" \
                "fthresh=${k}" \
                "selectTrials=${selectTrials[$i]}")
            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
        done
    done
done



# # use solqp, redos
# for i in 2; do 
#     for k in 0.05; do
#         for j in 2 1 0; do
#             for l in 0.05 0.1; do
#                 for n in 1; do
#                     for p in 0; do
#                         tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#                         args=("nThreads=${nThreads}" \
#                             "nYr=${nYr}" \
#                             "simName=$tmpSimName" \
#                             "projName=${projName}" \
#                             "reps=${reps}" \
#                             "nFounderPops=${nFounderPops}" \
#                             "selectRGSC=0.2" \
#                             "useTruth=${j}" \
#                             "nFam=${nFam[$i]}" \
#                             "famSize=${famSize[$i]}"\
#                             "limitN=${n}" \
#                             "pullCycle=${p}" \
#                             "selFuncOut=(solqpOut)" \
#                             "selFuncIn=(solqp)" \
#                             "fthresh=${k}" \
#                             "fthreshOut=${l}" \
#                             "selectTrials=${selectTrials[$i]}")
#                         Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#                     done
#                 done
#             done
#         done
#     done
# done


