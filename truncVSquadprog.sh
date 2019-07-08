#!/bin/bash

projName=${projectName}

# rgscint=(0.05 0.1 0.2 0.3 0.4 0.5)
nFam=(4 10 20 40)
famSize=(25 50 100 200)

reps=10
nFounderPops=5
nThreads=25
nYr=20

# Removed the low heritability trial
selectTrials=("c(0.5, 0.5, 0.4, 0.5)" \
              "c(0.5, 0.4, 0.3, 1/3)" \
              "c(0.25, 0.2, 0.4, 0.25)" \
              "c(0.25, 0.25, 0.2, 0.2)" )

# h2=("c(0.1, 0.3, 0.3, 0.3, 0.3)" \
#     "c(0.1, 0.3, 0.3, 0.3, 0.3)" \
#     "c(0.1, 0.3, 0.3, 0.3, 0.3)" \
#     "c(0.1, 0.3, 0.3, 0.3, 0.3)" )

# trialReps=("c(1, 1, 2, 3, 3)" \
#            "c(1, 1, 2, 3, 3)" \
#            "c(1, 1, 2, 3, 3)" \
#            "c(1, 1, 2, 3, 3)")

# trialLocs=("c(1, 1, 2, 5, 5)" \
#            "c(1, 1, 2, 5, 5)" \
#            "c(1, 1, 2, 5, 5)" \
#            "c(1, 1, 2, 5, 5)")

# returnVDPtoRGSC=("c(0, 0, 0, 0, 0, 0)" \
#                  "c(0, 0, 0, 0, 0, 0)" \
#                  "c(0, 0, 0, 0, 0, 0)" \
#                  "c(0, 0, 0, 0, 0, 0)")


# use solqp, full factorial... DONT RUN!
# for i in {0..3}; do 
#     for k in 0.05 0.01 0.005; do
#         for j in 2 1 0; do
#             for l in 0.1 0.2 0.3; do
#                 for n in 0 1 2; do
#                     for p in {0..3}; do
#                         tmpSimName="truncVSquadprog30yr1000QTL_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#                         args=("nThreads=${nThreads}" \
#                             "nYr=${nYr}" \
#                             "simName=$tmpSimName" \
#                             "projName=${projectName}" \
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
#                             "withinFamInt=1" \
#                             "fthresh=${k}" \
#                             "fthreshOut=${l}" \
#                             "selectTrials=${selectTrials[$i]}")
#                         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#                     done
#                 done
#             done
#         done
#     done
# done



# use solqp, test, 
for i in 1; do 
    for k in 0.01; do
        for j in 0; do
            for l in 0.2; do
                for n in 1; do
                    for p in 0; do
                        tmpSimName="TESTtruncVSquadprog30yr1000QTL_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
                        args=("nThreads=${nThreads}" \
                            "nYr=${nYr}" \
                            "simName=$tmpSimName" \
                            "projName=TEST" \
                            "reps=2" \
                            "nFounderPops=2" \
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
                        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                    done
                done
            done
        done
    done
done


# use solqp, vary fthreshout, 
for i in 1; do 
    for k in 0.01; do
        for j in 2 1 0; do
            for l in 0.1 0.2 0.3; do
                for n in 1; do
                    for p in 0; do
                        tmpSimName="truncVSquadprog30yr1000QTL_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
                        args=("nThreads=${nThreads}" \
                            "nYr=${nYr}" \
                            "simName=$tmpSimName" \
                            "projName=${projectName}" \
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
                        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
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
            for l in 0.1; do
                for n in 0 1 2; do
                    for p in 0; do
                        tmpSimName="truncVSquadprog30yr1000QTL_quadprog_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
                        args=("nThreads=${nThreads}" \
                            "nYr=${nYr}" \
                            "simName=$tmpSimName" \
                            "projName=${projectName}" \
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
                        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                    done
                done
            done
        done
    done
done


# use solqp, dont branch....
for i in 1; do 
    for k in 0.01; do
        for j in 2 1 0; do
            tmpSimName="truncVSquadprog30yr1000QTL_quadprog_fin${k}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
                args=("nThreads=${nThreads}" \
                    "nYr=${nYr}" \
                    "simName=$tmpSimName" \
                    "projName=${projectName}" \
                    "reps=${reps}" \
                    "nFounderPops=${nFounderPops}" \
                    "selectRGSC=0.2" \
                    "useTruth=${j}" \
                    "nFam=${nFam[$i]}" \
                    "famSize=${famSize[$i]}"\
                    "selFuncOut=(solqpOut)" \
                    "fthresh=${k}" \
                    "selectTrials=${selectTrials[$i]}")
                Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
            done
        done
    done
done


# truncation selection

for i in {0..3}; do 
    for j in 2 1 0; do
        tmpSimName="truncVSquadprog30yr1000QTL_trunc_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
        args=("nThreads=${nThreads}" \
            "nYr=${nYr}" \
            "simName=$tmpSimName" \
            "projName=${projectName}" \
            "reps=${reps}" \
            "nFounderPops=${nFounderPops}" \
            "useTruth=${j}" \
            "nFam=${nFam[$i]}" \
            "famSize=${famSize[$i]}"\
            "selectTrials=${selectTrials[$i]}")
        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
     done
done


