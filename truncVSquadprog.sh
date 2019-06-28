#!/bin/bash

# rgscint=(0.05 0.1 0.2 0.3 0.4 0.5)
nFam=(4 10 20 40)
famSize=(25 50 100 200)

# 100, 
selectTrials=("c(1, 0.5, 0.5, 0.4, 0.5)" \
              "c(0.5, 0.5, 0.4, 0.3, 1/3)" \
              "c(0.5, 0.25, 0.2, 0.4, 0.25)" \
              "c(0.25, 0.25, 0.25, 0.2, 0.2)" )

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



# use true values solqp
for i in {0..3}; do 
    for k in 0.01 0.005; do
        for j in 2 1 0; do
          # i=0
          # k=0.01
          # j=2
            tmpSimName="truncVSquadprog30yr1000QTL_quadprog_fin${k}_fout0.2_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
            args=("nThreads=5" \
                "nYr=30" \
                "simName=$tmpSimName" \
                "projName=truncVSquadprog" \
                "reps=2" \
                "nFounderPops=2" \
                "selectRGSC=0.2" \
                "useTruth=${j}" \
                "nFam=${nFam[$i]}" \
                "famSize=${famSize[$i]}"\
                # "h2=${h2[$i]}" \
                "pullCycle=1"
                "selFuncOut=(solqpOut)" \
                "selFuncIn=(solqp)" \
                "inbreedFunc=(withinFamSel)" \
                "withinFamInt=1" \
                "fthresh=${k}" \
                "fthreshOut=0.2" \
                "selectTrials=${selectTrials[$i]}")
                # "trialReps=${trialReps[$i]}" \
                # "trialLocs=${trialLocs[$i]}" \
                # "returnVDPtoRGSC=${returnVDPtoRGSC[$i]}")
            # echo ${args[@]}
            # printf '%s\n' "${args[@]}"
            Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
            # Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}"
         done
    done
done


# # use true values solqp
# for i in {0..3}; do 
#     for k in 0.01 0.005; do
#         for j in 2 1 0; do
#             tmpSimName="truncVSquadprog30yr1000QTL_quadprog_fin${k}_fout0.2_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#             args=("nThreads=30" \
#                 "nYr=30" \
#                 "simName=$tmpSimName" \
#                 "projName=truncVSquadprog" \
#                 "selectRGSC=0.2" \
#                 "useTruth=$j" \
#                 "nFam=${nFam[$i]}" \
#                 "famSize=${famSize[$i]}"\
#                 "h2=${h2[$i]}" \
#                 "pullCycle=1"
#                 "selFuncOut=solqpOut" \
#                 "selFuncIn=solqp" \
#                 "inbreedFunc=withinFamSel" \
#                 "withinFamInt=1" \
#                 "fthresh=$k" \
#                 "fthreshOut=0.2" \
#                 "selectTrials=${selectTrials[$i]}" \
#                 "trialReps=${trialReps[$i]}" \
#                 "trialLocs=${trialLocs[$i]}" \
#                 "returnVDPtoRGSC=${returnVDPtoRGSC[$i]}")
#             Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#          done
#     done
# done

# # use true values truncation
# for i in {0..3}; do 
#     for j in 2 1 0; do
#         tmpSimName="truncVSquadprog30yr1000QTL_trunc_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#         args=("nThreads=30" \
#             "nYr=30" \
#             "simName=$tmpSimName" \
#             "projName=truncVSquadprog" \
#             "selectRGSC=0.2" \
#             "useTruth=$j" \
#             "nFam=${nFam[$i]}" \
#             "famSize=${famSize[$i]}"\
#             "h2=${h2[$i]}" \
#             "selFuncOut=NULL" \
#             "selFuncIn=NULL" \
#             "inbreedFunc=withinFamSel" \
#             "withinFamInt=1" \
#             "selectTrials=${selectTrials[$i]}" \
#             "trialReps=${trialReps[$i]}" \
#             "trialLocs=${trialLocs[$i]}" \
#             "returnVDPtoRGSC=${returnVDPtoRGSC[$i]}")
#         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#     done
# done

