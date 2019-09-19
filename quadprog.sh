#!/bin/bash

projName=quadprog

nFam=(10 20 40)
famSize=(50 75 100)


reps=1
nFounderPops=100
nThreads=35
nYr=30
QTL=1000QTL
lgen=4


selectTrials=("c(0.5, 0.2, 0.5, 0.4)" \
              "c(0.5, 0.10, 0.4, 1/3)" \
              "c(0.5, 0.05, 0.4, 0.25)" )

# # test
# for i in 1; do 
#     for k in 0.01; do
#         for j in 0; do
#             for l in 0.2; do
#                 for n in 1; do
#                     for p in 0; do
#                         tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#                         args=("nThreads=${nThreads}" \
#                             "nYr=${nYr}" \
#                             "lgen=$lgen" \
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


# use solqp, dont branch....
for i in {0..2}; do 
# for i in 2; do 
   for k in 0.02 0.05 0.001; do
       for j in 0; do
           tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_pull3_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
               "selFuncIn=(solqp)" \
               "fthresh=${k}" \
               "selectTrials=${selectTrials[$i]}")
           Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
       done
   done
done

# use solqp, dont branch....
for i in {0..2}; do 
# for i in 2; do 
   for k in 0.01 0.005; do
       for j in 0; do
           tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_pull3_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
               "selFuncIn=(solqp)" \
               "fthresh=${k}" \
               "selectTrials=${selectTrials[$i]}")
           Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
       done
   done
done



# use solqp, vary fthreshout, 
for i in {0..2}; do 
    for k in 0.005 0.01; do
        for j in 0; do
            for l in 0.05 0.1 0.2; do
                for n in 1; do
                    for p in 0 1 2; do
                        tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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


###############################
# rm below ?
#####################

# # use solqp, vary limitN
# for i in 1; do 
#     for k in 0.01; do
#         for j in 0; do
#             for l in 0.1 0.2; do
#                 for n in 0 1 2; do
#                     for p in 0; do
#                         tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#                         args=("nThreads=${nThreads}" \
#                             "nYr=${nYr}" \
#                             "lgen=$lgen" \
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


# # use solqp, pull generation
# for i in {0..2}; do 
#     for k in 0.05 0.01 0.005; do
#         for j in 2 1 0; do
#             for l in 0.1 0.2; do
#                 for n in 1; do
#                     for p in {0..2}; do
#                         tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#                         args=("nThreads=${nThreads}" \
#                             "nYr=${nYr}" \
#                             "lgen=$lgen" \
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


# # use solqp, redos
# for i in 2; do 
#     for k in 0.05; do
#         for j in 2 1 0; do
#             for l in 0.05 0.1; do
#                 for n in 1; do
#                     for p in 0; do
#                         tmpSimName="${projName}${nYr}yr${QTL}_fin${k}_fout${l}_N${n}_pull${p}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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


