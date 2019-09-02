#!/bin/bash

projName=phenoRGSC

nFam=(10 20 40)
famSize=(50 75 100)


reps=1
nFounderPops=100
nThreads=50
nYr=30
QTL=1000QTL
lgen=4

selectTrials=("c(0.5, 0.2, 0.5, 0.4)" \
              "c(0.5, 0.10, 0.4, 1/3)" \
              "c(0.5, 0.05, 0.4, 0.25)" )

# use solqp, TEST
# for i in 1; do 
#     for k in 0.01; do
#         for j in 0; do
#             for l in 0.2; do
#                 for n in 1; do
#                     for p in 0; do
#                         for q in 0.4; do
#                             for r in 0 1; do
#                                 tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_fout${l}_N${n}_pull${p}_phRS${q}_sepTrn${r}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#                                 args=("nThreads=${nThreads}" \
#                                     "nYr=${nYr}" \
#                                     "lgen=${lgen}" \
#                                     "simName=$tmpSimName" \
#                                     "projName=${projName}" \
#                                     "reps=${reps}" \
#                                     "nFounderPops=${nFounderPops}" \
#                                     "selectRGSC=0.2" \
#                                     "useTruth=${j}" \
#                                     "nFam=${nFam[$i]}" \
#                                     "famSize=${famSize[$i]}"\
#                                     "limitN=${n}" \
#                                     "pullCycle=${p}" \
#                                     "selFuncOut=(solqpOut)" \
#                                     "selFuncIn=(solqp)" \
#                                     "fthresh=${k}" \
#                                     "fthreshOut=${l}" \
#                                     "phenoRGSC=${q}", \
#                                     "separateTrain=${r}", \
#                                     "selectTrials=${selectTrials[$i]}")
#                                 Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#                             done
#                         done
#                     done
#                 done
#             done
#         done
#     done
# done


# branch, update model with RGSC phenotypes
for i in 1 0 2; do 
    for k in 0.005 0.01; do
        for j in 0; do
            # for l in 0.2; do # already ran with 0.2
            for l in 0.1; do # need to run with 0.1!!!!!!!!!!!!!!!!!!
                for n in 1; do
                    for p in 1; do
                        for q in 0.2 0.4 0.6; do
                            for r in 0; do
                                tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_fout${l}_N${n}_pull${p}_phRS${q}_sepTrn${r}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
                                args=("nThreads=${nThreads}" \
                                    "nYr=${nYr}" \
                                    "lgen=${lgen}" \
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
                                    "phenoRGSC=${q}", \
                                    "separateTrain=${r}", \
                                    "selectTrials=${selectTrials[$i]}")
                                Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
                            done
                        done
                    done
                done
            done
        done
    done
done


# truncation selection update with RGSC phenotypes
for i in 1 0 2; do 
	for j in 0; do
        for r in 0.3; do 
        	for q in 0.2 0.4 0.6; do
	            tmpSimName="${projName}${nYr}yr${QTL}_trunc_phRS${q}_truth${j}_rgsc${r}_vdp${nFam[$i]}x${famSize[$i]}"
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
	                "separateTrain=0", \
	                "selectTrials=${selectTrials[$i]}")
	            Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
        	done
        done
    done
done




# use solqp, dont branch, update with RGSC phenotypes
for i in 1 0 2; do 
   for k in 0.005 0.01; do
       for j in 0; do
        	for q in 0.2 0.4 0.6; do
	           tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_pull3_phRS${q}_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
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
	               "separateTrain=0", \
	               "selectTrials=${selectTrials[$i]}")
	           Rscript $(pwd)/singleTraitProgram.R "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
            done
        done
    done
done



# # use solqp, vary limitN
# for i in 1; do 
#     for k in 0.01; do
#         for j in 2 1 0; do
#             for l in 0.1 0.2; do
#                 for n in 0 1 2; do
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
#                         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#                     done
#                 done
#             done
#         done
#     done
# done



# # use solqp, pull generation
# for i in {0..2}; do 
#     for k in 0.05 0.005; do
#         for j in 2 1 0; do
#             for l in 0.1 0.2; do
#                 for n in 1; do
#                     for p in {0..2}; do
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
#                         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#                     done
#                 done
#             done
#         done
#     done
# done


# # use solqp, dont branch....
# for i in {0..2}; do 
#     for k in 0.05 0.01 0.005; do
#         for j in 2 0 1; do
#             tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin${k}_pull3_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#             args=("nThreads=${nThreads}" \
#                 "nYr=${nYr}" \
#                 "simName=$tmpSimName" \
#                 "projName=${projName}" \
#                 "reps=${reps}" \
#                 "nFounderPops=${nFounderPops}" \
#                 "selectRGSC=0.2" \
#                 "useTruth=${j}" \
#                 "nFam=${nFam[$i]}" \
#                 "famSize=${famSize[$i]}"\
#                 "selFuncIn=(solqp)" \
#                 "fthresh=${k}" \
#                 "selectTrials=${selectTrials[$i]}")
#             Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#         done
#     done
# done


# # truncation selection

# for i in {0..2}; do 
# # for i in 1; do 
#     for j in 2 1 0; do
#         tmpSimName="${projName}${nYr}yr${QTL}_trunc_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#         args=("nThreads=${nThreads}" \
#             "nYr=${nYr}" \
#             "simName=$tmpSimName" \
#             "projName=${projName}" \
#             "reps=${reps}" \
#             "nFounderPops=${nFounderPops}" \
#             "useTruth=${j}" \
#             "nFam=${nFam[$i]}" \
#             "famSize=${famSize[$i]}"\
#             "selectTrials=${selectTrials[$i]}")
#         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#      done
# done




# # truncation selection
# for i in 1; do 
#     for j in 2 1 0; do
#         tmpSimName="${projName}${nYr}yr${QTL}_quadprog_fin_fout_N_pull_phRS_sepTrn_truth${j}_rgsc0.2_vdp${nFam[$i]}x${famSize[$i]}"
#         args=("nThreads=${nThreads}" \
#             "nYr=${nYr}" \
#             "simName=$tmpSimName" \
#             "projName=${projName}" \
#             "reps=${reps}" \
#             "nFounderPops=${nFounderPops}" \
#             "selectRGSC=0.2" \
#             "useTruth=${j}" \
#             "nFam=${nFam[$i]}" \
#             "famSize=${famSize[$i]}"\
#             "selectTrials=${selectTrials[$i]}")
#         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#     done
# done






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
#                         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt 
#                     done
#                 done
#             done
#         done
#     done
# done


