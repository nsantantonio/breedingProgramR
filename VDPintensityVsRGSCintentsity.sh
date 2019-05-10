#! /bin/bash
rgscint=(0.05 0.1 0.2 0.3 0.4 0.5)
nFam=(4 10 20 40)
famSize=(25 50 100 200)

# 100, 
selectTrials=("c(1, 0.5, 0.5, 0.4, 0.5)" \
              "c(0.5, 0.5, 0.4, 0.3, 1/3)" \
              "c(0.5, 0.25, 0.2, 0.4, 0.25)" \
              "c(0.25, 0.25, 0.25, 0.2, 0.2)" )

h2=("c(0.1, 0.3, 0.3, 0.3, 0.3)" \
    "c(0.1, 0.3, 0.3, 0.3, 0.3)" \
    "c(0.1, 0.3, 0.3, 0.3, 0.3)" \
    "c(0.1, 0.3, 0.3, 0.3, 0.3)" )

trialReps=("c(1, 1, 2, 3, 3)" \
           "c(1, 1, 2, 3, 3)" \
           "c(1, 1, 2, 3, 3)" \
           "c(1, 1, 2, 3, 3)")

trialLocs=("c(1, 1, 2, 5, 5)" \
           "c(1, 1, 2, 5, 5)" \
           "c(1, 1, 2, 5, 5)" \
           "c(1, 1, 2, 5, 5)")

returnVDPtoRGSC=("c(0, 0, 0, 0, 0, 0)" \
                 "c(0, 0, 0, 0, 0, 0)" \
                 "c(0, 0, 0, 0, 0, 0)" \
                 "c(0, 0, 0, 0, 0, 0)")


for i in {0..3}; do 
    for j in ${rgscint[@]}; do
        tmpSimName="VDPvsRGSCintensity_rgsc${j}_vdp${nFam[$i]}x${famSize[$i]}"
        args=("nThreads=2" \
            "simName=$tmpSimName" \
            "selectRGSC=${j}" \
            "nFam=${nFam[$i]}" \
            "famSize=${famSize[$i]}"\
            "h2=${h2[$i]}" \
            "selectTrials=${selectTrials[$i]}" \
            "trialReps=${trialReps[$i]}" \
            "trialLocs=${trialLocs[$i]}" \
            "returnVDPtoRGSC=${returnVDPtoRGSC[$i]}")
        # Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt &
    done
        tmpSimName="VDPvsRGSCintensity_trad_vdp${nFam[$i]}x${famSize[$i]}"
        args[1]="simName=$tmpSimName"
        args+=("traditional=TRUE")
        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt &
done




for i in {0..3}; do 
    for j in ${rgscint[@]}; do
        tmpSimName="VDPvsRGSCintensityACquant_rgsc${j}_vdp${nFam[$i]}x${famSize[$i]}"
        args=("nThreads=2" \
            "simName=$tmpSimName" \
            "selectRGSC=${j}" \
            "selFuncIn=ACquant" \
            "nFam=${nFam[$i]}" \
            "famSize=${famSize[$i]}"\
            "h2=${h2[$i]}" \
            "selectTrials=${selectTrials[$i]}" \
            "trialReps=${trialReps[$i]}" \
            "trialLocs=${trialLocs[$i]}" \
            "returnVDPtoRGSC=${returnVDPtoRGSC[$i]}")
        # Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt &
    done
        tmpSimName="VDPvsRGSCintensityACquant_trad_vdp${nFam[$i]}x${famSize[$i]}"
        args[1]="simName=$tmpSimName"
        args+=("traditional=TRUE")
        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt &
done


# fix indexing!!!!!!!!!!!!!!

# for i in ${rgscint[@]}; do
#   for j in {0..3}; do 
#         tmpSimName="VDPvsRGSCintensity100QTL_rgsc${i}_vdp${nFam[$j]}x${famSize[$j]}"
#         args=("nThreads=2" \
#             "simName=$tmpSimName" \
#             "selectRGSC=${i}" \
#             "nFam=${nFam[$j]}" \
#             "famSize=${famSize[$j]}"\
#             "h2=${h2[$j]}" \
#             "selectTrials=${selectTrials[$j]}" \
#             "trialReps=${trialReps[$j]}" \
#             "trialLocs=${trialLocs[$j]}" \
#             "returnVDPtoRGSC=${returnVDPtoRGSC[$j]}" \
#             "nQTL=10")
#         # echo ${tmpSimName}
#         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt &
#     done
# done

# for i in ${rgscint[@]}; do
#   for j in {0..3}; do 
#         tmpSimName="VDPvsRGSCintensity10QTL_rgsc${i}_vdp${nFam[$j]}x${famSize[$j]}"
#         args=("nThreads=2" \
#             "simName=$tmpSimName" \
#             "selectRGSC=${i}" \
#             "nFam=${nFam[$j]}" \
#             "famSize=${famSize[$j]}"\
#             "h2=${h2[$j]}" \
#             "selectTrials=${selectTrials[$j]}" \
#             "trialReps=${trialReps[$j]}" \
#             "trialLocs=${trialLocs[$j]}" \
#             "returnVDPtoRGSC=${returnVDPtoRGSC[$j]}" \
#             "nQTL=1")
#         # echo ${tmpSimName}
#         Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" &> $(pwd)/logs/log_${tmpSimName}.txt &
#     done
# done
