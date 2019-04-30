#! /bin/bash
# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings    
h2=("c(0.1)" \
    "c(0.1, 0.3)" \
    "c(0.1, 0.3, 0.3)" \
    "c(0.1, 0.3, 0.3, 0.3)" \
    "c(0.1, 0.3, 0.3, 0.3, 0.3)")

selectTrials=("c(0.01)" \
              "c(0.5, 0.02)" \
              "c(0.5, 0.5, 0.04)" \
              "c(0.5, 0.5, 0.4, 0.1)" \
              "c(0.5, 0.5, 0.4, 0.3, 1/3)")

trialReps=("c(1)" \
           "c(1, 1)" \
           "c(1, 1, 2)" \
           "c(1, 1, 2, 3)" \
           "c(1, 1, 2, 3, 3)")

trialLocs=("c(1)" \
           "c(1, 1)" \
           "c(1, 1, 2)" \
           "c(1, 1, 2, 5)" \
           "c(1, 1, 2, 5, 5)")

returnVDPtoRGSC=("c(0, 0)" \
                 "c(0, 0, 0)" \
                 "c(0, 0, 0, 0)" \
                 "c(0, 0, 0, 0, 0)" \
                 "c(0, 0, 0, 0, 0, 0)")

VDPnames=("trial1_0.01" \
           "trial2_0.02" \
           "trial3_0.04" \
           "trial4_0.1" \
           "trial5_0.33")

for i in {0..4}; do
	for j in truncCross expDistPairs maxVar; do 
        args=("nThreads=2" \
            "simName=VDPintensity_${VDPnames[$i]}_${j}" \
            "selFuncIn=(${j})" \
            "h2=${h2[$i]}" \
            "selectTrials=${selectTrials[$i]}" \
            "trialReps=${trialReps[$i]}" \
            "trialLocs=${trialLocs[$i]}" \
            "returnVDPtoRGSC=${returnVDPtoRGSC[$i]}" \
            )
        # echo "${args[@]}"
        Rscript $(pwd)/singleTraitProgram.R --args "${args[@]}" > $(pwd)/logs/log_VDPintensity_${VDPnames[$i]}_${j}.txt &
    done
done
