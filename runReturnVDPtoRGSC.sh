#! /bin/sh
pardir=$(pwd) # just in case $PWD doesnt exist on your system

# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings
    
VDPtoRGSC=("c(0.1, 0, 0, 0, 0, 0)" \
           "c(0, 0.1, 0, 0, 0, 0)" \
           "c(0, 0, 0.2, 0, 0, 0)" \
           "c(0, 0, 0, 0.50, 0, 0)" \
           "c(0, 0, 0, 0, 1, 0)" \
           "c(0, 0, 0, 0, 0, 1)" 
)

VDPnames=("trial1_0.1" \
           "trial2_0.1" \
           "trial3_0.2" \
           "trial4_0.5" \
           "trial5_1" \
           "variety_1"
)

for i in {0..5}; do
    args=("simName=returnVDPtoRGSC_${VDPnames[$i]}" \
        "returnVDPtoRGSC=${VDPtoRGSC[$i]}"
    )
    # echo ${args[1]} ${args[2]}
    Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" &
done
