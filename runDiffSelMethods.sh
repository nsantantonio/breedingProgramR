#! /bin/bash
pardir=$(pwd) # just in case $PWD doesnt exist on your system

args=("seed=12345" \
    "simName=default")
Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/default.txt 


# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings
for i in truncSel expDist simDHdist; do
	for j in truncCross expDistPairs simDHdistPairs maxVar; do 
        args=("seed=12345" \
            "simName=${i}_${j}" \
            "selFuncOut=(${i})" \
            "selFuncIn=(${j})")
        # echo ${args[1]} ${args[2]}
        Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_${i}_${j}.txt &
    done
done

for i in truncSel expDist simDHdist; do
    for j in expDistPairs simDHdistPairs maxVar; do 
        args=("seed=12345" \
            "simName=${i}_${j}_noTrunc" \
            "selFuncOut=(${i})" \
            "selFuncIn=(${j})"\
            "selectRGSC=1")
        # echo ${args[1]} ${args[2]}
        Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_${i}_${j}.txt &
    done
done

# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings
# for i in truncSel expDist simDHdist; do
# 	for j in truncCross expDistPairs simDHdistPairs maxVar; do 
#         args=("seed=12345" \
#             "simName=${i}_${j}_useTRUE" \
#             "selFuncOut=(${i})" \
#             "selFuncIn=(${j})" \
#             "useTrue=TRUE")
#         # echo ${args[1]} ${args[2]}
#         Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_${i}_${j}_useTRUE.txt
#     done
# done

# for i in truncSel expDist simDHdist; do
# 	for j in truncCross expDistPairs simDHdistPairs; do 
#         args=("seed=12345" \
#             "simName=${i}_${j}_w0.25" \
#             "selFuncOut=(${i})" \
#             "selFuncIn=(${j})" \
#             "weight=0.25"
#         )
#         # echo ${args[1]} ${args[2]}
#         Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_${i}_${j}_w0.25.txt
#     done
# done

# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings
# for i in truncSel expDist simDHdist; do
#     for j in truncCross expDistPairs simDHdistPairs; do 
#         args=("seed=12345" \
#             "simName=${i}_${j}_famSize20" \
#             "selFuncOut=(${i})" \
#             "selFuncIn=(${j})" \
#             "famSize=20"
#         )
#         # echo ${args[1]} ${args[2]}
#         Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_${i}_${j}_famSize20.txt
#     done
# done



# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings
# for i in truncSel expDist; do
#     for j in truncCross expDistPairs; do 
#         args=("seed=12345" \
#             "simName=${i}_${j}_selectRGSC0.1" \
#             "selFuncOut=(${i})" \
#             "selFuncIn=(${j})" \
#             "selectRGSC=0.1" \
#         )
#         # echo ${args[1]} ${args[2]}
#         Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_${i}_${j}.txt &
#     done
# done


# args=("seed=12345" \
#     "simName=truncSel_maxVar" \
#     "selFuncOut=(truncSel)" \
#     "selFuncIn=(maxVar)" 
# )
# Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}" > ${pardir}/logs/log_truncSel_maxVar.txt 
