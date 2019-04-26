 #! /bin/sh
pardir=$(pwd) # just in case $PWD doesnt exist on your system

# Note: functions, NULL or anything else that needs to be evaluated must be put in parentheses so they are evaluated properly, otherwise they will be treated as strings

args=("seed=12345" \
    "nThreads=20" \
    "simName='tradtional'" \
    "simFunc='fitFuncSingleVariate.R'" \
    "nYr=20" \
    "reps=10" \
    "lgen=5" \
    "useTrue=FALSE" \
    "traditional=FALSE" \
# founder parameters
    "founderRData='founderPop/testAlphaSimR1000SegSite.RData'" \
    "founderh2=0.3" \
    "simpleFounder=FALSE" \
    "founderBurnIn=3" \
# selection parameters
    "selectRGSC=0.2" \
    "RGSCprogenyPerCross=1" \
    "selectIn='ebv'" \
    "selectOut='ebv'" \
    "selectVDP='pheno'" \
    "returnVDPcrit='pheno'" \
    "selFuncOut=(truncSel)" \
    "selFuncIn=(truncCross)" \
    "withinFamInt=1" \
    "setXint=(NULL)" \
    "skip=(NULL)" \
# family parameters
    "nFounder=100" \
    "nNuclear=100" \
    "nFam=10" \
    "famSize=50" \
    "ssd=FALSE" \
    "selF2=FALSE" \
    "nF2=1" \
# genetic parameters
    "Vg=1" \
    "updateVg=FALSE" \
    "h2=c(0.1, 0.3, 0.3, 0.3, 0.3)" \
    "selectTrials=c(0.5, 0.5, 0.4, 0.3, 1/3)" \
    "trialReps=c(1, 1, 2, 3, 3)" \
    "trialLocs=c(1, 1, 2, 5, 5)" \
    "cyclePerYr=2" \
    "returnVDPtoRGSC=c(0, 0, 0, 0, 0, 0)" \
# chromosome parameters
    "nChrom=10" \
    "nLoci=100" \
    "nM=100" \
    "nQTL=100"
    )

Rscript ${pardir}/singleTraitProgram.R --args "${args[@]}"