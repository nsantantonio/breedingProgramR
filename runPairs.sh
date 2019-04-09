pardir=$(pwd) # just in case $PWD doesnt exist on your system 

args=("seed=12345" \
	"nThreads=20" \
	"simName='simNoTitle'" \
	"simFunc='fitFuncSingleVariate.R'" \
	"nGen=20" \
	"reps=10" \
	"lgen=5" \
# founder parameters
	"founderRData='founderPop/testAlphaSimR1000SegSite.RData'" \
	"founderh2=0.3" \
	"simpleFounder=FALSE" \
# selection parameters
	"RGSCintensity=0.2" \
	"selectInRGSC='ebv'" \
	"selectOutRGSC='ebv'" \
	"selectVDP='pheno'" \
	"returnVDPcrit='pheno'" \
	"selFunc=getExpDist" \
	"withinFamInt=1" \
	"useQuantile=0.9" \
	"skip=NULL" \
# family parameters
	"nFounder=10" \
	"nNuclear=50" \
	"nFam=10" \
	"famSize=10" \
	"ssd=FALSE" \
	"selF2=FALSE" \
	"nF2=10" \
# genetic parameters
	# kinship="SNP" \
	"Vg=1" \
	"updateVg=FALSE" \
	"h2=c(0.1,0.3,0.3,0.3,0.3)" \
	"nYr=20" \
	"selectTrials=c(0.5,0.5,0.4,0.3,1/3)" \
	"trialReps=c(1,1,2,3,3)" \
	"trialLocs=c(1,1,2,5,5)" \
	"GScylcePerYr=2" \
	"returnVDPtoRGSC=c(0,0,0,0,0,0)" \
# chromosome parameters
	"nChrom=10" \
	"nLoci=100" \
	"nM=100" \
	"nQTL=100" \
	)


# for i in ${args[@]}; do
#   echo $i
# done

Rscript ${pardir}/sparseSkipVDPsimple.R --args "${args[@]}"