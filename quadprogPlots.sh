projName=optCont

if [ ! -d figures/${projName} ]; then
	mkdir -p figures/${projName}
else
	echo "directory exists!"
fi


for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
        'results/quadprog/quadprog30yr1000QTL_fin0.001_pull3_truth0_rgsc0.2_vdp${v}.RData', 
        'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp${v}.RData', 
        'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp${v}.RData',
        'results/quadprog/quadprog30yr1000QTL_fin0.02_pull3_truth0_rgsc0.2_vdp${v}.RData',
        'results/quadprog/quadprog30yr1000QTL_fin0.05_pull3_truth0_rgsc0.2_vdp${v}.RData')" \
        "figDir=figures/${projName}" \
        "figName=quadprogPull3_${v}.pdf" \
        "labels=c('trad', 'fin0.001', 'fin0.005', 'fin0.01', 'fin0.02', 'fin0.05')"
done
# fin of 0.005 to 0.01 seems best to balance short term with long term gains. 

Rscript plotbyfileStandardizeByTrad.R --args "vdp=c('10x50','20x75','40x100')" \
    "tradfilenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp')" \
    "filenames=c('results/quadprog/quadprog30yr1000QTL_fin0.001_pull3_truth0_rgsc0.2_vdp', 
        'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp', 
        'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp',
        'results/quadprog/quadprog30yr1000QTL_fin0.02_pull3_truth0_rgsc0.2_vdp',
        'results/quadprog/quadprog30yr1000QTL_fin0.05_pull3_truth0_rgsc0.2_vdp')" \
    "figDir=figures/${projName}" \
    "figName=quadprog.pdf" \
    "labels=c('fin0.001', 'fin0.005', 'fin0.01', 'fin0.02', 'fin0.05')"
