projName=truncation

if [ ! -d figures/${projName} ]; then
	mkdir -p figures/${projName}
else
	echo "directory exists!"
fi

for v in 10x50 20x75 40x100; do
    for r in 0.3 0.2; do
        Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
        'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc${r}_vdp${v}.RData',
        'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.2_truth0_rgsc${r}_vdp${v}.RData',
        'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc${r}_vdp${v}.RData',
        'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.6_truth0_rgsc${r}_vdp${v}.RData')" \
            "figDir=figures/${projName}" \
            "figName=trunc_rgsc${r}_phRS_${v}.pdf" \
            "labels=c('trad', 'rgsc${r}', 'rgsc${r}_phRS0.2', 'rgsc${r}_phRS0.4', 'rgsc${r}_phRS0.6')"
    done
done
# phRS0.3 seems 

for r in 0.3 0.2; do
    Rscript plotbyfileStandardizeByTrad.R --args "vdp=c('10x50','20x75','40x100')" \
        "tradfilenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp')" \
        "filenames=c('results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc${r}_vdp',
        'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.2_truth0_rgsc${r}_vdp',
            'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc${r}_vdp',
            'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.6_truth0_rgsc${r}_vdp')" \
        "figDir=figures/${projName}" \
        "figName=truncationPhRS_rgsc${r}.pdf" \
        "labels=c( 'rgsc${r}', 'rgsc${r}_phRS0.2', 'rgsc${r}_phRS0.4', 'rgsc${r}_phRS0.6')"
done