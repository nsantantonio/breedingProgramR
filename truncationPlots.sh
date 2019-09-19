projName=truncation

if [ ! -d figures/${projName} ]; then
	mkdir -p figures/${projName}
else
	echo "directory exists!"
fi


for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
    'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.1_vdp${v}.RData',
    'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.2_vdp${v}.RData',
    'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp${v}.RData',
    'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.4_vdp${v}.RData',
    'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.5_vdp${v}.RData')" \
        "figDir=figures/${projName}" \
        "figName=trunc_${v}.pdf" \
        "labels=c('trad', 'rgsc0.1', 'rgsc0.2', 'rgsc0.3', 'rgsc0.4', 'rgsc0.5')"
done
# rgsc0.3 seems the best to balance short term with long term gains. 


Rscript plotbyfileStandardizeByTrad.R --args "vdp=c('10x50','20x75','40x100')" \
    "tradfilenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp')" \
    "filenames=c('results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.1_vdp',
        'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.2_vdp',
        'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp',
        'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.4_vdp',
        'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.5_vdp')" \
    "figDir=figures/${projName}" \
    "figName=truncation.pdf" \
    "labels=c('rgsc0.1', 'rgsc0.2', 'rgsc0.3', 'rgsc0.4', 'rgsc0.5')"
