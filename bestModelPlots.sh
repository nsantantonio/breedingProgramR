projName=bestModels

if [ ! -d figures/${projName} ]; then
    mkdir -p figures/${projName}
else
    echo "directory exists!"
fi


# by phRS
for p in 0 ; do
    for v in 10x50 20x75 40x100; do
        for fin in 0.005 0.01; do
            for fout in 0.1; do
                for ph in 0.4; do
                    Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
                        'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp${v}.RData',
                        'results/quadprog/quadprog30yr1000QTL_fin${fin}_pull3_truth0_rgsc0.2_vdp${v}.RData', 
                        'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout${fout}_N1_pull${p}_truth0_rgsc0.2_vdp${v}.RData',
                        'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin${fin}_fout${fout}_N1_pull${p}_phRS${ph}_sepTrn0_truth0_rgsc0.2_vdp${v}.RData')" \
                        "figDir=figures/${projName}" \
                        "figName=quadprogBranchPhenoRGSC_fin${fin}_pull${p}_byFout_${v}.pdf" \
                        "labels=c('traditional', 'truncation', 'optimalContribution', 'optContBranch', 'optContBranch_PhenoRS${p}')"
                done
            done
        done
    done
done


Rscript plotbyfileStandardizeByTrad.R --args "vdp=c('10x50','20x75','40x100')" \
    "tradfilenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp')" \
    "filenames=c( 'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp',
                        'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp', 
                        'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp',
                        'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull0_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp')" \
    "figDir=figures/${projName}" \
    "figName=best.pdf" \
    "labels=c('truncation', 'optimalContribution', 'optContBranch', 'optContBranch_PhenoRS')"
