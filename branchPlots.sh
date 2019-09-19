projName=optContBranch

if [ ! -d figures/${projName} ]; then
    mkdir -p figures/${projName}
else
    echo "directory exists!"
fi


for fin in 0.005 0.01; do
Rscript plotbyfileStandardizeByTrad.R --args "vdp=c('10x50','20x75','40x100')" \
    "tradfilenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp')" \
    "filenames=c('results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp',
    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull2_truth0_rgsc0.2_vdp',
    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull1_truth0_rgsc0.2_vdp',
    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull0_truth0_rgsc0.2_vdp')" \
    "figDir=figures/${projName}" \
    "figName=branch_fin${fin}_byPull.pdf" \
    "labels=c('fin${fin}_pull3', 'fin${fin}_pull2', 'fin${fin}_pull1', 'fin${fin}_pull0')"
done



# by fout
for p in 0 1 2; do
    for v in 10x50 20x75 40x100; do
        for fin in 0.005 0.01; do
            Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
                'results/quadprog/quadprog30yr1000QTL_fin${fin}_pull3_truth0_rgsc0.2_vdp${v}.RData', 
                'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.05_N1_pull${p}_truth0_rgsc0.2_vdp${v}.RData',
                'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull${p}_truth0_rgsc0.2_vdp${v}.RData',
                'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.2_N1_pull${p}_truth0_rgsc0.2_vdp${v}.RData')" \
                "figDir=figures/${projName}" \
                "figName=quadprogBranch_fin${fin}_pull${p}_byFout_${v}.pdf" \
                "labels=c('trad', 'fin${fin}_pull3', 'fout0.05_pull${p}', 'fout0.1_pull${p}', 'fout0.2_pull${p}')"
        done
    done
done


# by fin
for p in 0 1 2; do
    for v in 10x50 20x75 40x100; do
        for fout in 0.05 0.1 0.2; do
            Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
                'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp${v}.RData', 
                'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp${v}.RData', 
                'results/quadprog/quadprog30yr1000QTL_fin0.005_fout${fout}_N1_pull${p}_truth0_rgsc0.2_vdp${v}.RData',
                'results/quadprog/quadprog30yr1000QTL_fin0.01_fout${fout}_N1_pull${p}_truth0_rgsc0.2_vdp${v}.RData')" \
                "figDir=figures/${projName}" \
                "figName=quadprogBranch_fout${fout}_pull${p}_byFin_${v}.pdf" \
                "labels=c('trad', 'fin0.005_pull3', 'fin0.01_pull3', 'fin0.005_pull${p}', 'fin0.01_pull${p}')"
        done
    done
done



# by pull
for v in 10x50 20x75 40x100; do
    for fin in 0.005 0.01; do
        for fout in 0.05 0.1 0.2; do
                Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
                    'results/quadprog/quadprog30yr1000QTL_fin${fin}_pull3_truth0_rgsc0.2_vdp${v}.RData', 
                    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout${fout}_N1_pull2_truth0_rgsc0.2_vdp${v}.RData',
                    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout${fout}_N1_pull1_truth0_rgsc0.2_vdp${v}.RData',
                    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout${fout}_N1_pull0_truth0_rgsc0.2_vdp${v}.RData')" \
                    "figDir=figures/${projName}" \
                    "figName=quadprogBranch_fin${fin}_fout${fout}_byPull_${v}.pdf" \
                    "labels=c('trad', 'pull3', 'pull2', 'pull1', 'pull0')"
            done
        done
    done
done


for fin in 0.005 0.01; do
Rscript plotbyfileStandardizeByTrad.R --args "vdp=c('10x50','20x75','40x100')" \
    "tradfilenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp')" \
    "filenames=c('results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp',
    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull2_truth0_rgsc0.2_vdp',
    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull1_truth0_rgsc0.2_vdp',
    'results/quadprog/quadprog30yr1000QTL_fin${fin}_fout0.1_N1_pull0_truth0_rgsc0.2_vdp')" \
    "figDir=figures/${projName}" \
    "figName=branch_fin${fin}_byPull.pdf" \
    "labels=c('fin${fin}_pull3', 'fin${fin}_pull2', 'fin${fin}_pull1', 'fin${fin}_pull0')"
done