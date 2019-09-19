projName=optCont

if [ ! -d figures/${projName} ]; then
	mkdir -p figures/${projName}
else
	echo "directory exists!"
fi


# remake the 40x100 plots! They had failed runs. Rerunning today (Sept 9 2019).  
for v in 10x50 20x75 40x100; do
    for fin in 0.005 0.01; do
        Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData',
                'results/quadprog/quadprog30yr1000QTL_fin${fin}_pull3_truth0_rgsc0.2_vdp${v}.RData', 
                'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin${fin}_pull3_phRS0.2_truth0_rgsc0.2_vdp${v}.RData',
                'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin${fin}_pull3_phRS0.4_truth0_rgsc0.2_vdp${v}.RData',
                'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin${fin}_pull3_phRS0.6_truth0_rgsc0.2_vdp${v}.RData')" \
                "figDir=figures/${projName}" \
                "figName=quadprog_fin${fin}_phRS_${v}.pdf" \
                "labels=c('trad', 'fin${fin}', 'fin${fin}_phRS0.2', 'fin${fin}_phRS0.4', 'fin${fin}_phRS0.6')"
            done
        done
    done
done
