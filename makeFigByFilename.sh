# traditional
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad3_intWithin1_intAcross1_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=trad2vstrad3.pdf" \
    "labels=c('trad2', 'trad3')"


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=acrossInt.pdf" \
    "labels=c('across1', 'across0.5')"

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin0.5_intAcross0.5_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt.pdf" \
    "labels=c('within1', 'within0.5', 'within0.2')"

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt0.2.pdf" \
    "labels=c('trad2', 'trad3')"


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin0.1_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt0.2vs0.1.pdf" \
    "labels=c('trad2_0.2', 'trad2_0.1')"

# Best, Limit I
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditionalLimitI/traditionalLimitI30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=tradVsTradBestVsTradLimitI.pdf" \
    "labels=c('trad2', 'trad2_Best', 'trad2_LimitI')"



## quad prog, pull 3

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/quadprog" \
    "figName=tradVsQquadprog20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog3fin0.01')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/quadprog" \
    "figName=tradVsQuadprog10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog3fin0.01')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/quadprog" \
    "figName=trad3VsQuadprog20x75.pdf" \
    "labels=c('trad3', 'quadprog3fin0.005', 'quadprog3fin0.01')"    



# quadprog branch


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.005fout0.2.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog2fin0.005', 'quadprog1fin0.005', 'quadprog0fin0.005')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.01fout0.2.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fin0.01', 'quadprog1fin0.01', 'quadprog0fin0.01')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.005fout0.1.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog2fin0.005', 'quadprog1fin0.005', 'quadprog0fin0.005')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.01fout0.1.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fin0.01', 'quadprog1fin0.01', 'quadprog0fin0.01')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fout0.1', 'quadprog2fout0.2', 'quadprog0fout0.1', 'quadprog0fou0.2')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull2_truth0_rgsc0.2_vdp10x50.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull2_truth0_rgsc0.2_vdp10x50.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fout0.1', 'quadprog2fout0.2', 'quadprog0fout0.1', 'quadprog0fou0.2')"    


# truncation

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.1_vdp20x75.RData',
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.2_vdp20x75.RData',
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp20x75.RData')" \
    "figDir=figures/trunc" \
    "figName=trunc20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'rgsc0.1', 'rgsc0.2', 'rgsc0.3')"


# quadprog branch with RGSC pheno

#


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.005fout0.1.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog2fin0.005', 'quadprog1fin0.005', 'quadprog0fin0.005')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.3_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull0_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPull0Gen0.005fout0.1_phRS.pdf" \
    "labels=c('trad2', 'trunc0.2', 'trunc0.3', 'quadprog3fin0.005', 'quadprog0', 'quadprog0phRS0.4')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.3_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull1_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPull1Gen0.005fout0.1_phRS.pdf" \
    "labels=c('trad2', 'trunc0.2', 'trunc0.3', 'quadprog3fin0.005', 'quadprog0', 'quadprog0phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.3_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPull2Gen0.005fout0.1_phRS.pdf" \
    "labels=c('trad2', 'trunc0.2', 'trunc0.3', 'quadprog3fin0.005', 'quadprog0', 'quadprog0phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull1_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPull1Gen0.005fout0.1_phRS.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0', 'quadprog0phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPull2Gen0.005fout0.1_phRS.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0', 'quadprog0phRS0.4')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull0_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull1_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=pullGen20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog1fout0.1phRS0.4', 'quadprog2fout0.1phRS0.4')"    


## pull

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull0_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull1_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=pullGen20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog1fout0.1phRS0.4', 'quadprog2fout0.1phRS0.4')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull0_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull1_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=pullGen10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog1fout0.1phRS0.4', 'quadprog2fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp40x100.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull0_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull1_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData'
)" \
    "figDir=figures/branch" \
    "figName=pullGen40x100.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog1fout0.1phRS0.4', 'quadprog2fout0.1phRS0.4')"    

#########

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.01.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.01phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.01phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    


# 10 x 50
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.01_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005phRS0.4_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    

# 40 x 100

# something is causeing this to fail?

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.3_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005phRS0.4_40x100.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'trunc0.2phRS0.4', 'trunc0.3phRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4_40x100.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    


# truncation with RGSC pheno
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.2_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.6_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=truncPhenoRGSC20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'rgsc0.2', 'phenoRGSC0.2', 'phenoRGSC0.4', 'phenoRGSC0.6')"

# truncation with RGSC pheno
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.2_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.6_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=trunc0.3PhenoRGSC20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'rgsc0.3', 'phenoRGSC0.2', 'phenoRGSC0.4', 'phenoRGSC0.6')"


# compare VDP size


Rscript plotbyfilename.R --args "filenames=c('results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData', 
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4_40x100.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    

'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData',

########################################################################################

# traditional
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad3_intWithin1_intAcross1_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=trad2vstrad3.pdf" \
    "labels=c('trad2', 'trad3')"


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=acrossInt.pdf" \
    "labels=c('across1', 'across0.5')"

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin0.5_intAcross0.5_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt.pdf" \
    "labels=c('within1', 'within0.5', 'within0.2')"

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt0.2.pdf" \
    "labels=c('trad2', 'trad3')"


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin0.1_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt0.2vs0.1.pdf" \
    "labels=c('trad2_0.2', 'trad2_0.1')"

## quad prog, pull 3

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/quadprog" \
    "figName=tradVsQquadprog20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog3fin0.01')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/quadprog" \
    "figName=tradVsQuadprog10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog3fin0.01')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/quadprog" \
    "figName=trad3VsQuadprog20x75.pdf" \
    "labels=c('trad3', 'quadprog3fin0.005', 'quadprog3fin0.01')"    



# quadprog branch


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.005fout0.2.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog2fin0.005', 'quadprog1fin0.005', 'quadprog0fin0.005')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.01fout0.2.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fin0.01', 'quadprog1fin0.01', 'quadprog0fin0.01')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.005fout0.1.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog2fin0.005', 'quadprog1fin0.005', 'quadprog0fin0.005')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull1_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchPullGen0.01fout0.1.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fin0.01', 'quadprog1fin0.01', 'quadprog0fin0.01')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull2_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fout0.1', 'quadprog2fout0.2', 'quadprog0fout0.1', 'quadprog0fou0.2')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull2_truth0_rgsc0.2_vdp10x50.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull2_truth0_rgsc0.2_vdp10x50.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog2fout0.1', 'quadprog2fout0.2', 'quadprog0fout0.1', 'quadprog0fou0.2')"    


# truncation

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.1_vdp20x75.RData',
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.2_vdp20x75.RData',
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp20x75.RData')" \
    "figDir=figures/trunc" \
    "figName=trunc20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'rgsc0.1', 'rgsc0.2', 'rgsc0.3')"


# quadprog branch with RGSC pheno

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.01.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.01phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.01phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    


# 10 x 50
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.01_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.01_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.2_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.2', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.1phRS0.6')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005phRS0.4_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4_10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    

# 40 x 100

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.1phenoRGSCfin0.005phRS0.4_40x100.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.1', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4')"    


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp40x100.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.2_N1_pull0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.2_N1_pull2_phRS0.4_sepTrn0_truth0_rgsc0.2_vdp40x100.RData')" \
    "figDir=figures/branch" \
    "figName=branchFout0.2phenoRGSCfin0.005phRS0.4_40x100.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog0fout0.2', 'truncphRS0.4', 'quadprog0fout0.1phRS0.4', 'quadprog0fout0.2phRS0.4')"    


# truncation with RGSC pheno
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.2_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.6_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=truncPhenoRGSC20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'rgsc0.2', 'phenoRGSC0.2', 'phenoRGSC0.4', 'phenoRGSC0.6')"

# truncation with RGSC pheno
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 
'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData', 
'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.2_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.4_truth0_rgsc0.2_vdp20x75.RData',
'results/phenoRGSC/phenoRGSC30yr1000QTL_trunc_phRS0.6_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/branch" \
    "figName=trunc0.3PhenoRGSC20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.01', 'rgsc0.3', 'phenoRGSC0.2', 'phenoRGSC0.4', 'phenoRGSC0.6')"
