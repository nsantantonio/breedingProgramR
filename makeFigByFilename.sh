# test
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad3_intWithin1_intAcross1_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=trad2vstrad3.pdf" \
    "labels=c('trad2', 'trad3')"


Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=acrossInt.pdf" \
    "labels=c('across1', 'across0.5')"

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad2_intWithin0.5_intAcross0.5_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt.pdf" \
    "labels=c('within1', 'within0.5', 'within0.2')"

Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData')" \
    "figDir=figures/trad" \
    "figName=withinInt0.2.pdf" \
    "labels=c('trad2', 'trad3')"




Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData', 'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp20x75.RData', 'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp20x75.RData')" \
    "figDir=figures/tradVsQuadprog" \
    "figName=quadprog20x75.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog3fin0.01')"    



Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData', 'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData', 'results/quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData')" \
    "figDir=figures/tradVsQuadprog" \
    "figName=tradVsQuadprog10x50.pdf" \
    "labels=c('trad2', 'quadprog3fin0.005', 'quadprog3fin0.01')"    


traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData  
traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp20x75.RData
traditional30yr1000QTL_trad2_intWithin0.2_intAcross1_truth0_vdp20x75.RData    
traditional30yr1000QTL_trad3_intWithin0.2_intAcross1_truth0_vdp20x75.RData
traditional30yr1000QTL_trad2_intWithin0.5_intAcross0.5_truth0_vdp20x75.RData  
traditional30yr1000QTL_trad3_intWithin0.5_intAcross0.5_truth0_vdp20x75.RData
traditional30yr1000QTL_trad2_intWithin0.5_intAcross1_truth0_vdp20x75.RData    
traditional30yr1000QTL_trad3_intWithin0.5_intAcross1_truth0_vdp20x75.RData
traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp20x75.RData    
traditional30yr1000QTL_trad3_intWithin1_intAcross0.5_truth0_vdp20x75.RData
traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData      
traditional30yr1000QTL_trad3_intWithin1_intAcross1_truth0_vdp20x75.RData


quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp10x50.RData
quadprog/quadprog30yr1000QTL_fin0.01_pull3_truth0_rgsc0.2_vdp10x50.RData
