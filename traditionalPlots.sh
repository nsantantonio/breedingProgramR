projName=traditional

if [ ! -d figures/${projName} ]; then
	mkdir -p figures/${projName}
else
	echo "directory exists!"
fi

for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData', 
        'results/traditional/traditional30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=trad2_trad3_${v}.pdf" \
    "labels=c('trad2', 'trad3')"
done
# trad 2 is clearly better

for v in 10x50 20x75 40x100; do
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp${v}.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=acrossInt_${v}.pdf" \
    "labels=c('across1', 'across0.5')"
done
# across 0.5 is clearly better for med, large VDP


for v in 10x50 20x75 40x100; do
Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp${v}.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin0.5_intAcross0.5_truth0_vdp${v}.RData', 
'results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=withinInt_${v}.pdf" \
    "labels=c('within1', 'within0.5', 'within0.2')"
done
# Within 0.2 is clearly better

# only ran for 20x75
for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData', 
        'results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData', 
        'results/traditionalLimitI/traditionalLimitI30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=trad2VsTradBestVsTradLimitI_${v}.pdf" \
    "labels=c('trad2', 'trad2_Best', 'trad2_LimitI')"
done
# Best is clearly better for the traditional program. 



# just check that Best is arbitrary for acrossInt
for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross1_truth0_vdp${v}.RData', 
        'results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=tradBest_acrossInt_${v}.pdf" \
    "labels=c('across1', 'across0.5')"
done


# just check that Best is arbitrary for withinInt
for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin1_intAcross0.5_truth0_vdp${v}.RData', 
        'results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.5_intAcross0.5_truth0_vdp${v}.RData', 
        'results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=tradBest_withinInt_${v}.pdf" \
    "labels=c('within1', 'within0.5', 'within0.2')"
done
# not sure I understand why withinInt = 0.2 is better... Selects other within generation families if only one. 

for v in 10x50 20x75 40x100; do
    Rscript plotbyfilename.R --args "filenames=c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData', 
        'results/traditionalBest/traditionalBest30yr1000QTL_trad3_intWithin0.2_intAcross0.5_truth0_vdp${v}.RData')" \
    "figDir=figures/${projName}" \
    "figName=trad2Best_trad3Best_${v}.pdf" \
    "labels=c('trad2', 'trad3')"
done