
for j in traditionalBest truncation quadprog phenoRGSC; do 
    resDir=results/${j}/ 
    sumDir=resultSummary/${j}/ 
    fnames=$(diff ${resDir} ${sumDir}| grep "^Only in ${resDir}"| sed 's@^.*: @@')

    for i in ${fnames[@]}; do
    echo "results/${j}/${i}"
    Rscript processResults.R Rscript plotbyfilename.R --args "filenames=${resDir}${i}"
    done
done
