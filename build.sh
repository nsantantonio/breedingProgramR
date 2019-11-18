pkg=breedingProgramR
R -e 'library(devtools);document()'
cd ../
R CMD build $pkg
R CMD INSTALL $pkg
cd $pkg