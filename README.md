# GSImperfectLD   by CE Rabier and S Grusea 

## R code for the manuscript "Prediction in high dimensional linear models and application to genomic selection under imperfect linkage disequilibrium"



For simulated data, it requires the R package hypred of Franck Technow. In order to install hypred, download it from https://cran.r-project.org/src/contrib/Archive/hypred/,
and then you need to recompile the source in order to make it work with last version of R. So, write in your terminal: R CMD build hypred.

For analyzing real data, you first need to download Spindel's data, available from dryad at https://doi.org/10.5061/dryad.7369p 
Next, you need to consider the folder random subsets once you have downloaded Spindel' s data
