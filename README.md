# Leakiness-sex-expression-and-reproductive-effort-Mercurialis-annua-range-expansion
Supporting data for the manuscript "Geographical gradients in leaky sex expression and reproductive effort in a dioecious plant are consistent with selection during range expansion"

DATA files:
1. GPS26withdistancesMINFERRY.cvs GPS locations of the sites and the calculated distance metrics.
2. Seed phenotype.xlsx Seed weight measured as bulk of each populations, and measured of pooled 50 seeds per population.
3. phenobiomassfull.xlsx phenotyping data of all plants.
4. seeds pictures folder contains pictures and results from ImageJ for seed sizes of each population.

ANALYSIS files:
1. FunNiche-climateOFFICIAL.R Analysis of climate data downloaded from CHELSA. Generate PC1 and PC2 scores for all sites.
2. PrepareData.R Take input from the data files and climate file generated after running FunNiche-climateOFFICIAL.R, clean and prepare dataframes for analysis.
3. fulldataanalysis08OFFICIAL.R Performed GLMMs, variation partitioning, other stat test on all the prepared dataframes generated after running previous two analysis files.
