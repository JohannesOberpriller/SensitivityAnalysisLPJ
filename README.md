# SensitivityAnalysisLPJ

This repository contains the Code used to prepare the analysis of uncertainty and sensitivity of LPJ-GUESS for the manuscript "Sensitivity and uncertainty analysis of the LPJ – GUESS vegetation model to climate and parameter uncertainty across European forests" published in XXXX 

The repository is structured, such that the complete case study can be reanalysed.
To do so the following has to be done (course description, more description in the individual folders for each steps):
1. Get the development version of the [rLPJGUESS package] (https://github.com/biometry/rLPJGUESS/blob/develop/README.md) 
2. Get the climate data from the IPSL-CM5 Earth System Model (https://link.springer.com/article/10.1007/s00382-012-1636-1)
3. Create all files required to run the model using the scripts in ParameterMetaData 
4. Run the created R scripts on a cluster (shell scripts for submitting the data are provided, else results for runs are provided here (https://zenodo.org/record/4670295#.YKIkI-tCRqs))
5. Analyse the data for mean effects using the files in Calculation_effect_sizes
6. Make the plots using the scripts in Plotting 
