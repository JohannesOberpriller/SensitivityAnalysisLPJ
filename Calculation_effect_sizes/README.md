# Calculation of effect sizes 

In this folder, you find the functions to calculate the effect sizes <br>
To do so, we first compare the results with the results of the default parametrisation (files: comparison) and then calculate Linear regressions on the differences <br>

The following scripts are contained in the folder: <br>
**Combining_mixed_results.R**: Combines the runs of the different comparisons for mixed simulations. Required to save RAM <br>
**comparison_Fag_syl.R**: Comparison of runs with default parametrization for Fagus sylvatica <br>
**comparison_pin_syl.R**: Comparison of runs with default parametrization for Pinus sylvestris <br>
**comparison_pic_abi.R**:Comparison of runs with default parametrization for Picea abies <br>
**comparison_mixed.R**: Comparison of runs with default parametrization for mixed results <br>
**comparison_complete.R**: Sources all the different comparisons <br>
**Linear_Regressions.R**: Calculates the linear regressions and random forests for all the monoculture runs <br> 
**Linear_Regressions_mixed.R** Calculates the linear regressions and random forests for all the mixed runs <br>
