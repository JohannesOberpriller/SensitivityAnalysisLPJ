library(rLPJGUESS)
library(zoo)
print(getwd())
#define the ouputs which have been investigated before
typeList <- c("cmass", "lai","agpp","cpool","anpp","cflux","fpc","speciesdiam","dens")
scaleLPJ_PFT <- "normal"
#set the main directory to LPJrunTest
mainDir <- file.path(getwd(), "LPJrunTest")

#read in the data
sites <- readRDS("./EnvironmentalData/sites_data.rds")

results <- try(readRDS(paste0("./Runs_results/results_Pic_abi.rds")))
# for each site we have to run the model with its default configuration
# then we have to compare the difference between the default and the sampled runs
for(i in 1:nrow(sites)){
  # get default parameters from script
  defaultparameters <- InferParameterAndDesignList(list(main = paste0(mainDir,"/main_temp_comparison.ins")),
                                                   NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                   vectorvaluedparams = c("rootdist","eps_mon",
                                                                          "storfrac_mon","photo",
                                                                          "fertdates","fertrate"))

  # get the correct gridlist filename and parameter filename
  gridList_name = paste0(sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".txt")
  parameter_name = paste0("Pic_abi_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds")

  # make the templates runable
  AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                  defaultlist = defaultparameters$defaultlist,
                  MainTemplateDestination = "./LPJrunTest/main_new.ins",
                  PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                  NameMainFile = "main.ins", NamePftFile = "pft.ins")

  # produce runable parameters
  parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Pic_abi"))

  LPJsetup <- mainDir
  #define setup of the run
  LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                          file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                          file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                          file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                          file.temp = file.path(mainDir,"Data", "tasAdjust_rcp_45_1860_2100_mon_mean_bias_corrected.nc"),
                          file.prec = file.path(mainDir,"Data", "prAdjust_rcp_45_1860_2100_mon_mean_bias_corrected.nc"),
                          file.insol = file.path(mainDir,"Data", "rsdsAdjust_rcp_45_1860_2100_mon_mean_bias_corrected.nc"),
                          file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                          file.wetdays = file.path(mainDir,"Data","wetdays_rcp_45_1860_2100_mon_sum_bias_corrected6.nc"),
                          template1 = "pft_new.ins",
                          template2 = "main_new.ins",
                          variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                          variable.ndep = "RCP45",
                          variable.prec = "prAdjust", delete = F, save = F, processing = T,
                          plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                          gridList = gridList_name, parallel = "parameters",
                          defaultlist = parameters$defaultfiles)



  LPJsettings_PFT$design <- parameters$design
  # run model with default
  results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                        typeList = typeList, settings = LPJsettings_PFT)




  #load the data from the sensitivity runs


  #check if it actually there
  if("try-error" %in% class(results)){
    next
  }
  #if not there go to next site, else calculate difference
  else{

    deviations_anpp_climate_change = matrix(nrow = 2,ncol = length(results))
    deviations_anpp_steady_climate =  matrix(nrow = 2,ncol = length(results))
    deviations_anpp_complete =  matrix(nrow = 2,ncol = length(results))
    deviations_agpp_climate_change =  matrix(nrow = 2,ncol = length(results))
    deviations_agpp_steady_climate = matrix(nrow = 2,ncol = length(results))
    deviations_agpp_complete =  matrix(nrow = 2,ncol = length(results))
    deviations_cmass_climate_change =  matrix(nrow = 2,ncol = length(results))
    deviations_cmass_steady_climate =  matrix(nrow = 2,ncol = length(results))
    deviations_cmass_complete =  matrix(nrow = 2,ncol = length(results))
    deviations_cpool_climate_change =  matrix(nrow = 2,ncol = length(results))
    deviations_cpool_steady_climate =  matrix(nrow = 2,ncol = length(results))
    deviations_cpool_complete =  matrix(nrow = 2,ncol = length(results))
    deviations_cflux_climate_change = matrix(nrow = 2,ncol = length(results))
    deviations_cflux_steady_climate = matrix(nrow = 2,ncol = length(results))
    deviations_cflux_complete = matrix(nrow = 2,ncol = length(results))
    parameter_lists = vector(mode = "list", length = length(results))
    for(j in 1:length(results)){


      ## calculate deviations for gpp
      deviations_agpp_climate_change[1,j] = sum(results[[i]]$agpp[rownames(results[[i]]$agpp) %in% 2000:2099,j] -
                                                  results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2000:2099,"Total"])
      deviations_agpp_steady_climate[1,j] = sum(results[[i]]$agpp[rownames(results[[i]]$agpp) %in% 2100:2199,j] -
                                                  results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2100:2199,"Total"])
      deviations_agpp_complete[1,j] = sum(results[[i]]$agpp[rownames(results[[i]]$agpp) %in% 2000:2199,j] -
                                            results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2000:2199,"Total"])

      deviations_agpp_climate_change[2,j] =  sum(results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2001:2100,"Total"])
      deviations_agpp_steady_climate[2,j] = sum(results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2101:2200,"Total"])
      deviations_agpp_complete[2,j] =  sum(results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2001:2200,"Total"])


      ## calculate deviations for cpool
      deviations_cpool_climate_change[1,j] = sum(results[[i]]$cpool[rownames(results[[i]]$cpool) %in% 2000:2099,j] -
                                                   results_Que@dataTypes$cpool[results_Que@dataTypes$cpool[,'Year'] %in% 2000:2099,"Total"])
      deviations_cpool_steady_climate[1,j] = sum(results[[i]]$cpool[rownames(results[[i]]$cpool) %in% 2100:2199,j] -
                                                   results_Que@dataTypes$cpool[results_Que@dataTypes$cpool[,'Year'] %in% 2100:2199,"Total"])
      deviations_cpool_complete[1,j] = sum(results[[i]]$cpool[rownames(results[[i]]$cpool) %in% 2000:2199,j] -
                                             results_Que@dataTypes$cpool[results_Que@dataTypes$cpool[,'Year'] %in% 2000:2199,"Total"])

      deviations_cpool_climate_change[2,j] =  sum(results_Que@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2001:2100,"Total"] )
      deviations_cpool_steady_climate[2,j] = sum(results_Que@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2101:2200,"Total"] )
      deviations_cpool_complete[2,j] = sum(results_Que@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2001:2200,"Total"] )

      ## calculate deviations for cflux

      deviations_cflux_climate_change[1,j] = sum(results[[i]]$cflux[rownames(results[[i]]$cflux) %in% 2000:2099,j] -
                                                   results_Que@dataTypes$cflux[results_Que@dataTypes$cflux[,'Year'] %in% 2000:2099,"NEE"])
      deviations_cflux_steady_climate[1,j] = sum(results[[i]]$cflux[rownames(results[[i]]$cflux) %in% 2100:2199,j] -
                                                   results_Que@dataTypes$cflux[results_Que@dataTypes$cflux[,'Year'] %in% 2100:2199,"NEE"])
      deviations_cflux_complete[1,j] = sum(results[[i]]$cflux[rownames(results[[i]]$cflux) %in% 2000:2199,j] -
                                             results_Que@dataTypes$cflux[results_Que@dataTypes$cflux[,'Year'] %in% 2000:2199,"NEE"])

      deviations_cflux_climate_change[2,j] = sum(results_Que@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2001:2100,"NEE"] )
      deviations_cflux_steady_climate[2,j] = sum(results_Que@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2101:2200,"NEE"] )
      deviations_cflux_complete[2,j] = sum(results_Que@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2001:2200,"NEE"] )


    }

    #save the results
    results_list = list("agpp" = list("steady_climate" = deviations_agpp_steady_climate,
                                      "climate_change" = deviations_agpp_climate_change,
                                      "complete" = deviations_agpp_complete),
                        "cpool" = list("steady_climate" = deviations_cpool_steady_climate,
                                       "climate_change" = deviations_cpool_climate_change,
                                       "complete" = deviations_cpool_complete),
                        'cflux'= list('steady_climate' = deviations_cflux_steady_climate,
                                      'climate_change' = deviations_cflux_climate_change,
                                      'complete' = deviations_cflux_complete),

    )



    saveRDS(results_list,paste0("./LPJrunTest/Results2/Lin_",parameter_name), version = 2)
    gc()
  }
}
