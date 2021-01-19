library(rLPJGUESS)
library(zoo)
#define the ouputs which have been investigated before
typeList <- c("cmass", "lai","agpp","cpool","anpp","cflux","fpc","speciesdiam","dens")
scaleLPJ_PFT <- "normal"
#set the main directory to LPJrunTest
mainDir <- file.path(getwd(), "LPJrunTest")

#read in the data
sites <- readRDS("sites_data.rds")

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
  parameter_name = paste0("Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds")

  # make the templates runable
  AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                  defaultlist = defaultparameters$defaultlist,
                  MainTemplateDestination = "./LPJrunTest/main_new.ins",
                  PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                  NameMainFile = "main.ins", NamePftFile = "pft.ins")

  # produce runable parameters
  parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Fag_syl","Pic_abi","Pin_syl"))

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

  results <- try(readRDS(paste0("./../Results_sensi/",parameter_name)))
  if("try-error" %in% class(results)){
    print(parameter_name)
    next
  }
  #if not there go to next site, else calculate difference
  else{

    base_line_Fag_syl_climate_change = matrix(nrow =2, ncol =1)
    base_line_Fag_syl_steady_climate = matrix(nrow =2, ncol =1)
    base_line_Fag_syl_complete = matrix(nrow =2, ncol =1)

    base_line_Fag_syl_climate_change[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Fag_syl"])/
                              sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"]))
    base_line_Fag_syl_climate_change[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"])

    base_line_Fag_syl_steady_climate[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Fag_syl"])/
                                               sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"]))
    base_line_Fag_syl_steady_climate[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"])

    base_line_Fag_syl_complete[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Fag_syl"])/
                                               sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"]))
    base_line_Fag_syl_complete[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"])


    base_line_Pin_syl_climate_change = matrix(nrow =2, ncol =1)
    base_line_Pin_syl_steady_climate = matrix(nrow =2, ncol =1)
    base_line_Pin_syl_complete = matrix(nrow =2, ncol =1)

    base_line_Pin_syl_climate_change[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Pin_syl"])/
                                               sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"]))
    base_line_Pin_syl_climate_change[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"])

    base_line_Pin_syl_steady_climate[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Pin_syl"])/
                                               sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"]))
    base_line_Pin_syl_steady_climate[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"])

    base_line_Pin_syl_complete[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Pin_syl"])/
                                         sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"]))
    base_line_Pin_syl_complete[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"])

    base_line_Pic_abi_climate_change = matrix(nrow =2, ncol =1)
    base_line_Pic_abi_steady_climate = matrix(nrow =2, ncol =1)
    base_line_Pic_abi_complete = matrix(nrow =2, ncol =1)

    base_line_Pic_abi_climate_change[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Pic_abi"])/
                                               sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"]))
    base_line_Pic_abi_climate_change[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"])

    base_line_Pic_abi_steady_climate[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Pic_abi"])/
                                               sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"]))
    base_line_Pic_abi_steady_climate[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"])

    base_line_Pic_abi_complete[1,1] = (sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Pic_abi"])/
                                         sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"]))
    base_line_Pic_abi_complete[2,1] =  sum(results_Que@dataTypes$lai[results_Que@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"])

    deviations_Fag_syl_lai_climate_change = matrix(nrow = 1,ncol = length(results))
    deviations_Fag_syl_lai_steady_climate =  matrix(nrow = 1,ncol = length(results))
    deviations_Fag_syl_lai_complete =  matrix(nrow = 1,ncol = length(results))
    deviations_Pin_syl_lai_climate_change =  matrix(nrow = 1,ncol = length(results))
    deviations_Pin_syl_lai_steady_climate = matrix(nrow = 1,ncol = length(results))
    deviations_Pin_syl_lai_complete =  matrix(nrow = 1,ncol = length(results))
    deviations_Pic_abi_lai_climate_change =  matrix(nrow = 1,ncol = length(results))
    deviations_Pic_abi_lai_steady_climate =  matrix(nrow = 1,ncol = length(results))
    deviations_Pic_abi_lai_complete =  matrix(nrow = 1,ncol = length(results))

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
      ## calculate deviations for npp
      deviations_Fag_syl_lai_climate_change[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2100,"Fag_syl"])/
                                                   sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"]))

      deviations_Fag_syl_lai_steady_climate[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2101:2200,"Fag_syl"])/
                                                      sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"]))

      deviations_Fag_syl_lai_complete[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2200,"Fag_syl"])/
                                                      sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"]))



      ## calculate deviations for npp
      deviations_Pic_abi_lai_climate_change[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2100,"Pic_abi"])/
                                                      sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"]))

      deviations_Pic_abi_lai_steady_climate[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2101:2200,"Pic_abi"])/
                                                      sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"]))

      deviations_Pic_abi_lai_complete[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2200,"Pic_abi"])/
                                                sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"]))

      ## calculate deviations for cmass

      deviations_Pin_syl_lai_climate_change[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2100,"Pin_syl"])/
                                                      sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2100,"Total"]))

      deviations_Pin_syl_lai_steady_climate[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2101:2200,"Pin_syl"])/
                                                      sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2101:2200,"Total"]))

      deviations_Pin_syl_lai_complete[1,j] = (sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2200,"Pin_syl"])/
                                                sum(results[[j]]@dataTypes$lai[results[[j]]@dataTypes$lai[,'Year'] %in% 2001:2200,"Total"]))


      deviations_anpp_climate_change[1,j] = sum(results[[j]]@dataTypes$anpp[results[[j]]@dataTypes$anpp[,'Year'] %in% 2001:2100,"Total"] -
                                                   results_Que@dataTypes$anpp[results_Que@dataTypes$anpp[,'Year'] %in% 2001:2100,"Total"])
      deviations_anpp_steady_climate[1,j] = sum(results[[j]]@dataTypes$anpp[results[[j]]@dataTypes$anpp[,'Year'] %in% 2101:2200,"Total"] -
                                                   results_Que@dataTypes$anpp[results_Que@dataTypes$anpp[,'Year'] %in% 2101:2200,"Total"])
      deviations_anpp_complete[1,j] = sum(results[[j]]@dataTypes$anpp[results[[j]]@dataTypes$anpp[,'Year'] %in% 2001:2200,"Total"] -
                                             results_Que@dataTypes$anpp[results_Que@dataTypes$anpp[,'Year'] %in% 2001:2200,"Total"])

      deviations_anpp_climate_change[2,j] = sum(results_Que@dataTypes$anpp[results_Que@dataTypes$anpp[,'Year']  %in% 2001:2100,"Total"])
      deviations_anpp_steady_climate[2,j] = sum(results_Que@dataTypes$anpp[results_Que@dataTypes$anpp[,'Year']  %in% 2101:2200,"Total"])
      deviations_anpp_complete[2,j] = sum(results_Que@dataTypes$anpp[results_Que@dataTypes$anpp[,'Year']  %in% 2001:2200,"Total"])


      deviations_agpp_climate_change[1,j] = sum(results[[j]]@dataTypes$agpp[results[[j]]@dataTypes$agpp[,'Year'] %in% 2001:2100,"Total"] -
                                                   results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2001:2100,"Total"])
      deviations_agpp_steady_climate[1,j] = sum(results[[j]]@dataTypes$agpp[results[[j]]@dataTypes$agpp[,'Year']  %in% 2101:2200,"Total"] -
                                                   results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2101:2200,"Total"])
      deviations_agpp_complete[1,j] = sum(results[[j]]@dataTypes$agpp[results[[j]]@dataTypes$agpp[,'Year']  %in% 2001:2200,"Total"] -
                                             results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2001:2200,"Total"])

      deviations_agpp_climate_change[2,j] =  sum(results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2001:2100,"Total"])
      deviations_agpp_steady_climate[2,j] = sum(results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2101:2200,"Total"])
      deviations_agpp_complete[2,j] =  sum(results_Que@dataTypes$agpp[results_Que@dataTypes$agpp[,'Year'] %in% 2001:2200,"Total"])

      ## calculate deviations for cmass
      deviations_cmass_climate_change[1,j] = sum(results[[j]]@dataTypes$cmass[results[[j]]@dataTypes$cmass[,'Year'] %in% 2001:2100,"Total"] -
                                                    results_Que@dataTypes$cmass[results_Que@dataTypes$cmass[,'Year'] %in% 2001:2100,"Total"])
      deviations_cmass_steady_climate[1,j] = sum(results[[j]]@dataTypes$cmass[results[[j]]@dataTypes$cmass[,'Year'] %in% 2101:2200,"Total"] -
                                                    results_Que@dataTypes$cmass[results_Que@dataTypes$cmass[,'Year'] %in% 2101:2200,"Total"])
      deviations_cmass_complete[1,j] = sum(results[[j]]@dataTypes$cmass[results[[j]]@dataTypes$cmass[,'Year'] %in% 2001:2200,"Total"] -
                                              results_Que@dataTypes$cmass[results_Que@dataTypes$cmass[,'Year'] %in% 2001:2200,"Total"])

      deviations_cmass_climate_change[2,j] = sum(results_Que@dataTypes$cmass[results[[j]]@dataTypes$cmass[,'Year'] %in% 2001:2100,"Total"] )
      deviations_cmass_steady_climate[2,j] = sum(results_Que@dataTypes$cmass[results[[j]]@dataTypes$cmass[,'Year'] %in% 2101:2200,"Total"] )
      deviations_cmass_complete[2,j] = sum(results_Que@dataTypes$cmass[results[[j]]@dataTypes$cmass[,'Year'] %in% 2001:2200,"Total"] )

      ## calculate deviations for cpool
      deviations_cpool_climate_change[1,j] = sum(results[[j]]@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2001:2100,"Total"] -
                                                    results_Que@dataTypes$cpool[results_Que@dataTypes$cpool[,'Year'] %in% 2001:2100,"Total"])
      deviations_cpool_steady_climate[1,j] = sum(results[[j]]@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2101:2200,"Total"] -
                                                    results_Que@dataTypes$cpool[results_Que@dataTypes$cpool[,'Year'] %in% 2101:2200,"Total"])
      deviations_cpool_complete[1,j] = sum(results[[j]]@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2001:2200,"Total"] -
                                              results_Que@dataTypes$cpool[results_Que@dataTypes$cpool[,'Year'] %in% 2001:2200,"Total"])

      deviations_cpool_climate_change[2,j] =  sum(results_Que@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2001:2100,"Total"] )
      deviations_cpool_steady_climate[2,j] = sum(results_Que@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2101:2200,"Total"] )
      deviations_cpool_complete[2,j] = sum(results_Que@dataTypes$cpool[results[[j]]@dataTypes$cpool[,'Year'] %in% 2001:2200,"Total"] )

      ## calculate deviations for cflux

      deviations_cflux_climate_change[1,j] = sum(results[[j]]@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2001:2100,"NEE"] -
                                                    results_Que@dataTypes$cflux[results_Que@dataTypes$cflux[,'Year'] %in% 2001:2100,"NEE"])
      deviations_cflux_steady_climate[1,j] = sum(results[[j]]@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2101:2200,"NEE"] -
                                                    results_Que@dataTypes$cflux[results_Que@dataTypes$cflux[,'Year'] %in% 2101:2200,"NEE"])
      deviations_cflux_complete[1,j] = sum(results[[j]]@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2001:2200,"NEE"] -
                                              results_Que@dataTypes$cflux[results_Que@dataTypes$cflux[,'Year'] %in% 2001:2200,"NEE"])


      deviations_cflux_climate_change[2,j] = sum(results_Que@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2001:2100,"NEE"] )
      deviations_cflux_steady_climate[2,j] = sum(results_Que@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2101:2200,"NEE"] )
      deviations_cflux_complete[2,j] = sum(results_Que@dataTypes$cflux[results[[j]]@dataTypes$cflux[,'Year'] %in% 2001:2200,"NEE"] )


      parameters = readRDS("Parameter_mixed.rds")
      drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
      parameternames = c(as.character(parameters$NameRLPJ), drivernames)



      parameter_lists[[j]] = as.numeric(results[[j]]@runInfo$parameterList[which(names(results[[j]]@runInfo$parameterList) %in% parameternames)])



    }


    results_list = list("Fag_syl" = list("steady_climate" = deviations_Fag_syl_lai_steady_climate,
                                      "climate_change" = deviations_Fag_syl_lai_climate_change,
                                      "complete" = deviations_Fag_syl_lai_complete),
                        "Pic_abi" = list("steady_climate" = deviations_Pic_abi_lai_steady_climate,
                                         "climate_change" = deviations_Pic_abi_lai_climate_change,
                                         "complete" = deviations_Pic_abi_lai_complete),
                        "Pin_syl" = list("steady_climate" = deviations_Pin_syl_lai_steady_climate,
                                         "climate_change" = deviations_Pin_syl_lai_climate_change,
                                         "complete" = deviations_Pin_syl_lai_complete),
                        "base_line_Fag_syl" = list("steady_climate" = base_line_Fag_syl_steady_climate ,
                                               "climate_change" = base_line_Fag_syl_climate_change ,
                                               "complete" = base_line_Fag_syl_complete),
                        "base_line_Pic_abi" = list("steady_climate" = base_line_Pic_abi_steady_climate ,
                                               "climate_change" = base_line_Pic_abi_climate_change ,
                                               "complete" = base_line_Pic_abi_complete),
                        "base_line_Pin_syl" = list("steady_climate" = base_line_Pin_syl_steady_climate ,
                                               "climate_change" = base_line_Pin_syl_climate_change ,
                                               "complete" = base_line_Pin_syl_complete),
                        "anpp" = list("steady_climate" = deviations_anpp_steady_climate,
                                      "climate_change" = deviations_anpp_climate_change,
                                      "complete" = deviations_anpp_complete),
                        "agpp" = list("steady_climate" = deviations_agpp_steady_climate,
                                      "climate_change" = deviations_agpp_climate_change,
                                      "complete" = deviations_agpp_complete),
                        "cmass" = list("steady_climate" = deviations_cmass_steady_climate,
                                       "climate_change" = deviations_cmass_climate_change,
                                       "complete" = deviations_cmass_complete),
                        "cpool" = list("steady_climate" = deviations_cpool_steady_climate,
                                       "climate_change" = deviations_cpool_climate_change,
                                       "complete" = deviations_cpool_complete),
                        'cflux'= list('steady_climate' = deviations_cflux_steady_climate,
                                      'climate_change' = deviations_cflux_climate_change,
                                      'complete' = deviations_cflux_complete),
                        'parameters' = parameter_lists)


    saveRDS(results_list,paste0("./LPJrunTest/Results/Lin_",parameter_name), version = 2)
    rm(results_list, results, parameter_lists)
    gc()
  }
}
