## This script combines the different results to an overall result
## This has to be done, because of RAM constraints


## set up of different things to fill later

lower_bounds = seq(from = 1, to = 191, by = 10)
upper_bounds = seq(from = 10, to = 200, by = 10)
upper_bounds_lai = upper_bounds
upper_bounds_lai[which(upper_bounds_lai == 110)] = 100

effects_cflux_climate_change = matrix(ncol =74, nrow = 200)
effects_cflux_complete = matrix(ncol =74, nrow = 200)
effects_cflux_steady_climate = matrix(ncol =74, nrow = 200)

effects_cpool_climate_change = matrix(ncol =74, nrow = 200)
effects_cpool_complete = matrix(ncol =74, nrow = 200)
effects_cpool_steady_climate = matrix(ncol =74, nrow = 200)

effects_agpp_climate_change = matrix(ncol =74, nrow = 200)
effects_agpp_complete = matrix(ncol =74, nrow = 200)
effects_agpp_steady_climate = matrix(ncol =74, nrow = 200)

effects_cflux_rf_climate_change = matrix(ncol =74, nrow = 200)
effects_cflux_rf_complete = matrix(ncol =74, nrow = 200)
effects_cflux_rf_steady_climate = matrix(ncol =74, nrow = 200)

effects_cpool_rf_climate_change = matrix(ncol =74, nrow = 200)
effects_cpool_rf_complete = matrix(ncol =74, nrow = 200)
effects_cpool_rf_steady_climate = matrix(ncol =74, nrow = 200)

effects_agpp_rf_climate_change = matrix(ncol =74, nrow = 200)
effects_agpp_rf_complete = matrix(ncol =74, nrow = 200)
effects_agpp_rf_steady_climate = matrix(ncol =74, nrow = 200)



interactions_cflux_steady_climate = list()
interactions_cflux_complete = list()
interactions_cflux_climate_change = list()

interactions_cpool_steady_climate = list()
interactions_cpool_complete = list()
interactions_cpool_climate_change = list()

interactions_agpp_steady_climate = list()
interactions_agpp_complete = list()
interactions_agpp_climate_change = list()


## read the data in 10 site steps
for(i in 1:length(lower_bounds)){
  ## extract the stuff of the data which is important for later analysis
  file_carbon_results = paste0("./LPJrunTest/Results/Lin_Mixed_effects_",lower_bounds[i],"to",upper_bounds[i],".rds")
  file_lai_results = paste0("./LPJrunTest/Results/Lin_Mixed_effects_",lower_bounds[i],"to",upper_bounds_lai[i],"_lai.rds")
  results_carbon = readRDS(file_carbon_results)
  results_lai = readRDS(file_lai_results)
  for(j in 1:length(results_carbon)){
    effects_cflux_climate_change[lower_bounds[i]+j-1,] = results_carbon[[j]][["cflux"]][["climate_change"]]
    effects_cflux_steady_climate[lower_bounds[i]+j-1,] = results_carbon[[j]][["cflux"]][["steady_climate"]]
    effects_cflux_complete[lower_bounds[i]+j-1,] = results_carbon[[j]][["cflux"]][["complete"]]

    effects_cflux_rf_climate_change[lower_bounds[i]+j-1,] = results_carbon[[j]][["cflux_rf"]][["climate_change"]]
    effects_cflux_rf_steady_climate[lower_bounds[i]+j-1,] = results_carbon[[j]][["cflux_rf"]][["steady_climate"]]
    effects_cflux_rf_complete[lower_bounds[i]+j-1,] = results_carbon[[j]][["cflux_rf"]][["complete"]]

    interactions_cflux_climate_change[[lower_bounds[i]+j-1]] = results_carbon[[j]][["cflux_interactions"]][['climate_change']]
    interactions_cflux_complete[[lower_bounds[i]+j-1]] = results_carbon[[j]][["cflux_interactions"]][['complete']]
    interactions_cflux_steady_climate[[lower_bounds[i]+j-1]] = results_carbon[[j]][["cflux_interactions"]][['steady_climate']]

    effects_cpool_climate_change[lower_bounds[i]+j-1,] = results_carbon[[j]][["cpool"]][["climate_change"]]
    effects_cpool_steady_climate[lower_bounds[i]+j-1,] = results_carbon[[j]][["cpool"]][["steady_climate"]]
    effects_cpool_complete[lower_bounds[i]+j-1,] = results_carbon[[j]][["cpool"]][["complete"]]

    effects_cpool_rf_climate_change[lower_bounds[i]+j-1,] = results_carbon[[j]][["cpool_rf"]][["climate_change"]]
    effects_cpool_rf_steady_climate[lower_bounds[i]+j-1,] = results_carbon[[j]][["cpool_rf"]][["steady_climate"]]
    effects_cpool_rf_complete[lower_bounds[i]+j-1,] = results_carbon[[j]][["cpool_rf"]][["complete"]]

    interactions_cpool_climate_change[[lower_bounds[i]+j-1]] = results_carbon[[j]][["cpool_interactions"]][['climate_change']]
    interactions_cpool_complete[[lower_bounds[i]+j-1]] = results_carbon[[j]][["cpool_interactions"]][['complete']]
    interactions_cpool_steady_climate[[lower_bounds[i]+j-1]] = results_carbon[[j]][["cpool_interactions"]][['steady_climate']]

    effects_agpp_climate_change[lower_bounds[i]+j-1,] = results_carbon[[j]][["agpp"]][["climate_change"]]
    effects_agpp_steady_climate[lower_bounds[i]+j-1,] = results_carbon[[j]][["agpp"]][["steady_climate"]]
    effects_agpp_complete[lower_bounds[i]+j-1,] = results_carbon[[j]][["agpp"]][["complete"]]

    effects_agpp_rf_climate_change[lower_bounds[i]+j-1,] = results_carbon[[j]][["agpp_rf"]][["climate_change"]]
    effects_agpp_rf_steady_climate[lower_bounds[i]+j-1,] = results_carbon[[j]][["agpp_rf"]][["steady_climate"]]
    effects_agpp_rf_complete[lower_bounds[i]+j-1,] = results_carbon[[j]][["agpp_rf"]][["complete"]]

    interactions_agpp_climate_change[[lower_bounds[i]+j-1]] = results_carbon[[j]][["agpp_interactions"]][['climate_change']]
    interactions_agpp_complete[[lower_bounds[i]+j-1]] = results_carbon[[j]][["agpp_interactions"]][['complete']]
    interactions_agpp_steady_climate[[lower_bounds[i]+j-1]] = results_carbon[[j]][["agpp_interactions"]][['steady_climate']]

    # if(length(results_lai[[j]][["lai"]][["steady_climate"]])==1){
    #   effects_lai_steady_climate[lower_bounds[i]+j-1,] = NA
    #   interactions_lai_steady_climate[[lower_bounds[i]+j-1]] = NA
    # }
    # else{
    #   effects_lai_steady_climate[lower_bounds[i]+j-1,] = results_lai[[j]][["lai"]][["steady_climate"]]
    #   effects_lai_rf_steady_climate[lower_bounds[i]+j-1,] = results_lai[[j]][["lai_rf"]][["steady_climate"]]
    #   interactions_lai_steady_climate[[lower_bounds[i]+j-1]] = results_lai[[j]][["lai_interactions"]][['steady_climate']]
    # }
    #
    # if(length(results_lai[[j]][["lai"]][["climate_change"]])==1){
    #   effects_lai_climate_change[lower_bounds[i]+j-1,] = NA
    #   interactions_lai_climate_change[[lower_bounds[i]+j-1]] = NA
    # }
    # else{
    #   effects_lai_climate_change[lower_bounds[i]+j-1,] = results_lai[[j]][["lai"]][["climate_change"]]
    #   effects_lai_rf_climate_change[lower_bounds[i]+j-1,] = results_lai[[j]][["lai_rf"]][["climate_change"]]
    #   interactions_lai_climate_change[[lower_bounds[i]+j-1]] = results_lai[[j]][["lai_interactions"]][['climate_change']]
    # }
    #
    # if(length(results_lai[[j]][["lai"]][["complete"]])==1){
    #   effects_lai_complete[lower_bounds[i]+j-1,] = NA
    #   interactions_lai_complete[[lower_bounds[i]+j-1]] = NA
    # }
    # else{
    #   effects_lai_complete[lower_bounds[i]+j-1,] = results_lai[[j]][["lai"]][["complete"]]
    #   effects_lai_rf_complete[lower_bounds[i]+j-1,] = results_lai[[j]][["lai_rf"]][["complete"]]
    #   interactions_lai_complete[[lower_bounds[i]+j-1]] = results_lai[[j]][["lai_interactions"]][['complete']]
    # }
  }
}


result_list = list("agpp" = list("steady_climate" = effects_agpp_steady_climate,
                                "climate_change" = effects_agpp_climate_change,
                                "complete" = effects_agpp_complete),
                  "agpp_rf" = list("steady_climate" = effects_agpp_rf_steady_climate,
                                "climate_change" = effects_agpp_rf_climate_change,
                                "complete" = effects_agpp_rf_complete),
                  "agpp_interactions" = list("steady_climate" = interactions_agpp_steady_climate,
                                            "climate_change" = interactions_agpp_climate_change,
                                            "complete" = interactions_agpp_complete),
                  "cpool" = list("steady_climate" = effects_cpool_steady_climate,
                                 "climate_change" = effects_cpool_climate_change,
                                 "complete" = effects_cpool_complete),
                  "cpool_rf" = list("steady_climate" = effects_cpool_rf_steady_climate,
                                   "climate_change" = effects_cpool_rf_climate_change,
                                   "complete" = effects_cpool_rf_complete),
                    "cpool_interactions" = list("steady_climate" = interactions_cpool_steady_climate,
                                               "climate_change" = interactions_cpool_climate_change,
                                               "complete" = interactions_cpool_complete),
                  'cflux'= list("steady_climate" = effects_cflux_steady_climate,
                                "climate_change" = effects_cflux_climate_change,
                                "complete" = effects_cflux_complete),
                  "cflux_rf" = list("steady_climate" = effects_cflux_rf_steady_climate,
                                   "climate_change" = effects_cflux_rf_climate_change,
                                   "complete" = effects_cflux_rf_complete),
                  "cflux_interactions" = list("steady_climate" = interactions_cflux_steady_climate,
                                              "climate_change" = interactions_cflux_climate_change,
                                                "complete" = interactions_cflux_complete)
)

saveRDS(result_list,"LPJrunTest/Results/Mixed_effects_lin.rds")

## read in data to derive names
getwd()
mixed_results =  readRDS(paste0("./../Results_sensi/Mixed_results/Mixed_0.25_42.25.rds"))

# compare the names with the previous names to make sure, we have the right order
# get how often parameters are used in mixed and give them weights based on how often they are there
# in the end we have a map with names, weights and positions
parameters = readRDS("ParameterMetaData/Parameter_mixed.rds")
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameternames = c(as.character(parameters$NameRLPJ), drivernames)
parameter_names = names(mixed_results[[j]]@runInfo$parameterList[which(names(mixed_results[[j]]@runInfo$parameterList) %in% parameternames)])

## define prefixes, which could be annoying in when comparing results
prefixes = c("common_","Pic_abi_","Fag_syl_","Pin_syl_",
             "tree_","shade_intolerant_","run_","shade_tolerant_",
             "boreal_","temperate_","MixedC_", "needle_leaved_",
             "broad_leaved_")

# delete these prefixes
for(name in prefixes){
  parameter_names = gsub(name,"",parameter_names)
}
mapping_matrix = list()
for(i in 1:length(unique(parameter_names))){
  mapping_matrix[[i]] = which(parameter_names == unique(parameter_names)[i])
}
weight_matrix = list()
for(i in 1:length(mapping_matrix)){
  if(length(mapping_matrix[[i]]) == 1){
    weight_matrix[[i]] = 1
  }
  else if(length(mapping_matrix[[i]]) == 2){
    weight_matrix[[i]] = c(2/3,1/3)
    print(i)
  }
  else{
    weight_matrix[[i]] = c(1/3,1/3,1/3)
  }
}

mapping_list = list("parameternames" = unique(parameter_names),
                    "position_mapping" = mapping_matrix,
                    "weight_mapping" = weight_matrix)

saveRDS(mapping_list,"./LPJrunTest/Results/mapping_mixed_lin.rds")









