# loading the sites data
sites <- readRDS("EnvironmentalData/sites_data.rds")

# function to remap effects for mixed simulations
remap_paramters <- function(parameters, weighting_scheme, position_scheme){
  mapped_parameters = matrix(ncol = length(weighting_scheme),
                             nrow = nrow(parameters))
  for(site in 1:nrow(parameters)){
    for(parameter in 1:length(position_scheme)){
      mapped_parameters[site,parameter] = sum(abs(parameters[site,position_scheme[[parameter]]])*
                                                weighting_scheme[[parameter]])
    }
  }
  return(mapped_parameters)
}


remap_paramters_rel <- function(parameters, weighting_scheme, position_scheme){
  mapped_parameters = matrix(ncol = length(weighting_scheme),
                             nrow = nrow(parameters))
  for(site in 1:nrow(parameters)){
    for(parameter in 1:length(position_scheme)){
      mapped_parameters[site,parameter] = sum(parameters[site,position_scheme[[parameter]]]*
                                                weighting_scheme[[parameter]])
    }
  }
  return(mapped_parameters)
}

## getting the ordered parameter names ##

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameters_pic_abi = c(as.character(parameters$Pic_abi$NameRLPJ), drivernames)
parameters_fag_syl = c(as.character(parameters$Fag_syl$NameRLPJ), drivernames)
parameters_pin_syl = c(as.character(parameters$Pin_syl$NameRLPJ), drivernames)
mixed_results = readRDS("./LPJrunTest/Results/mapping_mixed.rds")


results =  readRDS(paste0("./LPJrunTest/Results/Pic_abi_0.25_42.25.rds"))
parameternames_ordered_pic_abi = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pic_abi)])

results =  readRDS(paste0("./LPJrunTest/Results/Fag_syl_0.25_42.25.rds"))
parameternames_ordered_fag_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_fag_syl)])

results =  readRDS(paste0("./LPJrunTest/Results/Pin_syl_0.25_42.25.rds"))
parameternames_ordered_pin_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pin_syl)])

parameternames_ordered_pic_abi2 = substring(parameternames_ordered_pic_abi,regexpr("_", parameternames_ordered_pic_abi) + 1)
parameternames_ordered_fag_syl2 = substring(parameternames_ordered_fag_syl,regexpr("_", parameternames_ordered_fag_syl) + 1)
parameternames_ordered_pin_syl2 = substring(parameternames_ordered_pin_syl,regexpr("_", parameternames_ordered_pin_syl) + 1)
parameternames_ordered_mixed2 = mixed_results$parameternames
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  parameternames_ordered_pic_abi2 = gsub(noisy_things[i],"",parameternames_ordered_pic_abi2)
  parameternames_ordered_fag_syl2 = gsub(noisy_things[i],"",parameternames_ordered_fag_syl2)
  parameternames_ordered_pin_syl2 = gsub(noisy_things[i],"",parameternames_ordered_pin_syl2)
  parameternames_ordered_mixed2 = gsub(noisy_things[i],"",parameternames_ordered_mixed2)
}

ordering_fag_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_fag_syl2)
ordering_pin_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_pin_syl2)
ordering_mixed = match(parameternames_ordered_pic_abi2, parameternames_ordered_mixed2)

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
drivers = c("ph","co2","prec","temp","insol","ndep")

variablenames = c(as.character(parameters$Pic_abi[,"NameRLPJ"]),drivers)
variablenames2 = substring(variablenames,regexpr("_", variablenames) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  variablenames2 = gsub(noisy_things[i],"",variablenames2)
}

grouping = read.csv("ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")
grouping = c(as.character(grouping[,1]), rep("Drivers",6))

order_grouping = match(parameternames_ordered_pic_abi2, variablenames2)
grouping2 = grouping[order_grouping]



### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects.rds")

effects_mixed = readRDS("LPJrunTest/Results/mixed_effects.rds")

###### cpool ######

#### climate change ###


mean_effects_Pic_abi_cpool_climate_change = effects_Pic_abi[["cpool"]][["climate_change"]]
growth_sites_Pic_abi_climate_change = which(mean_effects_Pic_abi_cpool_climate_change[,39]/100 > 2)
main_effects_Pic_abi_cpool_climate_change = mean_effects_Pic_abi_cpool_climate_change[growth_sites_Pic_abi_climate_change,]

main_effects_Pic_abi_cpool_weight_climate_change_abs = apply(abs(main_effects_Pic_abi_cpool_climate_change[,1:38]),2,mean)
main_effects_Pic_abi_cpool_weight_climate_change_abs = main_effects_Pic_abi_cpool_weight_climate_change_abs/sum(abs(main_effects_Pic_abi_cpool_weight_climate_change_abs))


main_effects_Pic_abi_cpool_weight_climate_change_rel = apply(main_effects_Pic_abi_cpool_climate_change[,1:38],2,mean)
main_effects_Pic_abi_cpool_weight_climate_change_rel = main_effects_Pic_abi_cpool_weight_climate_change_rel/sum(abs(main_effects_Pic_abi_cpool_weight_climate_change_rel))



mean_effects_Fag_syl_cpool_climate_change = effects_Fag_syl[["cpool"]][["climate_change"]]
growth_sites_Fag_syl_climate_change = which(mean_effects_Fag_syl_cpool_climate_change[,39]/100 > 2)
main_effects_Fag_syl_cpool_climate_change = mean_effects_Fag_syl_cpool_climate_change[growth_sites_Fag_syl_climate_change,]

main_effects_Fag_syl_cpool_weight_climate_change_abs = apply(abs(main_effects_Fag_syl_cpool_climate_change[,1:38]),2,mean)
main_effects_Fag_syl_cpool_weight_climate_change_abs = main_effects_Fag_syl_cpool_weight_climate_change_abs/sum(abs(main_effects_Fag_syl_cpool_weight_climate_change_abs))


main_effects_Fag_syl_cpool_weight_climate_change_rel = apply(main_effects_Fag_syl_cpool_climate_change[,1:38],2,mean)
main_effects_Fag_syl_cpool_weight_climate_change_rel = main_effects_Fag_syl_cpool_weight_climate_change_rel/sum(abs(main_effects_Fag_syl_cpool_weight_climate_change_rel))


main_effects_Fag_syl_cpool_weight_climate_change_abs = main_effects_Fag_syl_cpool_weight_climate_change_abs[ordering_fag_syl]
main_effects_Fag_syl_cpool_weight_climate_change_rel = main_effects_Fag_syl_cpool_weight_climate_change_rel[ordering_fag_syl]



mean_effects_Pin_syl_cpool_climate_change = effects_Pin_syl[["cpool"]][["climate_change"]]
growth_sites_Pin_syl_climate_change = which(mean_effects_Pin_syl_cpool_climate_change[,39]/100 > 2)
main_effects_Pin_syl_cpool_climate_change = mean_effects_Pin_syl_cpool_climate_change[growth_sites_Pin_syl_climate_change,]

main_effects_Pin_syl_cpool_weight_climate_change_abs = apply(abs(main_effects_Pin_syl_cpool_climate_change[,1:38]),2,mean)
main_effects_Pin_syl_cpool_weight_climate_change_abs = main_effects_Pin_syl_cpool_weight_climate_change_abs/sum(abs(main_effects_Pin_syl_cpool_weight_climate_change_abs))


main_effects_Pin_syl_cpool_weight_climate_change_rel = apply(main_effects_Pin_syl_cpool_climate_change[,1:38],2,mean)
main_effects_Pin_syl_cpool_weight_climate_change_rel = main_effects_Pin_syl_cpool_weight_climate_change_rel/sum(abs(main_effects_Pin_syl_cpool_weight_climate_change_rel))


main_effects_Pin_syl_cpool_weight_climate_change_abs = main_effects_Pin_syl_cpool_weight_climate_change_abs[ordering_pin_syl]
main_effects_Pin_syl_cpool_weight_climate_change_rel = main_effects_Pin_syl_cpool_weight_climate_change_rel[ordering_pin_syl]


mean_effects_mixed_cpool_climate_change = effects_mixed[["cpool"]][["climate_change"]]
growth_sites_mixed_climate_change = which(mean_effects_mixed_cpool_climate_change[,74]/100 > 2)
mean_effects_mixed_cpool_climate_change = mean_effects_mixed_cpool_climate_change[growth_sites_mixed_climate_change,]

remaped_effects_cpool_climate_change_abs = remap_paramters(mean_effects_mixed_cpool_climate_change, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool_weight_climate_change_abs = apply(abs(remaped_effects_cpool_climate_change_abs),2,mean)

main_effects_mixed_cpool_weight_climate_change_abs = main_effects_mixed_cpool_weight_climate_change_abs/sum(abs(main_effects_mixed_cpool_weight_climate_change_abs))

remaped_effects_cpool_climate_change_rel = remap_paramters_rel(mean_effects_mixed_cpool_climate_change, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool_weight_climate_change_rel = apply(remaped_effects_cpool_climate_change_rel,2,mean)

main_effects_mixed_cpool_weight_climate_change_rel = main_effects_mixed_cpool_weight_climate_change_rel/sum(abs(main_effects_mixed_cpool_weight_climate_change_rel))



main_effects_mixed_cpool_weight_climate_change_abs = main_effects_mixed_cpool_weight_climate_change_abs[ordering_mixed]
main_effects_mixed_cpool_weight_climate_change_rel = main_effects_mixed_cpool_weight_climate_change_rel[ordering_mixed]



### steady_climate



mean_effects_Pic_abi_cpool_steady_climate = effects_Pic_abi[["cpool"]][["steady_climate"]]
growth_sites_Pic_abi_steady_climate = which(mean_effects_Pic_abi_cpool_steady_climate[,39]/100 > 2)
main_effects_Pic_abi_cpool_steady_climate = mean_effects_Pic_abi_cpool_steady_climate[growth_sites_Pic_abi_steady_climate,]

main_effects_Pic_abi_cpool_weight_steady_climate_abs = apply(abs(main_effects_Pic_abi_cpool_steady_climate[,1:38]),2,mean)
main_effects_Pic_abi_cpool_weight_steady_climate_abs = main_effects_Pic_abi_cpool_weight_steady_climate_abs/sum(abs(main_effects_Pic_abi_cpool_weight_steady_climate_abs))


main_effects_Pic_abi_cpool_weight_steady_climate_rel = apply(main_effects_Pic_abi_cpool_steady_climate[,1:38],2,mean)
main_effects_Pic_abi_cpool_weight_steady_climate_rel = main_effects_Pic_abi_cpool_weight_steady_climate_rel/sum(abs(main_effects_Pic_abi_cpool_weight_steady_climate_rel))



mean_effects_Fag_syl_cpool_steady_climate = effects_Fag_syl[["cpool"]][["steady_climate"]]
growth_sites_Fag_syl_steady_climate = which(mean_effects_Fag_syl_cpool_steady_climate[,39]/100 > 2)
main_effects_Fag_syl_cpool_steady_climate = mean_effects_Fag_syl_cpool_steady_climate[growth_sites_Fag_syl_steady_climate,]

main_effects_Fag_syl_cpool_weight_steady_climate_abs = apply(abs(main_effects_Fag_syl_cpool_steady_climate[,1:38]),2,mean)
main_effects_Fag_syl_cpool_weight_steady_climate_abs = main_effects_Fag_syl_cpool_weight_steady_climate_abs/sum(abs(main_effects_Fag_syl_cpool_weight_steady_climate_abs))


main_effects_Fag_syl_cpool_weight_steady_climate_rel = apply(main_effects_Fag_syl_cpool_steady_climate[,1:38],2,mean)
main_effects_Fag_syl_cpool_weight_steady_climate_rel = main_effects_Fag_syl_cpool_weight_steady_climate_rel/sum(abs(main_effects_Fag_syl_cpool_weight_steady_climate_rel))


main_effects_Fag_syl_cpool_weight_steady_climate_abs = main_effects_Fag_syl_cpool_weight_steady_climate_abs[ordering_fag_syl]
main_effects_Fag_syl_cpool_weight_steady_climate_rel = main_effects_Fag_syl_cpool_weight_steady_climate_rel[ordering_fag_syl]



mean_effects_Pin_syl_cpool_steady_climate = effects_Pin_syl[["cpool"]][["steady_climate"]]
growth_sites_Pin_syl_steady_climate = which(mean_effects_Pin_syl_cpool_steady_climate[,39]/100 > 2)
main_effects_Pin_syl_cpool_steady_climate = mean_effects_Pin_syl_cpool_steady_climate[growth_sites_Pin_syl_steady_climate,]

main_effects_Pin_syl_cpool_weight_steady_climate_abs = apply(abs(main_effects_Pin_syl_cpool_steady_climate[,1:38]),2,mean)
main_effects_Pin_syl_cpool_weight_steady_climate_abs = main_effects_Pin_syl_cpool_weight_steady_climate_abs/sum(abs(main_effects_Pin_syl_cpool_weight_steady_climate_abs))


main_effects_Pin_syl_cpool_weight_steady_climate_rel = apply(main_effects_Pin_syl_cpool_steady_climate[,1:38],2,mean)
main_effects_Pin_syl_cpool_weight_steady_climate_rel = main_effects_Pin_syl_cpool_weight_steady_climate_rel/sum(abs(main_effects_Pin_syl_cpool_weight_steady_climate_rel))


main_effects_Pin_syl_cpool_weight_steady_climate_abs = main_effects_Pin_syl_cpool_weight_steady_climate_abs[ordering_pin_syl]
main_effects_Pin_syl_cpool_weight_steady_climate_rel = main_effects_Pin_syl_cpool_weight_steady_climate_rel[ordering_pin_syl]


mean_effects_mixed_cpool_steady_climate = effects_mixed[["cpool"]][["steady_climate"]]
growth_sites_mixed_steady_climate = which(mean_effects_mixed_cpool_steady_climate[,74]/100 > 2)
mean_effects_mixed_cpool_steady_climate = mean_effects_mixed_cpool_steady_climate[growth_sites_mixed_steady_climate,]

remaped_effects_cpool_steady_climate_abs = remap_paramters(mean_effects_mixed_cpool_steady_climate, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool_weight_steady_climate_abs = apply(abs(remaped_effects_cpool_steady_climate_abs),2,mean)

main_effects_mixed_cpool_weight_steady_climate_abs = main_effects_mixed_cpool_weight_steady_climate_abs/sum(abs(main_effects_mixed_cpool_weight_steady_climate_abs))

remaped_effects_cpool_steady_climate_rel = remap_paramters_rel(mean_effects_mixed_cpool_steady_climate, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool_weight_steady_climate_rel = apply(remaped_effects_cpool_steady_climate_rel,2,mean)

main_effects_mixed_cpool_weight_steady_climate_rel = main_effects_mixed_cpool_weight_steady_climate_rel/sum(abs(main_effects_mixed_cpool_weight_steady_climate_rel))



main_effects_mixed_cpool_weight_steady_climate_abs = main_effects_mixed_cpool_weight_steady_climate_abs[ordering_mixed]
main_effects_mixed_cpool_weight_steady_climate_rel = main_effects_mixed_cpool_weight_steady_climate_rel[ordering_mixed]

####
df_effects_cpool_climate_change_abs = rbind(main_effects_Fag_syl_cpool_weight_climate_change_abs,
                                        main_effects_Pin_syl_cpool_weight_climate_change_abs,
                                        main_effects_Pic_abi_cpool_weight_climate_change_abs,
                                        main_effects_mixed_cpool_weight_climate_change_abs)

df_effects_cpool_climate_change_rel = rbind(main_effects_Fag_syl_cpool_weight_climate_change_rel,
                                            main_effects_Pin_syl_cpool_weight_climate_change_rel,
                                            main_effects_Pic_abi_cpool_weight_climate_change_rel,
                                            main_effects_mixed_cpool_weight_climate_change_rel)

effects_cpool_climate_change_abs = apply(df_effects_cpool_climate_change_abs, 2, mean)
effects_cpool_climate_change_rel = apply(df_effects_cpool_climate_change_rel, 2, mean)

df_effects_cpool_steady_climate_abs = rbind(main_effects_Fag_syl_cpool_weight_steady_climate_abs,
                                            main_effects_Pin_syl_cpool_weight_steady_climate_abs,
                                            main_effects_Pic_abi_cpool_weight_steady_climate_abs,
                                            main_effects_mixed_cpool_weight_steady_climate_abs)

df_effects_cpool_steady_climate_rel = rbind(main_effects_Fag_syl_cpool_weight_steady_climate_rel,
                                            main_effects_Pin_syl_cpool_weight_steady_climate_rel,
                                            main_effects_Pic_abi_cpool_weight_steady_climate_rel,
                                            main_effects_mixed_cpool_weight_steady_climate_rel)

effects_cpool_steady_climate_abs = apply(df_effects_cpool_steady_climate_abs, 2, mean)
effects_cpool_steady_climate_rel = apply(df_effects_cpool_steady_climate_rel, 2, mean)

##

df_effect_cpool_abs = rbind(effects_cpool_climate_change_abs,
                   effects_cpool_steady_climate_abs,
                   grouping2)

df_effect_cpool_rel = rbind(effects_cpool_climate_change_rel,
                            effects_cpool_steady_climate_rel,
                            grouping2)



###### cflux ######

#### climate change ###



mean_effects_Pic_abi_cflux_climate_change = effects_Pic_abi[["cflux"]][["climate_change"]]
main_effects_Pic_abi_cflux_climate_change = mean_effects_Pic_abi_cflux_climate_change[growth_sites_Pic_abi_climate_change,]

main_effects_Pic_abi_cflux_weight_climate_change_abs = apply(abs(main_effects_Pic_abi_cflux_climate_change[,1:38]),2,mean)
main_effects_Pic_abi_cflux_weight_climate_change_abs = main_effects_Pic_abi_cflux_weight_climate_change_abs/sum(abs(main_effects_Pic_abi_cflux_weight_climate_change_abs))


main_effects_Pic_abi_cflux_weight_climate_change_rel = apply(main_effects_Pic_abi_cflux_climate_change[,1:38],2,mean)
main_effects_Pic_abi_cflux_weight_climate_change_rel = main_effects_Pic_abi_cflux_weight_climate_change_rel/sum(abs(main_effects_Pic_abi_cflux_weight_climate_change_rel))



mean_effects_Fag_syl_cflux_climate_change = effects_Fag_syl[["cflux"]][["climate_change"]]
main_effects_Fag_syl_cflux_climate_change = mean_effects_Fag_syl_cflux_climate_change[growth_sites_Fag_syl_climate_change,]

main_effects_Fag_syl_cflux_weight_climate_change_abs = apply(abs(main_effects_Fag_syl_cflux_climate_change[,1:38]),2,mean)
main_effects_Fag_syl_cflux_weight_climate_change_abs = main_effects_Fag_syl_cflux_weight_climate_change_abs/sum(abs(main_effects_Fag_syl_cflux_weight_climate_change_abs))


main_effects_Fag_syl_cflux_weight_climate_change_rel = apply(main_effects_Fag_syl_cflux_climate_change[,1:38],2,mean)
main_effects_Fag_syl_cflux_weight_climate_change_rel = main_effects_Fag_syl_cflux_weight_climate_change_rel/sum(abs(main_effects_Fag_syl_cflux_weight_climate_change_rel))


main_effects_Fag_syl_cflux_weight_climate_change_abs = main_effects_Fag_syl_cflux_weight_climate_change_abs[ordering_fag_syl]
main_effects_Fag_syl_cflux_weight_climate_change_rel = main_effects_Fag_syl_cflux_weight_climate_change_rel[ordering_fag_syl]



mean_effects_Pin_syl_cflux_climate_change = effects_Pin_syl[["cflux"]][["climate_change"]]
main_effects_Pin_syl_cflux_climate_change = mean_effects_Pin_syl_cflux_climate_change[growth_sites_Pin_syl_climate_change,]

main_effects_Pin_syl_cflux_weight_climate_change_abs = apply(abs(main_effects_Pin_syl_cflux_climate_change[,1:38]),2,mean)
main_effects_Pin_syl_cflux_weight_climate_change_abs = main_effects_Pin_syl_cflux_weight_climate_change_abs/sum(abs(main_effects_Pin_syl_cflux_weight_climate_change_abs))


main_effects_Pin_syl_cflux_weight_climate_change_rel = apply(main_effects_Pin_syl_cflux_climate_change[,1:38],2,mean)
main_effects_Pin_syl_cflux_weight_climate_change_rel = main_effects_Pin_syl_cflux_weight_climate_change_rel/sum(abs(main_effects_Pin_syl_cflux_weight_climate_change_rel))


main_effects_Pin_syl_cflux_weight_climate_change_abs = main_effects_Pin_syl_cflux_weight_climate_change_abs[ordering_pin_syl]
main_effects_Pin_syl_cflux_weight_climate_change_rel = main_effects_Pin_syl_cflux_weight_climate_change_rel[ordering_pin_syl]


mean_effects_mixed_cflux_climate_change = effects_mixed[["cflux"]][["climate_change"]]
mean_effects_mixed_cflux_climate_change = mean_effects_mixed_cflux_climate_change[growth_sites_mixed_climate_change,]

remaped_effects_cflux_climate_change_abs = remap_paramters(mean_effects_mixed_cflux_climate_change, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux_weight_climate_change_abs = apply(abs(remaped_effects_cflux_climate_change_abs),2,mean)

main_effects_mixed_cflux_weight_climate_change_abs = main_effects_mixed_cflux_weight_climate_change_abs/sum(abs(main_effects_mixed_cflux_weight_climate_change_abs))

remaped_effects_cflux_climate_change_rel = remap_paramters_rel(mean_effects_mixed_cflux_climate_change, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux_weight_climate_change_rel = apply(remaped_effects_cflux_climate_change_rel,2,mean)

main_effects_mixed_cflux_weight_climate_change_rel = main_effects_mixed_cflux_weight_climate_change_rel/sum(abs(main_effects_mixed_cflux_weight_climate_change_rel))



main_effects_mixed_cflux_weight_climate_change_abs = main_effects_mixed_cflux_weight_climate_change_abs[ordering_mixed]
main_effects_mixed_cflux_weight_climate_change_rel = main_effects_mixed_cflux_weight_climate_change_rel[ordering_mixed]



### steady_climate



mean_effects_Pic_abi_cflux_steady_climate = effects_Pic_abi[["cflux"]][["steady_climate"]]
main_effects_Pic_abi_cflux_steady_climate = mean_effects_Pic_abi_cflux_steady_climate[growth_sites_Pic_abi_steady_climate,]

main_effects_Pic_abi_cflux_weight_steady_climate_abs = apply(abs(main_effects_Pic_abi_cflux_steady_climate[,1:38]),2,mean)
main_effects_Pic_abi_cflux_weight_steady_climate_abs = main_effects_Pic_abi_cflux_weight_steady_climate_abs/sum(abs(main_effects_Pic_abi_cflux_weight_steady_climate_abs))


main_effects_Pic_abi_cflux_weight_steady_climate_rel = apply(main_effects_Pic_abi_cflux_steady_climate[,1:38],2,mean)
main_effects_Pic_abi_cflux_weight_steady_climate_rel = main_effects_Pic_abi_cflux_weight_steady_climate_rel/sum(abs(main_effects_Pic_abi_cflux_weight_steady_climate_rel))



mean_effects_Fag_syl_cflux_steady_climate = effects_Fag_syl[["cflux"]][["steady_climate"]]
main_effects_Fag_syl_cflux_steady_climate = mean_effects_Fag_syl_cflux_steady_climate[growth_sites_Fag_syl_steady_climate,]

main_effects_Fag_syl_cflux_weight_steady_climate_abs = apply(abs(main_effects_Fag_syl_cflux_steady_climate[,1:38]),2,mean)
main_effects_Fag_syl_cflux_weight_steady_climate_abs = main_effects_Fag_syl_cflux_weight_steady_climate_abs/sum(abs(main_effects_Fag_syl_cflux_weight_steady_climate_abs))


main_effects_Fag_syl_cflux_weight_steady_climate_rel = apply(main_effects_Fag_syl_cflux_steady_climate[,1:38],2,mean)
main_effects_Fag_syl_cflux_weight_steady_climate_rel = main_effects_Fag_syl_cflux_weight_steady_climate_rel/sum(abs(main_effects_Fag_syl_cflux_weight_steady_climate_rel))


main_effects_Fag_syl_cflux_weight_steady_climate_abs = main_effects_Fag_syl_cflux_weight_steady_climate_abs[ordering_fag_syl]
main_effects_Fag_syl_cflux_weight_steady_climate_rel = main_effects_Fag_syl_cflux_weight_steady_climate_rel[ordering_fag_syl]



mean_effects_Pin_syl_cflux_steady_climate = effects_Pin_syl[["cflux"]][["steady_climate"]]
main_effects_Pin_syl_cflux_steady_climate = mean_effects_Pin_syl_cflux_steady_climate[growth_sites_Pin_syl_steady_climate,]

main_effects_Pin_syl_cflux_weight_steady_climate_abs = apply(abs(main_effects_Pin_syl_cflux_steady_climate[,1:38]),2,mean)
main_effects_Pin_syl_cflux_weight_steady_climate_abs = main_effects_Pin_syl_cflux_weight_steady_climate_abs/sum(abs(main_effects_Pin_syl_cflux_weight_steady_climate_abs))


main_effects_Pin_syl_cflux_weight_steady_climate_rel = apply(main_effects_Pin_syl_cflux_steady_climate[,1:38],2,mean)
main_effects_Pin_syl_cflux_weight_steady_climate_rel = main_effects_Pin_syl_cflux_weight_steady_climate_rel/sum(abs(main_effects_Pin_syl_cflux_weight_steady_climate_rel))


main_effects_Pin_syl_cflux_weight_steady_climate_abs = main_effects_Pin_syl_cflux_weight_steady_climate_abs[ordering_pin_syl]
main_effects_Pin_syl_cflux_weight_steady_climate_rel = main_effects_Pin_syl_cflux_weight_steady_climate_rel[ordering_pin_syl]


mean_effects_mixed_cflux_steady_climate = effects_mixed[["cflux"]][["steady_climate"]]
mean_effects_mixed_cflux_steady_climate = mean_effects_mixed_cflux_steady_climate[growth_sites_mixed_steady_climate,]

remaped_effects_cflux_steady_climate_abs = remap_paramters(mean_effects_mixed_cflux_steady_climate, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux_weight_steady_climate_abs = apply(abs(remaped_effects_cflux_steady_climate_abs),2,mean)

main_effects_mixed_cflux_weight_steady_climate_abs = main_effects_mixed_cflux_weight_steady_climate_abs/sum(abs(main_effects_mixed_cflux_weight_steady_climate_abs))

remaped_effects_cflux_steady_climate_rel = remap_paramters_rel(mean_effects_mixed_cflux_steady_climate, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux_weight_steady_climate_rel = apply(remaped_effects_cflux_steady_climate_rel,2,mean)

main_effects_mixed_cflux_weight_steady_climate_rel = main_effects_mixed_cflux_weight_steady_climate_rel/sum(abs(main_effects_mixed_cflux_weight_steady_climate_rel))



main_effects_mixed_cflux_weight_steady_climate_abs = main_effects_mixed_cflux_weight_steady_climate_abs[ordering_mixed]
main_effects_mixed_cflux_weight_steady_climate_rel = main_effects_mixed_cflux_weight_steady_climate_rel[ordering_mixed]

####
df_effects_cflux_climate_change_abs = rbind(main_effects_Fag_syl_cflux_weight_climate_change_abs,
                                            main_effects_Pin_syl_cflux_weight_climate_change_abs,
                                            main_effects_Pic_abi_cflux_weight_climate_change_abs,
                                            main_effects_mixed_cflux_weight_climate_change_abs)

df_effects_cflux_climate_change_rel = rbind(main_effects_Fag_syl_cflux_weight_climate_change_rel,
                                            main_effects_Pin_syl_cflux_weight_climate_change_rel,
                                            main_effects_Pic_abi_cflux_weight_climate_change_rel,
                                            main_effects_mixed_cflux_weight_climate_change_rel)

effects_cflux_climate_change_abs = apply(df_effects_cflux_climate_change_abs, 2, mean)
effects_cflux_climate_change_rel = apply(df_effects_cflux_climate_change_rel, 2, mean)

df_effects_cflux_steady_climate_abs = rbind(main_effects_Fag_syl_cflux_weight_steady_climate_abs,
                                            main_effects_Pin_syl_cflux_weight_steady_climate_abs,
                                            main_effects_Pic_abi_cflux_weight_steady_climate_abs,
                                            main_effects_mixed_cflux_weight_steady_climate_abs)

df_effects_cflux_steady_climate_rel = rbind(main_effects_Fag_syl_cflux_weight_steady_climate_rel,
                                            main_effects_Pin_syl_cflux_weight_steady_climate_rel,
                                            main_effects_Pic_abi_cflux_weight_steady_climate_rel,
                                            main_effects_mixed_cflux_weight_steady_climate_rel)

effects_cflux_steady_climate_abs = apply(df_effects_cflux_steady_climate_abs, 2, mean)
effects_cflux_steady_climate_rel = apply(df_effects_cflux_steady_climate_rel, 2, mean)

##

df_effect_cflux_abs = rbind(effects_cflux_climate_change_abs,
                            effects_cflux_steady_climate_abs,
                            grouping2)

df_effect_cflux_rel = rbind(effects_cflux_climate_change_rel,
                            effects_cflux_steady_climate_rel,
                            grouping2)


###### agpp ######

#### climate change ###



mean_effects_Pic_abi_agpp_climate_change = effects_Pic_abi[["agpp"]][["climate_change"]]
main_effects_Pic_abi_agpp_climate_change = mean_effects_Pic_abi_agpp_climate_change[growth_sites_Pic_abi_climate_change,]

main_effects_Pic_abi_agpp_weight_climate_change_abs = apply(abs(main_effects_Pic_abi_agpp_climate_change[,1:38]),2,mean)
main_effects_Pic_abi_agpp_weight_climate_change_abs = main_effects_Pic_abi_agpp_weight_climate_change_abs/sum(abs(main_effects_Pic_abi_agpp_weight_climate_change_abs))


main_effects_Pic_abi_agpp_weight_climate_change_rel = apply(main_effects_Pic_abi_agpp_climate_change[,1:38],2,mean)
main_effects_Pic_abi_agpp_weight_climate_change_rel = main_effects_Pic_abi_agpp_weight_climate_change_rel/sum(abs(main_effects_Pic_abi_agpp_weight_climate_change_rel))



mean_effects_Fag_syl_agpp_climate_change = effects_Fag_syl[["agpp"]][["climate_change"]]
main_effects_Fag_syl_agpp_climate_change = mean_effects_Fag_syl_agpp_climate_change[growth_sites_Fag_syl_climate_change,]

main_effects_Fag_syl_agpp_weight_climate_change_abs = apply(abs(main_effects_Fag_syl_agpp_climate_change[,1:38]),2,mean)
main_effects_Fag_syl_agpp_weight_climate_change_abs = main_effects_Fag_syl_agpp_weight_climate_change_abs/sum(abs(main_effects_Fag_syl_agpp_weight_climate_change_abs))


main_effects_Fag_syl_agpp_weight_climate_change_rel = apply(main_effects_Fag_syl_agpp_climate_change[,1:38],2,mean)
main_effects_Fag_syl_agpp_weight_climate_change_rel = main_effects_Fag_syl_agpp_weight_climate_change_rel/sum(abs(main_effects_Fag_syl_agpp_weight_climate_change_rel))


main_effects_Fag_syl_agpp_weight_climate_change_abs = main_effects_Fag_syl_agpp_weight_climate_change_abs[ordering_fag_syl]
main_effects_Fag_syl_agpp_weight_climate_change_rel = main_effects_Fag_syl_agpp_weight_climate_change_rel[ordering_fag_syl]



mean_effects_Pin_syl_agpp_climate_change = effects_Pin_syl[["agpp"]][["climate_change"]]
main_effects_Pin_syl_agpp_climate_change = mean_effects_Pin_syl_agpp_climate_change[growth_sites_Pin_syl_climate_change,]

main_effects_Pin_syl_agpp_weight_climate_change_abs = apply(abs(main_effects_Pin_syl_agpp_climate_change[,1:38]),2,mean)
main_effects_Pin_syl_agpp_weight_climate_change_abs = main_effects_Pin_syl_agpp_weight_climate_change_abs/sum(abs(main_effects_Pin_syl_agpp_weight_climate_change_abs))


main_effects_Pin_syl_agpp_weight_climate_change_rel = apply(main_effects_Pin_syl_agpp_climate_change[,1:38],2,mean)
main_effects_Pin_syl_agpp_weight_climate_change_rel = main_effects_Pin_syl_agpp_weight_climate_change_rel/sum(abs(main_effects_Pin_syl_agpp_weight_climate_change_rel))


main_effects_Pin_syl_agpp_weight_climate_change_abs = main_effects_Pin_syl_agpp_weight_climate_change_abs[ordering_pin_syl]
main_effects_Pin_syl_agpp_weight_climate_change_rel = main_effects_Pin_syl_agpp_weight_climate_change_rel[ordering_pin_syl]


mean_effects_mixed_agpp_climate_change = effects_mixed[["agpp"]][["climate_change"]]
mean_effects_mixed_agpp_climate_change = mean_effects_mixed_agpp_climate_change[growth_sites_mixed_climate_change,]

remaped_effects_agpp_climate_change_abs = remap_paramters(mean_effects_mixed_agpp_climate_change, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp_weight_climate_change_abs = apply(abs(remaped_effects_agpp_climate_change_abs),2,mean)

main_effects_mixed_agpp_weight_climate_change_abs = main_effects_mixed_agpp_weight_climate_change_abs/sum(abs(main_effects_mixed_agpp_weight_climate_change_abs))

remaped_effects_agpp_climate_change_rel = remap_paramters_rel(mean_effects_mixed_agpp_climate_change, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp_weight_climate_change_rel = apply(remaped_effects_agpp_climate_change_rel,2,mean)

main_effects_mixed_agpp_weight_climate_change_rel = main_effects_mixed_agpp_weight_climate_change_rel/sum(abs(main_effects_mixed_agpp_weight_climate_change_rel))



main_effects_mixed_agpp_weight_climate_change_abs = main_effects_mixed_agpp_weight_climate_change_abs[ordering_mixed]
main_effects_mixed_agpp_weight_climate_change_rel = main_effects_mixed_agpp_weight_climate_change_rel[ordering_mixed]



### steady_climate



mean_effects_Pic_abi_agpp_steady_climate = effects_Pic_abi[["agpp"]][["steady_climate"]]
main_effects_Pic_abi_agpp_steady_climate = mean_effects_Pic_abi_agpp_steady_climate[growth_sites_Pic_abi_steady_climate,]

main_effects_Pic_abi_agpp_weight_steady_climate_abs = apply(abs(main_effects_Pic_abi_agpp_steady_climate[,1:38]),2,mean)
main_effects_Pic_abi_agpp_weight_steady_climate_abs = main_effects_Pic_abi_agpp_weight_steady_climate_abs/sum(abs(main_effects_Pic_abi_agpp_weight_steady_climate_abs))


main_effects_Pic_abi_agpp_weight_steady_climate_rel = apply(main_effects_Pic_abi_agpp_steady_climate[,1:38],2,mean)
main_effects_Pic_abi_agpp_weight_steady_climate_rel = main_effects_Pic_abi_agpp_weight_steady_climate_rel/sum(abs(main_effects_Pic_abi_agpp_weight_steady_climate_rel))



mean_effects_Fag_syl_agpp_steady_climate = effects_Fag_syl[["agpp"]][["steady_climate"]]
main_effects_Fag_syl_agpp_steady_climate = mean_effects_Fag_syl_agpp_steady_climate[growth_sites_Fag_syl_steady_climate,]

main_effects_Fag_syl_agpp_weight_steady_climate_abs = apply(abs(main_effects_Fag_syl_agpp_steady_climate[,1:38]),2,mean)
main_effects_Fag_syl_agpp_weight_steady_climate_abs = main_effects_Fag_syl_agpp_weight_steady_climate_abs/sum(abs(main_effects_Fag_syl_agpp_weight_steady_climate_abs))


main_effects_Fag_syl_agpp_weight_steady_climate_rel = apply(main_effects_Fag_syl_agpp_steady_climate[,1:38],2,mean)
main_effects_Fag_syl_agpp_weight_steady_climate_rel = main_effects_Fag_syl_agpp_weight_steady_climate_rel/sum(abs(main_effects_Fag_syl_agpp_weight_steady_climate_rel))


main_effects_Fag_syl_agpp_weight_steady_climate_abs = main_effects_Fag_syl_agpp_weight_steady_climate_abs[ordering_fag_syl]
main_effects_Fag_syl_agpp_weight_steady_climate_rel = main_effects_Fag_syl_agpp_weight_steady_climate_rel[ordering_fag_syl]



mean_effects_Pin_syl_agpp_steady_climate = effects_Pin_syl[["agpp"]][["steady_climate"]]
main_effects_Pin_syl_agpp_steady_climate = mean_effects_Pin_syl_agpp_steady_climate[growth_sites_Pin_syl_steady_climate,]

main_effects_Pin_syl_agpp_weight_steady_climate_abs = apply(abs(main_effects_Pin_syl_agpp_steady_climate[,1:38]),2,mean)
main_effects_Pin_syl_agpp_weight_steady_climate_abs = main_effects_Pin_syl_agpp_weight_steady_climate_abs/sum(abs(main_effects_Pin_syl_agpp_weight_steady_climate_abs))


main_effects_Pin_syl_agpp_weight_steady_climate_rel = apply(main_effects_Pin_syl_agpp_steady_climate[,1:38],2,mean)
main_effects_Pin_syl_agpp_weight_steady_climate_rel = main_effects_Pin_syl_agpp_weight_steady_climate_rel/sum(abs(main_effects_Pin_syl_agpp_weight_steady_climate_rel))


main_effects_Pin_syl_agpp_weight_steady_climate_abs = main_effects_Pin_syl_agpp_weight_steady_climate_abs[ordering_pin_syl]
main_effects_Pin_syl_agpp_weight_steady_climate_rel = main_effects_Pin_syl_agpp_weight_steady_climate_rel[ordering_pin_syl]


mean_effects_mixed_agpp_steady_climate = effects_mixed[["agpp"]][["steady_climate"]]
mean_effects_mixed_agpp_steady_climate = mean_effects_mixed_agpp_steady_climate[growth_sites_mixed_steady_climate,]

remaped_effects_agpp_steady_climate_abs = remap_paramters(mean_effects_mixed_agpp_steady_climate, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp_weight_steady_climate_abs = apply(abs(remaped_effects_agpp_steady_climate_abs),2,mean)

main_effects_mixed_agpp_weight_steady_climate_abs = main_effects_mixed_agpp_weight_steady_climate_abs/sum(abs(main_effects_mixed_agpp_weight_steady_climate_abs))

remaped_effects_agpp_steady_climate_rel = remap_paramters_rel(mean_effects_mixed_agpp_steady_climate, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp_weight_steady_climate_rel = apply(remaped_effects_agpp_steady_climate_rel,2,mean)

main_effects_mixed_agpp_weight_steady_climate_rel = main_effects_mixed_agpp_weight_steady_climate_rel/sum(abs(main_effects_mixed_agpp_weight_steady_climate_rel))



main_effects_mixed_agpp_weight_steady_climate_abs = main_effects_mixed_agpp_weight_steady_climate_abs[ordering_mixed]
main_effects_mixed_agpp_weight_steady_climate_rel = main_effects_mixed_agpp_weight_steady_climate_rel[ordering_mixed]

####
df_effects_agpp_climate_change_abs = rbind(main_effects_Fag_syl_agpp_weight_climate_change_abs,
                                            main_effects_Pin_syl_agpp_weight_climate_change_abs,
                                            main_effects_Pic_abi_agpp_weight_climate_change_abs,
                                            main_effects_mixed_agpp_weight_climate_change_abs)

df_effects_agpp_climate_change_rel = rbind(main_effects_Fag_syl_agpp_weight_climate_change_rel,
                                            main_effects_Pin_syl_agpp_weight_climate_change_rel,
                                            main_effects_Pic_abi_agpp_weight_climate_change_rel,
                                            main_effects_mixed_agpp_weight_climate_change_rel)

effects_agpp_climate_change_abs = apply(df_effects_agpp_climate_change_abs, 2, mean)
effects_agpp_climate_change_rel = apply(df_effects_agpp_climate_change_rel, 2, mean)

df_effects_agpp_steady_climate_abs = rbind(main_effects_Fag_syl_agpp_weight_steady_climate_abs,
                                            main_effects_Pin_syl_agpp_weight_steady_climate_abs,
                                            main_effects_Pic_abi_agpp_weight_steady_climate_abs,
                                            main_effects_mixed_agpp_weight_steady_climate_abs)

df_effects_agpp_steady_climate_rel = rbind(main_effects_Fag_syl_agpp_weight_steady_climate_rel,
                                            main_effects_Pin_syl_agpp_weight_steady_climate_rel,
                                            main_effects_Pic_abi_agpp_weight_steady_climate_rel,
                                            main_effects_mixed_agpp_weight_steady_climate_rel)

effects_agpp_steady_climate_abs = apply(df_effects_agpp_steady_climate_abs, 2, mean)
effects_agpp_steady_climate_rel = apply(df_effects_agpp_steady_climate_rel, 2, mean)

##

df_effect_agpp_abs = rbind(effects_agpp_climate_change_abs,
                            effects_agpp_steady_climate_abs,
                            grouping2)

df_effect_agpp_rel = rbind(effects_agpp_climate_change_rel,
                            effects_agpp_steady_climate_rel,
                            grouping2)



#### Analysis for Pic abi ####


mean_effects_Pic_abi_cmass = effects_Pic_abi[["cmass"]][["complete"]]

mean_effects_Pic_abi_cmass_scaled = t(apply(mean_effects_Pic_abi_cmass[,1:38],1,rescaling_abs))

main_effects_Pic_abi_cmass_weight = apply(mean_effects_Pic_abi_cmass_scaled,2,weighted.mean,
                                          w = mean_effects_Pic_abi_cmass[,39])



mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
mean_effects_Pic_abi_cpool_scaled = t(apply(mean_effects_Pic_abi_cpool[,1:38],1,rescaling_abs))

main_effects_Pic_abi_cpool_weight = apply(mean_effects_Pic_abi_cpool_scaled,2,weighted.mean,
                                          w = abs(mean_effects_Pic_abi_cpool[,39]))

mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]]
mean_effects_Pic_abi_cflux_scaled = t(apply(mean_effects_Pic_abi_cflux[,1:38],1,rescaling_abs))

main_effects_Pic_abi_cflux_weight = apply(mean_effects_Pic_abi_cflux_scaled,2,weighted.mean,
                                          w = abs(mean_effects_Pic_abi_cflux[,39]))

mean_effects_Pic_abi_agpp = effects_Pic_abi[["agpp"]][["complete"]]
mean_effects_Pic_abi_agpp_scaled = t(apply(mean_effects_Pic_abi_agpp[,1:38],1,rescaling_abs))

main_effects_Pic_abi_agpp_weight = apply(mean_effects_Pic_abi_agpp_scaled,2,weighted.mean,
                                         w = mean_effects_Pic_abi_agpp[,39])


#### Analysis for Fag syl ####

mean_effects_Fag_syl_cmass = effects_Fag_syl[["cmass"]][["complete"]]
mean_effects_Fag_syl_cmass_scaled = t(apply(mean_effects_Fag_syl_cmass[,1:38],1,rescaling_abs))

main_effects_Fag_syl_cmass_weight = apply(mean_effects_Fag_syl_cmass_scaled,2,weighted.mean,
                                          w = mean_effects_Fag_syl_cmass[,39])



mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
mean_effects_Fag_syl_cpool_scaled = t(apply(mean_effects_Fag_syl_cpool[,1:38],1,rescaling_abs))

main_effects_Fag_syl_cpool_weight = apply(mean_effects_Fag_syl_cpool_scaled,2,weighted.mean,
                                          w = abs(mean_effects_Fag_syl_cpool[,39]))

mean_effects_Fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]]
mean_effects_Fag_syl_cflux_scaled = t(apply(mean_effects_Fag_syl_cflux[,1:38],1,rescaling_abs))

main_effects_Fag_syl_cflux_weight = apply(mean_effects_Fag_syl_cflux_scaled,2,weighted.mean,
                                          w = abs(mean_effects_Fag_syl_cflux[,39]))

mean_effects_Fag_syl_agpp = effects_Fag_syl[["agpp"]][["complete"]]
mean_effects_Fag_syl_agpp_scaled = t(apply(mean_effects_Fag_syl_agpp[,1:38],1,rescaling_abs))

main_effects_Fag_syl_agpp_weight = apply(mean_effects_Fag_syl_agpp_scaled,2,weighted.mean,
                                         w = mean_effects_Fag_syl_agpp[,39])


#### Analysis for Pic abi ####


mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
mean_effects_Pic_abi_cpool = mean_effects_Pic_abi_cpool[growth_sites_Pic_abi,]

main_effects_Pic_abi_cpool = apply(abs(mean_effects_Pic_abi_cpool[,1:38]),2,mean)
main_effects_Pic_abi_cpool_weight_abs = main_effects_Pic_abi_cpool/sum(abs(main_effects_Pic_abi_cpool))


main_effects_Pic_abi_cpool = apply(mean_effects_Pic_abi_cpool[,1:38],2,mean)
main_effects_Pic_abi_cpool_weight = main_effects_Pic_abi_cpool/sum(abs(main_effects_Pic_abi_cpool))


mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]]
mean_effects_Pic_abi_cflux = mean_effects_Pic_abi_cflux[growth_sites_Pic_abi,]


main_effects_Pic_abi_cflux = apply(abs(mean_effects_Pic_abi_cflux[,1:38]),2,mean)
main_effects_Pic_abi_cflux_weight_abs = main_effects_Pic_abi_cflux/sum(abs(main_effects_Pic_abi_cflux))
sd_Pic_abic_cmass_abs = apply(abs(mean_effects_Pic_abi_cmass[,1:38]),2,sd)

main_effects_Pic_abi_cflux = apply(mean_effects_Pic_abi_cflux[,1:38],2,mean)
main_effects_Pic_abi_cflux_weight = main_effects_Pic_abi_cflux/sum(abs(main_effects_Pic_abi_cflux))


mean_effects_Pic_abi_agpp = effects_Pic_abi[["agpp"]][["complete"]]
mean_effects_Pic_abi_agpp = mean_effects_Pic_abi_agpp[growth_sites_Pic_abi,]

main_effects_Pic_abi_agpp = apply(abs(mean_effects_Pic_abi_agpp[,1:38]),2,mean)
main_effects_Pic_abi_agpp_weight_abs = main_effects_Pic_abi_agpp/sum(abs(main_effects_Pic_abi_agpp))

main_effects_Pic_abi_agpp = apply(mean_effects_Pic_abi_agpp[,1:38],2,mean)
main_effects_Pic_abi_agpp_weight = main_effects_Pic_abi_agpp/sum(abs(main_effects_Pic_abi_agpp))

## sd estimates ##

sd_Pic_abi_cpool_abs = apply((mean_effects_Pic_abi_cpool[,1:38]/mean_effects_Pic_abi_cpool[,39]),2,sd,na.rm = T)
sd_Pic_abi_cflux_abs = apply((mean_effects_Pic_abi_cflux[,1:38]/mean_effects_Pic_abi_cflux[,39]),2,sd, na.rm = T)
sd_Pic_abi_agpp_abs = apply((mean_effects_Pic_abi_agpp[,1:38]/mean_effects_Pic_abi_agpp[,39]),2,sd, na.rm = T)

#### Analysis for Fag syl ####



mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
growth_sites_Fag_syl = which(mean_effects_Fag_syl_cpool[,39]/200 > 2)
mean_effects_Fag_syl_cpool = mean_effects_Fag_syl_cpool[growth_sites_Fag_syl,]

main_effects_Fag_syl_cpool = apply(abs(mean_effects_Fag_syl_cpool[,1:38]),2,mean)
main_effects_Fag_syl_cpool_weight_abs = main_effects_Fag_syl_cpool/sum(abs(main_effects_Fag_syl_cpool))

main_effects_Fag_syl_cpool = apply(mean_effects_Fag_syl_cpool[,1:38],2,mean)
main_effects_Fag_syl_cpool_weight = main_effects_Fag_syl_cpool/sum(abs(main_effects_Fag_syl_cpool))

mean_effects_Fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]]
mean_effects_Fag_syl_cflux = mean_effects_Fag_syl_cflux[growth_sites_Fag_syl,]

main_effects_Fag_syl_cflux = apply(abs(mean_effects_Fag_syl_cflux[,1:38]),2,mean)
main_effects_Fag_syl_cflux_weight_abs = main_effects_Fag_syl_cflux/sum(abs(main_effects_Fag_syl_cflux))

main_effects_Fag_syl_cflux = apply(mean_effects_Fag_syl_cflux[,1:38],2,mean)
main_effects_Fag_syl_cflux_weight = main_effects_Fag_syl_cflux/sum(abs(main_effects_Fag_syl_cflux))


mean_effects_Fag_syl_agpp = effects_Fag_syl[["agpp"]][["complete"]]
mean_effects_Fag_syl_agpp = mean_effects_Fag_syl_agpp[growth_sites_Fag_syl,]

main_effects_Fag_syl_agpp = apply(abs(mean_effects_Fag_syl_agpp[,1:38]),2,mean)
main_effects_Fag_syl_agpp_weight_abs = main_effects_Fag_syl_agpp/sum(abs(main_effects_Fag_syl_agpp))

main_effects_Fag_syl_agpp = apply(mean_effects_Fag_syl_agpp[,1:38],2,mean)
main_effects_Fag_syl_agpp_weight = main_effects_Fag_syl_agpp/sum(abs(main_effects_Fag_syl_agpp))

## sd estimates ##

sd_Fag_syl_cpool_abs = apply(abs(mean_effects_Fag_syl_cpool[,1:38]/mean_effects_Fag_syl_cpool[,39]),2,sd, na.rm = T)
sd_Fag_syl_cflux_abs = apply(abs(mean_effects_Fag_syl_cflux[,1:38]/mean_effects_Fag_syl_cflux[,39]),2,sd, na.rm = T)
sd_Fag_syl_agpp_abs = apply(abs(mean_effects_Fag_syl_agpp[,1:38]/mean_effects_Fag_syl_agpp[,39]),2,sd, na.rm = T)

#### Analysis for Pin syl ####



mean_effects_Pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
growth_sites_Pin_syl = which(mean_effects_Pin_syl_cpool[,39]/200 > 2)
mean_effects_Pin_syl_cpool = mean_effects_Pin_syl_cpool[growth_sites_Pin_syl,]

main_effects_Pin_syl_cpool = apply(abs(mean_effects_Pin_syl_cpool[,1:38]),2,mean)
main_effects_Pin_syl_cpool_weight_abs = main_effects_Pin_syl_cpool/sum(abs(main_effects_Pin_syl_cpool))

main_effects_Pin_syl_cpool = apply(mean_effects_Pin_syl_cpool[,1:38],2,mean)
main_effects_Pin_syl_cpool_weight = main_effects_Pin_syl_cpool/sum(abs(main_effects_Pin_syl_cpool))

mean_effects_Pin_syl_cflux = effects_Pin_syl[["cflux"]][["complete"]]
mean_effects_Pin_syl_cflux = mean_effects_Pin_syl_cflux[growth_sites_Pin_syl,]

main_effects_Pin_syl_cflux = apply(abs(mean_effects_Pin_syl_cflux[,1:38]),2,mean)
main_effects_Pin_syl_cflux_weight_abs = main_effects_Pin_syl_cflux/sum(abs(main_effects_Pin_syl_cflux))

main_effects_Pin_syl_cflux = apply(mean_effects_Pin_syl_cflux[,1:38],2,mean)
main_effects_Pin_syl_cflux_weight = main_effects_Pin_syl_cflux/sum(abs(main_effects_Pin_syl_cflux))


mean_effects_Pin_syl_agpp = effects_Pin_syl[["agpp"]][["complete"]]
mean_effects_Pin_syl_agpp = mean_effects_Pin_syl_agpp[growth_sites_Pin_syl,]


main_effects_Pin_syl_agpp = apply(abs(mean_effects_Pin_syl_agpp[,1:38]),2,mean)
main_effects_Pin_syl_agpp_weight_abs = main_effects_Pin_syl_agpp/sum(abs(main_effects_Pin_syl_agpp))

main_effects_Pin_syl_agpp = apply(mean_effects_Pin_syl_agpp[,1:38],2,mean)
main_effects_Pin_syl_agpp_weight = main_effects_Pin_syl_agpp/sum(abs(main_effects_Pin_syl_agpp))

## sd estimates ##

sd_Pin_syl_cpool_abs = apply((mean_effects_Pin_syl_cpool[,1:38]/mean_effects_Pin_syl_cpool[,39]),2,sd, na.rm = T)
sd_Pin_syl_cflux_abs = apply((mean_effects_Pin_syl_cflux[,1:38]/mean_effects_Pin_syl_cflux[,39]),2,sd, na.rm = T)
sd_Pin_syl_agpp_abs = apply((mean_effects_Pin_syl_agpp[,1:38]/mean_effects_Pin_syl_agpp[,39]),2,sd, na.rm = T)

#### Analysis for mixed ####


mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
mean_effects_mixed_cpool = mean_effects_mixed_cpool[growth_sites_mixed,]

remaped_effects_cpool_abs = remap_paramters(mean_effects_mixed_cpool, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool = apply(abs(remaped_effects_cpool_abs),2,mean)

main_effects_mixed_cpool_weight_abs = main_effects_mixed_cpool/sum(abs(main_effects_mixed_cpool))

remaped_effects_cpool = remap_paramters_rel(mean_effects_mixed_cpool, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool = apply(remaped_effects_cpool,2,mean)

main_effects_mixed_cpool_weight = main_effects_mixed_cpool/sum(abs(main_effects_mixed_cpool))


mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]]
mean_effects_mixed_cflux = mean_effects_mixed_cflux[growth_sites_mixed,]

remaped_effects_cflux_abs = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux = apply(abs(remaped_effects_cflux_abs),2,mean)

main_effects_mixed_cflux_weight_abs = main_effects_mixed_cflux/sum(abs(main_effects_mixed_cflux))

remaped_effects_cflux = remap_paramters_rel(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux = apply(remaped_effects_cflux,2,mean)

main_effects_mixed_cflux_weight = main_effects_mixed_cflux/sum(abs(main_effects_mixed_cflux))


mean_effects_mixed_agpp = effects_mixed[["agpp"]][["complete"]]
mean_effects_mixed_agpp = mean_effects_mixed_agpp[growth_sites_mixed,]

remaped_effects_agpp_abs = remap_paramters(mean_effects_mixed_agpp, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp = apply(abs(remaped_effects_agpp_abs),2,mean)

main_effects_mixed_agpp_weight_abs = main_effects_mixed_agpp/sum(abs(main_effects_mixed_agpp))

remaped_effects_agpp = remap_paramters_rel(mean_effects_mixed_agpp, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp = apply(remaped_effects_agpp,2,mean)

main_effects_mixed_agpp_weight = main_effects_mixed_agpp/sum(abs(main_effects_mixed_agpp))


#### reordering of effects such that all are the same ####

## reordering relative effects

main_effects_Pin_syl_cpool_weight = main_effects_Pin_syl_cpool_weight[ordering_pin_syl]
main_effects_Pin_syl_cflux_weight = main_effects_Pin_syl_cflux_weight[ordering_pin_syl]
main_effects_Pin_syl_agpp_weight = main_effects_Pin_syl_agpp_weight[ordering_pin_syl]

main_effects_Fag_syl_cpool_weight = main_effects_Fag_syl_cpool_weight[ordering_fag_syl]
main_effects_Fag_syl_cflux_weight = main_effects_Fag_syl_cflux_weight[ordering_fag_syl]
main_effects_Fag_syl_agpp_weight = main_effects_Fag_syl_agpp_weight[ordering_fag_syl]

main_effects_mixed_cpool_weight = main_effects_mixed_cpool_weight[ordering_mixed]
main_effects_mixed_cflux_weight = main_effects_mixed_cflux_weight[ordering_mixed]
main_effects_mixed_agpp_weight = main_effects_mixed_agpp_weight[ordering_mixed]

## reordering absolute effects

main_effects_Pin_syl_cpool_weight_abs = main_effects_Pin_syl_cpool_weight_abs[ordering_pin_syl]
main_effects_Pin_syl_cflux_weight_abs = main_effects_Pin_syl_cflux_weight_abs[ordering_pin_syl]
main_effects_Pin_syl_agpp_weight_abs = main_effects_Pin_syl_agpp_weight_abs[ordering_pin_syl]

main_effects_Fag_syl_cpool_weight_abs = main_effects_Fag_syl_cpool_weight_abs[ordering_fag_syl]
main_effects_Fag_syl_cflux_weight_abs = main_effects_Fag_syl_cflux_weight_abs[ordering_fag_syl]
main_effects_Fag_syl_agpp_weight_abs = main_effects_Fag_syl_agpp_weight_abs[ordering_fag_syl]

main_effects_mixed_cpool_weight_abs = main_effects_mixed_cpool_weight_abs[ordering_mixed]
main_effects_mixed_cflux_weight_abs = main_effects_mixed_cflux_weight_abs[ordering_mixed]
main_effects_mixed_agpp_weight_abs = main_effects_mixed_agpp_weight_abs[ordering_mixed]

df_mono_agpp_abs = rbind(main_effects_Fag_syl_agpp_weight_abs,
                          main_effects_Pin_syl_agpp_weight_abs,
                          main_effects_Pic_abi_agpp_weight_abs)

df_mono_cpool_abs = rbind(main_effects_Fag_syl_cpool_weight_abs,
                          main_effects_Pin_syl_cpool_weight_abs,
                          main_effects_Pic_abi_cpool_weight_abs)

df_mono_cflux_abs = rbind(main_effects_Fag_syl_cflux_weight_abs,
                      main_effects_Pin_syl_cflux_weight_abs,
                      main_effects_Pic_abi_cflux_weight_abs)

df_mono_agpp_rel = rbind(main_effects_Fag_syl_agpp_weight,
                         main_effects_Pin_syl_agpp_weight,
                         main_effects_Pic_abi_agpp_weight)

df_mono_cpool_rel = rbind(main_effects_Fag_syl_cpool_weight,
                          main_effects_Pin_syl_cpool_weight,
                          main_effects_Pic_abi_cpool_weight)

df_mono_cflux_rel = rbind(main_effects_Fag_syl_cflux_weight,
                          main_effects_Pin_syl_cflux_weight,
                          main_effects_Pic_abi_cflux_weight)


effects_mono_agpp_abs = apply(df_mono_agpp_abs, 2, mean)
effects_mono_cflux_abs = apply(df_mono_cflux_abs, 2, mean)
effects_mono_cpool_abs = apply(df_mono_cpool_abs, 2, mean)



effects_mono_agpp_rel = apply(df_mono_agpp_rel, 2, mean)
effects_mono_cflux_rel = apply(df_mono_cflux_rel, 2, mean)
effects_mono_cpool_rel = apply(df_mono_cpool_rel, 2, mean)

df_effect_comp_agpp_abs = rbind(effects_mono_agpp_abs,
                       main_effects_mixed_agpp_weight_abs,
                       grouping2)

df_effect_comp_cpool_abs = rbind(effects_mono_cpool_abs,
                            main_effects_mixed_cpool_weight_abs,
                            grouping2)

df_effect_comp_cflux_abs = rbind(effects_mono_cflux_abs,
                            main_effects_mixed_cflux_weight_abs,
                            grouping2)

df_effect_comp_agpp_rel = rbind(effects_mono_agpp_rel,
                                main_effects_mixed_agpp_weight,
                                grouping2)

df_effect_comp_cpool_rel = rbind(effects_mono_cpool_rel,
                                 main_effects_mixed_cpool_weight,
                                 grouping2)

df_effect_comp_cflux_rel = rbind(effects_mono_cflux_rel,
                             main_effects_mixed_cflux_weight,
                             grouping2)

#### Plotting the stuff all together ####


library(RColorBrewer)


graphics.off()
par("mar")
par(mar=c(1,1,1,1))
dev.off()
pdf("./Figures/Comparison_future_abs.pdf",width=14,height=10)
layout(matrix(c(1,3,3,
                2,3,3), nrow = 2,ncol = 3, byrow = TRUE),
       widths = c(1,1,1,0.32),
       heights =  c(1,1,1,1))

colors_me = c("gold","burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")

#plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
# text(paste("I) Comparison climate change - stable climate"),
#      cex = 1.6, col = "black", font = 2, srt = 90, x = 0.5,y = 0.5)
plot(x = df_effect_cpool_abs[1,], y = df_effect_cpool_abs[2,],log = "xy" ,
     type = "p", col = colors_me[as.numeric(as.factor(df_effect_cpool_abs[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(0.01, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(0.01, 0.2), las = 1, cex = 1.5)
axis(side = 2, at = c(0.01,seq(0.1,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
axis(side = 1, at = c(0.01,seq(0.1,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
segments(0.01,0.01,0.2,0.2)
title(main = "a)",  adj = 0.,cex.main =2)
title(xlab = "Uncertainty contributions in a changing climate", line = 2, cex.lab = 1.3)
title(ylab = "Uncertainty contributions in a stable future climate", line = 2, cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(df_effect_cpool_abs[1,]> 0.08)],
                        parameternames_ordered_pic_abi2[which(df_effect_cpool_abs[2,]>0.08)])
par(xpd =T)
text(x = as.numeric(df_effect_cpool_abs[1,][union(which(df_effect_cpool_abs[1,]> 0.08),which(df_effect_cpool_abs[2,]> 0.08))]),
     y=  as.numeric(df_effect_cpool_abs[2,][union(which(df_effect_cpool_abs[1,]> 0.08),which(df_effect_cpool_abs[2,]> 0.08))]),
     labels=names_important, cex = 1.5, pos=3, offset = 1.5)
par(xpd =F)

plot(x = df_effect_agpp_abs[1,], y = df_effect_agpp_abs[2,], type = "p", log = "xy",
     col = colors_me[as.numeric(as.factor(df_effect_agpp_abs[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(0.01, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(0.01, 0.2), las = 1, cex = 1.5)
axis(side = 2, at = c(0.01,seq(0.05,0.2, by = 0.05)), labels = T, pos = 0.01, las = 1)
axis(side = 1, at = c(0.01,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
segments(0.01,0.01,0.2,0.2)

title(main = "b) ",adj = 0., cex.main =2)
title(xlab =  "Uncertainty contributions in a changing climate", line = 2,cex.lab = 1.3)
title(ylab =  "Uncertainty contributions in a stable future climate", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(df_effect_agpp_abs[1,]> 0.1)],
                        parameternames_ordered_pic_abi2[which(df_effect_agpp_abs[2,]>0.1)])
par(xpd =T)
text(x = as.numeric(df_effect_agpp_abs[1,][union(which(df_effect_agpp_abs[1,]> 0.1),which(df_effect_agpp_abs[2,]> 0.1))]),
     y=  as.numeric(df_effect_agpp_abs[2,][union(which(df_effect_agpp_abs[1,]> 0.1),which(df_effect_agpp_abs[2,]> 0.1))]),
     labels=names_important, cex = 1.5, pos=3, offset = 1.5)
par(xpd =F)

plot(x = df_effect_cflux_abs[1,], y = df_effect_cflux_abs[2,], type = "p", log = "xy",
     col = colors_me[as.numeric(as.factor(df_effect_cflux_abs[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(0.01, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(0.01, 0.2), las = 1, cex = 2)
axis(side = 2, at = c(0.01,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
axis(side = 1, at = c(0.01,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
segments(0.01,0.01,0.2,0.2)

title(main = "c) ", adj = 0., cex.main =2)
title(xlab = " Uncertainty contributions in a changing climate", line = 2,cex.lab = 2.)
title(ylab = " Uncertainty contributions in a stable future climate", line = 2,cex.lab = 2.)
names_important = union(parameternames_ordered_pic_abi2[which(df_effect_cflux_abs[1,]> 0.05)],
                        parameternames_ordered_pic_abi2[which(df_effect_cflux_abs[2,]>0.05)])
par(xpd =T)
text(x = as.numeric(df_effect_cflux_abs[1,][union(which(df_effect_cflux_abs[1,]> 0.05),which(df_effect_cflux_abs[2,]> 0.05))]),
     y=  as.numeric(df_effect_cflux_abs[2,][union(which(df_effect_cflux_abs[1,]> 0.05),which(df_effect_cflux_abs[2,]> 0.05))]),
     labels=names_important, cex =2, pos=3, offset = 1.5)
par(xpd =F)
par(xpd = T)
legend(x = 0.12, y = 0.1,legend = unique(df_effect_comp_agpp_abs[3,]), col = unique(colors_me[as.numeric(as.factor(df_effect_comp_agpp_abs[3,]))]), pch = 19,
       bty = "n", title = as.expression(bquote(bold("Processes"))),
       cex = 1.98)
par(xpd = F)

dev.off()

pdf("./Figures/Comparison_mono_abs.pdf",width=14,height=10)
layout(matrix(c(1,3,3,
                2,3,3), nrow = 2,ncol = 3, byrow = TRUE),
       widths = c(1,1,1,0.32),
       heights =  c(1,1,1,1))

plot(x = df_effect_comp_cpool_abs[1,], y = df_effect_comp_cpool_abs[2,], type = "p", log = "xy",
     col = colors_me[as.numeric(as.factor(df_effect_comp_cpool_abs[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(0.01, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(0.01, 0.2), las = 1, cex = 1.5)
axis(side = 2, at = c(0.01,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
axis(side = 1, at = c(0.01,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.01, las = 1)
segments(0.01,0.01,0.2,0.2)
title(main = "a) ",adj = 0., cex.main = 2)
title(xlab = " Uncertainty contributions in mono stands", line = 2,cex.lab = 1.3)
title(ylab = " Uncertainty contributions in mixed stands", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(df_effect_comp_cpool_abs[1,]> 0.08)],
                        parameternames_ordered_pic_abi2[which(df_effect_comp_cpool_abs[2,]>0.08)])
par(xpd =T)
text(x = as.numeric(df_effect_comp_cpool_abs[1,][union(which(df_effect_comp_cpool_abs[1,]> 0.08),which(df_effect_comp_cpool_abs[2,]> 0.08))]),
     y=  as.numeric(df_effect_comp_cpool_abs[2,][union(which(df_effect_comp_cpool_abs[1,]> 0.08),which(df_effect_comp_cpool_abs[2,]> 0.08))]),
     labels=names_important, cex = 1.5, pos=3, offset = 1.5)
par(xpd =F)


plot(x = df_effect_comp_agpp_abs[1,], y = df_effect_comp_agpp_abs[2,], type = "p", log = "xy",
     col = colors_me[as.numeric(as.factor(df_effect_comp_agpp_abs[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(0.005, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(0.005, 0.2), las = 1, cex = 1.5)
axis(side = 2, at = c(0.005,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.005, las = 1)
axis(side = 1, at = c(0.005,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.005, las = 1)
segments(0.005,0.005,0.2,0.2)
title(main = "b) ",adj = 0., cex.main = 2)
title(xlab = " Uncertainty contributions in mono stands", line = 2,cex.lab = 1.5)
title(ylab = " Uncertainty contributions in mixed stands", line = 2,cex.lab = 1.5)
names_important = union(parameternames_ordered_pic_abi2[which(df_effect_comp_agpp_abs[1,]> 0.08)],
                        parameternames_ordered_pic_abi2[which(df_effect_comp_agpp_abs[2,]>0.08)])
par(xpd =T)
text(x = as.numeric(df_effect_comp_agpp_abs[1,][union(which(df_effect_comp_agpp_abs[1,]> 0.08),which(df_effect_comp_agpp_abs[2,]> 0.08))]),
     y=  as.numeric(df_effect_comp_agpp_abs[2,][union(which(df_effect_comp_agpp_abs[1,]> 0.08),which(df_effect_comp_agpp_abs[2,]> 0.08))]),
     labels=names_important, cex = 1.5, pos=3, offset = 1.5)
par(xpd =F)


plot(x = df_effect_comp_cflux_abs[1,], y = df_effect_comp_cflux_abs[2,], type = "p", log = "xy",
     col = colors_me[as.numeric(as.factor(df_effect_comp_cflux_abs[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(0.005, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(0.005, 0.2), las = 1, cex = 2)
axis(side = 2, at = c(0.005,seq(0.05,.2 ,by = 0.05)), labels = T, pos = 0.005, las = 1)
axis(side = 1,at = c(0.005,seq(0.05,0.2 ,by = 0.05)), labels = T, pos = 0.005, las = 1)
segments(0.005,0.005,0.2,0.2)
title(main = "c) ",adj = 0., cex.main = 2)
title(xlab = "Uncertainty contributions in mono stands", line = 2,cex.lab =2.)
title(ylab = "Uncertainty contributions in mixed stands", line = 2,cex.lab = 2.)
names_important = union(parameternames_ordered_pic_abi2[which(df_effect_comp_cflux_abs[1,]> 0.05)],
                        parameternames_ordered_pic_abi2[which(df_effect_comp_cflux_abs[2,]>0.05)])
par(xpd =T)
text(x = as.numeric(df_effect_comp_cflux_abs[1,][union(which(df_effect_comp_cflux_abs[1,]> 0.05),which(df_effect_comp_cflux_abs[2,]> 0.05))]),
                    y=  as.numeric(df_effect_comp_cflux_abs[2,][union(which(df_effect_comp_cflux_abs[1,]> 0.05),which(df_effect_comp_cflux_abs[2,]> 0.05))]),
     labels=names_important, cex= 2, pos=3, offset = 1.5)
legend(x = 0.1, y = 0.08,legend = unique(df_effect_comp_agpp_abs[3,]), col = unique(colors_me[as.numeric(as.factor(df_effect_comp_agpp_abs[3,]))]), pch = 19,
       bty = "n", title = as.expression(bquote(bold("Processes"))),
       cex = 1.98)
par(xpd = F)
dev.off()


pdf("./Figures/Comparison_future_now.pdf",width=20,height=20)
layout(matrix(c(1,2,4,4,5,
                1,3,4,4,5), nrow = 2,ncol = 5, byrow = TRUE),
       widths = c(0.32,1,1,1,0.32),
       heights =  c(1,1,1,1))

colors_me = c("gold","burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
plot(x = df_effect_cpool_rel[1,], y = df_effect_cpool_rel[2,],
     type = "p", col = colors_me[as.numeric(as.factor(df_effect_cpool_rel[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(-0.2, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(-0.2, 0.2), las = 1)
axis(side = 2, at = seq(-0.2,0.2 ,by = 0.05), labels = T, pos = 0., las = 1)
axis(side = 1, at = seq(-0.2,0.2 ,by = 0.05), labels = T, pos = 0.,las = 1)
segments(-0.2,-.2,0.2,0.2)
title(main = "a)",  adj = 0., cex = 2)
title(xlab = " Uncertainty contributions in a changing climate", line = 2, cex.lab = 1.3)
title(ylab = "Uncertainty contributions in a stable future climate", line = 2, cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_cpool_rel[1,]))> 0.1)],
                        parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_cpool_rel[2,]))>0.1)])
par(xpd =T)
text(x = as.numeric(df_effect_cpool_rel[1,][union(which(abs(as.numeric(df_effect_cpool_rel[1,]))> 0.1),which(abs(as.numeric(df_effect_cpool_rel[2,]))> 0.1))]),
     y=  as.numeric(df_effect_cpool_rel[2,][union(which(abs(as.numeric(df_effect_cpool_rel[1,]))> 0.1),which(abs(as.numeric(df_effect_cpool_rel[2,]))> 0.1))]),
     labels=names_important, cex = 1.5, pos=3, offset = 1.5)
par(xpd =F)

par(xpd =T)
plot(x = df_effect_agpp_rel[1,], y = df_effect_agpp_rel[2,], type = "p",
     col = colors_me[as.numeric(as.factor(df_effect_agpp_rel[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(-0.2, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(-0.2, 0.2), las = 1)
axis(side = 2, at = seq(-0.2,0.2, by = 0.05), labels = T,pos = 0., las = 1)
axis(side = 1, at = seq(-0.2,0.2 ,by = 0.05), labels = T,pos = 0.,  las = 1)
segments(-0.2,-0.2,0.2,0.2)

title(main = "b) ",adj = 0., cex = 2)
title(xlab =  "Uncertainty contributions in a changing climate", line = 2,cex.lab = 1.3)
title(ylab =  "Uncertainty contributions in a stable future climate", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_agpp_rel[1,]))> 0.1)],
                        parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_agpp_rel[2,]))>0.1)])

text(x = as.numeric(df_effect_agpp_rel[1,][union(which(abs(as.numeric(df_effect_agpp_rel[1,]))> 0.1),which(abs(as.numeric(df_effect_agpp_rel[2,]))> 0.1))]),
     y = as.numeric(df_effect_agpp_rel[2,][union(which(abs(as.numeric(df_effect_agpp_rel[1,]))> 0.1),which(abs(as.numeric(df_effect_agpp_rel[2,]))> 0.1))]),
     labels=names_important, cex = 1.5, pos=3, offset = 1.5)
par(xpd =F)

plot(x = df_effect_cflux_rel[1,], y = df_effect_cflux_rel[2,], type = "p",
     col = colors_me[as.numeric(as.factor(df_effect_cflux_rel[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(-0.2, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(-0.2, 0.2), las = 1)
axis(side = 2, at = seq(-0.2,0.2 ,by = 0.05), labels = T, pos = 0., las = 1)
axis(side = 1, at = seq(-0.2,0.2 ,by = 0.05), labels = T, pos = 0., las = 1)
segments(-0.2,-0.2,0.2,0.2)

title(main = "c) ", adj = 0., cex = 2)
title(xlab = " Uncertainty contributions in a changing climate", line = 2,cex.lab = 1.3)
title(ylab = " Uncertainty contributions in a stable future climate", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_cflux_rel[1,]))> 0.05)],
                        parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_cflux_rel[2,]))>0.05)])
par(xpd =T)
text(x = as.numeric(df_effect_cflux_rel[1,][union(which(abs(as.numeric(df_effect_cflux_rel[1,]))> 0.05),which(abs(as.numeric(df_effect_cflux_rel[2,]))> 0.05))]),
     y=  as.numeric(df_effect_cflux_rel[2,][union(which(abs(as.numeric(df_effect_cflux_rel[1,]))> 0.05),which(abs(as.numeric(df_effect_cflux_rel[2,]))> 0.05))]),
     labels=names_important, cex = 2., pos=3, offset = 1.5)
par(xpd =F)

dev.off()

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text( paste("II) Comparison mono and mixed stands"),
      cex = 1.6, col = "black", font =2, srt = 90, x = 0.5, y = 0.5)

plot(x = df_effect_comp_cpool_rel[1,], y = df_effect_comp_cpool_rel[2,], type = "p",
     col = colors_me[as.numeric(as.factor(df_effect_comp_cpool_rel[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(-0.2, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(-0.2, 0.2), las = 1)
axis(side = 2, at = seq(-0.2,0.2 ,by = 0.05), labels = T,pos = 0., las = 1)
axis(side = 1, at = seq(-0.2,0.2 ,by = 0.05), labels = T,pos = 0., las = 1)
segments(-0.2,-0.2,0.2,0.2)
title(main = "d) ",adj = 0.)
title(xlab = " Uncertainty contributions in mono stands", line = 2,cex.lab = 1.3)
title(ylab = " Uncertainty contributions in mixed stands", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_comp_cpool_rel[1,]))> 0.1)],
                        parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_comp_cpool_rel[2,]))>0.1)])
par(xpd =T)
text(x = as.numeric(df_effect_comp_cpool_rel[1,][union(which(abs(as.numeric(df_effect_comp_cpool_rel[1,]))> 0.1),which(abs(as.numeric(df_effect_comp_cpool_rel[2,]))> 0.1))]),
     y=  as.numeric(df_effect_comp_cpool_rel[2,][union(which(abs(as.numeric(df_effect_comp_cpool_rel[1,]))> 0.1),which(abs(as.numeric(df_effect_comp_cpool_rel[2,]))> 0.1))]),
     labels=names_important, cex = 1.3, pos=3, offset = 1.5)
par(xpd =F)


plot(x = df_effect_comp_agpp_rel[1,], y = df_effect_comp_agpp_rel[2,], type = "p",
     col = colors_me[as.numeric(as.factor(df_effect_comp_agpp_rel[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(-0.2, 0.2),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(-0.2, 0.2), las = 1)
axis(side = 2, at = seq(-0.2,0.2 ,by = 0.05), labels = T,pos = 0., las = 1)
axis(side = 1, at = seq(-0.2,0.2 ,by = 0.05), labels = T,pos = 0., las = 1)
segments(-0.2,-0.2,0.2,0.2)
title(main = "e) ",adj = 0.)
title(xlab = " Uncertainty contributions in mono stands", line = 2,cex.lab = 1.3)
title(ylab = " Uncertainty contributions in mixed stands", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_comp_agpp_rel[1,]))> 0.08)],
                        parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_comp_agpp_rel[2,]))>0.08)])
par(xpd =T)
text(x = as.numeric(df_effect_comp_agpp_rel[1,][union(which(abs(as.numeric(df_effect_comp_agpp_rel[1,]))> 0.08),which(abs(as.numeric(df_effect_comp_agpp_rel[2,]))> 0.08))]),
     y=  as.numeric(df_effect_comp_agpp_rel[2,][union(which(abs(as.numeric(df_effect_comp_agpp_rel[1,]))> 0.08),which(abs(as.numeric(df_effect_comp_agpp_rel[2,]))> 0.08))]),
     labels=names_important, cex = 1.3, pos=3, offset = 1.5)
par(xpd =F)


plot(x = df_effect_comp_cflux_rel[1,], y = df_effect_comp_cflux_rel[2,], type = "p",
     col = colors_me[as.numeric(as.factor(df_effect_comp_cflux_rel[3,]))], pch = 19,
     xlab = "", ylab = "", xlim = c(-0.1, 0.1),bty = "n", xaxt ="n", yaxt = "n",
     ylim = c(-0.1, 0.1), las = 1)
axis(side = 2, at = seq(-.1,.1 ,by = 0.05), labels = T,pos = 0.,  las = 1)
axis(side = 1, at = seq(-.1,0.1 ,by = 0.05), labels = T,pos = 0.,  las = 1)
segments(-0.1,-0.1,0.1,0.1)
title(main = "f) ",adj = 0.)
title(xlab = "Uncertainty contributions in mono stands", line = 2,cex.lab = 1.3)
title(ylab = "Uncertainty contributions in mixed stands", line = 2,cex.lab = 1.3)
names_important = union(parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_comp_cflux_rel[1,]))> 0.08)],
                        parameternames_ordered_pic_abi2[which(abs(as.numeric(df_effect_comp_cflux_rel[2,]))>0.08)])
par(xpd =T)
text(x = as.numeric(df_effect_comp_cflux_rel[1,][union(which(abs(as.numeric(df_effect_comp_cflux_rel[1,]))> 0.08),which(abs(as.numeric(df_effect_comp_cflux_rel[2,]))> 0.08))]),
     y=  as.numeric(df_effect_comp_cflux_rel[2,][union(which(abs(as.numeric(df_effect_comp_cflux_rel[1,]))> 0.08),which(abs(as.numeric(df_effect_comp_cflux_rel[2,]))> 0.08))]),
     labels=names_important, cex= 1.3, pos=3, offset = 1.5)
par(xpd =F)


plot.new()
par(xpd = T)
legend(x = -0.65, y = 0.7,legend = unique(df_effect_comp_agpp_rel[3,]), col = unique(colors_me[as.numeric(as.factor(df_effect_comp_agpp_rel[3,]))]), pch = 19,
       bty = "n", title = as.expression(bquote(bold("Processes"))),
       cex = 1.5)
par(xpd = F)
dev.off()


