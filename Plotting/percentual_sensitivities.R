# function to rescale effects from zero to one
rescaling <- function(x) {(x-min(x))/(max(x)-min(x))}

# function to rescale effects from 0 to 1
rescaling_abs <- function(x) {abs(x)/max(abs(x))}



## getting the ordered parameter names ##

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
for(i in 1:length(parameters)){
  parameters[[i]][,2] = sub(",",".",parameters[[i]][,2])
  parameters[[i]][,3] = sub(",",".",parameters[[i]][,3])
  parameters[[i]][,4] = sub(",",".",parameters[[i]][,4])
}

parameters_only_pic_abi = parameters$Pic_abi$NameRLPJ
parameters_only_fag_syl = parameters$Fag_syl$NameRLPJ
parameters_only_pin_syl = parameters$Pin_syl$NameRLPJ

drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameters_pic_abi = c(as.character(parameters$Pic_abi$NameRLPJ), drivernames)
parameters_fag_syl = c(as.character(parameters$Fag_syl$NameRLPJ), drivernames)
parameters_pin_syl = c(as.character(parameters$Pin_syl$NameRLPJ), drivernames)



ordering_param_only_fag_syl = match(gsub("Pic_abi","",parameters$Pic_abi$NameRLPJ), gsub("Fag_syl","",parameters$Fag_syl$NameRLPJ))
ordering_param_only_pin_syl = match(gsub("Pic_abi","",parameters$Pic_abi$NameRLPJ), gsub("Pin_syl","",parameters$Pin_syl$NameRLPJ))

# loading the results to get the names in the prespecified order

results =  readRDS(paste0("./LPJrunTest/Results/Pic_abi_0.25_42.25.rds"))
parameternames_ordered_pic_abi = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pic_abi)])

results =  readRDS(paste0("./LPJrunTest/Results/Fag_syl_0.25_42.25.rds"))
parameternames_ordered_fag_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_fag_syl)])

results =  readRDS(paste0("./LPJrunTest/Results/Pin_syl_0.25_42.25.rds"))
parameternames_ordered_pin_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pin_syl)])

mixed_results = readRDS("./LPJrunTest/Results/mapping_mixed.rds")

driver_components = which(mixed_results$parameternames %in% c("temp_change","ph_change","prec_change",
                                                              "insol_change","ndep_change","co2_change"))

## get the order into the parameters by a split

parameternames_ordered_pic_abi = substring(parameternames_ordered_pic_abi,regexpr("_", parameternames_ordered_pic_abi) + 1)
parameternames_ordered_fag_syl = substring(parameternames_ordered_fag_syl,regexpr("_", parameternames_ordered_fag_syl) + 1)
parameternames_ordered_pin_syl = substring(parameternames_ordered_pin_syl,regexpr("_", parameternames_ordered_pin_syl) + 1)
parameternames_ordered_mixed = mixed_results$parameternames



parameters_only_pic_abi = substring(parameters_only_pic_abi,regexpr("_", parameters_only_pic_abi) + 1)
parameters_only_fag_syl = substring(parameters_only_fag_syl,regexpr("_", parameters_only_fag_syl) + 1)
parameters_only_pin_syl = substring(parameters_only_pin_syl,regexpr("_", parameters_only_pin_syl) + 1)
parameters_only_mixed = mixed_results$parameternames[-driver_components]

## eliminate the noisy things in the names, that do not allow matching the data

noisy_things = c("intolerant","tolerant","tree","abi","_","syl","leaved")
for(i in 1:length(noisy_things)){
  parameternames_ordered_pic_abi = gsub(noisy_things[i],"",parameternames_ordered_pic_abi)
  parameternames_ordered_fag_syl = gsub(noisy_things[i],"",parameternames_ordered_fag_syl)
  parameternames_ordered_pin_syl = gsub(noisy_things[i],"",parameternames_ordered_pin_syl)
  parameters_only_fag_syl = gsub(noisy_things[i],"",parameters_only_fag_syl)
  parameters_only_pin_syl = gsub(noisy_things[i],"",parameters_only_pin_syl)
  parameters_only_pic_abi = gsub(noisy_things[i],"",parameters_only_pic_abi)
  parameters_only_mixed = gsub(noisy_things[i],"",parameters_only_mixed)
}

ordering_fag_syl = match(parameternames_ordered_pic_abi,parameternames_ordered_fag_syl)
ordering_pin_syl = match(parameternames_ordered_pic_abi,parameternames_ordered_pin_syl)


ordering_param_only_fag_syl = match(parameters_only_pic_abi, parameters_only_fag_syl)
ordering_param_only_pin_syl = match(parameters_only_pic_abi, parameters_only_pin_syl)
ordering_param_only_mixed = match(parameters_only_pic_abi, parameters_only_mixed)




### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects_lin.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects_lin.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects_lin.rds")

effects_mixed = readRDS("LPJrunTest/Results/mixed_effects_lin.rds")

drivernames = c("ph_change","temp_change","co2_change","insol_change","prec_change","ndep_change")

for(i in 1:length(drivernames)){
  parameternames_ordered_pic_abi = gsub(drivernames[i],"NA",parameternames_ordered_pic_abi)
  parameternames_ordered_fag_syl = gsub(drivernames[i],"NA",parameternames_ordered_fag_syl)
  parameternames_ordered_pin_syl = gsub(drivernames[i],"NA",parameternames_ordered_pin_syl)
}

#### Analysis for Pic abi ####

# function to get the sensitivities by dividing to the range and default value
get_sensi = function(effect_sizes,
                     maximum,
                     minimum,
                     default){

  sensitivities = vector("numeric", length(effect_sizes))
  for(i in 1:length(sensitivities)){
    sensitivities[i] = abs(effect_sizes[i])/(abs(maximum[i]-minimum[i])*abs(default[i]))
  }
  return(sensitivities)
}




parameters$Pic_abi$Maximum.Value = gsub(",",".",parameters$Pic_abi$Maximum.Value)
parameters$Pic_abi$Minimum.Value = gsub(",",".",parameters$Pic_abi$Minimum.Value)
parameters$Pic_abi$Default = gsub(",",".",parameters$Pic_abi$Default)


mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
mean_effects_Pic_abi_cpool = mean_effects_Pic_abi_cpool[growth_sites_Pic_abi,]
main_effects_Pic_abi_cpool = mean_effects_Pic_abi_cpool/mean(abs(mean_effects_Pic_abi_cpool[,39]))
main_effects_Pic_abi_cpool = apply(main_effects_Pic_abi_cpool[,1:38],2,mean,na.rm =T)
main_effects_Pic_abi_cpool = main_effects_Pic_abi_cpool[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Pic_abi_cpool = get_sensi(main_effects_Pic_abi_cpool,
                                     as.numeric(parameters$Pic_abi$Maximum.Value),
                                     as.numeric(parameters$Pic_abi$Minimum.Value),
                                     as.numeric(parameters$Pic_abi$Default))


mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]]
mean_effects_Pic_abi_cflux = mean_effects_Pic_abi_cflux[growth_sites_Pic_abi,]
main_effects_Pic_abi_cflux = mean_effects_Pic_abi_cflux/mean(abs(mean_effects_Pic_abi_cflux[,39]))
main_effects_Pic_abi_cflux = apply(main_effects_Pic_abi_cflux[,1:38],2,mean,na.rm =T)
main_effects_Pic_abi_cflux = main_effects_Pic_abi_cflux[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Pic_abi_cflux = get_sensi(main_effects_Pic_abi_cflux,
                                     as.numeric(parameters$Pic_abi$Maximum.Value),
                                     as.numeric(parameters$Pic_abi$Minimum.Value),
                                     as.numeric(parameters$Pic_abi$Default))


mean_effects_Pic_abi_agpp = effects_Pic_abi[["agpp"]][["complete"]]
mean_effects_Pic_abi_agpp = mean_effects_Pic_abi_agpp[growth_sites_Pic_abi,]
main_effects_Pic_abi_agpp = mean_effects_Pic_abi_agpp/mean(abs(mean_effects_Pic_abi_agpp[,39]))
main_effects_Pic_abi_agpp = apply(main_effects_Pic_abi_agpp[,1:38],2,mean,na.rm =T)
main_effects_Pic_abi_agpp = main_effects_Pic_abi_agpp[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Pic_abi_agpp = get_sensi(main_effects_Pic_abi_agpp,
                                     as.numeric(parameters$Pic_abi$Maximum.Value),
                                     as.numeric(parameters$Pic_abi$Minimum.Value),
                                     as.numeric(parameters$Pic_abi$Default))



#### Analysis for Fag syl ####


parameters$Fag_syl$Maximum.Value = gsub(",",".",parameters$Fag_syl$Maximum.Value)
parameters$Fag_syl$Minimum.Value = gsub(",",".",parameters$Fag_syl$Minimum.Value)
parameters$Fag_syl$Default = gsub(",",".",parameters$Fag_syl$Default)


mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
growth_sites_Fag_syl = which(mean_effects_Fag_syl_cpool[,39]/200 > 2)
mean_effects_Fag_syl_cpool = mean_effects_Fag_syl_cpool[growth_sites_Fag_syl,]
main_effects_Fag_syl_cpool = mean_effects_Fag_syl_cpool/mean(abs(mean_effects_Fag_syl_cpool[,39]))
main_effects_Fag_syl_cpool = apply(main_effects_Fag_syl_cpool[,1:38],2,mean,na.rm =T)
main_effects_Fag_syl_cpool = main_effects_Fag_syl_cpool[ordering_fag_syl]
main_effects_Fag_syl_cpool = main_effects_Fag_syl_cpool[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Fag_syl_cpool = get_sensi(main_effects_Fag_syl_cpool,
                                     as.numeric(parameters$Fag_syl$Maximum.Value[ordering_param_only_fag_syl]),
                                     as.numeric(parameters$Fag_syl$Minimum.Value[ordering_param_only_fag_syl]),
                                     as.numeric(parameters$Fag_syl$Default[ordering_param_only_fag_syl]))


mean_effects_Fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]]
mean_effects_Fag_syl_cflux = mean_effects_Fag_syl_cflux[growth_sites_Fag_syl,]
main_effects_Fag_syl_cflux = mean_effects_Fag_syl_cflux/mean(abs(mean_effects_Fag_syl_cflux[,39]))
main_effects_Fag_syl_cflux = apply(main_effects_Fag_syl_cflux[,1:38],2,mean,na.rm =T)
main_effects_Fag_syl_cflux = main_effects_Fag_syl_cflux[ordering_fag_syl]
main_effects_Fag_syl_cflux = main_effects_Fag_syl_cflux[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Fag_syl_cflux = get_sensi(main_effects_Fag_syl_cflux,
                                     as.numeric(parameters$Fag_syl$Maximum.Value[ordering_param_only_fag_syl]),
                                     as.numeric(parameters$Fag_syl$Minimum.Value[ordering_param_only_fag_syl]),
                                     as.numeric(parameters$Fag_syl$Default[ordering_param_only_fag_syl]))


mean_effects_Fag_syl_agpp = effects_Fag_syl[["agpp"]][["complete"]]
mean_effects_Fag_syl_agpp = mean_effects_Fag_syl_agpp[growth_sites_Fag_syl,]
main_effects_Fag_syl_agpp = mean_effects_Fag_syl_agpp/mean(abs(mean_effects_Fag_syl_agpp[,39]))
main_effects_Fag_syl_agpp = apply(main_effects_Fag_syl_agpp[,1:38],2,mean)
main_effects_Fag_syl_agpp = main_effects_Fag_syl_agpp[ordering_fag_syl]
main_effects_Fag_syl_agpp = main_effects_Fag_syl_agpp[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Fag_syl_agpp = get_sensi(main_effects_Fag_syl_agpp,
                                     as.numeric(parameters$Fag_syl$Maximum.Value[ordering_param_only_fag_syl]),
                                     as.numeric(parameters$Fag_syl$Minimum.Value[ordering_param_only_fag_syl]),
                                     as.numeric(parameters$Fag_syl$Default[ordering_param_only_fag_syl]))



#### Analysis for Pin syl ####

parameters$Pin_syl$Maximum.Value = gsub(",",".",parameters$Pin_syl$Maximum.Value)
parameters$Pin_syl$Minimum.Value = gsub(",",".",parameters$Pin_syl$Minimum.Value)
parameters$Pin_syl$Default = gsub(",",".",parameters$Pin_syl$Default)


mean_effects_Pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
growth_sites_Pin_syl = which(mean_effects_Pin_syl_cpool[,39]/200 > 2)
mean_effects_Pin_syl_cpool = mean_effects_Pin_syl_cpool[growth_sites_Pin_syl,]
main_effects_Pin_syl_cpool = mean_effects_Pin_syl_cpool/mean(abs(mean_effects_Pin_syl_cpool[,39]))
main_effects_Pin_syl_cpool = apply(main_effects_Pin_syl_cpool[,1:38],2,mean, na.rm =T)
main_effects_Pin_syl_cpool = main_effects_Pin_syl_cpool[ordering_pin_syl]
main_effects_Pin_syl_cpool = main_effects_Pin_syl_cpool[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Pin_syl_cpool = get_sensi(main_effects_Pin_syl_cpool,
                                     as.numeric(parameters$Pin_syl$Maximum.Value[ordering_param_only_pin_syl]),
                                     as.numeric(parameters$Pin_syl$Minimum.Value[ordering_param_only_pin_syl]),
                                     as.numeric(parameters$Pin_syl$Default[ordering_param_only_pin_syl]))

mean_effects_Pin_syl_cflux = effects_Pin_syl[["cflux"]][["complete"]]
mean_effects_Pin_syl_cflux = mean_effects_Pin_syl_cflux[growth_sites_Pin_syl,]
main_effects_Pin_syl_cflux = mean_effects_Pin_syl_cflux/mean(abs(mean_effects_Pin_syl_cflux[,39]))
main_effects_Pin_syl_cflux = apply(main_effects_Pin_syl_cflux[,1:38],2,mean,na.rm =T)
main_effects_Pin_syl_cflux = main_effects_Pin_syl_cflux[ordering_pin_syl]
main_effects_Pin_syl_cflux = main_effects_Pin_syl_cflux[which(parameternames_ordered_pic_abi != "NA")]
main_sensi_Pin_syl_cflux = get_sensi(main_effects_Pin_syl_cflux,
                                     as.numeric(parameters$Pin_syl$Maximum.Value[ordering_param_only_pin_syl]),
                                     as.numeric(parameters$Pin_syl$Minimum.Value[ordering_param_only_pin_syl]),
                                     as.numeric(parameters$Pin_syl$Default[ordering_param_only_pin_syl]))


mean_effects_Pin_syl_agpp = effects_Pin_syl[["agpp"]][["complete"]]
mean_effects_Pin_syl_agpp = mean_effects_Pin_syl_agpp[growth_sites_Pin_syl,]
main_effects_Pin_syl_agpp = mean_effects_Pin_syl_agpp/mean(abs(mean_effects_Pin_syl_agpp[,39]))
main_effects_Pin_syl_agpp = apply(main_effects_Pin_syl_agpp[,1:38],2,mean,na.rm =T)
main_effects_Pin_syl_agpp = main_effects_Pin_syl_agpp[ordering_pin_syl]
main_effects_Pin_syl_agpp = main_effects_Pin_syl_agpp[which(parameternames_ordered_pic_abi != "NA")]

main_sensi_Pin_syl_agpp = get_sensi(main_effects_Pin_syl_agpp,
                                    as.numeric(parameters$Pin_syl$Maximum.Value[ordering_param_only_pin_syl]),
                                    as.numeric(parameters$Pin_syl$Minimum.Value[ordering_param_only_pin_syl]),
                                    as.numeric(parameters$Pin_syl$Default[ordering_param_only_pin_syl]))
#### Analysis for mixed ####

results =  readRDS(paste0("./../Results_sensi/Mixed_0.25_42.25.rds"))
parametermixed = readRDS("ParameterMetaData/Parameter_mixed.rds")
parameternames_ordered_mixed = names(unlist(results[[1]]@runInfo$parameterList)[which(names(unlist(results[[1]]@runInfo$parameterList)) %in% parametermixed$NameRLPJ)])

parameternames_ordered_mixed_try =  substring(parameternames_ordered_mixed,regexpr("_", parameternames_ordered_mixed) + 1)
new_noisy = c("syl_","abi_","leaved_","tolerant_","tree_","intolerant_")
for(i in 1:length(new_noisy)){
  parameternames_ordered_mixed_try  = gsub(new_noisy[i],"",parameternames_ordered_mixed_try )
}
mixed_parameter_only = c()
for(i in 1:length(mixed_results$parameternames)){
  if(any(parameternames_ordered_mixed_try == mixed_results$parameternames[i])){
    positions = mixed_results$position_mapping[[i]]
    mixed_parameter_only = c(mixed_parameter_only, positions)
  }
  else{
    next
  }
}
mixed_parameter_only = sort(mixed_parameter_only)



# function for the mixed simulations to caucluate sensitivties by dividing to the prespecified range


calculate_sensi <- function(uncertainty, nameofparameters, nameslist,
                            minimum, maximum,default, reference){
  sensitivities = matrix(ncol = ncol(uncertainty),
                         nrow = nrow(uncertainty))
  default = as.numeric(gsub(",",".",default))
  minimum = as.numeric(gsub(",",".",minimum))
  maximum = as.numeric(gsub(",",".",maximum))
  for(parameter in 1:ncol(uncertainty)){
    name = which(nameslist[parameter] == nameofparameters)
    sensitivities[,parameter] = abs(uncertainty[,parameter])/
      (abs(maximum[name]-minimum[name])*abs(default[name])*mean(abs(reference)))
  }
  return(sensitivities)
}


remap_sensi <- function(sensitivities, order_pic_abi,
                        weighting_scheme, position_scheme,
                        driver_components){

  # get the positions of drivers in the sensi

  driver_positions = unlist(position_scheme[driver_components])

  ## find the driver components in the position and weighting scheme
  ## and eliminate them from the scheme

  position_scheme = position_scheme[-driver_components]
  weighting_scheme = weighting_scheme[-driver_components]


  ## get the right positions in the new matrix

  for(i in 1:length(position_scheme)){
    for(j in 1:length(position_scheme[[i]])){
      if(position_scheme[[i]][j] < min(driver_positions)){
        position_scheme[[i]][j] = position_scheme[[i]][j]
      }
      else if(position_scheme[[i]][j] > max(driver_positions)){
        position_scheme[[i]][j] = position_scheme[[i]][j] - length(driver_positions)
      }
      else{
        position_scheme[[i]][j] = position_scheme[[i]][j]
      }
    }
  }


  ## remap the parameters to parameters only

  mapped_parameters = matrix(ncol = length(weighting_scheme),
                             nrow = nrow(sensitivities))
  for(site in 1:nrow(sensitivities)){
    for(parameter in 1:length(position_scheme)){
      print(position_scheme[[parameter]])
      mapped_parameters[site,parameter] = sum(abs(sensitivities[site,position_scheme[[parameter]]])*
                                                weighting_scheme[[parameter]])
    }
  }
  mean_mapped_parameters = apply(FUN = mean, X = abs(mapped_parameters), MARGIN = 2,na.rm= T)

  mean_mapped_parameters = mean_mapped_parameters[order_pic_abi]

  return(mean_mapped_parameters)
}

## Applying the analysis for the mixed case for all setting which have been run

mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
mean_effects_mixed_cpool = mean_effects_mixed_cpool[growth_sites_mixed,]
mean_effects_mixed_cpool = mean_effects_mixed_cpool[,c(mixed_parameter_only,74)]

mean_sensi_mixed_cpool = calculate_sensi(mean_effects_mixed_cpool[,1:(ncol(mean_effects_mixed_cpool)-1)],
                                         nameofparameters = parametermixed$NameRLPJ, nameslist = parameternames_ordered_mixed ,
                minimum = parametermixed$Minimum.Value,maximum = parametermixed$Maximum.Value , default = parametermixed$Default,
                reference = mean_effects_mixed_cpool[,ncol(mean_effects_mixed_cpool)])


main_sensi_mixed_cpool = remap_sensi(mean_sensi_mixed_cpool, order_pic_abi = ordering_param_only_mixed,
            weighting_scheme = mixed_results$weight_mapping,
            position_scheme = mixed_results$position_mapping,
            driver_components = driver_components)


mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]]
mean_effects_mixed_cflux = mean_effects_mixed_cflux[growth_sites_mixed,]
mean_effects_mixed_cflux = mean_effects_mixed_cflux[-which(abs(mean_effects_mixed_cflux[,74])<2),]
mean_effects_mixed_cflux = mean_effects_mixed_cflux[,c(mixed_parameter_only,74)]

mean_sensi_mixed_cflux = calculate_sensi(mean_effects_mixed_cflux[,1:(ncol(mean_effects_mixed_cflux)-1)], nameofparameters = parametermixed$NameRLPJ, nameslist = parameternames_ordered_mixed ,
                                         minimum = parametermixed$Minimum.Value,maximum = parametermixed$Maximum.Value , default = parametermixed$Default,
                                         reference = mean_effects_mixed_cflux[,ncol(mean_effects_mixed_cflux)])


main_sensi_mixed_cflux =  remap_sensi(mean_sensi_mixed_cflux, order_pic_abi = ordering_param_only_mixed,
                                      weighting_scheme = mixed_results$weight_mapping,
                                      position_scheme = mixed_results$position_mapping,
                                      driver_components = driver_components)


mean_effects_mixed_agpp = effects_mixed[["agpp"]][["complete"]]
mean_effects_mixed_agpp = mean_effects_mixed_agpp[growth_sites_mixed,]
mean_effects_mixed_agpp = mean_effects_mixed_agpp[,c(mixed_parameter_only,74)]

mean_sensi_mixed_agpp = calculate_sensi(mean_effects_mixed_agpp[,1:(ncol(mean_effects_mixed_agpp)-1)], nameofparameters = parametermixed$NameRLPJ, nameslist = parameternames_ordered_mixed ,
                                         minimum = parametermixed$Minimum.Value,maximum = parametermixed$Maximum.Value , default = parametermixed$Default,
                                         reference = mean_effects_mixed_agpp[,ncol(mean_effects_mixed_agpp)])


main_sensi_mixed_agpp =  remap_sensi(mean_sensi_mixed_agpp, order_pic_abi = ordering_param_only_mixed,
                                      weighting_scheme = mixed_results$weight_mapping,
                                      position_scheme = mixed_results$position_mapping,
                                      driver_components = driver_components)



parameters = readRDS("ParameterMetaData/Parameter_list.rds")
variablenames = c(as.character(parameters$Pic_abi[,"NameRLPJ"]))
variablenames2 = substring(variablenames,regexpr("_", variablenames) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  variablenames2 = gsub(noisy_things[i],"",variablenames2)
}

grouping = read.csv("ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")$Group

order_grouping = match(parameters_only_pic_abi,variablenames2)
grouping2 = grouping[order_grouping]




mean_sensi_abs = data.frame("Group" = grouping2,
                              "cpool_fag_syl" = main_sensi_Fag_syl_cpool[1:32],
                              "cflux_fag_syl" = main_sensi_Fag_syl_cflux[1:32],
                              "agpp_fag_syl" = main_sensi_Fag_syl_agpp[1:32],
                              "cpool_pin_syl" = main_sensi_Pin_syl_cpool[1:32],
                              "cflux_pin_syl" = main_sensi_Pin_syl_cflux[1:32],
                              "agpp_pin_syl" = main_sensi_Pin_syl_agpp[1:32],
                              "cpool_pic_abi" = main_sensi_Pic_abi_cpool[1:32],
                              "cflux_pic_abi" = main_sensi_Pic_abi_cflux[1:32],
                              "agpp_pic_abi" = main_sensi_Pic_abi_agpp[1:32],
                              "cpool_mixed" = main_sensi_mixed_cpool[1:32],
                              "cflux_mixed" = main_sensi_mixed_cflux[1:32],
                              "agpp_mixed" = main_sensi_mixed_agpp[1:32],
                              "Names" = parameters_only_pic_abi)


mean_sensi_abs = mean_sensi_abs[order(mean_sensi_abs$Group),]
mean_sensi_abs$Names
summary(mean_sensi_abs)

library(plotrix)
library(scales)

t_col <- function(color, percent = 20, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

  ## Get RGB values for named color
  rgb.val <- col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

  ## Save the color
  invisible(t.col)
}

log10Tck <- function(side, type){
  lim <- switch(side,
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceil(lim[2])
  return(switch(type,
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

colors_parameter = c("burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")
spaces = table(mean_sensi_abs$Group)
colors_full = vector()
for(i in 1:6){
  colors_full = c(colors_full,rep(colors_parameter[i],spaces[i]))
}

pdf("./Figures/Mean_sensi.pdf", width = 18.0, height = 12)

layout(matrix(c(1,2,3,4),byrow =T,nrow =4), heights = c(2,2,2,1))
old_mar = c(5.1,4.1,4.1,2.1)
par(mar = c(0.6,4.1,9.1,12.1))
plot(y = rep(0,33),x=1:33,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n',
     ylim = c(0,18))
mtext("Percentage sensitivity [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                             mean_sensi_abs$cpool_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32), log ="y")
points(y = mean_sensi_abs$cpool_fag_syl , x= 1:32, pch = 15,
       ylim = c(0,18), col = "darkgreen",
       cex = 1.2)
points(y = mean_sensi_abs$cpool_pic_abi,
       x = 1:nrow(mean_sensi_abs), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_abs$cpool_pin_syl,
       x = 1:nrow(mean_sensi_abs), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_abs$cpool_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_sensi_abs))
points(y = apply(rbind(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl),2,mean),
       x = 1:nrow(mean_sensi_abs), col = "black", pch = "-", cex=3)

axis(1, at = 1:32, labels = rep("",33),srt = 45, las =2 , pos = 0)
axis(2, at = c(0,seq(2,18,2)), las = 2, line = -0.5)
ablineclip(v = cumsum(table(mean_sensi_abs$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
spacings = (cumsum(table(mean_sensi_abs$Group)) + c(0,cumsum(table(mean_sensi_abs$Group))[-length(cumsum(table(mean_sensi_abs$Group)))]))/2 +0.5
grouping_variables = names(table(mean_sensi_abs$Group))
par(xpd =  T)
text(c(names(table(mean_sensi_abs$Group))[which(names(table(mean_sensi_abs$Group)) != "Drivers")]),
     x = c(spacings) +c(0.5,0,1.5,2.5,0,0.5),  y = 22 + c(1,0,0,1,0,-1), font = 2, cex =1.7, srt = 30)
par(xpd = F)
ablineclip(h = 18, x1 =0, x2 = 33.5)
title(main = "a)", line = 2, adj = 0.01, cex.main = 2)


par(mar = c(0.6,4.1,5.6,12.1))
plot(y = rep(0,33),x=1:33,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n',
     ylim = c(0,220))
mtext("Percentage sensitivity [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                             mean_sensi_abs$cflux_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_sensi_abs$cflux_fag_syl , x= 1:32, pch = 15,
       ylim = c(0,220), col = "darkgreen",
       cex = 1.2)
points(y = mean_sensi_abs$cflux_pic_abi,
       x = 1:nrow(mean_sensi_abs), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_abs$cflux_pin_syl,
       x = 1:nrow(mean_sensi_abs), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_abs$cflux_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_sensi_abs))
points(y = apply(rbind(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl),2,mean),
       x = 1:nrow(mean_sensi_abs), col = "black", pch = "-", cex=3)

axis(1, at = 1:32, labels = rep("",33),srt = 45, las =2 , pos = 0)
axis(2, at = seq(0,220,length.out = 23), las = 2, line = -0.5)
ablineclip(v = cumsum(table(mean_sensi_abs$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
spacings = (cumsum(table(mean_sensi_abs$Group)) + c(0,cumsum(table(mean_sensi_abs$Group))[-length(cumsum(table(mean_sensi_abs$Group)))]))/2 +0.5
grouping_variables = names(table(mean_sensi_abs$Group))
par(xpd=T)
legend(x = 34, y = 180,
       legend = c("Fag. syl.","Pic. abi.",'Pin. syl.',"Mixed",'Mean Mono'),
       col = c('darkgreen', 'darkblue','brown', "purple",'black'), bty = 'n',
       title = as.expression(bquote(italic(bold("Species")))),
       pch = c(15,19,17,8,NA), cex = 1.8, lty = c(NA,NA,NA,NA,1),
       lwd = 2)
par(xpd = F)
ablineclip(h = 220, x1 =0, x2 = 33.5)
title(main = "b)", line = 2, adj = 0.01, cex.main = 2)

par(mar = c(1.1,4.1,6.1,12.1))
plot(y = rep(0.0001,33),x=1:33,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n',
     ylim = c(0.,40))
mtext("Percentage sensitivity [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                             mean_sensi_abs$agpp_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_sensi_abs$agpp_fag_syl , x= 1:32, pch = 15,
       ylim = c(0.0,40), col = "darkgreen",
       cex = 1.2)
points(y = mean_sensi_abs$agpp_pic_abi,
       x = 1:nrow(mean_sensi_abs), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_abs$agpp_pin_syl,
       x = 1:nrow(mean_sensi_abs), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_abs$agpp_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_sensi_abs))
points(y = apply(rbind(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl),2,mean),
       x = 1:nrow(mean_sensi_abs), col = "black", pch = "-", cex=3)
axis(1, at = 1:32, labels = rep("",32),
     srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0.,130,length.out = 14),2), las = 2, line = -0.5)
ablineclip(v = cumsum(table(mean_sensi_abs$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
spacings = (cumsum(table(mean_sensi_abs$Group)) + c(0,cumsum(table(mean_sensi_abs$Group))[-length(cumsum(table(mean_sensi_abs$Group)))]))/2 +0.5
grouping_variables = names(table(mean_sensi_abs$Group))
ablineclip(h = 40, x1 =0, x2 = 33.5)
ablineclip(v = c(39.5), y1=0)
title(main = "c)", line = 1.7, adj = 0.01, cex.main = 2)

par(mar = c(1.1,4.1,0.1,12.1))
plot(y = rep(0,33),x=1:33,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",33)),xaxt='n', bty="n", yaxt = 'n', cex.axis =2,
     ylim = c(-0.01,0.01), cex.main = 1.5)
par(xpd = T)
text(x = 1:33, y = 0.01, c(as.character(mean_sensi_abs$Names),""),srt = 270,
     adj = c(0.0,0.), cex =1.5)
par(xpd = F)

dev.off()


