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


drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameters_pic_abi = c(as.character(parameters$Pic_abi$NameRLPJ))
parameters_fag_syl = c(as.character(parameters$Fag_syl$NameRLPJ))
parameters_pin_syl = c(as.character(parameters$Pin_syl$NameRLPJ))

#
#
#
# ordering_param_only_fag_syl = match(gsub("Pic_abi","",parameters$Pic_abi$NameRLPJ), gsub("Fag_syl","",parameters$Fag_syl$NameRLPJ))
# ordering_param_only_pin_syl = match(gsub("Pic_abi","",parameters$Pic_abi$NameRLPJ), gsub("Pin_syl","",parameters$Pin_syl$NameRLPJ))

# loading the results to get the names in the prespecified order

results =  readRDS(paste0("./LPJrunTest/Results/Pic_abi_0.25_42.25.rds"))
parameternames_ordered_pic_abi = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% c(parameters_pic_abi,drivernames))])

results =  readRDS(paste0("./LPJrunTest/Results/Fag_syl_0.25_42.25.rds"))
parameternames_ordered_fag_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% c(parameters_fag_syl,drivernames))])

results =  readRDS(paste0("./LPJrunTest/Results/Pin_syl_0.25_42.25.rds"))
parameternames_ordered_pin_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% c(parameters_pin_syl,drivernames))])

mixed_results = readRDS("./LPJrunTest/Results/mapping_mixed.rds")

driver_components = which(mixed_results$parameternames %in% c("temp_change","ph_change","prec_change",
                                                              "insol_change","ndep_change","co2_change"))

## get the order into the parameters by a split

parameternames_ordered_pic_abi = substring(parameternames_ordered_pic_abi,regexpr("_", parameternames_ordered_pic_abi) + 1)
parameternames_ordered_fag_syl = substring(parameternames_ordered_fag_syl,regexpr("_", parameternames_ordered_fag_syl) + 1)
parameternames_ordered_pin_syl = substring(parameternames_ordered_pin_syl,regexpr("_", parameternames_ordered_pin_syl) + 1)
parameternames_ordered_mixed = mixed_results$parameternames

parameters_pic_abi = substring(parameters_pic_abi,regexpr("_", parameters_pic_abi) + 1)
parameters_fag_syl = substring(parameters_fag_syl,regexpr("_", parameters_fag_syl) + 1)
parameters_pin_syl = substring(parameters_pin_syl,regexpr("_", parameters_pin_syl) + 1)



## eliminate the noisy things in the names, that do not allow matching the data

noisy_things = c("intolerant","tolerant","tree","abi","_","syl","leaved","common")
for(i in 1:length(noisy_things)){
  parameternames_ordered_pic_abi = gsub(noisy_things[i],"",parameternames_ordered_pic_abi)
  parameternames_ordered_fag_syl = gsub(noisy_things[i],"",parameternames_ordered_fag_syl)
  parameternames_ordered_pin_syl = gsub(noisy_things[i],"",parameternames_ordered_pin_syl)
  parameters_fag_syl = gsub(noisy_things[i],"",parameters_fag_syl)
  parameters_pin_syl = gsub(noisy_things[i],"",parameters_pin_syl)
  parameters_pic_abi = gsub(noisy_things[i],"",parameters_pic_abi)
  #parameters_only_mixed = gsub(noisy_things[i],"",parameters_only_mixed)
}

order_to_ranges_pic_abi = vector(length = length(parameternames_ordered_pic_abi))
order_to_ranges_fag_syl = vector(length = length(parameternames_ordered_pic_abi))
order_to_ranges_pin_syl = vector(length = length(parameternames_ordered_pic_abi))

for(i in 1:length(parameternames_ordered_fag_syl)){
  position = try(which(parameternames_ordered_pic_abi[i] == parameters_pic_abi))
  if(length(position) ==0){
    order_to_ranges_pic_abi[i] = NA
  }
  else{
    order_to_ranges_pic_abi[i] = position
  }
  position = try(which(parameternames_ordered_fag_syl[i] == parameters_fag_syl))
  if(length(position) ==0){
    order_to_ranges_fag_syl[i] = NA
  }
  else{
    order_to_ranges_fag_syl[i] = position
  }
  position = try(which(parameternames_ordered_pin_syl[i] == parameters_pin_syl))
  if(length(position) ==0){
    order_to_ranges_pin_syl[i] = NA
  }
  else{
    order_to_ranges_pin_syl[i] = position
  }
}

order_to_ranges_pic_abi= order_to_ranges_pic_abi[which(!is.na(order_to_ranges_pic_abi))]
order_to_ranges_pin_syl= order_to_ranges_pin_syl[which(!is.na(order_to_ranges_pin_syl))]
order_to_ranges_fag_syl= order_to_ranges_fag_syl[which(!is.na(order_to_ranges_fag_syl))]



### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects_lin.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects_lin.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects_lin.rds")

effects_mixed = readRDS("LPJrunTest/Results/mixed_effects_lin.rds")

drivernames = c("phchange","tempchange","co2change","insolchange","precchange","ndepchange")

# for(i in 1:length(drivernames)){
#   parameternames_ordered_pic_abi = gsub(drivernames[i],"NA",parameternames_ordered_pic_abi)
#   parameternames_ordered_fag_syl = gsub(drivernames[i],"NA",parameternames_ordered_fag_syl)
#   parameternames_ordered_pin_syl = gsub(drivernames[i],"NA",parameternames_ordered_pin_syl)
# }

#### Analysis for Pic abi ####

climaticranges <- readRDS("./EnvironmentalData/climateranges_list.rds")
climaticranges$temp[1,] = climaticranges$temp[1,] -273.15 # conversion of Kelvin to degree Celsius

mean_climatic_ranges = cbind(climaticranges$ph[1,],
                             climaticranges$co2[1,],
                             climaticranges$prec[1,],
                             climaticranges$temp[1,],
                             climaticranges$insol[1,],
                             climaticranges$ndep[1,])

low_climatic_ranges = cbind(climaticranges$ph[2,],
                            climaticranges$co2[2,],
                            climaticranges$prec[2,],
                            climaticranges$temp[2,],
                            climaticranges$insol[2,],
                            climaticranges$ndep[2,])

high_climatic_ranges = cbind(climaticranges$ph[3,],
                             climaticranges$co2[3,],
                             climaticranges$prec[3,],
                             climaticranges$temp[3,],
                             climaticranges$insol[3,],
                             climaticranges$ndep[3,])

drivernames_short = c("phchange","co2change","precchange","tempchange","insolchange","ndepchange")



parameters$Pic_abi$Maximum.Value = gsub(",",".",parameters$Pic_abi$Maximum.Value)
parameters$Pic_abi$Minimum.Value = gsub(",",".",parameters$Pic_abi$Minimum.Value)
parameters$Pic_abi$Default = gsub(",",".",parameters$Pic_abi$Default)

# function to get the sensitivities by dividing to the range and default value
get_sensi = function(effect_sizes,
                     maximum,
                     minimum,
                     default){

  sensitivities = vector("numeric", length(effect_sizes))
  for(i in 1:length(sensitivities)){
    sensitivities[i] = effect_sizes[i]*abs(mean(c(maximum[i],minimum[i])))/(abs(maximum[i]-minimum[i]))
  }
  return(sensitivities)
}

get_sensi_drivers = function(site_uncertainties, site_low_range, site_high_range, site_mean){
  sensitivity = matrix(nrow = nrow(site_uncertainties), ncol = ncol(site_uncertainties))
  for(site in 1:nrow(site_uncertainties)){
    for(driver in 1:ncol(site_uncertainties)){
      sensitivity[site,driver] = (site_uncertainties[site, driver]*abs(mean(site_low_range[site,driver],site_high_range[site,driver])))/
        (abs(site_high_range[site,driver]-site_low_range[site,driver]))
    }
  }
  return(sensitivity)
}




get_sensi_flux = function(uncertainties, climate_low,climate_mean, climate_high,
                          cpool_flux,parameternames_ordered,drivernames,
                          parameter_max, parameter_min, parameter_default){

  growth_sites = which(cpool_flux[,ncol(cpool_flux)]/200 > 2)
  uncertainties = uncertainties[growth_sites,]
  uncertainties = uncertainties/mean(abs(uncertainties[,ncol(uncertainties)]))
  climate_low = climate_low[growth_sites,]
  climate_high = climate_high[growth_sites,]
  climate_mean = climate_mean[growth_sites,]
  parameter_uncertainties = uncertainties[,which(!(parameternames_ordered %in% drivernames))]
  driver_uncertainties = uncertainties[,which(parameternames_ordered %in% drivernames)]

  averaged_param_uncertainties = apply(FUN = mean, MARGIN =2, X = parameter_uncertainties, na.rm =T)

  averaged_param_sensitivities = get_sensi(averaged_param_uncertainties,
                                           maximum = parameter_max,
                                           minimum = parameter_min,
                                           default = parameter_default)

  driver_sensitivities = get_sensi_drivers(driver_uncertainties,
                                           site_low_range = climate_low,
                                           site_high_range = climate_high ,
                                           site_mean = climate_mean)


  averaged_driver_sensitivities = apply(X = driver_sensitivities,
                                        MARGIN = 2, FUN = mean, na.rm =T)

  averaged_sensitivties = c(averaged_param_sensitivities, averaged_driver_sensitivities)
  names(averaged_sensitivties) = c(parameternames_ordered[which(!(parameternames_ordered %in% drivernames))], drivernames)
  return(averaged_sensitivties)
}



main_sensi_Pic_abi_cpool = get_sensi_flux(uncertainties = effects_Pic_abi[["cpool"]][["complete"]],
               climate_low = low_climatic_ranges,
               climate_mean = mean_climatic_ranges,
               climate_high = high_climatic_ranges,
               cpool_flux = effects_Pic_abi[["cpool"]][["complete"]],
               parameternames_ordered = parameternames_ordered_pic_abi,
               drivernames = drivernames_short,
               parameter_max = as.numeric(parameters$Pic_abi$Maximum.Value[order_to_ranges_pic_abi]),
               parameter_min = as.numeric(parameters$Pic_abi$Minimum.Value[order_to_ranges_pic_abi]),
               parameter_default = as.numeric(parameters$Pic_abi$Default[order_to_ranges_pic_abi]))


main_sensi_Pic_abi_cflux = get_sensi_flux(uncertainties = effects_Pic_abi[["cflux"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Pic_abi[["cpool"]][["complete"]],
                                          parameternames_ordered = parameternames_ordered_pic_abi,
                                          drivernames = drivernames_short,
                                          parameter_max = as.numeric(parameters$Pic_abi$Maximum.Value[order_to_ranges_pic_abi]),
                                          parameter_min = as.numeric(parameters$Pic_abi$Minimum.Value[order_to_ranges_pic_abi]),
                                          parameter_default = as.numeric(parameters$Pic_abi$Default[order_to_ranges_pic_abi]))



main_sensi_Pic_abi_agpp = get_sensi_flux(uncertainties = effects_Pic_abi[["agpp"]][["complete"]],
                                         climate_low = low_climatic_ranges,
                                         climate_mean = mean_climatic_ranges,
                                         climate_high = high_climatic_ranges,
                                         cpool_flux = effects_Pic_abi[["cpool"]][["complete"]],
                                         parameternames_ordered = parameternames_ordered_pic_abi,
                                         drivernames = drivernames_short,
                                         parameter_max = as.numeric(parameters$Pic_abi$Maximum.Value[order_to_ranges_pic_abi]),
                                         parameter_min = as.numeric(parameters$Pic_abi$Minimum.Value[order_to_ranges_pic_abi]),
                                         parameter_default = as.numeric(parameters$Pic_abi$Default[order_to_ranges_pic_abi]))


#### Analysis for Fag syl ####


parameters$Fag_syl$Maximum.Value = gsub(",",".",parameters$Fag_syl$Maximum.Value)
parameters$Fag_syl$Minimum.Value = gsub(",",".",parameters$Fag_syl$Minimum.Value)
parameters$Fag_syl$Default = gsub(",",".",parameters$Fag_syl$Default)



main_sensi_Fag_syl_cpool = get_sensi_flux(uncertainties = effects_Fag_syl[["cpool"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Fag_syl[["cpool"]][["complete"]],
                                          parameternames_ordered = parameternames_ordered_fag_syl,
                                          drivernames = drivernames_short,
                                          parameter_max = as.numeric(parameters$Fag_syl$Maximum.Value[order_to_ranges_fag_syl]),
                                          parameter_min = as.numeric(parameters$Fag_syl$Minimum.Value[order_to_ranges_fag_syl]),
                                          parameter_default = as.numeric(parameters$Fag_syl$Default[order_to_ranges_fag_syl]))



main_sensi_Fag_syl_cflux = get_sensi_flux(uncertainties = effects_Fag_syl[["cflux"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Fag_syl[["cpool"]][["complete"]],
                                          parameternames_ordered = parameternames_ordered_fag_syl,
                                          drivernames = drivernames_short,
                                          parameter_max = as.numeric(parameters$Fag_syl$Maximum.Value[order_to_ranges_fag_syl]),
                                          parameter_min = as.numeric(parameters$Fag_syl$Minimum.Value[order_to_ranges_fag_syl]),
                                          parameter_default = as.numeric(parameters$Fag_syl$Default[order_to_ranges_fag_syl]))


main_sensi_Fag_syl_agpp = get_sensi_flux(uncertainties = effects_Fag_syl[["agpp"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Fag_syl[["cpool"]][["complete"]],
                                         parameternames_ordered = parameternames_ordered_fag_syl,
                                         drivernames = drivernames_short,
                                         parameter_max = as.numeric(parameters$Fag_syl$Maximum.Value[order_to_ranges_fag_syl]),
                                         parameter_min = as.numeric(parameters$Fag_syl$Minimum.Value[order_to_ranges_fag_syl]),
                                         parameter_default = as.numeric(parameters$Fag_syl$Default[order_to_ranges_fag_syl]))




#### Analysis for Pin syl ####

parameters$Pin_syl$Maximum.Value = gsub(",",".",parameters$Pin_syl$Maximum.Value)
parameters$Pin_syl$Minimum.Value = gsub(",",".",parameters$Pin_syl$Minimum.Value)
parameters$Pin_syl$Default = gsub(",",".",parameters$Pin_syl$Default)


main_sensi_Pin_syl_cpool = get_sensi_flux(uncertainties = effects_Pin_syl[["cpool"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Pin_syl[["cpool"]][["complete"]],
                                          parameternames_ordered = parameternames_ordered_pin_syl,
                                          drivernames = drivernames_short,
                                          parameter_max = as.numeric(parameters$Pin_syl$Maximum.Value[order_to_ranges_pin_syl]),
                                          parameter_min = as.numeric(parameters$Pin_syl$Minimum.Value[order_to_ranges_pin_syl]),
                                          parameter_default = as.numeric(parameters$Pin_syl$Default[order_to_ranges_pin_syl]))

main_sensi_Pin_syl_cflux = get_sensi_flux(uncertainties = effects_Pin_syl[["cflux"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Pin_syl[["cpool"]][["complete"]],
                                          parameternames_ordered = parameternames_ordered_pin_syl,
                                          drivernames = drivernames_short,
                                          parameter_max = as.numeric(parameters$Pin_syl$Maximum.Value[order_to_ranges_pin_syl]),
                                          parameter_min = as.numeric(parameters$Pin_syl$Minimum.Value[order_to_ranges_pin_syl]),
                                          parameter_default = as.numeric(parameters$Pin_syl$Default[order_to_ranges_pin_syl]))



main_sensi_Pin_syl_agpp = get_sensi_flux(uncertainties = effects_Pin_syl[["agpp"]][["complete"]],
                                          climate_low = low_climatic_ranges,
                                          climate_mean = mean_climatic_ranges,
                                          climate_high = high_climatic_ranges,
                                          cpool_flux = effects_Pin_syl[["cpool"]][["complete"]],
                                         parameternames_ordered = parameternames_ordered_pin_syl,
                                         drivernames = drivernames_short,
                                         parameter_max = as.numeric(parameters$Pin_syl$Maximum.Value[order_to_ranges_pin_syl]),
                                         parameter_min = as.numeric(parameters$Pin_syl$Minimum.Value[order_to_ranges_pin_syl]),
                                         parameter_default = as.numeric(parameters$Pin_syl$Default[order_to_ranges_pin_syl]))

#### Analysis for mixed ####
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
results =  readRDS(paste0("./../Results_sensi/Mixed_0.25_42.25.rds"))
parametermixed = readRDS("ParameterMetaData/Parameter_mixed.rds")
parameternames_ordered_mixed = names(unlist(results[[1]]@runInfo$parameterList)[which(names(unlist(results[[1]]@runInfo$parameterList)) %in% c(as.character(parametermixed$NameRLPJ),
                                                                                                                                               drivernames))])
rm(results)
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
    sensitivities[,parameter] = uncertainty[,parameter]*abs(mean(c(maximum[name], minimum[name])))/
      (abs(maximum[name]-minimum[name]))
  }
  return(sensitivities)
}


remap_sensi <- function(sensitivities,
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
      mapped_parameters[site,parameter] = sum(sensitivities[site,position_scheme[[parameter]]]*
                                                weighting_scheme[[parameter]])
    }
  }
  mean_mapped_parameters = apply(FUN = mean, X = mapped_parameters, MARGIN = 2,na.rm= T)

  mean_mapped_parameters = mean_mapped_parameters

  return(mean_mapped_parameters)
}



calculate_sensi_mixed = function(uncertainties, cpool_flux, drivernames_short,nameofparameters, nameslist,
                                 parameter_min,parameter_max,parameter_default, climate_low,
                                 climate_high, climate_mean, mixed_results,
                                 drivernames){

  growth_sites = which(cpool_flux[,ncol(cpool_flux)]/200 > 2)
  uncertainties = uncertainties[growth_sites,]
  uncertainties = uncertainties/mean(abs(uncertainties[,ncol(uncertainties)]))
  climate_low = climate_low[growth_sites,]
  climate_high = climate_high[growth_sites,]
  climate_mean = climate_mean[growth_sites,]

  parameter_uncertainties = uncertainties[,which(!(nameslist %in% drivernames))]
  driver_uncertainties = uncertainties[,which(nameslist %in% drivernames)]

  nameslist_without_drivers = nameslist[which(!(nameslist %in% drivernames))]

  parameter_sensitivities = calculate_sensi(parameter_uncertainties,
                                            nameofparameters= nameofparameters ,
                                            nameslist = nameslist_without_drivers,
                                            minimum= parameter_min,
                                            maximum = parameter_max,
                                            default = parameter_default,
                                            reference = uncertainties[,ncol(uncertainties)])


  driver_components = which(gsub("_","",mixed_results$parameternames) %in% drivernames_short)

  remapped_parameter_sensitivities = remap_sensi(sensitivities = parameter_sensitivities,
                                                 weighting_scheme = mixed_results$weight_mapping,
                                                 position_scheme = mixed_results$position_mapping,
                                                 driver_components = driver_components)

  names_parameter = gsub("_","",mixed_results$parameternames)[which(!( gsub("_","",mixed_results$parameternames) %in% drivernames_short))]

  driver_uncertainties = driver_uncertainties/mean(abs(uncertainties[,ncol(uncertainties)]))

  driver_sensitivities = get_sensi_drivers(driver_uncertainties,
                                           site_low_range = climate_low,
                                           site_high_range = climate_high ,
                                           site_mean = climate_mean)

  averaged_driver_sensitivities = apply(X = driver_sensitivities,
                                        MARGIN = 2, FUN = mean, na.rm =T)

  sensitivities = c(remapped_parameter_sensitivities, averaged_driver_sensitivities)
  names(sensitivities) = c(names_parameter, drivernames_short)

  return(sensitivities)

}


main_sensi_mixed_cpool = calculate_sensi_mixed(uncertainties = effects_mixed[["cpool"]][["complete"]],
                                    cpool_flux =  effects_mixed[["cpool"]][["complete"]],
                                    drivernames_short = drivernames_short,
                                    nameofparameters = parametermixed$NameRLPJ,
                                    nameslist = parameternames_ordered_mixed,
                                    parameter_min = parametermixed$Minimum.Value,
                                    parameter_max = parametermixed$Maximum.Value,
                                    parameter_default = parametermixed$Default,
                                    climate_low = low_climatic_ranges,
                                    climate_high = high_climatic_ranges,
                                    climate_mean = mean_climatic_ranges,
                                    mixed_results = mixed_results,
                                    drivernames = drivernames)

main_sensi_mixed_cflux = calculate_sensi_mixed(uncertainties = effects_mixed[["cflux"]][["complete"]],
                                    cpool_flux =  effects_mixed[["cpool"]][["complete"]],
                                    drivernames_short = drivernames_short,
                                    nameofparameters = parametermixed$NameRLPJ,
                                    nameslist = parameternames_ordered_mixed,
                                    parameter_min = parametermixed$Minimum.Value,
                                    parameter_max = parametermixed$Maximum.Value,
                                    parameter_default = parametermixed$Default,
                                    climate_low = low_climatic_ranges,
                                    climate_high = high_climatic_ranges,
                                    climate_mean = mean_climatic_ranges,
                                    mixed_results = mixed_results,
                                    drivernames = drivernames)

main_sensi_mixed_agpp = calculate_sensi_mixed(uncertainties = effects_mixed[["agpp"]][["complete"]],
                                    cpool_flux =  effects_mixed[["cpool"]][["complete"]],
                                    drivernames_short = drivernames_short,
                                    nameofparameters = parametermixed$NameRLPJ,
                                    nameslist = parameternames_ordered_mixed,
                                    parameter_min = parametermixed$Minimum.Value,
                                    parameter_max = parametermixed$Maximum.Value,
                                    parameter_default = parametermixed$Default,
                                    climate_low = low_climatic_ranges,
                                    climate_high = high_climatic_ranges,
                                    climate_mean = mean_climatic_ranges,
                                    mixed_results = mixed_results,
                                    drivernames = drivernames)




grouping = c(as.character(read.csv("ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")$Group), rep("Drivers",6))
names = c(as.character(read.csv("ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")$Name), gsub("change","",drivernames_short))
names = as.character(gsub("_","",names))
grouping2 = grouping[order(grouping)]
names_ordered = gsub(" ","",names[order(grouping)])


order_fag_syl = match(names_ordered,gsub("change","",names(main_sensi_Fag_syl_agpp)))
order_pic_abi = match(names_ordered,gsub("change","",names(main_sensi_Pic_abi_agpp)))
order_pin_syl = match(names_ordered,gsub("change","",names(main_sensi_Pin_syl_agpp)))
order_mixed = match(names_ordered,gsub("change","",names(main_sensi_mixed_agpp)))


mean_sensi_abs = data.frame("Group" = grouping2,
                              "cpool_fag_syl" = main_sensi_Fag_syl_cpool[order_fag_syl],
                              "cflux_fag_syl" = main_sensi_Fag_syl_cflux[order_fag_syl],
                              "agpp_fag_syl" = main_sensi_Fag_syl_agpp[order_fag_syl],
                              "cpool_pin_syl" = main_sensi_Pin_syl_cpool[order_pin_syl],
                              "cflux_pin_syl" = main_sensi_Pin_syl_cflux[order_pin_syl],
                              "agpp_pin_syl" = main_sensi_Pin_syl_agpp[order_pin_syl],
                              "cpool_pic_abi" = main_sensi_Pic_abi_cpool[order_pic_abi],
                              "cflux_pic_abi" = main_sensi_Pic_abi_cflux[order_pic_abi],
                              "agpp_pic_abi" = main_sensi_Pic_abi_agpp[order_pic_abi],
                              "cpool_mixed" = main_sensi_mixed_cpool[order_mixed],
                              "cflux_mixed" = main_sensi_mixed_cflux[order_mixed],
                              "agpp_mixed" = main_sensi_mixed_agpp[order_mixed],
                              "Names" = names_ordered)

mean_sensi_drivers = mean_sensi_abs[mean_sensi_abs$Group == "Drivers",]

mean_sensi_parameters = mean_sensi_abs[mean_sensi_abs$Group != "Drivers",]
mean_sensi_parameters$Group = as.character(mean_sensi_parameters$Group)
mean_sensi_parameters$Group = as.factor(mean_sensi_parameters$Group)

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
spaces = table(mean_sensi_parameters$Group)
colors_full = vector()
for(i in 1:6){
  colors_full = c(colors_full,rep(colors_parameter[i],spaces[i]))
}

pdf("./Figures/Mean_sensi.pdf", width = 18.0, height = 12)


layout(matrix(c(1,2,3,4),byrow =T,nrow =4), heights = c(2,2,2,1))
old_mar = c(5.1,4.1,4.1,2.1)
par(mar = c(0.6,4.1,9.1,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n',
     ylim = c(min(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                    mean_sensi_abs$cpool_mixed))*1.1,max(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                                                           mean_sensi_abs$cpool_mixed))*1.1), cex.main = 1.5)
mtext("Factorial sensitivities [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_sensi_parameters$cpool_fag_syl,mean_sensi_parameters$cpool_pic_abi,mean_sensi_parameters$cpool_pin_syl,
                             mean_sensi_parameters$cpool_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_sensi_parameters$cpool_fag_syl , x= 1:32, pch = 15, col = "darkgreen",
       cex = 1.2)
points(y = mean_sensi_parameters$cpool_pic_abi,
       x = 1:nrow(mean_sensi_parameters), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_parameters$cpool_pin_syl,
       x = 1:nrow(mean_sensi_parameters), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_parameters$cpool_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_sensi_parameters))
points(y = apply(rbind(mean_sensi_parameters$cpool_fag_syl,mean_sensi_parameters$cpool_pic_abi,mean_sensi_parameters$cpool_pin_syl),2,mean),
       x = 1:nrow(mean_sensi_parameters), col = "black", pch = "-", cex=3)

axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(min(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                             mean_sensi_abs$cpool_mixed))*1.1,max(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                               mean_sensi_abs$cpool_mixed))*1.1,length.out = 10),2), las = 2, line = -0.5)
barplot(height = apply(rbind(mean_sensi_drivers$cpool_fag_syl,mean_sensi_drivers$cpool_pic_abi,mean_sensi_drivers$cpool_pin_syl,
                             mean_sensi_drivers$cpool_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= rep("gold2",6), percent = 65)),
        add = T,axes = F, width = 0.85, space = c(33.5/0.85,rep(0.15/0.85,5)))
points(y =mean_sensi_drivers$cpool_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,.15),cex = 1.2)
ablineclip(v = cumsum(table(mean_sensi_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = min(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                      mean_sensi_abs$cpool_mixed))*1.1)
points(y = mean_sensi_drivers$cpool_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_drivers$cpool_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_drivers$cpool_mixed, col = "purple",pch = 8,
       cex = 1.2,x = 34:39)
points(y = apply(rbind(mean_sensi_drivers$cpool_fag_syl,mean_sensi_drivers$cpool_pic_abi,mean_sensi_drivers$cpool_pin_syl),2,mean),
       x = 34:39, col = "black", pch = "-", cex=3)
ablineclip(v = 33,col = alpha("black",0.7),y1 = min(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                                                      mean_sensi_abs$cpool_mixed))*1.1)
spacings = (cumsum(table(mean_sensi_parameters$Group)) + c(0,cumsum(table(mean_sensi_parameters$Group))[-length(cumsum(table(mean_sensi_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(mean_sensi_parameters$Group))
par(xpd =  T)
text(c(names(table(mean_sensi_parameters$Group))[which(names(table(mean_sensi_parameters$Group)) != "Drivers")],"Drivers"),
     x = c(spacings,36) ,  y = max(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                                                          mean_sensi_abs$cpool_mixed))*1.2, font = 2, cex =1.8, srt = 30, pos = 4, offset =  c(0,20))
par(xpd = F)
ablineclip(h = max(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                     mean_sensi_abs$cpool_mixed))*1.1, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=min(c(mean_sensi_abs$cpool_fag_syl,mean_sensi_abs$cpool_pic_abi,mean_sensi_abs$cpool_pin_syl,
                                 mean_sensi_abs$cpool_mixed))*1.1)
title(main = "a)", line = 6, adj = 0.01, cex.main = 2)

par(mar = c(0.6,4.1,5.6,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n',
     ylim = c(min(c(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                    mean_sensi_abs$cflux_mixed))*1.1,max(c(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                                                           mean_sensi_abs$cflux_mixed))*1.1))
mtext("Factorial sensitivities [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_sensi_parameters$cflux_fag_syl,mean_sensi_parameters$cflux_pic_abi,mean_sensi_parameters$cflux_pin_syl,
                             mean_sensi_parameters$cflux_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_sensi_parameters$cflux_fag_syl , x= 1:32, pch = 15, col = "darkgreen",
       cex = 1.2)
points(y = mean_sensi_parameters$cflux_pic_abi,
       x = 1:nrow(mean_sensi_parameters), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_parameters$cflux_pin_syl,
       x = 1:nrow(mean_sensi_parameters), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_parameters$cflux_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_sensi_parameters))
points(y = apply(rbind(mean_sensi_parameters$cflux_fag_syl,mean_sensi_parameters$cflux_pic_abi,mean_sensi_parameters$cflux_pin_syl),2,mean),
       x = 1:nrow(mean_sensi_parameters), col = "black", pch = "-", cex=3)

axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(min(c(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                             mean_sensi_abs$cflux_mixed))*1.1,max(c(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                                                                    mean_sensi_abs$cflux_mixed))*1.1,length.out = 10),2), las = 2, line = -0.5)
barplot(height = apply(rbind(mean_sensi_drivers$cflux_fag_syl,mean_sensi_drivers$cflux_pic_abi,mean_sensi_drivers$cflux_pin_syl,
                             mean_sensi_drivers$cflux_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= rep("gold2",6), percent = 65)),
        add = T,axes = F, width = 0.85, space = c(33.5/0.85,rep(0.15/0.85,5)))
points(y =mean_sensi_drivers$cflux_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,.15),cex = 1.2)
ablineclip(v = cumsum(table(mean_sensi_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = min(c(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                      mean_sensi_abs$cflux_mixed))*1.1)
points(y = mean_sensi_drivers$cflux_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_drivers$cflux_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_drivers$cflux_mixed, col = "purple",pch = 8,
       cex = 1.2,x = 34:39)
points(y = apply(rbind(mean_sensi_drivers$cflux_fag_syl,mean_sensi_drivers$cflux_pic_abi,mean_sensi_drivers$cflux_pin_syl),2,mean),
       x = 34:39, col = "black", pch = "-", cex=3)
ablineclip(v = 33,col = alpha("black",0.7),y1 = min(c(mean_sensi_abs$cflux_fag_syl,mean_sensi_abs$cflux_pic_abi,mean_sensi_abs$cflux_pin_syl,
                                                      mean_sensi_abs$cflux_mixed))*1.1)
par(xpd =T)
legend(x = 40, y = 1.1,
       legend = c("Fag. syl.","Pic. abi.",'Pin. syl.',"Mixed",'Mean Mono'),
       col = c('darkgreen', 'darkblue','brown', "purple",'black'), bty = 'n',
       title = as.expression(bquote(italic(bold("Species")))),
       pch = c(15,19,17,8,NA), cex = 1.8, lty = c(NA,NA,NA,NA,1),
       lwd = 2)
par(xpd = F)
ablineclip(h = 220, x1 =0, x2 = 33.5)
title(main = "b)", line = 2, adj = 0.01, cex.main = 2)

par(mar = c(1.1,4.1,6.1,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n',
     ylim = c(min(c(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                    mean_sensi_abs$agpp_mixed))*1.1,max(c(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                                                           mean_sensi_abs$agpp_mixed))*1.1))
mtext("Factorial sensitivities [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_sensi_parameters$agpp_fag_syl,mean_sensi_parameters$agpp_pic_abi,mean_sensi_parameters$agpp_pin_syl,
                             mean_sensi_parameters$agpp_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_sensi_parameters$agpp_fag_syl , x= 1:32, pch = 15, col = "darkgreen",
       cex = 1.2)
points(y = mean_sensi_parameters$agpp_pic_abi,
       x = 1:nrow(mean_sensi_parameters), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_parameters$agpp_pin_syl,
       x = 1:nrow(mean_sensi_parameters), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_parameters$agpp_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_sensi_parameters))
points(y = apply(rbind(mean_sensi_parameters$agpp_fag_syl,mean_sensi_parameters$agpp_pic_abi,mean_sensi_parameters$agpp_pin_syl),2,mean),
       x = 1:nrow(mean_sensi_parameters), col = "black", pch = "-", cex=3)

axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(min(c(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                             mean_sensi_abs$agpp_mixed))*1.1,max(c(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                                                                    mean_sensi_abs$agpp_mixed))*1.1,length.out = 10),2), las = 2, line = -0.5)
barplot(height = apply(rbind(mean_sensi_drivers$agpp_fag_syl,mean_sensi_drivers$agpp_pic_abi,mean_sensi_drivers$agpp_pin_syl,
                             mean_sensi_drivers$agpp_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= rep("gold2",6), percent = 65)),
        add = T,axes = F, width = 0.85, space = c(33.5/0.85,rep(0.15/0.85,5)))
points(y =mean_sensi_drivers$agpp_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,.15),cex = 1.2)
ablineclip(v = cumsum(table(mean_sensi_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = min(c(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                      mean_sensi_abs$agpp_mixed))*1.1)
points(y = mean_sensi_drivers$agpp_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_sensi_drivers$agpp_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_sensi_drivers$agpp_mixed, col = "purple",pch = 8,
       cex = 1.2,x = 34:39)
points(y = apply(rbind(mean_sensi_drivers$agpp_fag_syl,mean_sensi_drivers$agpp_pic_abi,mean_sensi_drivers$agpp_pin_syl),2,mean),
       x = 34:39, col = "black", pch = "-", cex=3)
ablineclip(v = 33,col = alpha("black",0.7),y1 = min(c(mean_sensi_abs$agpp_fag_syl,mean_sensi_abs$agpp_pic_abi,mean_sensi_abs$agpp_pin_syl,
                                                      mean_sensi_abs$agpp_mixed))*1.1)
title(main = "c)", line = 1.7, adj = 0.01, cex.main = 2)

par(mar = c(1.1,4.1,0.1,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",33)),xaxt='n', bty="n", yaxt = 'n', cex.axis =2,
     ylim = c(-0.01,0.01), cex.main = 1.5)
par(xpd = T)
text(x = 1:39-0.5, y = 0.01, c(as.character(mean_sensi_parameters$Names),"",as.character(mean_sensi_drivers$Names)),srt = 270,
     adj = c(0.0,0.), cex =1.5)
par(xpd = F)

dev.off()


