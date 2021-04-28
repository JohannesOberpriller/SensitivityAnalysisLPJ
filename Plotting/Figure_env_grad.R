#loading the site data information
sites <- readRDS("EnvironmentalData/sites_data.rds")

#function to rescale parameters from 0 to 1
rescaling <- function(x) {(x-min(x))/(max(x)-min(x))}
#function to rescale effects from 0 to 1
rescaling_abs <- function(x) {abs(x)/max(abs(x))}

#function to remap effects of mixed simulations to parameters used for monocultures
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
ordering_mixed = match(parameternames_ordered_pic_abi2,parameternames_ordered_mixed2)

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
drivers = c("ph","co2","prec","temp","insol","ndep")

variablenames = c(as.character(parameters$Pic_abi[,"NameRLPJ"]),drivers)
variablenames2 = substring(variablenames,regexpr("_", variablenames) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  variablenames2 = gsub(noisy_things[i],"",variablenames2)
}

grouping = read.csv(".ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")
grouping = c(as.character(grouping[,1]), rep("Drivers",6))

order_grouping = match(parameternames_ordered_pic_abi2, variablenames2)
grouping2 = grouping[order_grouping]



### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects.rds")

effects_mixed = readRDS("LPJrunTest/Results/Mixed_effects.rds")


library(tidyverse)

## going over all settings and calculating the per process aggregated contributed uncertainty ##


# get the sites, where the species grows
mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
mean_effects_Pic_abi_cpool_scaled = mean_effects_Pic_abi_cpool

# loop over sites, aggregate the effects per process

aggregated_effects_Pic_abi_cpool_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pic_abi){
  aux_df_Pic_abi_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pic_abi_cpool_scaled[,1:38]))[,i])
  aggregated_effects_Pic_abi_cpool_scaled[i,] = (aggregate(effect ~ group, aux_df_Pic_abi_cpool_agg, sum)[,"effect"]
                                                 /sum(abs(aggregate(effect ~ group, aux_df_Pic_abi_cpool_agg, sum)[,"effect"])))
  }else{
    aggregated_effects_Pic_abi_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pic_abi_cpool_scaled) = aggregate(effect ~ group, aux_df_Pic_abi_cpool_agg, sum)[,"group"]


mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]]
mean_effects_Pic_abi_cflux_scaled = mean_effects_Pic_abi_cflux
aggregated_effects_Pic_abi_cflux_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pic_abi){
  aux_df_Pic_abi_cflux_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pic_abi_cflux_scaled[,1:38]))[,i])
  aggregated_effects_Pic_abi_cflux_scaled[i,] = (aggregate(effect ~ group, aux_df_Pic_abi_cflux_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_Pic_abi_cflux_agg, sum)[,"effect"]))
  }else{
    aggregated_effects_Pic_abi_cflux_scaled[i,] =  NA
  }
}
colnames(aggregated_effects_Pic_abi_cflux_scaled) = aggregate(effect ~ group, aux_df_Pic_abi_cflux_agg, sum)[,"group"]


mean_effects_Pic_abi_agpp = effects_Pic_abi[["agpp"]][["complete"]]
mean_effects_Pic_abi_agpp_scaled = mean_effects_Pic_abi_agpp
aggregated_effects_Pic_abi_agpp_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pic_abi){
  aux_df_Pic_abi_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pic_abi_agpp_scaled[,1:38]))[,i])
  aggregated_effects_Pic_abi_agpp_scaled[i,] = (aggregate(effect ~ group, aux_df_Pic_abi_agpp_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_Pic_abi_agpp_agg, sum)[,"effect"]))
  }else{
    aggregated_effects_Pic_abi_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pic_abi_agpp_scaled) = aggregate(effect ~ group, aux_df_Pic_abi_agpp_agg, sum)[,"group"]

#### Fag syl ####




mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
growth_sites_Fag_syl = which(mean_effects_Fag_syl_cpool[,39]/200 > 2)
mean_effects_Fag_syl_cpool_scaled = mean_effects_Fag_syl_cpool
aggregated_effects_Fag_syl_cpool_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Fag_syl){
  aux_df_Fag_syl_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Fag_syl_cpool_scaled[,1:38]))[ordering_fag_syl,i])
  aggregated_effects_Fag_syl_cpool_scaled[i,] = (aggregate(effect ~ group, aux_df_Fag_syl_cpool_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_Fag_syl_cpool_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_Fag_syl_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Fag_syl_cpool_scaled) = aggregate(effect ~ group, aux_df_Fag_syl_cpool_agg, sum)[,"group"]


mean_effects_Fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]]
mean_effects_Fag_syl_cflux_scaled = mean_effects_Fag_syl_cflux
aggregated_effects_Fag_syl_cflux_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Fag_syl){
  aux_df_Fag_syl_cflux_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Fag_syl_cflux_scaled[,1:38]))[ordering_fag_syl,i])
  aggregated_effects_Fag_syl_cflux_scaled[i,] = (aggregate(effect ~ group, aux_df_Fag_syl_cflux_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_Fag_syl_cflux_agg, sum)[,"effect"]))
  }else{
    aggregated_effects_Fag_syl_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Fag_syl_cflux_scaled) = aggregate(effect ~ group, aux_df_Fag_syl_cflux_agg, sum)[,"group"]


mean_effects_Fag_syl_agpp = effects_Fag_syl[["agpp"]][["complete"]]
mean_effects_Fag_syl_agpp_scaled = mean_effects_Fag_syl_agpp
aggregated_effects_Fag_syl_agpp_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Fag_syl){
  aux_df_Fag_syl_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Fag_syl_agpp_scaled[,1:38]))[ordering_fag_syl,i])
  aggregated_effects_Fag_syl_agpp_scaled[i,] = (aggregate(effect ~ group, aux_df_Fag_syl_agpp_agg, sum)[,"effect"]
                                                /sum(aggregate(effect ~ group, aux_df_Fag_syl_agpp_agg, sum)[,"effect"]))
  }else{
    aggregated_effects_Fag_syl_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Fag_syl_agpp_scaled) = aggregate(effect ~ group, aux_df_Fag_syl_agpp_agg, sum)[,"group"]


### Pin syl ###


mean_effects_Pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
growth_sites_Pin_syl = which(mean_effects_Pin_syl_cpool[,39]/200 > 2)
mean_effects_Pin_syl_cpool_scaled = mean_effects_Pin_syl_cpool
aggregated_effects_Pin_syl_cpool_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pin_syl){
  aux_df_Pin_syl_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pin_syl_cpool_scaled[,1:38]))[ordering_pin_syl,i])
  aggregated_effects_Pin_syl_cpool_scaled[i,] = (aggregate(effect ~ group, aux_df_Pin_syl_cpool_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_Pin_syl_cpool_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_Pin_syl_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pin_syl_cpool_scaled) = aggregate(effect ~ group, aux_df_Pin_syl_cpool_agg, sum)[,"group"]


mean_effects_Pin_syl_cflux = effects_Pin_syl[["cflux"]][["complete"]]
mean_effects_Pin_syl_cflux_scaled = mean_effects_Pin_syl_cflux
aggregated_effects_Pin_syl_cflux_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:length(growth_sites_Pin_syl)){
  if(i %in% growth_sites_Pin_syl){
  aux_df_Pin_syl_cflux_agg = data.frame("group" = grouping2, "effect" = t(abs(mean_effects_Pin_syl_cflux_scaled[,1:38]))[ordering_pin_syl,i])
  aggregated_effects_Pin_syl_cflux_scaled[i,] = (aggregate(effect ~ group, aux_df_Pin_syl_cflux_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_Pin_syl_cflux_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_Pin_syl_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pin_syl_cflux_scaled) = aggregate(effect ~ group, aux_df_Pin_syl_cflux_agg, sum)[,"group"]


mean_effects_Pin_syl_agpp = effects_Pin_syl[["agpp"]][["complete"]]
mean_effects_Pin_syl_agpp_scaled = mean_effects_Pin_syl_agpp
aggregated_effects_Pin_syl_agpp_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pin_syl){
  aux_df_Pin_syl_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pin_syl_agpp_scaled[,1:38]))[ordering_pin_syl,i])
  aggregated_effects_Pin_syl_agpp_scaled[i,] = (aggregate(effect ~ group, aux_df_Pin_syl_agpp_agg, sum)[,"effect"]
                                                /sum(aggregate(effect ~ group, aux_df_Pin_syl_agpp_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_Pin_syl_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pin_syl_agpp_scaled) = aggregate(effect ~ group, aux_df_Pin_syl_agpp_agg, sum)[,"group"]


### Mixed stands ###

mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
mean_effects_mixed_cpool = mean_effects_mixed_cpool
mean_effects_mixed_cpool_scaled = remap_paramters(mean_effects_mixed_cpool, mixed_results$weight_mapping, mixed_results$position_mapping)

aggregated_effects_mixed_cpool_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_mixed){
  aux_df_mixed_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_mixed_cpool_scaled[,1:38]))[ordering_mixed,i])
  aggregated_effects_mixed_cpool_scaled[i,] = (aggregate(effect ~ group, aux_df_mixed_cpool_agg, sum)[,"effect"]
                                                 /sum(aggregate(effect ~ group, aux_df_mixed_cpool_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_mixed_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_mixed_cpool_scaled) = aggregate(effect ~ group, aux_df_mixed_cpool_agg, sum)[,"group"]


mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]]
mean_effects_mixed_cflux = mean_effects_mixed_cflux
mean_effects_mixed_cflux_scaled = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)


aggregated_effects_mixed_cflux_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_mixed){
  aux_df_mixed_cflux_agg = data.frame("group" = grouping2, "effect" = t(abs(mean_effects_mixed_cflux_scaled[,1:38]))[ordering_mixed,i])
  aggregated_effects_mixed_cflux_scaled[i,] = (aggregate(effect ~ group, aux_df_mixed_cflux_agg, sum)[,"effect"]
                                               /sum(aggregate(effect ~ group, aux_df_mixed_cflux_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_mixed_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_mixed_cflux_scaled) = aggregate(effect ~ group, aux_df_mixed_cflux_agg, sum)[,"group"]


mean_effects_mixed_agpp = effects_mixed[["agpp"]][["complete"]]
mean_effects_mixed_agpp = mean_effects_mixed_agpp
mean_effects_mixed_agpp_scaled = remap_paramters(mean_effects_mixed_agpp, mixed_results$weight_mapping, mixed_results$position_mapping)

aggregated_effects_mixed_agpp_scaled = matrix(nrow = nrow(sites), ncol = 7)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_mixed){
  aux_df_mixed_agpp_agg = data.frame("group" = grouping2, "effect" = t(abs(mean_effects_mixed_agpp_scaled))[ordering_mixed,i])
  aggregated_effects_mixed_agpp_scaled[i,] = (aggregate(effect ~ group, aux_df_mixed_agpp_agg, sum)[,"effect"]
                                               /sum(aggregate(effect ~ group, aux_df_mixed_agpp_agg, sum)[,"effect"]))
  }
  else{
    aggregated_effects_mixed_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_mixed_agpp_scaled) = aggregate(effect ~ group, aux_df_mixed_agpp_agg, sum)[,"group"]

average_effects_cpool =   (aggregated_effects_Fag_syl_cpool_scaled +
                           aggregated_effects_Pin_syl_cpool_scaled +
                           aggregated_effects_Pic_abi_cpool_scaled +
                           aggregated_effects_mixed_cpool_scaled)/4

average_effects_cflux = (aggregated_effects_Fag_syl_cflux_scaled +
                           aggregated_effects_Pin_syl_cflux_scaled +
                           aggregated_effects_Pic_abi_cflux_scaled+
                           aggregated_effects_mixed_cflux_scaled)/4

average_effects_agpp = (aggregated_effects_Fag_syl_agpp_scaled +
                           aggregated_effects_Pin_syl_agpp_scaled +
                           aggregated_effects_Pic_abi_agpp_scaled +
                          aggregated_effects_mixed_agpp_scaled)/4


#### Plots for temperature gradient #####


### cpool

average_effects_cpool <- cbind(average_effects_cpool, sites[["Temperature [C°]"]])
colnames(average_effects_cpool) = c(colnames(aggregated_effects_Pic_abi_cpool_scaled), "Temperature")

df_effects_cpool_temp = matrix(ncol = 4, nrow = 7)
colnames(df_effects_cpool_temp) = c("minimum","maximum","mean","effectsize")
for(i in 1:7){
  df_effects_cpool_temp[i,1] = min(average_effects_cpool[,i], na.rm = T)
  df_effects_cpool_temp[i,2] = max(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp[i,3] = mean(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp_lm = lm(average_effects_cpool[,i] ~average_effects_cpool[,8],na.action = na.omit)
  df_effects_cpool_temp[i,4] = coef(df_effects_cpool_temp_lm)[2]
}
rownames(df_effects_cpool_temp)= colnames(average_effects_cpool)[1:7]


### cflux


average_effects_cflux <- cbind(average_effects_cflux, sites[["Temperature [C°]"]] )
colnames(average_effects_cflux) = c(colnames(aggregated_effects_Pic_abi_cflux_scaled), "Temperature")


df_effects_cflux_temp = matrix(ncol = 4, nrow = 7)
colnames(df_effects_cflux_temp) = c("minimum","maximum","mean","effectsize")
for(i in 1:7){
  print(i)
  df_effects_cflux_temp[i,1] = min(average_effects_cflux[,i],na.rm = T)
  df_effects_cflux_temp[i,2] = max(average_effects_cflux[,i],na.rm = T)
  df_effects_cflux_temp[i,3] = mean(average_effects_cflux[,i],na.rm = T)
  df_effects_cflux_temp_lm = lm(average_effects_cflux[,i] ~average_effects_cflux[,8])
  df_effects_cflux_temp[i,4] = coef(df_effects_cflux_temp_lm)[2]
}
rownames(df_effects_cflux_temp)= colnames(average_effects_cflux)[1:7]





### agpp

average_effects_agpp <- cbind(average_effects_agpp, sites[["Temperature [C°]"]] )
colnames(average_effects_agpp) = c(colnames(aggregated_effects_Pic_abi_agpp_scaled), "Temperature")

df_effects_agpp_temp = matrix(ncol = 4, nrow = 7)
colnames(df_effects_agpp_temp) = c("minimum","maximum","mean","effectsize")
for(i in 1:7){
  print(i)
  df_effects_agpp_temp[i,1] = min(average_effects_agpp[,i],na.rm = T)
  df_effects_agpp_temp[i,2] = max(average_effects_agpp[,i],na.rm = T)
  df_effects_agpp_temp[i,3] = mean(average_effects_agpp[,i],na.rm = T)
  df_effects_agpp_temp_lm = lm(average_effects_agpp[,i] ~average_effects_agpp[,8])
  df_effects_agpp_temp[i,4] = coef(df_effects_agpp_temp_lm)[2]
}
rownames(df_effects_agpp_temp)= colnames(average_effects_agpp)[1:7]


average_effects_cpool = (aggregated_effects_Fag_syl_cpool_scaled +
                           aggregated_effects_Pin_syl_cpool_scaled +
                           aggregated_effects_Pic_abi_cpool_scaled +
                           aggregated_effects_mixed_cpool_scaled)/4

average_effects_cflux = (aggregated_effects_Fag_syl_cflux_scaled +
                           aggregated_effects_Pin_syl_cflux_scaled +
                           aggregated_effects_Pic_abi_cflux_scaled+
                           aggregated_effects_mixed_cflux_scaled)/4

average_effects_agpp = (aggregated_effects_Fag_syl_agpp_scaled +
                          aggregated_effects_Pin_syl_agpp_scaled +
                          aggregated_effects_Pic_abi_agpp_scaled +
                          aggregated_effects_mixed_agpp_scaled)/4

# #### Precipitation ######
#
#
# ### cpool
#
# average_effects_cpool <- cbind(average_effects_cpool, sites[["Precipitation [l/m^2]"]])
# colnames(average_effects_cpool) = c(colnames(aggregated_effects_Pic_abi_cpool_scaled), "Precipitation")
#
# df_effects_cpool_prec = matrix(ncol = 4, nrow = 7)
# colnames(df_effects_cpool_prec) = c("minimum","maximum","mean","effectsize")
# for(i in 1:7){
#   print(i)
#   df_effects_cpool_prec[i,1] = min(average_effects_cpool[,i],na.rm = T)
#   df_effects_cpool_prec[i,2] = max(average_effects_cpool[,i],na.rm = T)
#   df_effects_cpool_prec[i,3] = mean(average_effects_cpool[,i],na.rm = T)
#   df_effects_cpool_prec_lm = lm(average_effects_cpool[,i] ~log(average_effects_cpool[,8]))
#   df_effects_cpool_prec[i,4] = coef(df_effects_cpool_prec_lm)[2]
# }
# rownames(df_effects_cpool_prec)= colnames(average_effects_cpool)[1:7]
#
#
#
# ### cflux
#
#
# average_effects_cflux <- cbind(average_effects_cflux, sites[["Precipitation [l/m^2]"]])
# colnames(average_effects_cflux) = c(colnames(aggregated_effects_Pic_abi_cflux_scaled), "Precipitation")
#
# df_effects_cflux_prec = matrix(ncol = 4, nrow = 7)
# colnames(df_effects_cflux_prec) = c("minimum","maximum","mean","effectsize")
# for(i in 1:7){
#   print(i)
#   df_effects_cflux_prec[i,1] = min(average_effects_cflux[,i],na.rm = T)
#   df_effects_cflux_prec[i,2] = max(average_effects_cflux[,i],na.rm = T)
#   df_effects_cflux_prec[i,3] = mean(average_effects_cflux[,i],na.rm = T)
#   df_effects_cflux_prec_lm = lm(average_effects_cflux[,i] ~log(average_effects_cflux[,8]))
#   df_effects_cflux_prec[i,4] = coef(df_effects_cflux_prec_lm)[2]
# }
# rownames(df_effects_cflux_prec)= colnames(average_effects_cflux)[1:7]
#
#
# ### agpp
#
# average_effects_agpp <- cbind(average_effects_agpp, sites[["Precipitation [l/m^2]"]])
# colnames(average_effects_agpp) = c(colnames(aggregated_effects_Pic_abi_agpp_scaled), "Precipitation")
#
# df_effects_agpp_prec = matrix(ncol = 4, nrow = 7)
# colnames(df_effects_agpp_prec) = c("minimum","maximum","mean","effectsize")
# for(i in 1:7){
#   print(i)
#   df_effects_agpp_prec[i,1] = min(average_effects_agpp[,i],na.rm = T)
#   df_effects_agpp_prec[i,2] = max(average_effects_agpp[,i],na.rm = T)
#   df_effects_agpp_prec[i,3] = mean(average_effects_agpp[,i],na.rm = T)
#   df_effects_agpp_prec_lm = lm(average_effects_agpp[,i] ~ log(average_effects_agpp[,8]))
#   df_effects_agpp_prec[i,4] = coef(df_effects_agpp_prec_lm)[2]
# }
# rownames(df_effects_agpp_prec)= colnames(average_effects_agpp)[1:7]

pdf("Figures/Changing_uncertainty_cflux.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_cflux_temp[,4]*100)*1.5, max(df_effects_cflux_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,7),bty="n", type = "n")
par(xpd =T)
text(x = min(df_effects_cflux_temp[,4]*100)*1.3, y = 1:7, rownames(df_effects_cflux_temp), font =2, cex = 1.5, pos =2)
par(xpd =F)
abline(v = 1.25*min(df_effects_cflux_temp[,4])*100)
abline(v =0)
for(i in 1:7){
  arrows(x0 = 0, x1 = df_effects_cflux_temp[i,4]*100, y0 = i, y1=i, lwd =3, length = 0.1)
}
axis(side = 1, at = seq(from = -0.8, to = 0.8,length.out = 5))
mtext("Change in uncertatinty contributions on a temperature gradient [%/°C]",side =1,line = 2.5, cex =1.5)

dev.off()


pdf("Figures/Changing_uncertainty_agpp.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_agpp_temp[,4]*100)*1.7, max(df_effects_agpp_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,7),bty="n", type = "n")
par(xpd =T)
text(x = -0.8, y = 1:7, rownames(df_effects_agpp_temp), font =2, cex = 1.5, pos =2)
par(xpd =F)
abline(v = -0.8)
abline(v =0)
for(i in 1:7){
  arrows(x0 = 0, x1 = df_effects_agpp_temp[i,4]*100, y0 = i, y1=i, lwd =3, length = 0.1)
}
axis(side = 1, at = seq(from = -0.8, to = 1.2,length.out = 6))
mtext("Change in uncertatinty contributions on a temperature gradient [%/°C]",side =1,line = 2.5, cex =1.5)

dev.off()



pdf("Figures/Changing_uncertainty_cpool.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_cpool_temp[,4]*100)*1.5, max(df_effects_cpool_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,7),bty="n", type = "n")
par(xpd =T)
text(x = min(df_effects_cpool_temp[,4]*100)*1.3, y = 1:7, rownames(df_effects_cpool_temp), font =2, cex = 1.5, pos =2)
par(xpd =F)
abline(v = -0.5)
abline(v =0)
for(i in 1:7){
  arrows(x0 = 0, x1 = df_effects_cpool_temp[i,4]*100, y0 = i, y1=i, lwd =3, length = 0.1)
}
axis(side = 1, at = seq(from = -0.5, to = 0.5,length.out = 3))
mtext("Change in uncertatinty contributions on a temperature gradient [%/°C]",side =1,line = 2.5, cex =1.5)

dev.off()




aggregated_effects_Pic_abi_cpool_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pic_abi){
    aux_df_Pic_abi_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pic_abi_cpool_scaled[,1:38]))[,i])
    aggregated_effects_Pic_abi_cpool_scaled[i,] = aux_df_Pic_abi_cpool_agg[which(aux_df_Pic_abi_cpool_agg$group =="Drivers"),"effect"]/
                                                   sum(aux_df_Pic_abi_cpool_agg[,"effect"])
  }else{
    aggregated_effects_Pic_abi_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pic_abi_cpool_scaled) = c("ph","co2","prec","temp","insol","ndep")


aggregated_effects_Pic_abi_cflux_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pic_abi){
    aux_df_Pic_abi_cflux_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pic_abi_cflux_scaled[,1:38]))[,i])
    aggregated_effects_Pic_abi_cflux_scaled[i,] = aux_df_Pic_abi_cflux_agg[which(aux_df_Pic_abi_cflux_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Pic_abi_cflux_agg[,"effect"])
  }else{
    aggregated_effects_Pic_abi_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pic_abi_cflux_scaled) = c("ph","co2","prec","temp","insol","ndep")


aggregated_effects_Pic_abi_agpp_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pic_abi){
    aux_df_Pic_abi_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pic_abi_agpp_scaled[,1:38]))[,i])
    aggregated_effects_Pic_abi_agpp_scaled[i,] = aux_df_Pic_abi_agpp_agg[which(aux_df_Pic_abi_agpp_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Pic_abi_agpp_agg[,"effect"])
  }else{
    aggregated_effects_Pic_abi_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pic_abi_agpp_scaled) = c("ph","co2","prec","temp","insol","ndep")

#### Fag syl ####



aggregated_effects_Fag_syl_cpool_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Fag_syl){
    aux_df_Fag_syl_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Fag_syl_cpool_scaled[,1:38]))[ordering_fag_syl,i])
    aggregated_effects_Fag_syl_cpool_scaled[i,] = aux_df_Fag_syl_cpool_agg[which(aux_df_Fag_syl_cpool_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Fag_syl_cpool_agg[,"effect"])
  }else{
    aggregated_effects_Fag_syl_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Fag_syl_cpool_scaled) = c("ph","co2","prec","temp","insol","ndep")


aggregated_effects_Fag_syl_cflux_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Fag_syl){
    aux_df_Fag_syl_cflux_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Fag_syl_cflux_scaled[,1:38]))[ordering_fag_syl,i])
    aggregated_effects_Fag_syl_cflux_scaled[i,] = aux_df_Fag_syl_cflux_agg[which(aux_df_Fag_syl_cflux_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Fag_syl_cflux_agg[,"effect"])
  }else{
    aggregated_effects_Fag_syl_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Fag_syl_cflux_scaled) = c("ph","co2","prec","temp","insol","ndep")


aggregated_effects_Fag_syl_agpp_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Fag_syl){
    aux_df_Fag_syl_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Fag_syl_agpp_scaled[,1:38]))[ordering_fag_syl,i])
    aggregated_effects_Fag_syl_agpp_scaled[i,] = aux_df_Fag_syl_agpp_agg[which(aux_df_Fag_syl_agpp_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Fag_syl_agpp_agg[,"effect"])
  }else{
    aggregated_effects_Fag_syl_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Fag_syl_agpp_scaled) = c("ph","co2","prec","temp","insol","ndep")

### Pin syl ###


aggregated_effects_Pin_syl_cpool_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pin_syl){
    aux_df_Pin_syl_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pin_syl_cpool_scaled[,1:38]))[ordering_pin_syl,i])
    aggregated_effects_Pin_syl_cpool_scaled[i,] = aux_df_Pin_syl_cpool_agg[which(aux_df_Pin_syl_cpool_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Pin_syl_cpool_agg[,"effect"])
  }else{
    aggregated_effects_Pin_syl_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pin_syl_cpool_scaled) = c("ph","co2","prec","temp","insol","ndep")


aggregated_effects_Pin_syl_cflux_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pin_syl){
    aux_df_Pin_syl_cflux_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pin_syl_cflux_scaled[,1:38]))[ordering_pin_syl,i])
    aggregated_effects_Pin_syl_cflux_scaled[i,] = aux_df_Pin_syl_cflux_agg[which(aux_df_Pin_syl_cflux_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Pin_syl_cflux_agg[,"effect"])
  }else{
    aggregated_effects_Pin_syl_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pin_syl_cflux_scaled) = c("ph","co2","prec","temp","insol","ndep")


aggregated_effects_Pin_syl_agpp_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_Pin_syl){
    aux_df_Pin_syl_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_Pin_syl_agpp_scaled[,1:38]))[ordering_pin_syl,i])
    aggregated_effects_Pin_syl_agpp_scaled[i,] = aux_df_Pin_syl_agpp_agg[which(aux_df_Pin_syl_agpp_agg$group =="Drivers"),"effect"]/
      sum(aux_df_Pin_syl_agpp_agg[,"effect"])
  }else{
    aggregated_effects_Pin_syl_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_Pin_syl_agpp_scaled) = c("ph","co2","prec","temp","insol","ndep")

### Mixed stands ###

mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
mean_effects_mixed_cpool = mean_effects_mixed_cpool
mean_effects_mixed_cpool_scaled = remap_paramters(mean_effects_mixed_cpool, mixed_results$weight_mapping, mixed_results$position_mapping)

aggregated_effects_mixed_cpool_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_mixed){
    aux_df_mixed_cpool_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_mixed_cpool_scaled[,1:38]))[ordering_mixed,i])
    aggregated_effects_mixed_cpool_scaled[i,] =  aux_df_mixed_cpool_agg$effect[which(aux_df_mixed_cpool_agg$group == "Drivers")]/
      sum(aux_df_mixed_cpool_agg$effect)
  }
  else{
    aggregated_effects_mixed_cpool_scaled[i,] = NA
  }
}
colnames(aggregated_effects_mixed_cpool_scaled) = c("ph","co2","prec","temp","insol","ndep")


mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]]
mean_effects_mixed_cflux = mean_effects_mixed_cflux
mean_effects_mixed_cflux_scaled = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)


aggregated_effects_mixed_cflux_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_mixed){
    aux_df_mixed_cflux_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_mixed_cflux_scaled[,1:38]))[ordering_mixed,i])
    aggregated_effects_mixed_cflux_scaled[i,] =  aux_df_mixed_cflux_agg$effect[which(aux_df_mixed_cflux_agg$group == "Drivers")]/
      sum(aux_df_mixed_cflux_agg$effect)
  }
  else{
    aggregated_effects_mixed_cflux_scaled[i,] = NA
  }
}
colnames(aggregated_effects_mixed_cflux_scaled) = c("ph","co2","prec","temp","insol","ndep")


mean_effects_mixed_agpp = effects_mixed[["agpp"]][["complete"]]
mean_effects_mixed_agpp = mean_effects_mixed_agpp
mean_effects_mixed_agpp_scaled = remap_paramters(mean_effects_mixed_agpp, mixed_results$weight_mapping, mixed_results$position_mapping)

aggregated_effects_mixed_agpp_scaled = matrix(nrow = nrow(sites), ncol = 6)
for(i in 1:nrow(sites)){
  if(i %in% growth_sites_mixed){
    aux_df_mixed_agpp_agg = data.frame("group" = grouping2, "effect" = abs(t(mean_effects_mixed_agpp_scaled[,1:38]))[ordering_mixed,i])
    aggregated_effects_mixed_agpp_scaled[i,] =  aux_df_mixed_agpp_agg$effect[which(aux_df_mixed_agpp_agg$group == "Drivers")]/
      sum(aux_df_mixed_agpp_agg$effect)
  }
  else{
    aggregated_effects_mixed_agpp_scaled[i,] = NA
  }
}
colnames(aggregated_effects_mixed_agpp_scaled) = c("ph","co2","prec","temp","insol","ndep")


average_effects_cpool =   (aggregated_effects_Fag_syl_cpool_scaled +
                             aggregated_effects_Pin_syl_cpool_scaled +
                             aggregated_effects_Pic_abi_cpool_scaled +
                             aggregated_effects_mixed_cpool_scaled)/4

average_effects_cflux = (aggregated_effects_Fag_syl_cflux_scaled +
                           aggregated_effects_Pin_syl_cflux_scaled +
                           aggregated_effects_Pic_abi_cflux_scaled+
                           aggregated_effects_mixed_cflux_scaled)/4

average_effects_agpp = (aggregated_effects_Fag_syl_agpp_scaled +
                          aggregated_effects_Pin_syl_agpp_scaled +
                          aggregated_effects_Pic_abi_agpp_scaled +
                          aggregated_effects_mixed_agpp_scaled)/4


#### Plots for temperature gradient #####


### cpool

average_effects_cpool <- cbind(average_effects_cpool, sites[["Temperature [C°]"]])
colnames(average_effects_cpool) = c(colnames(aggregated_effects_Pic_abi_cpool_scaled), "Temperature")

df_effects_cpool_temp = matrix(ncol = 4, nrow = 6)
colnames(df_effects_cpool_temp) = c("minimum","maximum","mean","effectsize")
for(i in 1:6){
  df_effects_cpool_temp[i,1] = min(average_effects_cpool[,i], na.rm = T)
  df_effects_cpool_temp[i,2] = max(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp[i,3] = mean(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp_lm = lm(average_effects_cpool[,i] ~average_effects_cpool[,7],na.action = na.omit)
  df_effects_cpool_temp[i,4] = coef(df_effects_cpool_temp_lm)[2]
}
rownames(df_effects_cpool_temp)= colnames(average_effects_cpool)[1:6]


### cflux


average_effects_cflux <- cbind(average_effects_cflux, sites[["Temperature [C°]"]] )
colnames(average_effects_cflux) = c(colnames(aggregated_effects_Pic_abi_cflux_scaled), "Temperature")


df_effects_cflux_temp = matrix(ncol = 4, nrow = 6)
colnames(df_effects_cflux_temp) = c("minimum","maximum","mean","effectsize")
for(i in 1:6){
  print(i)
  df_effects_cflux_temp[i,1] = min(average_effects_cflux[,i],na.rm = T)
  df_effects_cflux_temp[i,2] = max(average_effects_cflux[,i],na.rm = T)
  df_effects_cflux_temp[i,3] = mean(average_effects_cflux[,i],na.rm = T)
  df_effects_cflux_temp_lm = lm(average_effects_cflux[,i] ~average_effects_cflux[,7])
  df_effects_cflux_temp[i,4] = coef(df_effects_cflux_temp_lm)[2]
}
rownames(df_effects_cflux_temp)= colnames(average_effects_cflux)[1:6]





### agpp

average_effects_agpp <- cbind(average_effects_agpp, sites[["Temperature [C°]"]] )
colnames(average_effects_agpp) = c(colnames(aggregated_effects_Pic_abi_agpp_scaled), "Temperature")

df_effects_agpp_temp = matrix(ncol = 4, nrow = 6)
colnames(df_effects_agpp_temp) = c("minimum","maximum","mean","effectsize")
for(i in 1:6){
  print(i)
  df_effects_agpp_temp[i,1] = min(average_effects_agpp[,i],na.rm = T)
  df_effects_agpp_temp[i,2] = max(average_effects_agpp[,i],na.rm = T)
  df_effects_agpp_temp[i,3] = mean(average_effects_agpp[,i],na.rm = T)
  df_effects_agpp_temp_lm = lm(average_effects_agpp[,i] ~average_effects_agpp[,7])
  df_effects_agpp_temp[i,4] = coef(df_effects_agpp_temp_lm)[2]
}
rownames(df_effects_agpp_temp)= colnames(average_effects_agpp)[1:6]


pdf("Figures/Changing_uncertainty_cflux_driver.pdf", width = 10, height = 10)
plot(1,xlim=c(-0.3, max(df_effects_cflux_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,6),bty="n", type = "n")
par(xpd =T)
text(x = -0.2, y = 1:6, rownames(df_effects_cflux_temp), font =2, cex = 1.5, pos =2)
par(xpd =F)
abline(v = -0.2)
abline(v =0)
for(i in 1:6){
  arrows(x0 = 0, x1 = df_effects_cflux_temp[i,4]*100, y0 = i, y1=i, lwd =3, length = 0.1)
}
axis(side = 1, at = seq(from = -0.2, to = 0.6,length.out = 5))
mtext("Change in uncertatinty contributions on a temperature gradient [%/°C]",side =1,line = 2.5, cex =1.5)

dev.off()


pdf("Figures/Changing_uncertainty_agpp_driver.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_agpp_temp[,4]*100)*1.7, max(df_effects_agpp_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,6),bty="n", type = "n")
par(xpd =T)
text(x = -0.8, y = 1:6, rownames(df_effects_agpp_temp), font =2, cex = 1.5, pos =2)
par(xpd =F)
abline(v = -0.8)
abline(v =0)
for(i in 1:6){
  arrows(x0 = 0, x1 = df_effects_agpp_temp[i,4]*100, y0 = i, y1=i, lwd =3, length = 0.1)
}
axis(side = 1, at = seq(from = -0.8, to = 2.4,length.out = 5))
mtext("Change in uncertatinty contributions on a temperature gradient [%/°C]",side =1,line = 2.5, cex =1.5)

dev.off()



pdf("Figures/Changing_uncertainty_cpool_driver.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_cpool_temp[,4]*100)*1.5, max(df_effects_cpool_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,6),bty="n", type = "n")
par(xpd =T)
text(x = -1, y = 1:6, rownames(df_effects_cpool_temp), font =2, cex = 1.5, pos =2)
par(xpd =F)
abline(v = -1)
abline(v =0)
for(i in 1:6){
  arrows(x0 = 0, x1 = df_effects_cpool_temp[i,4]*100, y0 = i, y1=i, lwd =3, length = 0.1)
}
axis(side = 1, at = seq(from = -1, to = 2,length.out = 7))
mtext("Change in uncertatinty contributions on a temperature gradient [%/°C]",side =1,line = 2.5, cex =1.5)

dev.off()



library(RColorBrewer)
library(ggplot2)
library(grid)
library(plotrix)

colorings = colors_me = c("gold","burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")

dev.off()
pdf("./Figures/environmental_gradient_analysis.pdf", width = 14,
    height = 11)

layout(matrix(c(1,1,4,4,
                1,1,4,4,
                1,1,4,4,
                2,2,5,5,
                3,3,6,6
), nrow = 5,byrow =T))


par(mar=c(8.1,4.1,4.1,8.1))
ordering_for_plots = c(1,2,3,5,6,4,7)
df_effects_cflux_temp = df_effects_cflux_temp[ordering_for_plots,]

plot(c(1,7),range(df_effects_cflux_temp),type="n",xlab="",ylab="Uncertainty contributions [%]",xaxt="n",las=1, bty ="n", cex.lab =1.8)
title("a)", cex.main =3, adj = 0)
axis(1,1:7,rownames(df_effects_cflux_temp), las =2,pos =0,cex.axis = 1.5)
rect(xleft = seq(0.8,9.8,1),ybottom = df_effects_cflux_temp[,"minimum"],xright = seq(1.2,10.2,1),ytop = df_effects_cflux_temp[,"maximum"],
     col = colorings[ordering_for_plots])
points(1:7,df_effects_cflux_temp[,"mean"],bg="white",pch=21,type="p", cex = 2)
arrows(x0=1:7,y0=df_effects_cflux_temp[,"mean"],x1=1:7,y1=df_effects_cflux_temp[,"mean"]+df_effects_cflux_temp[,"effectsize"]*5,
       angle=45,code=2,length=0.1, col = "black", lwd =3)
ablineclip(v = 1.45, y1 = 0)
ablineclip(v = 1.55, y1 = 0)
par(mar=c(5.1,4.1,4.1,2.1))


par(mar=c(2.1,4.1,4.1,8.1))

df_effects_cpool_temp = df_effects_cpool_temp[ordering_for_plots,]

plot(c(1,7),range(df_effects_cpool_temp),type="n",xlab="",ylab="Uncertainty contributions [%]",xaxt="n",las=1, bty ="n",cex.lab = 1.2)
title("b)", cex.main =3, adj = 0., line =2)
rect(xleft = seq(0.8,9.8,1),ybottom = df_effects_cpool_temp[,"minimum"],xright = seq(1.2,10.2,1),ytop = df_effects_cpool_temp[,"maximum"],
     col = colorings[ordering_for_plots])
points(1:7,df_effects_cpool_temp[,"mean"],bg="white",pch=21,type="p", cex = 2)
ablineclip(v = 1.45, y1 = 0)
ablineclip(v = 1.55, y1 = 0)
arrows(x0=1:7,y0=df_effects_cpool_temp[,"mean"],x1=1:7,y1=df_effects_cpool_temp[,"mean"]+df_effects_cpool_temp[,"effectsize"]*5,angle=45,code=2,length=0.1, col = "black",
       lwd =3)

par(mar=c(5.1,4.1,4.1,8.1))

par(mar=c(2.1,4.1,4.1,8.1))
df_effects_agpp_temp = df_effects_agpp_temp[ordering_for_plots,]
plot(c(1,7),range(df_effects_agpp_temp),type="n",xlab="",ylab="Uncertainty contributions [%]",xaxt="n",las=1, bty ="n", cex.lab = 1.2)
title("c)", cex.main =3, adj = 0, line =2)
rect(xleft = seq(0.8,9.8,1),ybottom = df_effects_agpp_temp[,"minimum"],xright = seq(1.2,10.2,1),ytop = df_effects_agpp_temp[,"maximum"],
     col = colorings[ordering_for_plots])
points(1:7,df_effects_agpp_temp[,"mean"],bg="white",pch=21,type="p", cex = 2)
ablineclip(v = 1.45, y1 = 0)
ablineclip(v = 1.55, y1 = 0)
arrows(x0=1:7,y0=df_effects_agpp_temp[,"mean"],x1=1:7,y1=df_effects_agpp_temp[,"mean"]+df_effects_agpp_temp[,"effectsize"]*5,angle=45,code=2,length=0.1, col = "black",
       lwd =3)

par(mar=c(5.1,4.1,4.1,8.1))

par(mar=c(6.1,4.1,4.1,8.1))
###precipitation plot
df_effects_cflux_prec = df_effects_cflux_prec[ordering_for_plots,]
plot(c(1,7),range(df_effects_cflux_prec),type="n",xlab="",ylab="Uncertainty contributions [%]",xaxt="n",las=1, bty ="n", cex.lab =1.8)
title("d)", cex.main =3, adj = 0)
axis(1,1:7,rownames(df_effects_cflux_prec), las =2,pos =0, cex.axis =1.5)
rect(xleft = seq(0.8,9.8,1),ybottom = df_effects_cflux_prec[,"minimum"],xright = seq(1.2,10.2,1),ytop = df_effects_cflux_prec[,"maximum"],
     col = colorings[ordering_for_plots])
points(1:7,df_effects_cflux_prec[,"mean"],bg="white",pch=21,type="p", cex = 2)
ablineclip(v = 1.45, y1 = 0)
ablineclip(v = 1.55, y1 = 0)
arrows(x0=1:7,y0=df_effects_cflux_prec[,"mean"],x1=1:7,y1=df_effects_cflux_prec[,"mean"]+df_effects_cflux_prec[,"effectsize"]*5,angle=45,code=2,length=0.1, col = "black",
       lwd =3)
par(xpd=T)
rect(xleft= 7.4,xright = 7.5, ybottom = 0.25, ytop = 0.35,
     col = "grey" )
points(x = 7.45, y = 0.28 ,bg = "white", pch =21, cex =1)
arrows(x0 = 7.45, y0 =0.28, x1= 7.45, y1 =0.26, lwd =1, length = 0.1)
lines(x= c(7.6,7.7), y = c(0.25,0.25))
text(x = 8.1, y = 0.25, "minimum", cex = 1)
lines(x= c(7.6,7.7), y = c(0.35,0.35))
text(x = 8.1, y = 0.35, "maximum", cex = 1.)
lines(x= c(7.6,7.7), y = c(0.28,0.28))
text(x = 8., y = 0.28, "mean", cex = 1)
#arrows(x0 = 0.52, x1= 0.52, y0=0.5,y1=0.7, lwd = 3 )
text(x = 7.3,y = 0.27,'{', cex = 1.5,  lwd = 0.11)
text(x = 6.3, y = 0.275, "relative effect size of", cex = 1.1)
text(x = 6.3, y= 0.265, "change in explained sensitivity",cex = 1.1)
arrows(x0= 7.3, y0= 0.33,x1=7.3,y1=0.35, lwd = 1., length = 0.1)
text(x = 6.3, y = 0.345,"sensitivity increases with", cex =1.1)
text(x = 6.3, y = 0.335, "increasing driver",cex = 1.1)
arrows(x0= 7.3, y0= 0.32,x1=7.3,y1=0.30, lwd = 1,length = 0.1)
text(x = 6.3, y = 0.315,"sensitivity decreases with", cex =1.1)
text(x = 6.3, y = 0.305, "increasing driver", cex =1.1)
par(xpd =F)

par(mar=c(5.1,4.1,4.1,8.1))


par(mar=c(1.1,4.1,4.1,8.1))
df_effects_cpool_prec = df_effects_cpool_prec[ordering_for_plots,]
plot(c(1,7),range(df_effects_cpool_prec),type="n",xlab="",ylab="Uncertainty contributions [%]",xaxt="n",las=1, bty ="n", cex.lab = 1.2)
title("e)", cex.main =3, adj = 0, line =2)
rect(xleft = seq(0.8,9.8,1),ybottom = df_effects_cpool_prec[,"minimum"],xright = seq(1.2,10.2,1),ytop = df_effects_cpool_prec[,"maximum"],
     col = colorings[ordering_for_plots])
points(1:7,df_effects_cpool_prec[,"mean"],bg="white",pch=21,type="p", cex = 2)
ablineclip(v = 1.45, y1 = 0)
ablineclip(v = 1.55, y1 = 0)
arrows(x0=1:7,y0=df_effects_cpool_prec[,"mean"],x1=1:7,y1=df_effects_cpool_prec[,"mean"]+df_effects_cpool_prec[,"effectsize"]*5,angle=45,code=2,length=0.1, col = "black",
       lwd =3)

par(mar=c(5.1,4.1,4.1,8.1))

par(mar=c(2.1,4.1,4.1,8.1))
df_effects_agpp_prec = df_effects_agpp_temp[ordering_for_plots,]
plot(c(1,7),range(df_effects_agpp_prec),type="n",xlab="",ylab="Uncertainty contributions [%]",xaxt="n",las=1, bty ="n", cex.lab = 1.2)
title("f)", cex.main =3, adj = 0, line =2)
rect(xleft = seq(0.8,9.8,1),ybottom = df_effects_agpp_prec[,"minimum"],xright = seq(1.2,10.2,1),ytop = df_effects_agpp_prec[,"maximum"],
     col = colorings[ordering_for_plots])
points(1:7,df_effects_agpp_prec[,"mean"],bg="white",pch=21,type="p", cex = 2)
ablineclip(v = 1.45, y1 = 0)
ablineclip(v = 1.55, y1 = 0)
arrows(x0=1:7,y0=df_effects_agpp_prec[,"mean"],x1=1:7,y1=df_effects_agpp_prec[,"mean"]+df_effects_agpp_prec[,"effectsize"]*5,angle=45,code=2,length=0.1, col = "black",
       lwd =3)

par(mar=c(5.1,4.1,4.1,8.1))

dev.off()



