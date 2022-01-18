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

grouping = read.csv("ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")
grouping = c(as.character(grouping[,1]), rep("Drivers",6))

order_grouping = match(parameternames_ordered_pic_abi2, variablenames2)
grouping2 = grouping[order_grouping]



### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects_lin2.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects_lin2.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects_lin2.rds")

effects_mixed = readRDS("LPJrunTest/Results/Mixed_effects_lin.rds")


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




average_effects_cpool =   (aggregated_effects_Fag_syl_cpool_scaled +
                           aggregated_effects_Pin_syl_cpool_scaled +
                           aggregated_effects_Pic_abi_cpool_scaled +
                           aggregated_effects_mixed_cpool_scaled)/4



#### Plots for temperature gradient #####


### cpool

average_effects_cpool <- cbind(average_effects_cpool, sites[["Temperature [C°]"]])
colnames(average_effects_cpool) = c(colnames(aggregated_effects_Pic_abi_cpool_scaled), "Temperature")

df_effects_cpool_temp = matrix(ncol = 6, nrow = 7)
colnames(df_effects_cpool_temp) = c("minimum","maximum","mean","effectsize","p-value", "R^2")
for(i in 1:7){
  df_effects_cpool_temp[i,1] = min(average_effects_cpool[,i], na.rm = T)
  df_effects_cpool_temp[i,2] = max(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp[i,3] = mean(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp_lm = lm(average_effects_cpool[,i] ~ average_effects_cpool[,8],na.action = na.omit)
  df_effects_cpool_temp[i,4] = coef(df_effects_cpool_temp_lm)[2]
  df_effects_cpool_temp[i,5] = summary(df_effects_cpool_temp_lm)$coefficients[2,4]
  df_effects_cpool_temp[i,6] = summary(df_effects_cpool_temp_lm)$r.squared
}
rownames(df_effects_cpool_temp)= colnames(average_effects_cpool)[1:7]



dev.off()


pdf("Figures/Changing_uncertainty_cpool.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_cpool_temp[,4]*100)*1.5, max(df_effects_cpool_temp[,4])*100)*1.2,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,7),bty="n", type = "n",cex.axis = 1.5)
par(xpd =T)
text(x = min(df_effects_cpool_temp[,4]*100)*1.3, y = 1:7, rownames(df_effects_cpool_temp), font =2, cex = 2.5, pos =2)
par(xpd =F)
abline(v = -0.6)
abline(v =0)
for(i in 1:7){
  arrows(x0 = 0, x1 = df_effects_cpool_temp[i,4]*100, y0 = i, y1=i, lwd =5, length = 0.1)
}
axis(side = 1, at = round(seq(from = -0.6, to = 0.4,length.out = 6),1), cex.axis = 1.5)
mtext("           Change in uncertainty contributions [%/°C]",side =1,line = 3.5, cex =2)

dev.off()


### Same analysis for drivers only

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





average_effects_cpool =   (aggregated_effects_Fag_syl_cpool_scaled +
                             aggregated_effects_Pin_syl_cpool_scaled +
                             aggregated_effects_Pic_abi_cpool_scaled +
                             aggregated_effects_mixed_cpool_scaled)/4



#### Plots for temperature gradient #####


### cpool

average_effects_cpool <- cbind(average_effects_cpool, sites[["Temperature [C°]"]])
colnames(average_effects_cpool) = c(colnames(aggregated_effects_Pic_abi_cpool_scaled), "Temperature")

df_effects_cpool_temp = matrix(ncol = 6, nrow = 6)
colnames(df_effects_cpool_temp) = c("minimum","maximum","mean","effectsize","p-value", "R^2")
for(i in 1:6){
  df_effects_cpool_temp[i,1] = min(average_effects_cpool[,i], na.rm = T)
  df_effects_cpool_temp[i,2] = max(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp[i,3] = mean(average_effects_cpool[,i],na.rm = T)
  df_effects_cpool_temp_lm = lm(average_effects_cpool[,i] ~ average_effects_cpool[,7],na.action = na.omit)
  df_effects_cpool_temp[i,4] = coef(df_effects_cpool_temp_lm)[2]
  df_effects_cpool_temp[i,5] = summary(df_effects_cpool_temp_lm)$coefficients[2,4]
  df_effects_cpool_temp[i,6] = summary(df_effects_cpool_temp_lm)$r.squared
}
rownames(df_effects_cpool_temp)= colnames(average_effects_cpool)[1:6]






pdf("Figures/Changing_uncertainty_cpool_driver.pdf", width = 10, height = 10)
plot(1,xlim=c(min(df_effects_cpool_temp[,4]*100)*1.1, max(df_effects_cpool_temp[,4])*100)*1.1,
     ylab = "", xlab = "", xaxt = "n", yaxt = 'n',  ylim = c(1,6),bty="n", type = "n")
par(xpd =T)
text(x = -0.75, y = 1:6, rownames(df_effects_cpool_temp), font =2, cex = 2.5, pos =2)
par(xpd =F)
abline(v = -0.75)
abline(v =0)
for(i in 1:6){
  arrows(x0 = 0, x1 = df_effects_cpool_temp[i,4]*100, y0 = i, y1=i, lwd =5, length = 0.1)
}
axis(side = 1, at = seq(from = -0.75, to = 0.75,length.out = 7), cex.axis = 1.5)
mtext("Change in uncertainty contributions [%/°C]",side =1,line = 3.5, cex =2.5)

dev.off()







