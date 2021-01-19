# function to rescale the parameters to 0-1
rescaling <- function(x) {(x-min(x))/(max(x)-min(x))}
# function to rescale the results to 0-1
rescaling_abs <- function(x) {abs(x)/max(abs(x))}

# function to reamp the parameters according to a scheme with weights and postion
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

# remap of parameters based on relative effects
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

# function to make colors transparent
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
## getting the ordered parameter names ##

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameters_pic_abi = c(as.character(parameters$Pic_abi$NameRLPJ), drivernames)
parameters_fag_syl = c(as.character(parameters$Fag_syl$NameRLPJ), drivernames)
parameters_pin_syl = c(as.character(parameters$Pin_syl$NameRLPJ), drivernames)

results =  readRDS(paste0("./LPJrunTest/Results/Pic_abi_0.25_42.25.rds"))
parameternames_ordered_pic_abi = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pic_abi)])

results =  readRDS(paste0("./LPJrunTest/Results/Fag_syl_0.25_42.25.rds"))
parameternames_ordered_fag_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_fag_syl)])

results =  readRDS(paste0("./LPJrunTest/Results/Pin_syl_0.25_42.25.rds"))
parameternames_ordered_pin_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pin_syl)])

mixed_results = readRDS("./LPJrunTest/Results/mapping_mixed.rds")


## now we go over each output, time and setting and calcuclate
## the mean effects averaged over all sites
## where we make a cutoff for sites with an biomass increase smaller 2 tC/ha







### Loading the results of the random forest ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects_lin.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects_lin.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects_lin.rds")

effects_mixed = readRDS("LPJrunTest/Results/mixed_effects_lin.rds")



#### Analysis for Pic abi ####


mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool_rf"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
mean_effects_Pic_abi_cpool = mean_effects_Pic_abi_cpool[growth_sites_Pic_abi,]

main_effects_Pic_abi_cpool = apply(abs(mean_effects_Pic_abi_cpool[,1:38]),2,mean)
main_effects_Pic_abi_cpool_weight_abs = main_effects_Pic_abi_cpool/sum(abs(main_effects_Pic_abi_cpool))


main_effects_Pic_abi_cpool = apply(mean_effects_Pic_abi_cpool[,1:38],2,mean)
main_effects_Pic_abi_cpool_weight = main_effects_Pic_abi_cpool/sum(abs(main_effects_Pic_abi_cpool))


mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux_rf"]][["complete"]]
mean_effects_Pic_abi_cflux = mean_effects_Pic_abi_cflux[growth_sites_Pic_abi,]


main_effects_Pic_abi_cflux = apply(abs(mean_effects_Pic_abi_cflux[,1:38]),2,mean)
main_effects_Pic_abi_cflux_weight_abs = main_effects_Pic_abi_cflux/sum(abs(main_effects_Pic_abi_cflux))
#sd_Pic_abic_cmass_abs = apply(abs(mean_effects_Pic_abi_cmass[,1:38]),2,sd)

main_effects_Pic_abi_cflux = apply(mean_effects_Pic_abi_cflux[,1:38],2,mean)
main_effects_Pic_abi_cflux_weight = main_effects_Pic_abi_cflux/sum(abs(main_effects_Pic_abi_cflux))


mean_effects_Pic_abi_agpp = effects_Pic_abi[["agpp_rf"]][["complete"]]
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



mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool_rf"]][["complete"]]
growth_sites_Fag_syl = which(mean_effects_Fag_syl_cpool[,39]/200 > 2)
mean_effects_Fag_syl_cpool = mean_effects_Fag_syl_cpool[growth_sites_Fag_syl,]

main_effects_Fag_syl_cpool = apply(abs(mean_effects_Fag_syl_cpool[,1:38]),2,mean)
main_effects_Fag_syl_cpool_weight_abs = main_effects_Fag_syl_cpool/sum(abs(main_effects_Fag_syl_cpool))

main_effects_Fag_syl_cpool = apply(mean_effects_Fag_syl_cpool[,1:38],2,mean)
main_effects_Fag_syl_cpool_weight = main_effects_Fag_syl_cpool/sum(abs(main_effects_Fag_syl_cpool))

mean_effects_Fag_syl_cflux = effects_Fag_syl[["cflux_rf"]][["complete"]]
mean_effects_Fag_syl_cflux = mean_effects_Fag_syl_cflux[growth_sites_Fag_syl,]

main_effects_Fag_syl_cflux = apply(abs(mean_effects_Fag_syl_cflux[,1:38]),2,mean)
main_effects_Fag_syl_cflux_weight_abs = main_effects_Fag_syl_cflux/sum(abs(main_effects_Fag_syl_cflux))

main_effects_Fag_syl_cflux = apply(mean_effects_Fag_syl_cflux[,1:38],2,mean)
main_effects_Fag_syl_cflux_weight = main_effects_Fag_syl_cflux/sum(abs(main_effects_Fag_syl_cflux))


mean_effects_Fag_syl_agpp = effects_Fag_syl[["agpp_rf"]][["complete"]]
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



mean_effects_Pin_syl_cpool = effects_Pin_syl[["cpool_rf"]][["complete"]]
growth_sites_Pin_syl = which(mean_effects_Pin_syl_cpool[,39]/200 > 2)
mean_effects_Pin_syl_cpool = mean_effects_Pin_syl_cpool[growth_sites_Pin_syl,]

main_effects_Pin_syl_cpool = apply(abs(mean_effects_Pin_syl_cpool[,1:38]),2,mean)
main_effects_Pin_syl_cpool_weight_abs = main_effects_Pin_syl_cpool/sum(abs(main_effects_Pin_syl_cpool))

main_effects_Pin_syl_cpool = apply(mean_effects_Pin_syl_cpool[,1:38],2,mean)
main_effects_Pin_syl_cpool_weight = main_effects_Pin_syl_cpool/sum(abs(main_effects_Pin_syl_cpool))

mean_effects_Pin_syl_cflux = effects_Pin_syl[["cflux_rf"]][["complete"]]
mean_effects_Pin_syl_cflux = mean_effects_Pin_syl_cflux[growth_sites_Pin_syl,]

main_effects_Pin_syl_cflux = apply(abs(mean_effects_Pin_syl_cflux[,1:38]),2,mean)
main_effects_Pin_syl_cflux_weight_abs = main_effects_Pin_syl_cflux/sum(abs(main_effects_Pin_syl_cflux))

main_effects_Pin_syl_cflux = apply(mean_effects_Pin_syl_cflux[,1:38],2,mean)
main_effects_Pin_syl_cflux_weight = main_effects_Pin_syl_cflux/sum(abs(main_effects_Pin_syl_cflux))


mean_effects_Pin_syl_agpp = effects_Pin_syl[["agpp_rf"]][["complete"]]
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


mean_effects_mixed_cpool = effects_mixed[["cpool_rf"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
mean_effects_mixed_cpool = mean_effects_mixed_cpool[growth_sites_mixed,]

remaped_effects_cpool_abs = remap_paramters(mean_effects_mixed_cpool, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool = apply(abs(remaped_effects_cpool_abs),2,mean)

main_effects_mixed_cpool_weight_abs = main_effects_mixed_cpool/sum(abs(main_effects_mixed_cpool))

remaped_effects_cpool = remap_paramters_rel(mean_effects_mixed_cpool, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cpool = apply(remaped_effects_cpool,2,mean)

main_effects_mixed_cpool_weight = main_effects_mixed_cpool/sum(abs(main_effects_mixed_cpool))


mean_effects_mixed_cflux = effects_mixed[["cflux_rf"]][["complete"]]
mean_effects_mixed_cflux = mean_effects_mixed_cflux[growth_sites_mixed,]

remaped_effects_cflux_abs = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux = apply(abs(remaped_effects_cflux_abs),2,mean)

main_effects_mixed_cflux_weight_abs = main_effects_mixed_cflux/sum(abs(main_effects_mixed_cflux))

remaped_effects_cflux = remap_paramters_rel(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_cflux = apply(remaped_effects_cflux,2,mean)

main_effects_mixed_cflux_weight = main_effects_mixed_cflux/sum(abs(main_effects_mixed_cflux))


mean_effects_mixed_agpp = effects_mixed[["agpp_rf"]][["complete"]]
mean_effects_mixed_agpp = mean_effects_mixed_agpp[growth_sites_mixed,]

remaped_effects_agpp_abs = remap_paramters(mean_effects_mixed_agpp, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp = apply(abs(remaped_effects_agpp_abs),2,mean)

main_effects_mixed_agpp_weight_abs = main_effects_mixed_agpp/sum(abs(main_effects_mixed_agpp))

remaped_effects_agpp = remap_paramters_rel(mean_effects_mixed_agpp, mixed_results$weight_mapping, mixed_results$position_mapping)

main_effects_mixed_agpp = apply(remaped_effects_agpp,2,mean)

main_effects_mixed_agpp_weight = main_effects_mixed_agpp/sum(abs(main_effects_mixed_agpp))

## sd estimates ##

sd_mixed_cpool_abs = apply((remaped_effects_cpool_abs[,1:38]/mean_effects_mixed_cpool[,74]),2,sd, na.rm =T)
sd_mixed_cflux_abs = apply((remaped_effects_cflux_abs[,1:38]/mean_effects_mixed_cflux[,74]),2,sd, na.rm = T)
sd_mixed_agpp_abs = apply((remaped_effects_agpp_abs[,1:38]/mean_effects_mixed_agpp[,74]),2,sd, na.rm = T)

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

## reordering standard deviations

sd_Pin_syl_cpool_abs = sd_Pin_syl_cpool_abs[ordering_pin_syl]
sd_Pin_syl_cflux_abs = sd_Pin_syl_cflux_abs[ordering_pin_syl]
sd_Pin_syl_agpp_abs = sd_Pin_syl_agpp_abs[ordering_pin_syl]

sd_Fag_syl_cpool_abs = sd_Fag_syl_cpool_abs[ordering_fag_syl]
sd_Fag_syl_cflux_abs = sd_Fag_syl_cflux_abs[ordering_fag_syl]
sd_Fag_syl_agpp_abs = sd_Fag_syl_agpp_abs[ordering_fag_syl]

sd_mixed_cpool_abs = sd_mixed_cpool_abs[ordering_mixed]
sd_mixed_cflux_abs = sd_mixed_cflux_abs[ordering_mixed]
sd_mixed_agpp_abs = sd_mixed_agpp_abs[ordering_mixed]



parameters = readRDS("ParameterMetaData/Parameter_list.rds")
drivers = c("co2","ndep","insol","temp","ph","prec")

variablenames = c(as.character(parameters$Pic_abi[,"NameRLPJ"]),drivers)
variablenames2 = substring(variablenames,regexpr("_", variablenames) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  variablenames2 = gsub(noisy_things[i],"",variablenames2)
}

grouping = read.csv("ParameterMetaData/Grouping_Fag_syl.csv", header = T, sep = ";")
grouping = c(as.character(grouping[,1]), rep("Drivers",6))

order_grouping = match(parameternames_ordered_pic_abi2,variablenames2)
grouping2 = grouping[order_grouping]




mean_effects_abs = data.frame("Group" = grouping2,
                              "cpool_fag_syl" = main_effects_Fag_syl_cpool_weight_abs,
                              "cflux_fag_syl" = main_effects_Fag_syl_cflux_weight_abs,
                              "agpp_fag_syl" = main_effects_Fag_syl_agpp_weight_abs,
                              "cpool_pin_syl" = main_effects_Pin_syl_cpool_weight_abs,
                              "cflux_pin_syl" = main_effects_Pin_syl_cflux_weight_abs,
                              "agpp_pin_syl" = main_effects_Pin_syl_agpp_weight_abs,
                              "cpool_pic_abi" = main_effects_Pic_abi_cpool_weight_abs,
                              "cflux_pic_abi" = main_effects_Pic_abi_cflux_weight_abs,
                              "agpp_pic_abi" = main_effects_Pic_abi_agpp_weight_abs,
                              "cpool_mixed" = main_effects_mixed_cpool_weight_abs,
                              "cflux_mixed" = main_effects_mixed_cflux_weight_abs,
                              "agpp_mixed" = main_effects_mixed_agpp_weight_abs,
                              "Names" = parameternames_ordered_pic_abi2)

summary(mean_effects_abs)
mean_effects_abs = mean_effects_abs[order(mean_effects_abs$Group),]

mean_effects_drivers = mean_effects_abs[mean_effects_abs$Group == "Drivers",]

mean_effects_parameters = mean_effects_abs[mean_effects_abs$Group != "Drivers",]
mean_effects_parameters$Group = as.numeric(mean_effects_parameters$Group)
mean_effects_parameters$Group = as.factor(mean_effects_parameters$Group)

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
colors_parameter = c("burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")
spaces = table(mean_effects_parameters$Group)
colors_full = vector()
for(i in 1:6){
  colors_full = c(colors_full,rep(colors_parameter[i],spaces[i]))
}

pdf("./Figures/Mean_effects_abs_rf.pdf", width = 14.0, height = 9.)

#par(mfrow = c(3,1),mai = c(0.1, 0.1, 0.1, 0.1))
layout(matrix(c(1,2,3,4),byrow =T,nrow =4), heights = c(2,2,2,1))
old_mar = c(5.1,4.1,4.1,2.1)
par(mar = c(0.6,4.1,9.1,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n', cex.axis =2,
     ylim = c(0,0.39), cex.main = 1.5)
mtext("Uncertainty contributions [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_effects_parameters$cpool_fag_syl,
                             mean_effects_parameters$cpool_pic_abi,
                             mean_effects_parameters$cpool_pin_syl,
                             mean_effects_parameters$cpool_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_effects_parameters$cpool_fag_syl , x= 1:32, pch = 15,
       ylim = c(0,0.39), col = "darkgreen",
       cex = 1.2)
points(y = mean_effects_parameters$cpool_pic_abi,
       x = 1:nrow(mean_effects_parameters), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_effects_parameters$cpool_pin_syl,
       x = 1:nrow(mean_effects_parameters), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_effects_parameters$cpool_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_effects_parameters))
points(y = apply(rbind(mean_effects_parameters$cpool_fag_syl,
                       mean_effects_parameters$cpool_pic_abi,
                       mean_effects_parameters$cpool_pin_syl),2,mean),
       x = 1:nrow(mean_effects_parameters), col = "black", pch = "-", cex=3)
axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0,0.40,length.out = 11),2), las = 2, line = -0.5)
barplot(height = apply(rbind(mean_effects_drivers$cpool_fag_syl,
                             mean_effects_drivers$cpool_pic_abi,
                             mean_effects_drivers$cpool_pin_syl,
                             mean_effects_drivers$cpool_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= rep("gold2",6), percent = 65)),
        add = T,axes = F, width = 0.85, space = c(33.5/0.85,rep(0.15/0.85,5)))
points(y =mean_effects_drivers$cpool_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,.38),cex = 1.2)
ablineclip(v = cumsum(table(mean_effects_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
points(y = mean_effects_drivers$cpool_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_effects_drivers$cpool_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_effects_drivers$cpool_mixed, col = "purple",pch = 8,
       cex = 1.2,x = 34:39)
points(y = apply(rbind(mean_effects_drivers$cpool_fag_syl,
                       mean_effects_drivers$cpool_pic_abi,
                       mean_effects_drivers$cpool_pin_syl),2,mean),
       x = 34:39, col = "black", pch = "-", cex=3)
ablineclip(v = 33,col = alpha("black",0.7),y1 = -0.19 )
spacings = (cumsum(table(mean_effects_parameters$Group)) + c(0,cumsum(table(mean_effects_parameters$Group))[-length(cumsum(table(mean_effects_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(mean_effects_parameters$Group))
par(xpd =  T)
text(c(names(table(mean_effects_abs$Group))[which(names(table(mean_effects_abs$Group)) != "Drivers")],"Drivers"),
     x = c(spacings,36),  y = 0.43, font = 2, cex =1.7, srt = 30,
     adj = c(0,0))
par(xpd = F)
#ablineclip(h = 0.19, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=-0.19)
title(main = "a)", line = 6, adj = 0.01, cex.main = 2)


par(mar = c(0.6,4.1,5.6,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n', cex.axis =2,
     ylim = c(0,0.4), cex.main = 1.5)
mtext("Uncertainty contributions [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_effects_parameters$cflux_fag_syl,
                             mean_effects_parameters$cflux_pic_abi,
                             mean_effects_parameters$cflux_pin_syl,
                             mean_effects_parameters$cflux_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_effects_parameters$cflux_fag_syl , x= 1:32, pch = 15,
       ylim = c(0,0.4), col = "darkgreen",
       cex = 1.2)
points(y = mean_effects_parameters$cflux_pic_abi,
       x = 1:nrow(mean_effects_parameters), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_effects_parameters$cflux_pin_syl,
       x = 1:nrow(mean_effects_parameters), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_effects_parameters$cflux_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_effects_parameters))
points(y = apply(rbind(mean_effects_parameters$cflux_fag_syl,
                       mean_effects_parameters$cflux_pic_abi,
                       mean_effects_parameters$cflux_pin_syl)
                 ,2,mean),
       x = 1:nrow(mean_effects_parameters), col = "black", pch = "-", cex=3)
axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0,0.4,length.out = 11),2), las = 2, line = -0.5)
barplot(height = apply(rbind(mean_effects_drivers$cflux_fag_syl,
                             mean_effects_drivers$cflux_pic_abi,
                             mean_effects_drivers$cflux_pin_syl,
                             mean_effects_drivers$cflux_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= rep("gold2",6), percent = 65)),
        add = T,axes = F, width = 0.85, space = c(33.5/0.85,rep(0.15/0.85,5)))
points(y =mean_effects_drivers$cflux_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,.4),cex = 1.2)
ablineclip(v = cumsum(table(mean_effects_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = 0, y2 =0.4)
points(y = mean_effects_drivers$cflux_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_effects_drivers$cflux_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_effects_drivers$cflux_mixed, col = "purple",pch = 8,
       cex = 1.2,x = 34:39)
points(y = apply(rbind(mean_effects_drivers$cflux_fag_syl,
                       mean_effects_drivers$cflux_pic_abi,
                       mean_effects_drivers$cflux_pin_syl)
                 ,2,mean),
       x = 34:39, col = "black", pch = "-", cex=3)
ablineclip(v = 33,col = alpha("black",0.7),y1 = -0.16, y2 =0.4 )
spacings = (cumsum(table(mean_effects_parameters$Group)) + c(0,cumsum(table(mean_effects_parameters$Group))[-length(cumsum(table(mean_effects_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(mean_effects_parameters$Group))
par(xpd=T)
legend(x = 40, y = 0.45,
       legend = c("Fag. syl.","Pic. abi.",'Pin. syl.',"Mixed",'Mean Mono'),
       col = c('darkgreen', 'darkblue','brown', "purple",'black'), bty = 'n',
       title = as.expression(bquote(italic(bold("Species")))),
       pch = c(15,19,17,8,NA), cex = 1.8, lty = c(NA,NA,NA,NA,1),
       lwd = 2)
par(xpd = F)
#ablineclip(h = 0.12, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=-0.16, y2 =0.4)
title(main = "b)", line = 2, adj = 0.01, cex.main = 2)

par(mar = c(1.1,4.1,6.1,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n', cex.axis =2,
     ylim = c(0,0.4), cex.main = 1.5)
mtext("Uncertainty contributions [%]", side =2, cex =1., line = 2.5)
barplot(height = apply(rbind(mean_effects_parameters$agpp_fag_syl,
                             mean_effects_parameters$agpp_pic_abi,
                             mean_effects_parameters$agpp_pin_syl,
                             mean_effects_parameters$agpp_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= colors_full, percent = 65)), add = T,axes = F,
        width = 0.85, space = c(0.5/0.85,rep(0.15/0.85,31)), pos = rep(1,32))
points(y = mean_effects_parameters$agpp_fag_syl , x= 1:32, pch = 15,
       ylim = c(0,0.4), col = "darkgreen",
       cex = 1.2)
points(y = mean_effects_parameters$agpp_pic_abi,
       x = 1:nrow(mean_effects_parameters), col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_effects_parameters$agpp_pin_syl,
       x = 1:nrow(mean_effects_parameters), col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_effects_parameters$agpp_mixed, col = "purple",pch =8,
       cex = 1.2,x = 1:nrow(mean_effects_parameters))
points(y = apply(rbind(mean_effects_parameters$agpp_fag_syl,
                       mean_effects_parameters$agpp_pic_abi,
                       mean_effects_parameters$agpp_pin_syl),2,mean),
       x = 1:nrow(mean_effects_parameters), col = "black", pch = "-", cex=3)
axis(1, at = 1:39, labels = rep("",39),
     srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0,0.4,length.out = 11),2), las = 2, line = -0.5)
barplot(height = apply(rbind(mean_effects_drivers$agpp_fag_syl,
                             mean_effects_drivers$agpp_pic_abi,
                             mean_effects_drivers$agpp_pin_syl,
                             mean_effects_drivers$agpp_mixed),2,mean),
        col = as.vector(sapply(FUN = t_col, X= rep("gold2",6), percent = 65)),
        add = T,axes = F, width = 0.85, space = c(33.5/0.85,rep(0.15/0.85,5)))
points(y =mean_effects_drivers$agpp_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,.4),cex = 1.2)
ablineclip(v = cumsum(table(mean_effects_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = -0.25, y2= 0.4)
points(y = mean_effects_drivers$agpp_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 1.2)
points(y = mean_effects_drivers$agpp_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 1.2)
points(y = mean_effects_drivers$agpp_mixed, col = "purple",pch = 8,
       cex = 1.2,x = 34:39)
points(y = apply(rbind(mean_effects_drivers$agpp_fag_syl,
                       mean_effects_drivers$agpp_pic_abi,
                       mean_effects_drivers$agpp_pin_syl)
                 ,2,mean),
       x = 34:39, col = "black", pch = "-", cex=3)
ablineclip(v = 33,col = alpha("black",0.7),y1 =0, y2= 0.4 )
spacings = (cumsum(table(mean_effects_parameters$Group)) + c(0,cumsum(table(mean_effects_parameters$Group))[-length(cumsum(table(mean_effects_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(mean_effects_parameters$Group))
#ablineclip(h = 0.14, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=-0.25, y2= 0.4)
title(main = "c)", line = 1.7, adj = 0.01, cex.main = 2)

par(mar = c(1.1,4.1,0.1,12.1))
plot(y = rep(0,39),x=1:39,
     ylab = "", xlab = "", xaxt = "n",
     col = c(rep("white",39)),xaxt='n', bty="n", yaxt = 'n', cex.axis =2,
     ylim = c(-0.01,0.01), cex.main = 1.5)
par(xpd = T)
text(x = 1:39, y = 0.01, c(as.character(mean_effects_parameters$Names),"",
                           as.character(mean_effects_drivers$Names)),srt = 270,
     adj = c(0.0,0.), cex =1.48)
par(xpd = F)

dev.off()



sd_effects = data.frame("Group" = grouping2,
                        "cpool_fag_syl" = sd_Fag_syl_cpool_abs,
                        "cflux_fag_syl" = sd_Fag_syl_cflux_abs,
                        "agpp_fag_syl" = sd_Fag_syl_agpp_abs,
                        "cpool_pin_syl" = sd_Pin_syl_cpool_abs,
                        "cflux_pin_syl" = sd_Pin_syl_cflux_abs,
                        "agpp_pin_syl" = sd_Pin_syl_agpp_abs,
                        "cpool_pic_abi" = sd_Pic_abi_cpool_abs,
                        "cflux_pic_abi" = sd_Pic_abi_cflux_abs,
                        "agpp_pic_abi" = sd_Pic_abi_agpp_abs,
                        "cpool_mixed" = sd_mixed_cpool_abs,
                        "cflux_mixed" = sd_mixed_cflux_abs,
                        "agpp_mixed" = sd_mixed_agpp_abs,
                        "Names" = parameternames_ordered_pic_abi2)



summary(sd_effects)
sd_effects = sd_effects[order(sd_effects$Group),]

sd_effects_drivers = sd_effects[sd_effects$Group == "Drivers",]

sd_effects_parameters = sd_effects[sd_effects$Group != "Drivers",]
sd_effects_parameters$Group = as.numeric(sd_effects_parameters$Group)
sd_effects_parameters$Group = as.factor(sd_effects_parameters$Group)


pdf("./Figures/sd_effects_lin_rf.pdf", width = 18.0, height = 12)

par(mfrow = c(3,1),mai = c(0.1, 0.1, 0.1, 0.1))
old_mar = c(5.1,4.1,4.1,2.1)
par(mar = c(1.6,4.1,9.1,12.1))
plot(y = c(sd_effects_parameters$cpool_fag_syl,rep(0,7)) , x= 1:39,
     ylab = "Uncertainty contributions [%]", xlab = "", xaxt = "n",
     col = c(rep("darkgreen",32),rep("white",7)),xaxt='n', bty="n", yaxt = 'n', pch = 15,
     ylim = c(0,1.2),
     cex = 0.8)
points(y = sd_effects_parameters$cpool_pic_abi,
       x = 1:nrow(sd_effects_parameters), col = "darkblue", pch = 19,cex = 0.8)
points(y = sd_effects_parameters$cpool_pin_syl,
       x = 1:nrow(sd_effects_parameters), col = "brown", pch = 17,
       cex = 0.8)
points(y = sd_effects_parameters$cpool_mixed, col = "purple",pch =8,
       cex = 0.8,x = 1:nrow(sd_effects_parameters))
points(y = apply(rbind(sd_effects_parameters$cpool_fag_syl,sd_effects_parameters$cpool_pic_abi,sd_effects_parameters$cpool_pin_syl),2,mean),
       x = 1:nrow(sd_effects_parameters), col = "black", pch = "-", cex =2)
axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0,1.2,length.out = 7),2), las = 2, line = -0.5)
points(y =sd_effects_drivers$cpool_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,1.2),cex = 0.8)
ablineclip(v = cumsum(table(sd_effects_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
points(y = sd_effects_drivers$cpool_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 0.8)
points(y = sd_effects_drivers$cpool_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 0.8)
points(y = sd_effects_drivers$cpool_mixed, col = "purple",pch = 8,
       cex = 0.8,x = 34:39)
points(y = apply(rbind(sd_effects_drivers$cpool_fag_syl,sd_effects_drivers$cpool_pic_abi,sd_effects_drivers$cpool_pin_syl),2,mean),
       x = 34:39, col = "black", pch = "-", cex =2)
ablineclip(v = 33,col = alpha("black",0.7),y1 = 0 )
spacings = (cumsum(table(sd_effects_parameters$Group)) + c(0,cumsum(table(sd_effects_parameters$Group))[-length(cumsum(table(sd_effects_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(sd_effects_parameters$Group))
par(xpd =  T)
text(c(names(table(sd_effects$Group))[which(names(table(sd_effects$Group)) != "Drivers")],"Drivers"),
     x = c(spacings,36) +c(0,0,1.5,2,0,0.5,0),  y = 1.2, font = 2, cex =1, srt = 45, adj = c(0,0))
par(xpd = F)
ablineclip(h = 1.2, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=0)
title(main = "a)  Total standing biomass", line = 6, adj = 0.01, cex.main = 1.5)


par(mar = c(4.6,4.1,6.1,12.1))
plot(y = c(sd_effects_parameters$cflux_fag_syl,rep(0,7)) , x= 1:39,
     ylab = "Uncertainty contributions [%]", xlab = "", xaxt = "n",
     col = c(rep("darkgreen",32),rep("white",7)),xaxt='n', bty="n", yaxt = 'n', pch = 15,
     ylim = c(0,62),
     cex = 0.8)
points(y = sd_effects_parameters$cflux_pic_abi,
       x = 1:nrow(sd_effects_parameters), col = "darkblue", pch = 19,cex = 0.8)
points(y = sd_effects_parameters$cflux_pin_syl,
       x = 1:nrow(sd_effects_parameters), col = "brown", pch = 17,
       cex = 0.8)
points(y = sd_effects_parameters$cflux_mixed, col = "purple",pch =8,
       cex = 0.8,x = 1:nrow(sd_effects_parameters))
points(y = apply(rbind(sd_effects_parameters$cflux_fag_syl,sd_effects_parameters$cflux_pic_abi,sd_effects_parameters$cflux_pin_syl)
                 ,2,mean),
       x = 1:nrow(sd_effects_parameters), col = "black", pch = "-", cex =2)
axis(1, at = 1:39, labels = rep("",39),srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0,60,length.out = 7),2), las = 2, line = -0.5)
points(y =sd_effects_drivers$cflux_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,60),cex = 0.8)
ablineclip(v = cumsum(table(sd_effects_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
points(y = sd_effects_drivers$cflux_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 0.8)
points(y = sd_effects_drivers$cflux_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 0.8)
points(y = sd_effects_drivers$cflux_mixed, col = "purple",pch = 8,
       cex = 0.8,x = 34:39)
points(y = apply(rbind(sd_effects_drivers$cflux_fag_syl,sd_effects_drivers$cflux_pic_abi,sd_effects_drivers$cflux_pin_syl)
                 ,2,mean),
       x = 34:39, col = "black", pch = "-", cex =2)
ablineclip(v = 33,col = alpha("black",0.7),y1 = 0 )
spacings = (cumsum(table(sd_effects_parameters$Group)) + c(0,cumsum(table(sd_effects_parameters$Group))[-length(cumsum(table(sd_effects_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(sd_effects_parameters$Group))
par(xpd=T)
legend(x = 40, y = 50,
       legend = c("Fagus sylvatica","Picea abies",'Pinus sylvestris',"Mixed",'sd Mono'),
       col = c('darkgreen', 'darkblue','brown', "purple",'black'), bty = 'n',
       title = as.expression(bquote(italic(bold("Species")))),
       pch = c(15,19,17,8,NA), cex = 1.3, lty = c(NA,NA,NA,NA,1),
       lwd = 2)
par(xpd = F)
ablineclip(h =62, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=0)
title(main = "b)  Carbon balance", line = 2, adj = 0.01, cex.main = 1.5)

par(mar = c(8.1,4.1,2.6,12.1))
plot(y = c(sd_effects_parameters$agpp_fag_syl,rep(0,7)) , x= 1:39,
     ylab = "Uncertainty contributions [%]", xlab = "", xaxt = "n",
     col = c(rep("darkgreen",32),rep("white",7)),xaxt='n', bty="n", yaxt = 'n', pch = 15,
     ylim = c(0,8),
     cex = 0.8)
points(y = sd_effects_parameters$agpp_pic_abi,
       x = 1:nrow(sd_effects_parameters), col = "darkblue", pch = 19,cex = 0.8)
points(y = sd_effects_parameters$agpp_pin_syl,
       x = 1:nrow(sd_effects_parameters), col = "brown", pch = 17,
       cex = 0.8)
points(y = sd_effects_parameters$agpp_mixed, col = "purple",pch =8,
       cex = 0.8,x = 1:nrow(sd_effects_parameters))
points(y = apply(rbind(sd_effects_parameters$agpp_fag_syl,sd_effects_parameters$agpp_pic_abi,sd_effects_parameters$agpp_pin_syl),2,mean),
       x = 1:nrow(sd_effects_parameters), col = "black", pch = "-", cex =2)
axis(1, at = 1:39, labels = c(as.character(sd_effects_parameters$Names),"",
                              as.character(sd_effects_drivers$Names)),
     srt = 45, las =2 , pos = 0)
axis(2, at = round(seq(0,8,length.out =5),2), las = 2, line = -0.5)
points(y =sd_effects_drivers$agpp_fag_syl, x = 34:39,col = "darkgreen",xaxt='n', bty="n", yaxt = 'n', pch = 15,
       ylim = c(0,8),cex = 0.8)
ablineclip(v = cumsum(table(sd_effects_parameters$Group))+0.5, col = alpha("black",0.7),
           y1 = 0)
points(y = sd_effects_drivers$agpp_pic_abi,
       x = 34:39, col = "darkblue", pch = 19,cex = 0.8)
points(y = sd_effects_drivers$agpp_pin_syl,
       x = 34:39, col = "brown", pch = 17,
       cex = 0.8)
points(y = sd_effects_drivers$agpp_mixed, col = "purple",pch = 8,
       cex = 0.8,x = 34:39)
points(y = apply(rbind(sd_effects_drivers$agpp_fag_syl,sd_effects_drivers$agpp_pic_abi,sd_effects_drivers$agpp_pin_syl)
                 ,2,mean),
       x = 34:39, col = "black", pch = "-", cex =2)
ablineclip(v = 33,col = alpha("black",0.7),y1 = 0 )
spacings = (cumsum(table(sd_effects_parameters$Group)) + c(0,cumsum(table(sd_effects_parameters$Group))[-length(cumsum(table(sd_effects_parameters$Group)))]))/2 +0.5
grouping_variables = names(table(sd_effects_parameters$Group))
#mtext(c(letters[1:6],"Drivers"),at = c(spacings,36), side = 3, font = 2, cex =1)
ablineclip(h = 8, x1 =0, x2 = 40.5)
ablineclip(v = c(39.5), y1=0)
title(main = "c)  Gross primary production", line = 1.7, adj = 0.01, cex.main = 1.5)


dev.off()
