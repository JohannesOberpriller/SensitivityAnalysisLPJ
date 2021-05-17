library(ranger)
# function to rescale the parameters from -0.5 to 0.5
rescaling <- function(x) {(x-min(x))/(max(x)-min(x))-0.5}

# functions to extract the interactions matrix from the results of a linear regression


get_interaction_matrix <-function(coefficients, num_parameters){
  interaction_matrix <- matrix(ncol = num_parameters, nrow = num_parameters)
  interaction_matrix[lower.tri(interaction_matrix)] = coefficients[(num_parameters+1):length(coefficients)]
  for(j in 2:ncol(interaction_matrix)){
    for(k in 1:(j-1)){
      interaction_matrix[k,j] = interaction_matrix[j,k]
    }
  }
  return(interaction_matrix)
}

# function to get the main effects

get_main_effects_per_site <- function(results, parameternames){

  # sites: a list of sites
  # parameternames: the names of the parameters


  effects_cpool_complete = matrix(nrow = 1, ncol = 74)
  effects_cpool_climate_change = matrix(nrow = 1, ncol = 74)
  effects_cpool_steady_climate = matrix(nrow = 1, ncol = 74)

  effects_cflux_complete = matrix(nrow = 1, ncol = 74)
  effects_cflux_climate_change = matrix(nrow = 1, ncol = 74)
  effects_cflux_steady_climate = matrix(nrow = 1, ncol = 74)

  effects_agpp_complete = matrix(nrow = 1, ncol = 74)
  effects_agpp_climate_change = matrix(nrow = 1, ncol = 74)
  effects_agpp_steady_climate = matrix(nrow = 1, ncol = 74)

  effects_cpool_complete_rf = matrix(nrow = 1, ncol = 74)
  effects_cpool_climate_change_rf = matrix(nrow = 1, ncol = 74)
  effects_cpool_steady_climate_rf = matrix(nrow = 1, ncol = 74)

  effects_cflux_complete_rf = matrix(nrow = 1, ncol = 74)
  effects_cflux_climate_change_rf = matrix(nrow = 1, ncol = 74)
  effects_cflux_steady_climate_rf = matrix(nrow = 1, ncol = 74)

  effects_agpp_complete_rf = matrix(nrow = 1, ncol = 74)
  effects_agpp_climate_change_rf = matrix(nrow = 1, ncol = 74)
  effects_agpp_steady_climate_rf = matrix(nrow = 1, ncol = 74)


  interactions_cpool_complete = list()
  interactions_cpool_climate_change = list()
  interactions_cpool_steady_climate = list()

  interactions_cflux_complete = list()
  interactions_cflux_climate_change = list()
  interactions_cflux_steady_climate = list()


  interactions_agpp_complete = list()
  interactions_agpp_climate_change = list()
  interactions_agpp_steady_climate = list()

  parameter_matrix = matrix(ncol = length(unique(parameternames)), nrow = length(results$parameters))
  for(j in 1:length(results$parameters)){
    parameter_matrix[j,] = results$parameters[[j]]
  }
  # apply random forest and linear regression to get the main effects of the parameters and env. drivers for all outputs

    parameters = apply(parameter_matrix, 2, rescaling)

      ### cpool stuff ###

      ### complete ###

      overall_cpool_complete = results[['cpool']][['complete']][1,]

      data_cpool_complete <- cbind(overall_cpool_complete,parameters)
      fit <- lm(overall_cpool_complete ~ (.)^2 , data = as.data.frame(data_cpool_complete))
      coefficients = coef(fit)[-1]
      effects_cpool_complete[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cpool_complete[1,ncol(effects_cpool_complete)] = results[['cpool']][['complete']][2,1]

      interaction_matrix_cpool_complete = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_cpool_complete[[1]] = interaction_matrix_cpool_complete

      fit_rf <- ranger(overall_cpool_complete ~ . ,
                       data = as.data.frame(data_cpool_complete)[complete.cases(overall_cpool_complete),],
                       importance = 'impurity', num.threads = 2)
      effects_cpool_complete_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cpool_complete_rf[1,ncol(effects_cpool_complete)] = results[['cpool']][['complete']][2,1]


      ### climate change ###

      overall_cpool_climate_change = results[['cpool']][['climate_change']][1,]

      data_cpool_climate_change <- cbind(overall_cpool_climate_change,parameters)
      fit <- lm(overall_cpool_climate_change ~ (.)^2 , data = as.data.frame(data_cpool_climate_change))
      coefficients = coef(fit)[-1]
      effects_cpool_climate_change[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cpool_climate_change[1,ncol(effects_cpool_complete)] = results[['cpool']][['climate_change']][2,1]

      interaction_matrix_cpool_climate_change = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_cpool_climate_change[[1]] = interaction_matrix_cpool_climate_change

      fit_rf <- ranger(overall_cpool_climate_change ~ . ,
                       data = as.data.frame(data_cpool_climate_change)[complete.cases(overall_cpool_climate_change),],
                       importance = 'impurity', num.threads = 2)
      effects_cpool_climate_change_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cpool_climate_change_rf[1,ncol(effects_cpool_climate_change)] = results[['cpool']][['climate_change']][2,1]


      ### steady climate ###

      overall_cpool_steady_climate= results[['cpool']][['steady_climate']][1,]

      data_cpool_steady_climate <- cbind(overall_cpool_steady_climate,parameters)
      fit <- lm(overall_cpool_steady_climate ~ (.)^2 , data = as.data.frame(data_cpool_steady_climate))
      coefficients = coef(fit)[-1]
      effects_cpool_steady_climate[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cpool_steady_climate[1,ncol(effects_cpool_complete)] = results[['cpool']][['steady_climate']][2,1]

      interaction_matrix_cpool_steady_climate = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_cpool_steady_climate[[1]] = interaction_matrix_cpool_steady_climate

      fit_rf <- ranger(overall_cpool_steady_climate ~ . ,
                       data = as.data.frame(data_cpool_steady_climate)[complete.cases(overall_cpool_steady_climate),],
                       importance = 'impurity', num.threads = 2)
      effects_cpool_steady_climate_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cpool_steady_climate_rf[1,ncol(effects_cpool_steady_climate)] = results[['cpool']][['steady_climate']][2,1]


      ## cflux related stuff

      ### complete ###

      overall_cflux_complete = results[['cflux']][['complete']][1,]

      data_cflux_complete <- cbind(overall_cflux_complete,parameters)
      fit <- lm(overall_cflux_complete ~ (.)^2 , data = as.data.frame(data_cflux_complete))
      coefficients = coef(fit)[-1]
      effects_cflux_complete[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cflux_complete[1,ncol(effects_cflux_complete)] = results[['cflux']][['complete']][2,1]

      interaction_matrix_cflux_complete = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_cflux_complete[[1]] = interaction_matrix_cflux_complete

      fit_rf <- ranger(overall_cflux_complete ~ . ,
                       data = as.data.frame(data_cflux_complete)[complete.cases(overall_cflux_complete),],
                       importance = 'impurity', num.threads = 2)
      effects_cflux_complete_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cflux_complete_rf[1,ncol(effects_cflux_complete)] = results[['cflux']][['complete']][2,1]


      ### climate change ###

      overall_cflux_climate_change = results[['cflux']][['climate_change']][1,]

      data_cflux_climate_change <- cbind(overall_cflux_climate_change,parameters)
      fit <- lm(overall_cflux_climate_change ~ (.)^2 , data = as.data.frame(data_cflux_climate_change))
      coefficients = coef(fit)[-1]
      effects_cflux_climate_change[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cflux_climate_change[1,ncol(effects_cflux_complete)] = results[['cflux']][['climate_change']][2,1]

      interaction_matrix_cflux_climate_change = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_cflux_climate_change[[1]] = interaction_matrix_cflux_climate_change

      fit_rf <- ranger(overall_cflux_climate_change ~ . ,
                       data = as.data.frame(data_cflux_climate_change)[complete.cases(overall_cflux_climate_change),],
                       importance = 'impurity', num.threads = 2)
      effects_cflux_climate_change_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cflux_climate_change_rf[1,ncol(effects_cflux_climate_change)] = results[['cflux']][['climate_change']][2,1]

      ### steady climate ###

      overall_cflux_steady_climate= results[['cflux']][['steady_climate']][1,]

      data_cflux_steady_climate <- cbind(overall_cflux_steady_climate,parameters)
      fit <- lm(overall_cflux_steady_climate ~ (.)^2 , data = as.data.frame(data_cflux_steady_climate))
      coefficients = coef(fit)[-1]
      effects_cflux_steady_climate[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cflux_steady_climate[1,ncol(effects_cflux_complete)] = results[['cflux']][['steady_climate']][2,1]

      interaction_matrix_cflux_steady_climate = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_cflux_steady_climate[[1]] = interaction_matrix_cflux_steady_climate


      fit_rf <- ranger(overall_cflux_steady_climate ~ . ,
                       data = as.data.frame(data_cflux_steady_climate)[complete.cases(overall_cflux_steady_climate),],
                       importance = 'impurity', num.threads = 2)
      effects_cflux_steady_climate_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cflux_steady_climate_rf[1,ncol(effects_cflux_steady_climate)] = results[['cflux']][['steady_climate']][2,1]

       ## agpp related stuff

      ### complete ###

      overall_agpp_complete = results[['agpp']][['complete']][1,]

      data_agpp_complete <- cbind(overall_agpp_complete,parameters)
      fit <- lm(overall_agpp_complete ~ (.)^2 , data = as.data.frame(data_agpp_complete))
      coefficients = coef(fit)[-1]
      effects_agpp_complete[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_agpp_complete[1,ncol(effects_agpp_complete)] = results[['agpp']][['complete']][2,1]

      interaction_matrix_agpp_complete = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_agpp_complete[[1]] = interaction_matrix_agpp_complete

      fit_rf <- ranger(overall_agpp_complete ~ . ,
                       data = as.data.frame(data_agpp_complete)[complete.cases(overall_agpp_complete),],
                       importance = 'impurity', num.threads = 2)
      effects_agpp_complete_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_agpp_complete_rf[1,ncol(effects_agpp_complete)] = results[['agpp']][['complete']][2,1]

      ### climate change ###

      overall_agpp_climate_change = results[['agpp']][['climate_change']][1,]

      data_agpp_climate_change <- cbind(overall_agpp_climate_change,parameters)
      fit <- lm(overall_agpp_climate_change ~ (.)^2 , data = as.data.frame(data_agpp_climate_change))
      coefficients = coef(fit)[-1]
      effects_agpp_climate_change[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_agpp_climate_change[1,ncol(effects_agpp_complete)] = results[['agpp']][['climate_change']][2,1]

      interaction_matrix_agpp_climate_change = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_agpp_climate_change[[1]] = interaction_matrix_agpp_climate_change

      fit_rf <- ranger(overall_agpp_climate_change ~ . ,
                       data = as.data.frame(data_agpp_climate_change)[complete.cases(overall_agpp_climate_change),],
                       importance = 'impurity', num.threads = 2)
      effects_agpp_climate_change_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_agpp_climate_change_rf[1,ncol(effects_agpp_climate_change)] = results[['agpp']][['climate_change']][2,1]
      ### steady climate ###

      overall_agpp_steady_climate= results[['agpp']][['steady_climate']][1,]

      data_agpp_steady_climate <- cbind(overall_agpp_steady_climate,parameters)
      fit <- lm(overall_agpp_steady_climate ~ (.)^2 , data = as.data.frame(data_agpp_steady_climate))
      coefficients = coef(fit)[-1]
      effects_agpp_steady_climate[1,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_agpp_steady_climate[1,ncol(effects_agpp_complete)] = results[['agpp']][['steady_climate']][2,1]

      interaction_matrix_agpp_steady_climate = get_interaction_matrix(coefficients, num_parameters = ncol(parameters))
      interactions_agpp_steady_climate[[1]] = interaction_matrix_agpp_steady_climate

      fit_rf <- ranger(overall_agpp_steady_climate ~ . ,
                       data = as.data.frame(data_agpp_steady_climate)[complete.cases(overall_agpp_steady_climate),],
                       importance = 'impurity', num.threads = 2)
      effects_agpp_steady_climate_rf[1,1:ncol(parameters)] = fit_rf$variable.importance
      effects_agpp_steady_climate_rf[1,ncol(effects_agpp_steady_climate)] = results[['agpp']][['steady_climate']][2,1]


  result_list = list("agpp" = list("steady_climate" = effects_agpp_steady_climate,
                                   "climate_change" = effects_agpp_climate_change,
                                   "complete" = effects_agpp_complete),
                     "agpp_interactions" = list("steady_climate" = interactions_agpp_steady_climate,
                                                "climate_change" = interactions_agpp_climate_change,
                                                "complete" = interactions_agpp_complete),
                     "agpp_rf" = list("steady_climate" = effects_agpp_steady_climate_rf,
                                      "climate_change" = effects_agpp_climate_change_rf,
                                      "complete" = effects_agpp_complete_rf),
                     "cpool" = list("steady_climate" = effects_cpool_steady_climate,
                                    "climate_change" = effects_cpool_climate_change,
                                    "complete" = effects_cpool_complete),
                     "cpool_interactions" = list("steady_climate" = interactions_cpool_steady_climate,
                                                 "climate_change" = interactions_cpool_climate_change,
                                                 "complete" = interactions_cpool_complete),
                     "cpool_rf" = list("steady_climate" = effects_cpool_steady_climate_rf,
                                      "climate_change" = effects_cpool_climate_change_rf,
                                      "complete" = effects_cpool_complete_rf),
                     'cflux'= list("steady_climate" = effects_cflux_steady_climate,
                                   "climate_change" = effects_cflux_climate_change,
                                   "complete" = effects_cflux_complete),
                     "cflux_interactions" = list("steady_climate" = interactions_cflux_steady_climate,
                                                 "climate_change" = interactions_cflux_climate_change,
                                                 "complete" = interactions_cflux_complete),
                     "cflux_rf" = list("steady_climate" = effects_cflux_steady_climate_rf,
                                      "climate_change" = effects_cflux_climate_change_rf,
                                      "complete" = effects_cflux_complete_rf)
                     )
  return(result_list)
}

### Loading environmental data and parameterlists
sites <- readRDS("EnvironmentalData/sites_data.rds")

parameters = readRDS("ParameterMetaData/Parameter_mixed.rds")
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
species = "Mixed"
parameternames = c(as.character(parameters$NameRLPJ), drivernames)


library(parallel)

cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))

# results = list()
# j =1
# for(i in 1:10){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
#
# effect_list_mixed_1to10 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                         parameternames = parameternames)
#
# saveRDS(effect_list_mixed_1to10, file = "LPJrunTest/Results/Lin_Mixed_effects_1to10.rds")
# rm(results)
#
# print("done")
# results = list()
# j =1
# for(i in 11:20){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j =j +1
# }

# effect_list_mixed_11to20 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                parameternames = parameternames)
#
# saveRDS(effect_list_mixed_11to20, file = "LPJrunTest/Results/Lin_Mixed_effects_11to20.rds")
# rm(results)
# print("done")

# results = list()

# j =1
# for(i in 21:30){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_21to30 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                parameternames = parameternames)
#
# saveRDS(effect_list_mixed_21to30, file = "LPJrunTest/Results/Lin_Mixed_effects_21to30.rds")
# rm(results)
# print("done")
# results = list()
#
# j =1
# for(i in 31:40){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }

# effect_list_mixed_31to40 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                parameternames = parameternames)
#
# saveRDS(effect_list_mixed_31to40, file = "LPJrunTest/Results/Lin_Mixed_effects_31to40.rds")
# rm(results)
# print("done")

# results = list()
# j =1
# for(i in 41:50){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_41to50 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                parameternames = parameternames)
#
# saveRDS(effect_list_mixed_41to50, file = "LPJrunTest/Results/Lin_Mixed_effects_41to50.rds")
# rm(results)
# print("done")
# results = list()
#
# j = 1
# for(i in 51:60){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_51to60 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                parameternames = parameternames)
#
# saveRDS(effect_list_mixed_51to60, file = "LPJrunTest/Results/Lin_Mixed_effects_51to60.rds")
# rm(results)
# print("done")
# results = list()
#
# # j =1
# # for(i in 61:70){
# #   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
# #   j = j+1
# # }
# #
# # effect_list_mixed_61to70 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
# #                                                 parameternames = parameternames)
#
# saveRDS(effect_list_mixed_61to70, file = "LPJrunTest/Results/Lin_Mixed_effects_61to70.rds")
# rm(results)
# rm(effect_list_mixed_61to70)
# print("done")
# results = list()
# j=1
# for(i in 71:80){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_71to80 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                 parameternames = parameternames)
#
# stopCluster(cl)
# saveRDS(effect_list_mixed_71to80, file = "LPJrunTest/Results/Lin_Mixed_effects_71to80.rds")
# rm(results)
# rm(effect_list_mixed_71to80)
# print("done")
# cl <- makeCluster(11)
# clusterExport(cl,c("get_interaction_matrix","rescaling"),
#               envir=environment())
# clusterEvalQ(cl, library("ranger"))
# results = list()
# j=1
# for(i in 81:90){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_81to90 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                 parameternames = parameternames)
#
# stopCluster(cl)
# saveRDS(effect_list_mixed_81to90, file = "LPJrunTest/Results/Lin_Mixed_effects_81to90.rds")
# rm(results)
# rm(effect_list_mixed_81to90)
# print("done")
# cl <- makeCluster(11)
# clusterExport(cl,c("get_interaction_matrix","rescaling"),
#               envir=environment())
# clusterEvalQ(cl, library("ranger"))
# results = list()
# j=1
# for(i in 91:100){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_91to100 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                 parameternames = parameternames)
#
# stopCluster(cl)
# saveRDS(effect_list_mixed_91to100, file = "LPJrunTest/Results/Lin_Mixed_effects_91to100.rds")
# rm(results)
# rm(effect_list_mixed_91to100)
# print("done")
# cl <- makeCluster(11)
# clusterExport(cl,c("get_interaction_matrix","rescaling"),
#               envir=environment())
# clusterEvalQ(cl, library("ranger"))
# results = list()
# j=1
# for(i in 101:110){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_101to110 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                 parameternames = parameternames)
#
# stopCluster(cl)
# saveRDS(effect_list_mixed_101to110, file = "LPJrunTest/Results/Lin_Mixed_effects_101to110.rds")
# rm(results)
# rm(effect_list_mixed_101to110)
# print("done")
# cl <- makeCluster(11)
# clusterExport(cl,c("get_interaction_matrix","rescaling"),
#               envir=environment())
# clusterEvalQ(cl, library("ranger"))
# results = list()
# j=1
# for(i in 111:120){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_111to120 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                   parameternames = parameternames)
#
# stopCluster(cl)
# saveRDS(effect_list_mixed_111to120, file = "LPJrunTest/Results/Lin_Mixed_effects_111to120.rds")
# rm(results)
# rm(effect_list_mixed_111to120)
# print("done")
# cl <- makeCluster(11)
# clusterExport(cl,c("get_interaction_matrix","rescaling"),
#               envir=environment())
# clusterEvalQ(cl, library("ranger"))
# results = list()
# j =1
# for(i in 121:130){
#   results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
#   j = j+1
# }
#
# effect_list_mixed_121to130 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
#                                                   parameternames = parameternames)
#
# stopCluster(cl)
# saveRDS(effect_list_mixed_121to130, file = "LPJrunTest/Results/Lin_Mixed_effects_121to130.rds")
# rm(results)
# rm(effect_list_mixed_121to130)
# print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j =1
for(i in 131:140){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
  }

effect_list_mixed_131to140 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)

stopCluster(cl)
saveRDS(effect_list_mixed_131to140, file = "LPJrunTest/Results/Lin_Mixed_effects_131to140.rds")
rm(results)
rm(effect_list_mixed_131to140)
print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j=1
for(i in 141:150){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
}

effect_list_mixed_141to150 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)

stopCluster(cl)
saveRDS(effect_list_mixed_141to150, file = "LPJrunTest/Results/Lin_Mixed_effects_141to150.rds")
rm(results)
rm(effect_list_mixed_141to150)
print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j=1
for(i in 151:160){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
}

effect_list_mixed_151to160 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)

stopCluster(cl)
saveRDS(effect_list_mixed_151to160, file = "LPJrunTest/Results/Lin_Mixed_effects_151to160.rds")
rm(results)
rm(effect_list_mixed_151to160)
print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j =1
for(i in 161:170){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
  }

effect_list_mixed_161to170 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)

stopCluster(cl)
saveRDS(effect_list_mixed_161to170, file = "LPJrunTest/Results/Lin_Mixed_effects_161to170.rds")
rm(results)
rm(effect_list_mixed_161to170)
print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j =1
for(i in 171:180){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
  }

effect_list_mixed_171to180 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)

stopCluster(cl)
saveRDS(effect_list_mixed_171to180, file = "LPJrunTest/Results/Lin_Mixed_effects_171to180.rds")
rm(results)
rm(effect_list_mixed_171to180)
print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j =1
for(i in 181:190){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
  }

effect_list_mixed_181to190 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)


stopCluster(cl)
saveRDS(effect_list_mixed_181to190, file = "LPJrunTest/Results/Lin_Mixed_effects_181to190.rds")
rm(results)
rm(effect_list_mixed_181to190)
print("done")
cl <- makeCluster(11)
clusterExport(cl,c("get_interaction_matrix","rescaling"),
              envir=environment())
clusterEvalQ(cl, library("ranger"))
results = list()
j = 1
for(i in 191:200){
  results[[j]] = readRDS(paste0("./LPJrunTest/Results/Lin_Mixed_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds"))
  j = j+1
  }

effect_list_mixed_191to200 <- parallel::parLapply(cl=cl, X = results, fun= get_main_effects_per_site,
                                                  parameternames = parameternames)

saveRDS(effect_list_mixed_191to200, file = "LPJrunTest/Results/Lin_Mixed_effects_191to200.rds")
rm(results)
rm(effect_list_mixed_191to200)
print("done")


stopCluster(cl)



