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

get_main_effects <- function(sites, parameternames, species){

  # sites: a list of sites
  # parameternames: the names of the parameters
  # species: the species for which the things should be calculated


  ## setup empty matrixes to store the results
  effects_cmass_complete = matrix(nrow = nrow(sites), ncol = 39)
  effects_cmass_climate_change = matrix(nrow = nrow(sites), ncol = 39)
  effects_cmass_steady_climate = matrix(nrow = nrow(sites), ncol = 39)

  effects_cmass_complete_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_cmass_climate_change_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_cmass_steady_climate_rf = matrix(nrow = nrow(sites), ncol = 39)

  effects_cpool_complete = matrix(nrow = nrow(sites), ncol = 39)
  effects_cpool_climate_change = matrix(nrow = nrow(sites), ncol = 39)
  effects_cpool_steady_climate = matrix(nrow = nrow(sites), ncol = 39)

  effects_cpool_complete_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_cpool_climate_change_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_cpool_steady_climate_rf = matrix(nrow = nrow(sites), ncol = 39)


  effects_cflux_complete = matrix(nrow = nrow(sites), ncol = 39)
  effects_cflux_climate_change = matrix(nrow = nrow(sites), ncol = 39)
  effects_cflux_steady_climate = matrix(nrow = nrow(sites), ncol = 39)

  effects_cflux_complete_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_cflux_climate_change_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_cflux_steady_climate_rf = matrix(nrow = nrow(sites), ncol = 39)


  effects_anpp_complete = matrix(nrow = nrow(sites), ncol = 39)
  effects_anpp_climate_change = matrix(nrow = nrow(sites), ncol = 39)
  effects_anpp_steady_climate = matrix(nrow = nrow(sites), ncol = 39)

  effects_anpp_complete_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_anpp_climate_change_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_anpp_steady_climate_rf = matrix(nrow = nrow(sites), ncol = 39)


  effects_agpp_complete = matrix(nrow = nrow(sites), ncol = 39)
  effects_agpp_climate_change = matrix(nrow = nrow(sites), ncol = 39)
  effects_agpp_steady_climate = matrix(nrow = nrow(sites), ncol = 39)

  effects_agpp_complete_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_agpp_climate_change_rf = matrix(nrow = nrow(sites), ncol = 39)
  effects_agpp_steady_climate_rf = matrix(nrow = nrow(sites), ncol = 39)



  interactions_cmass_complete = list()
  interactions_cmass_climate_change = list()
  interactions_cmass_steady_climate = list()

  interactions_cpool_complete = list()
  interactions_cpool_climate_change = list()
  interactions_cpool_steady_climate = list()

  interactions_cflux_complete = list()
  interactions_cflux_climate_change = list()
  interactions_cflux_steady_climate = list()

  interactions_anpp_complete = list()
  interactions_anpp_climate_change = list()
  interactions_anpp_steady_climate = list()

  interactions_agpp_complete = list()
  interactions_agpp_climate_change = list()
  interactions_agpp_steady_climate = list()

  # iterate over the sites

  for(i in 1:nrow(sites)){
    # load the parameters
    parameter_name = paste0(species,"_",sites[i,"Longitudinal"],"_",sites[i,"Latitudinal"],".rds")
    # load the results for the site
    results =  try(readRDS(paste0("./LPJrunTest/Results/",parameter_name)))

    if("try-error" %in% class(results)){
      main_effects_Pic_abi = rep(NA,ncol(parameters))
      interaction_matrix = NA
    }
    else{
      parameter_matrix = matrix(ncol = 38, nrow = length(results$parameters))
      colnames(parameter_matrix) = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameternames)])
      for(j in 1:length(results$parameters)){
        parameter_matrix[j,] = as.numeric(results$parameters[[j]][which(names(results$parameters[[j]]) %in% parameternames)])
      }
      parameters = apply(parameter_matrix, 2, rescaling)

      ## fitting of the linear regression and random forest to check consistency ##
      ### cpool stuff ###

      ### complete ###

      overall_cpool_complete = results[['cpool']][['complete']][1,]

      data_cpool_complete <- cbind(overall_cpool_complete,parameters)
      fit <- lm(overall_cpool_complete ~ (.)^2 , data = as.data.frame(data_cpool_complete))
      coefficients = coef(fit)[-1]
      effects_cpool_complete[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cpool_complete[i,ncol(effects_cpool_complete)] = results[['cpool']][['complete']][2,1]

      interaction_matrix_cpool_complete = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_cpool_complete[[i]] = interaction_matrix_cpool_complete

      fit_rf <- ranger(overall_cpool_complete ~ ., data = as.data.frame(data_cpool_complete), importance = 'impurity', num.threads = 10)
      effects_cpool_complete_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cpool_complete_rf[i,ncol(effects_cpool_complete)] = results[['cpool']][['complete']][2,1]


      ### climate change ###

      overall_cpool_climate_change = results[['cpool']][['climate_change']][1,]

      data_cpool_climate_change <- cbind(overall_cpool_climate_change,parameters)
      fit <- lm(overall_cpool_climate_change ~ (.)^2 , data = as.data.frame(data_cpool_climate_change))
      coefficients = coef(fit)[-1]
      effects_cpool_climate_change[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cpool_climate_change[i,ncol(effects_cpool_complete)] = results[['cpool']][['climate_change']][2,1]

      interaction_matrix_cpool_climate_change = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_cpool_climate_change[[i]] = interaction_matrix_cpool_climate_change

      fit_rf <- ranger(overall_cpool_climate_change ~ . , data = as.data.frame(data_cpool_climate_change), importance = 'impurity',num.threads = 10)
      effects_cpool_climate_change_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cpool_climate_change_rf[i,ncol(effects_cpool_climate_change)] = results[['cpool']][['climate_change']][2,1]

      ### steady climate ###

      overall_cpool_steady_climate= results[['cpool']][['steady_climate']][1,]

      data_cpool_steady_climate <- cbind(overall_cpool_steady_climate,parameters)
      fit <- lm(overall_cpool_steady_climate ~ (.)^2 , data = as.data.frame(data_cpool_steady_climate))
      coefficients = coef(fit)[-1]
      effects_cpool_steady_climate[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cpool_steady_climate[i,ncol(effects_cpool_complete)] = results[['cpool']][['steady_climate']][2,1]

      interaction_matrix_cpool_steady_climate = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_cpool_steady_climate[[i]] = interaction_matrix_cpool_steady_climate

      fit_rf <- ranger(overall_cpool_steady_climate ~ . , data = as.data.frame(data_cpool_steady_climate), importance = 'impurity',num.threads = 10)
      effects_cpool_steady_climate_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cpool_steady_climate_rf[i,ncol(effects_cpool_steady_climate)] = results[['cpool']][['steady_climate']][2,1]

      ## cflux related stuff

      ### complete ###

      overall_cflux_complete = results[['cflux']][['complete']][1,]

      data_cflux_complete <- cbind(overall_cflux_complete,parameters)
      fit <- lm(overall_cflux_complete ~ (.)^2 , data = as.data.frame(data_cflux_complete))
      coefficients = coef(fit)[-1]
      effects_cflux_complete[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cflux_complete[i,ncol(effects_cflux_complete)] = results[['cflux']][['complete']][2,1]

      interaction_matrix_cflux_complete = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_cflux_complete[[i]] = interaction_matrix_cflux_complete

      fit_rf <- ranger(overall_cflux_complete ~ ., data = as.data.frame(data_cflux_complete), importance = 'impurity',num.threads = 10)
      effects_cflux_complete_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cflux_complete_rf[i,ncol(effects_cflux_complete)] = results[['cflux']][['complete']][2,1]
      ### climate change ###

      overall_cflux_climate_change = results[['cflux']][['climate_change']][1,]

      data_cflux_climate_change <- cbind(overall_cflux_climate_change,parameters)
      fit <- lm(overall_cflux_climate_change ~ (.)^2 , data = as.data.frame(data_cflux_climate_change))
      coefficients = coef(fit)[-1]
      effects_cflux_climate_change[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cflux_climate_change[i,ncol(effects_cflux_complete)] = results[['cflux']][['climate_change']][2,1]

      interaction_matrix_cflux_climate_change = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_cflux_climate_change[[i]] = interaction_matrix_cflux_climate_change

      fit_rf <- ranger(overall_cflux_climate_change ~ . , data = as.data.frame(data_cflux_climate_change), importance = 'impurity',num.threads = 10)
      effects_cflux_climate_change_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cflux_climate_change_rf[i,ncol(effects_cflux_climate_change)] = results[['cflux']][['climate_change']][2,1]
      ### steady climate ###

      overall_cflux_steady_climate= results[['cflux']][['steady_climate']][1,]

      data_cflux_steady_climate <- cbind(overall_cflux_steady_climate,parameters)
      fit <- lm(overall_cflux_steady_climate ~ (.)^2 , data = as.data.frame(data_cflux_steady_climate))
      coefficients = coef(fit)[-1]
      effects_cflux_steady_climate[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_cflux_steady_climate[i,ncol(effects_cflux_complete)] = results[['cflux']][['steady_climate']][2,1]

      interaction_matrix_cflux_steady_climate = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_cflux_steady_climate[[i]] = interaction_matrix_cflux_steady_climate

      fit_rf <- ranger(overall_cflux_steady_climate ~ . , data = as.data.frame(data_cflux_steady_climate), importance = 'impurity',num.threads = 10)
      effects_cflux_steady_climate_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_cflux_steady_climate_rf[i,ncol(effects_cflux_steady_climate)] = results[['cflux']][['steady_climate']][2,1]

      ## agpp related stuff

      ### complete ###

      overall_agpp_complete = results[['agpp']][['complete']][1,]

      data_agpp_complete <- cbind(overall_agpp_complete,parameters)
      fit <- lm(overall_agpp_complete ~ (.)^2 , data = as.data.frame(data_agpp_complete))
      coefficients = coef(fit)[-1]
      effects_agpp_complete[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_agpp_complete[i,ncol(effects_agpp_complete)] = results[['agpp']][['complete']][2,1]

      interaction_matrix_agpp_complete = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_agpp_complete[[i]] = interaction_matrix_agpp_complete

      fit_rf <- ranger(overall_agpp_complete ~ . , data = as.data.frame(data_agpp_complete), importance = 'impurity',num.threads = 10)
      effects_agpp_complete_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_agpp_complete_rf[i,ncol(effects_agpp_complete)] = results[['agpp']][['complete']][2,1]
      ### climate change ###

      overall_agpp_climate_change = results[['agpp']][['climate_change']][1,]

      data_agpp_climate_change <- cbind(overall_agpp_climate_change,parameters)
      fit <- lm(overall_agpp_climate_change ~ (.)^2 , data = as.data.frame(data_agpp_climate_change))
      coefficients = coef(fit)[-1]
      effects_agpp_climate_change[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_agpp_climate_change[i,ncol(effects_agpp_complete)] = results[['agpp']][['climate_change']][2,1]

      interaction_matrix_agpp_climate_change = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_agpp_climate_change[[i]] = interaction_matrix_agpp_climate_change

      fit_rf <- ranger(overall_agpp_climate_change ~ . , data = as.data.frame(data_agpp_climate_change), importance = 'impurity',num.threads = 10)
      effects_agpp_climate_change_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_agpp_climate_change_rf[i,ncol(effects_agpp_climate_change)] = results[['agpp']][['climate_change']][2,1]
      ### steady climate ###

      overall_agpp_steady_climate= results[['agpp']][['steady_climate']][1,]

      data_agpp_steady_climate <- cbind(overall_agpp_steady_climate,parameters)
      fit <- lm(overall_agpp_steady_climate ~ (.)^2 , data = as.data.frame(data_agpp_steady_climate))
      coefficients = coef(fit)[-1]
      effects_agpp_steady_climate[i,1:ncol(parameters)] =  coefficients[1:ncol(parameters)]
      effects_agpp_steady_climate[i,ncol(effects_agpp_complete)] = results[['agpp']][['steady_climate']][2,1]

      interaction_matrix_agpp_steady_climate = get_interaction_matrix(coefficients, num_parameters = 38)
      interactions_agpp_steady_climate[[i]] = interaction_matrix_agpp_steady_climate

      fit_rf <- ranger(overall_agpp_steady_climate ~ . , data = as.data.frame(data_agpp_steady_climate), importance = 'impurity',num.threads = 10)
      effects_agpp_steady_climate_rf[i,1:ncol(parameters)] = fit_rf$variable.importance
      effects_agpp_steady_climate_rf[i,ncol(effects_agpp_steady_climate)] = results[['agpp']][['complete']][2,1]



    }
  }


  result_list = list(
                                  "agpp" = list("steady_climate" = effects_agpp_steady_climate,
                                                "climate_change" = effects_agpp_climate_change,
                                                "complete" = effects_agpp_complete),
                     "agpp_rf" = list("steady_climate" = effects_agpp_steady_climate_rf,
                                      "climate_change" = effects_agpp_climate_change_rf,
                                      "complete" = effects_agpp_complete_rf),
                     "agpp_interactions" = list("steady_climate" = interactions_agpp_steady_climate,
                                                "climate_change" = interactions_agpp_climate_change,
                                                "complete" = interactions_agpp_complete),
                                  "cpool" = list("steady_climate" = effects_cpool_steady_climate,
                                                 "climate_change" = effects_cpool_climate_change,
                                                 "complete" = effects_cpool_complete),
                     "cpool_rf" = list("steady_climate" = effects_cpool_steady_climate_rf,
                                      "climate_change" = effects_cpool_climate_change_rf,
                                      "complete" = effects_cpool_complete_rf),
                     "cpool_interactions" = list("steady_climate" = interactions_cpool_steady_climate,
                                                "climate_change" = interactions_cpool_climate_change,
                                                "complete" = interactions_cpool_complete),
                                  'cflux'= list("steady_climate" = effects_cflux_steady_climate,
                                                "climate_change" = effects_cflux_climate_change,
                                                "complete" = effects_cflux_complete),
                     "cflux_rf" = list("steady_climate" = effects_cflux_steady_climate_rf,
                                      "climate_change" = effects_cflux_climate_change_rf,
                                      "complete" = effects_cflux_complete_rf),
                     "cflux_interactions" = list("steady_climate" = interactions_cflux_steady_climate,
                                                "climate_change" = interactions_cflux_climate_change,
                                                "complete" = interactions_cflux_complete))

  return(result_list)
}

### Loading environmental data and parameterlists
sites <- readRDS("EnvironmentalData/sites_data.rds")

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameters_pic_abi = c(as.character(parameters$Pic_abi$NameRLPJ), drivernames)
parameters_fag_syl = c(as.character(parameters$Fag_syl$NameRLPJ), drivernames)
parameters_pin_syl = c(as.character(parameters$Pin_syl$NameRLPJ), drivernames)

effect_list_Pic_abi <- get_main_effects(sites = sites, parameternames = parameters_pic_abi, species = "Lin_Pic_abi")

saveRDS(effect_list_Pic_abi, file = "LPJrunTest/Results/Pic_abi_effects_lin2.rds")

effect_list_Pin_syl <- get_main_effects(sites = sites, parameternames = parameters_pin_syl, species = "Lin_Pin_syl")

saveRDS(effect_list_Pin_syl, "LPJrunTest/Results/Pin_syl_effects_lin2.rds")

effect_list_Fag_syl <- get_main_effects(sites = sites, parameternames = parameters_fag_syl, species = "Lin_Fag_syl")

saveRDS(effect_list_Fag_syl, "LPJrunTest/Results/Fag_syl_effects_lin2.rds")


