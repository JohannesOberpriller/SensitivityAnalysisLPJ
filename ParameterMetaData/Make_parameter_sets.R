## Load the sites we selected for a sensitivitiy analysis

sites <- readRDS("EnvironmentalData/sites_data.rds")

## write the long and lat coordinates to a txt file for later use

gridlist = sites[,c("Longitudinal","Latitudinal")]
write.table(gridlist, file = "LPJrunTest/gridlist_sensi.txt",sep="\t",row.names=FALSE, col.names = F)

## determine a mastersheet to get an order for runnning
mastersheet = data.frame(matrix(nrow = 600, ncol = 6))

colnames(mastersheet) = c("site_id","lon","lat","Fag_syl","Pic_abi", "Pin_syl")

mastersheet[,c("site_id","lon","lat")] = sites[,c("Site","Longitudinal","Latitudinal")]

mastersheet[1:200,4:6] = matrix(rep(c(1,0,0),200),ncol = 3, nrow = 200, byrow = T )
mastersheet[201:400,4:6] = matrix(rep(c(0,1,0),200),ncol = 3, nrow = 200, byrow = T )
mastersheet[401:600,4:6] = matrix(rep(c(0,0,1),200),ncol = 3, nrow = 200, byrow = T )

saveRDS(mastersheet, "./Templates/mastersheet.rds")

## load parameters, order them and group them for the species

parameterranges = read.csv("./ParameterMetaData/SyntheseRangesWithNames.csv", header = T, sep = ";")

parameterranges_common = parameterranges[1:11,]

parameterranges_pinsyl = parameterranges[12:33,]

parameterranges_picabi = parameterranges[34:55,]
parameterranges_picabi = parameterranges_picabi[-which(parameterranges_picabi[,"Default"] == "siehe allgemein"),]

parameterranges_fagsyl = parameterranges[56:77,]
parameterranges_fagsyl = parameterranges_fagsyl[-which(parameterranges_fagsyl[,"Default"] == "siehe allgemein"),]


Parameters_FagSyl = rbind(parameterranges_common, parameterranges_fagsyl)
Parameters_PicAbi = rbind(parameterranges_common, parameterranges_picabi)
Parameters_PinSyl = rbind(parameterranges_common, parameterranges_pinsyl)
Parameters_PinSyl = Parameters_PinSyl[-which(Parameters_PinSyl[,"NameRLPJ"] == "tree_k_allom2"),]

## save the parameters with ranges into an rds file for later use

Parameter_list = list("Fag_syl" = Parameters_FagSyl,
                     "Pic_abi" = Parameters_PicAbi,
                     "Pin_syl" = Parameters_PinSyl)


saveRDS(Parameter_list,"./ParameterMetaData/Parameter_list.rds")

## read the mixed parameters and save also them

parameterranges_mixed = read.csv("./ParameterMetaData/SyntheseRangesMixed.csv", header = T, sep = ";")

saveRDS(parameterranges_mixed, "./ParameterMetaData/Parameter_mixed.rds")




## function to write the data for the grids

make_gridlists <- function(gridlist){

  for(i in 1:nrow(sites)){
    individual_site = paste0(gridlist[i,"Longitudinal"]," ",gridlist[i,"Latitudinal"])
    if(file.exists(paste0("grids/",gridlist[i, "Longitudinal"],"_",gridlist[i,"Latitudinal"],".txt"))){
      writeLines(text = individual_site,
      con = paste0("grids/",gridlist[i, "Longitudinal"],"_",gridlist[i,"Latitudinal"],".txt"))
    }
    else{
      file.create(paste0("grids/",gridlist[i, "Longitudinal"],"_",gridlist[i,"Latitudinal"],".txt"))
      writeLines(text = individual_site,
                 con = paste0("grids/",gridlist[i, "Longitudinal"],"_",gridlist[i,"Latitudinal"],".txt"))
    }
  }
}




## function to make parameters readable for later use

make_parameter_data_frames <- function(mastersheet, climaticranges, white_parameters, parameter_range_list){


  # get the site you want --> get climatic ranges

  site_id = mastersheet["site_id"]
  print(site_id)

  ## need this as an input later

  names_climatic = names(climaticranges)
  climatic_ranges_site = matrix(ncol = length(climaticranges), nrow = 2)
  rownames(climatic_ranges_site) = c("minimum", "maximum")
  for(i in 1:length(names_climatic)){
    climatic_ranges_site[,i] = climaticranges[[i]][2:3,site_id]
  }

  colnames(climatic_ranges_site) = paste0("run_", names(climaticranges),"_change")

  # get what is simulated --> therefor get the parameter ranges


  if(mastersheet["Fag_syl"] == 1){
    parameters = t(parameter_range_list[["Fag_syl"]][,c("Minimum.Value","Maximum.Value")])
    colnames(parameters) = parameter_range_list[["Fag_syl"]][,"NameRLPJ"]
  }
  else if(mastersheet["Pic_abi"] == 1){
    parameters = t(parameter_range_list[["Pic_abi"]][,c("Minimum.Value","Maximum.Value")])
    colnames(parameters) = parameter_range_list[["Pic_abi"]][,"NameRLPJ"]
  }
  else{
    parameters = t(parameter_range_list[["Pin_syl"]][,c("Minimum.Value","Maximum.Value")])
    colnames(parameters) = parameter_range_list[["Pin_syl"]][,"NameRLPJ"]
  }


  ## translate the "white" parameters to real parameters


  complete_parameters = cbind(parameters,climatic_ranges_site)

  real_parameters = get_parameters_from_ranges(white_parameters,
                                               complete_parameters)

  longitudinal_coordinate = mastersheet["lon"]
  latitudinal_coordinate = mastersheet["lat"]

  if(mastersheet["Fag_syl"] == 1){
    saveRDS(real_parameters,paste0("./parameters/Fag_syl_",longitudinal_coordinate,"_",
                                   latitudinal_coordinate,".rds"), version = 2)
  }
  else if(mastersheet["Pic_abi"] == 1){
    saveRDS(real_parameters,paste0("./parameters/Pic_abi_",longitudinal_coordinate,"_",
                                   latitudinal_coordinate,".rds"), version = 2)
  }
  else{
    saveRDS(real_parameters,paste0("./parameters/Pin_syl_",longitudinal_coordinate,"_",
                                   latitudinal_coordinate,".rds"), version = 2)
  }
}

get_parameters_from_ranges <- function(white_parameters, complete_parameters){
  real_parameters = matrix(ncol = ncol(white_parameters), nrow = nrow(white_parameters))
  names_variables = colnames(complete_parameters)
  complete_parameters <- gsub("," , "." , complete_parameters)
  numeric_parameters <- mapply(complete_parameters, FUN=as.numeric)
  complete_parameters = matrix(numeric_parameters, nrow = nrow(complete_parameters),
                               ncol = ncol(complete_parameters))
  for(i in 1:nrow(white_parameters)){
    real_parameters[i,] = complete_parameters[1,] +
      white_parameters[i,]*(complete_parameters[2,]
                            - complete_parameters[1,])
  }
  real_parameters =  matrix(as.character(real_parameters),
                            nrow = nrow(real_parameters),
                            ncol = ncol(real_parameters))
  colnames(real_parameters) = names_variables
  return(real_parameters)
}


### generate the gridlists

make_gridlists(gridlist)

### generate the parameter_blocks

climaticranges <- readRDS("./EnvironmentalData/climateranges_list.rds")
parameter_list <- readRDS("./ParameterMetaData/Parameter_list.rds")
mastersheet <- readRDS("./Templates/mastersheet.rds")


white_parameters_mono = readRDS("./Templates/white_parameters_mono.rds")

## generate and save parameters to files

apply(X = mastersheet, FUN = make_parameter_data_frames, climaticranges = climaticranges,
       white_parameters = white_parameters_mono, parameter_range_list = parameter_list, MARGIN = 1)


## function to generate the rscripts

generate_R_scripts_mono <- function(mastersheet){

  # parameters of the function:
  # mastersheet: containing all combinations of simulations
  # climaticranges: ranges of climate per site
  # white_parameters: parameters always the same which are drawn from 0 to 1
  # parameter_range_list: list with parameters for the different mono trees
  # and the mixed stands

  # get the site you want --> get climatic ranges

  site_id = mastersheet["site_id"]


  # make input files and for each side

  longitudinal_coordinate = mastersheet["lon"]
  latitudinal_coordinate = mastersheet["lat"]
  # pick the correct mastertempla
  if(mastersheet["Fag_syl"] == 1){
    plain_text = readLines("Templates/template_Fag_syl.R")
    gridlist = paste0(longitudinal_coordinate,"_",latitudinal_coordinate,".txt")
    plain_text = sub("path_gridlist",gridlist,plain_text)
    plain_text = sub("input_parameters",paste0("parameters/Fag_syl_",longitudinal_coordinate,"_",
                                               latitudinal_coordinate,".rds"), plain_text)
    plain_text = sub("_save_the_file_",paste0("Fag_syl_",longitudinal_coordinate,"_",
                                               latitudinal_coordinate),plain_text)
    writeLines(plain_text,paste0("R_scripts/Fag_syl_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".R"))
  }
  else if(mastersheet["Pic_abi"] == 1){
    plain_text = readLines("Templates/template_Pic_abi.R")
    gridlist = paste0(longitudinal_coordinate,"_",latitudinal_coordinate,".txt")
    plain_text = sub("path_gridlist",gridlist,plain_text)
    plain_text = sub("input_parameters",paste0("parameters/Pic_abi_",longitudinal_coordinate,"_",
                                               latitudinal_coordinate,".rds"), plain_text)
    plain_text = sub("_save_the_file_",paste0("Pic_abi_",longitudinal_coordinate,"_",
                                              latitudinal_coordinate),plain_text)
    writeLines(plain_text,paste0("R_scripts/Pic_abi_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".R"))
  }
  else{
    plain_text = readLines("Templates/template_Pin_syl.R")
    gridlist = paste0(longitudinal_coordinate,"_",latitudinal_coordinate,".txt")
    plain_text = sub("path_gridlist",gridlist,plain_text)
    plain_text = sub("input_parameters",paste0("parameters/Pin_syl_",longitudinal_coordinate,"_",
                                               latitudinal_coordinate,".rds"), plain_text)
    plain_text = sub("_save_the_file_",paste0("Pin_syl_",longitudinal_coordinate,"_",
                                              latitudinal_coordinate),plain_text)
    writeLines(plain_text,paste0("R_scripts/Pin_syl_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".R"))
  }


}

## generate r scripts and save them in order to run the model later

apply(X = mastersheet, FUN = generate_R_scripts_mono, MARGIN = 1)

generate_submit_scripts <- function(mastersheet){
  site_id = mastersheet["site_id"]


  # make input files and for each side

  longitudinal_coordinate = mastersheet["lon"]
  latitudinal_coordinate = mastersheet["lat"]
  # pick the correct mastertemplate
  if(mastersheet["Fag_syl"] == 1){
    plain_text = readLines("Templates/submit_template.sh")
    plain_text = sub("_script_to_execute_",noquote(paste0("Fag_syl_",longitudinal_coordinate,"_",
                                               latitudinal_coordinate)), plain_text)
    writeLines(plain_text,paste0("submit_scripts/Fag_syl_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".sh"))
  }
  else if(mastersheet["Pic_abi"] == 1){
    plain_text = readLines("Templates/submit_template.sh")
    plain_text = sub("_script_to_execute_",noquote(paste0("Pic_abi_",longitudinal_coordinate,"_",
                                                  latitudinal_coordinate)), plain_text)
    writeLines(plain_text,paste0("submit_scripts/Pic_abi_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".sh"))
  }
  else{
    plain_text = readLines("Templates/submit_template.sh")
    plain_text = sub("_script_to_execute_",noquote(paste0("Pin_syl_",longitudinal_coordinate,"_",
                                                  latitudinal_coordinate)), plain_text)
    writeLines(plain_text,paste0("submit_scripts/Pin_syl_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".sh"))
  }

}
## generate submit scripts in order to run them later

apply(X = mastersheet, FUN = generate_submit_scripts, MARGIN = 1)

##### make parameter sets for the mixed stands ######

make_parameter_data_frames_mixed <- function(mastersheet, climaticranges, white_parameters, parameter_range_list){


  # get the site you want --> get climatic ranges

  site_id = mastersheet["site_id"]

  ## need this as an input later

  names_climatic = names(climaticranges)
  climatic_ranges_site = matrix(ncol = length(climaticranges), nrow = 2)
  rownames(climatic_ranges_site) = c("minimum", "maximum")
  for(i in 1:length(names_climatic)){
    climatic_ranges_site[,i] = climaticranges[[i]][2:3,site_id]
  }

  colnames(climatic_ranges_site) = paste0("run_", names(climaticranges),"_change")

  # get what is simulated --> therefor get the parameter ranges


  parameters = t(parameter_range_list[,c("Minimum.Value","Maximum.Value")])
  colnames(parameters) = parameter_range_list[,"NameRLPJ"]



  ## translate the "white" parameters to real parameters


  complete_parameters = cbind(parameters,climatic_ranges_site)

  real_parameters = get_parameters_from_ranges(white_parameters,
                                               complete_parameters)

  longitudinal_coordinate = mastersheet["lon"]
  latitudinal_coordinate = mastersheet["lat"]

  saveRDS(real_parameters,paste0("./parameters/Mixed_",longitudinal_coordinate,"_",
                                   latitudinal_coordinate,".rds"), version = 2)

}

## read in the pure data for the mixed simulations
parameterranges_mixed <- readRDS("ParameterMetaData/Parameter_mixed.rds")
climaticranges <- readRDS("EnvironmentalData/climateranges_list.rds")
parameter_list <- readRDS("ParameterMetaData/Parameter_list.rds")
mastersheet <- readRDS("Templates/mastersheet.rds")

white_parameters_mixed = readRDS("Templates/white_parameters_mixed.rds")

## make parameter sets and save them

apply(X = mastersheet[1:200,], FUN = make_parameter_data_frames_mixed, climaticranges = climaticranges,
      white_parameters = white_parameters_mixed, parameter_range_list = parameterranges_mixed, MARGIN = 1)

### parameter ranges for the mixed case ###

get_parameters_from_ranges <- function(white_parameters, complete_parameters){
  real_parameters = matrix(ncol = ncol(white_parameters), nrow = nrow(white_parameters))
  names_variables = colnames(complete_parameters)
  complete_parameters <- gsub("," , "." , complete_parameters)
  numeric_parameters <- mapply(complete_parameters, FUN=as.numeric)
  complete_parameters = matrix(numeric_parameters, nrow = nrow(complete_parameters),
                               ncol = ncol(complete_parameters))
  for(i in 1:nrow(white_parameters)){
    real_parameters[i,] = complete_parameters[1,] +
      white_parameters[i,]*(complete_parameters[2,]
                            - complete_parameters[1,])
  }
  real_parameters =  matrix(as.character(real_parameters),
                            nrow = nrow(real_parameters),
                            ncol = ncol(real_parameters))
  colnames(real_parameters) = names_variables
  return(real_parameters)
}


#### write the scripts for the mixed stands



generate_R_scripts_mixed <- function(mastersheet_mixed){


  # get the site you want --> get climatic ranges

  site_id = mastersheet_mixed["site_id"]
  print(site_id)

  # make input files and for each side

  longitudinal_coordinate = mastersheet_mixed["lon"]
  latitudinal_coordinate = mastersheet_mixed["lat"]

  # pick the correct mastertemplate

  plain_text = readLines("Templates/template_mixed.R")
  gridlist = paste0(longitudinal_coordinate,"_",latitudinal_coordinate,".txt")
  plain_text = sub("path_gridlist",gridlist,plain_text)
  plain_text = sub("input_parameters",paste0("parameters/Mixed_",longitudinal_coordinate,"_",
                                               latitudinal_coordinate,".rds"), plain_text)
  plain_text = sub("_save_the_file_",paste0("Mixed_",longitudinal_coordinate,"_",
                                            latitudinal_coordinate),plain_text)
  writeLines(text = plain_text, con = paste0("R_scripts/Mixed_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".R"))



}


apply(X = mastersheet[1:200,], FUN = generate_R_scripts_mixed, MARGIN = 1)

## function to generate the submit scripts on the cluster in order for running and save them
generate_submit_scripts_mixed <- function(mastersheet){
  site_id = mastersheet["site_id"]


  # make input files and for each side

  longitudinal_coordinate = mastersheet["lon"]
  latitudinal_coordinate = mastersheet["lat"]
  # pick the correct mastertemplate

  plain_text = readLines("Templates/submit_template_mixed.sh")
  plain_text = gsub("_script_to_execute_",noquote(paste0("Mixed_",longitudinal_coordinate,"_",
                                                  latitudinal_coordinate)), plain_text)
  writeLines(plain_text,paste0("submit_scripts/Mixed_",longitudinal_coordinate,"_",
                                 latitudinal_coordinate,".sh"))


}
## acutally generate submit scripts
apply(X = mastersheet[1:200,], FUN = generate_submit_scripts_mixed, MARGIN = 1)


