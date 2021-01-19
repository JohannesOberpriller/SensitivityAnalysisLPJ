library(rLPJGUESS)
##specify outputs one wants to look at
typeList <- c("cmass", "lai","agpp","cpool","anpp","cflux","fpc","speciesdiam","dens")
scaleLPJ_PFT <- "normal"
# set main dir to LPJrunTest
mainDir <- file.path(getwd(), "LPJrunTest")

# infer the default parameters from main
defaultparameters <- InferParameterAndDesignList(list(main = paste0(mainDir,"/main_temp_mono.ins")),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate"))
# adjust the templates in order to be used
AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

# generate runable parameters for the species one wants to investigate
parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Fag_syl","Pic_abi","Pin_syl"))

gridList_run = "path_gridlist"

# specify set-up of the run (path to the data, and so on)

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tasAdjust_rcp_45_1860_2100_mon_mean_bias_corrected.nc"),
                        file.prec = file.path(mainDir,"Data", "prAdjust_rcp_45_1860_2100_mon_mean_bias_corrected.nc"),
                        file.insol = file.path(mainDir,"Data", "rsdsAdjust_rcp_45_1860_2100_mon_mean_bias_corrected.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        file.wetdays = file.path(mainDir,"Data","wetdays_rcp_45_1860_2100_mon_sum_bias_corrected6.nc"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = gridList_run, parallel = "parameters",
                        defaultlist = parameters$defaultfiles)


# specify parallel setup
numCores <- 24
LPJsetup <- setupLPJParallel(numCores = numCores, clusterType = "SOCK",
                             mainDir = mainDir)

LPJsettings_PFT$design <- parameters$design


## Parameter with name from the

sensitivity_parameters = readRDS("input_parameters")

drivers = c("ph","co2","prec","temp","insol","ndep")

designs = c("run_nfix_a","run_nfix_b", "run_distinterval","run_nrelocfrac")

drivernames = paste0("run_",drivers,"_change")

parameters$design = parameters$design[-which(names(parameters$design) %in% c(drivernames,designs))]

## reset some parameters

LPJsettings_PFT$design <- parameters$design

for(i in 1:length(drivernames)){
  parameters$runParameters[[drivernames[i]]] = "0"
}

for(i in 1:length(designs)){
  parameters$runParameters[[designs[i]]] = "0"
}

parameters_new = matrix(rep(unlist(parameters$runParameters),nrow(sensitivity_parameters)), nrow = nrow(sensitivity_parameters),
                        ncol = length(parameters$runParameters), byrow = T)

colnames(parameters_new) = names(parameters$runParameters)
names_sensitivity_parameters = colnames(sensitivity_parameters)


for(i in 1:length(names_sensitivity_parameters)){
  parameters_new[,names_sensitivity_parameters[i]] = sensitivity_parameters[,i]
}


colnames(parameters_new) = names(parameters$runParameters)

## run simulations
results_Que <- runLPJ(x = LPJsetup, parameterList = parameters_new,
                      typeList = typeList, settings = LPJsettings_PFT)

saveRDS(results_Que,"LPJrunTest/Results/_save_the_file_.rds")
