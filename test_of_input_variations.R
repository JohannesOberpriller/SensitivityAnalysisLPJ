library(rLPJGUESS)
typeList <- c("cmass", "lai","nflux")
scaleLPJ_PFT <- "normal"
setwd("/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode")
mainDir <- file.path(getwd(), "LPJrunTest")

defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cru_ncep",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)


######Get the NDEP values from RCP 8.5 and 6.5 ########### 


LPJsettings_PFT$variable.ndep = "RCP26"

results_Que2 <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)

LPJsettings_PFT$variable.ndep = "RCP60"

results_Que3 <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)

LPJsettings_PFT$variable.ndep = "RCP85"

results_Que4 <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)

#########Test if environmental changes work ###############
### Need to rebuild the lpj exectuable with the new flags ############
## variable names: 
# 1. ndep_change 
# 2. insol_change 
# 3. temp_change 
# 4. ph_change 
# 5. co2_change 
# 6. prec_change 

## next steps ### 
## test effects of each of the changes ###
## since i implemented all to be additive i can see if they work 
## by setting all except one to zero. ##

### No changes ###



defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp1.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate"))


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_hurtt_RNDM_midpoint.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)



###Temperature changes ### 

defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp2.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)


#### Precipitation changes #####

defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp3.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)

#### Insol changes #####


defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp4.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)
#### co2 changes ####


defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp5.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)

#### ndep changes #####


defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp6.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)


#### ph changes #####


defaultparameters <- InferParameterAndDesignList(list(main = "/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/LPJrunTest/main_temp7.ins"),
                                                 NameMainFile = "main.ins",NamePftFile = "pft.ins",
                                                 vectorvaluedparams = c("rootdist","eps_mon",
                                                                        "storfrac_mon","photo",
                                                                        "fertdates","fertrate") )


AdjustTemplates(defaultparameters = defaultparameters$defaultparameters,
                defaultlist = defaultparameters$defaultlist,
                MainTemplateDestination = "./LPJrunTest/main_new.ins",
                PftTemplateDestination = "./LPJrunTest/pft_new.ins",
                NameMainFile = "main.ins", NamePftFile = "pft.ins")

parameters <- GetRunAbleParameters(defaultparameters = defaultparameters, PFTs = c("Que_rob"))

LPJsetup <- mainDir

LPJsettings_PFT <- list(file.co2 = file.path(mainDir,"Data/cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"Data", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir,"Data", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"Data", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir,"Data", "tas_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.prec = file.path(mainDir,"Data", "pr_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.insol = file.path(mainDir,"Data", "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew.nc"),
                        file.soildata = file.path(mainDir,"Data", "soilmap_center_interpolated.dat"),
                        template1 = "pft_new.ins",
                        template2 = "main_new.ins",
                        variable.temp = "tasAdjust", variable.insol = "rsdsAdjust",
                        variable.ndep = "RCP45",
                        variable.prec = "prAdjust", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_sensi.txt", parallel = "parameters", 
                        defaultlist = parameters$defaultfiles)



LPJsettings_PFT$design <- parameters$design

results_Que <- runLPJ(x = mainDir, parameterList = parameters$runParameters,
                      typeList = typeList, settings = LPJsettings_PFT)
