## This is the sensitivity analysis for the next point to ratisbona
## in the data and model given by Anja and Andy

## Data generation ##


graphics.off()
rm(list = ls())
library(rLPJGUESS)

typeList <- c("cmass", "lai")
scaleLPJ_SPP <- "global"
scaleLPJ_PFT <- "global"

LPJparameters_SPP <- getParameterList(scale = scaleLPJ_SPP, list = F)
LPJparameters_PFT <- getParameterList(scale = scaleLPJ_PFT, list = F)

setwd("/home/johannes/Documents/PhD/Code/LPJ_Guess_Project/R_Code/rLPJGUESS/rLPJGUESS_Freising/Testcode")
mainDir <- file.path(getwd(), "LPJrunTest")


LPJsettings_SPP <- list(file.co2 = file.path(mainDir, "crudata", "cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"crudata", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir, "crudata", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"crudata", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir, "temp.nc"),
                        file.prec = file.path(mainDir, "prec.nc"),
                        file.insol = file.path(mainDir, "rad.nc"),
                        template1 = "global.ins",
                        template2 = "guess.ins",
                        variable.temp = "temp", variable.insol = "rad",
                        variable.prec = "prec", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_SPP, mode = "cf",
                        gridList = "gridlist_nextratisbona.txt")
#
#
designLPJ <- getDesign(scaleLPJ_SPP, list = T)
designLPJ$run_vegmode <- "cohort"
designLPJ$run_ifcentury <- 0
designLPJ$run_iffire <- 0
designLPJ$run_ifnlim <- 0
designLPJ$run_ifstochestab <- 0
designLPJ$run_ifstochmort <- 0
designLPJ$run_patcharea <- 25^2
designLPJ$run_npatch <- 1
designLPJ$run_ifdisturb <- 0
designLPJ$run_nyear_spinup <- 1
designLPJ$run_freenyears <- 0
designLPJ$run_save_state <- 0
designLPJ$run_restart <- 0
designLPJ$run_state_path <- mainDir
LPJsettings_SPP$design <- designLPJ


LPJsetup <- mainDir


LPJsettings_PFT <- list(file.co2 = file.path(mainDir, "crudata", "cmip5_rec_co2_rcp4p5_1901_2100.txt"),
                        file.cru = file.path(mainDir,"crudata", "Cruncep_1901_2015.bin"),
                        file.cru.misc = file.path(mainDir, "crudata", "Cruncep_1901_2015misc.bin"),
                        file.ndep = file.path(mainDir,"crudata", "GlobalNitrogenDeposition.bin"),
                        file.temp = file.path(mainDir, "temp.nc"),
                        file.prec = file.path(mainDir, "prec.nc"),
                        file.insol = file.path(mainDir, "rad.nc"),
                        template1 = "global.ins",
                        template2 = "guess.ins",
                        variable.temp = "temp", variable.insol = "rad",
                        variable.prec = "prec", delete = F, save = F, processing = T,
                        plot.data = F, save.plots = F, scale = scaleLPJ_PFT, mode = "cf",
                        gridList = "gridlist_nextratisbona.txt")

designLPJ <- getDesign(scaleLPJ_PFT, list = T)
designLPJ$run_vegmode <- "cohort"
designLPJ$run_ifcentury <- 0
designLPJ$run_iffire <- 0
designLPJ$run_ifnlim <- 0
designLPJ$run_ifstochestab <- 0
designLPJ$run_ifstochmort <- 0
designLPJ$run_patcharea <- 25^2
designLPJ$run_npatch <- 1
designLPJ$run_ifdisturb <- 0
designLPJ$run_nyear_spinup <- 1
designLPJ$run_freenyears <- 0
designLPJ$run_save_state <- 0
designLPJ$run_restart <- 0

designLPJ$run_state_path <- mainDir
LPJsettings_PFT$design <- designLPJ

spp <- as.matrix(LPJparameters_PFT[grep("_include", rownames(LPJparameters_PFT)), ])
PFTs <- c("TeNE", "TeBS", "IBS", "TeBE", "C3G")
PFTsRows <- c()
for (i in PFTs) PFTsRows <- c(PFTsRows, grep(i, rownames(spp)))
spp[-PFTsRows, 1] <- 0
parameterStandard_PFT <- as.matrix(LPJparameters_PFT[-grep("_include",
                                                           rownames(LPJparameters_PFT)), ])
parameterStandard_PFT <- rbind(parameterStandard_PFT, spp)

ptm <- proc.time()
results_PFT <- runLPJ(x = mainDir, parameterList = parameterStandard_PFT,
                      typeList = typeList, settings = LPJsettings_PFT)
print(proc.time() - ptm)

# Code actually runs until this point without any problems

spp <- as.matrix(LPJparameters_SPP[grep("_include", rownames(LPJparameters_SPP)),
                                   ])
species <- c("Abi_alb", "Bet_pen", "Car_bet", "Fag_syl", "Cor_ave",
             "Fra_exc", "Que_rob", "Til_cor", "C3_gr")
speciesRows <- c()
for (i in species) speciesRows <- c(speciesRows, grep(i, rownames(spp)))
spp[-speciesRows, 1] <- 0
parameterStandard <- as.matrix(LPJparameters_SPP[-grep("_include",
                                                       rownames(LPJparameters_SPP)), ])
parameterStandard_SPP <- rbind(parameterStandard, spp)

ptm <- proc.time()
results_SPP <- runLPJ(x = LPJsetup, parameterList = parameterStandard_SPP,
                      typeList = typeList, settings = LPJsettings_SPP)
print(proc.time() - ptm)

LPJsettings <- LPJsettings_SPP

spp <- as.matrix(LPJparameters_SPP[grep("_include", rownames(LPJparameters_SPP)),])
species1 <- c("Fag_syl")
speciesRows <- c()

for (i in species1) speciesRows <- c(speciesRows, grep(i, rownames(spp)))

spp[-speciesRows,1] <- 0
parameterStandard <- as.matrix(LPJparameters_SPP[grep("_include", rownames(LPJparameters_SPP)),])
parameterStandard <- rbind(parameterStandard, spp)

ptm <- proc.time()
results_Fag <- runLPJ(x = LPJsetup, parameterList = parameterStandard,
                      typeList = typeList, settings = LPJsettings)

print(proc.time() - ptm)



###### Running the sensitivity analysis#########

typeList <- c("cmass", "lai")
calibList <- typeList
scaleLPJ <- "global"
species <- c("Fag_syl")
referenceData <- results_Fag@dataTypes
for (i in calibList) referenceData[[i]] <- referenceData[[i]][,species]


LPJsettings$file.co2 <- file.path(mainDir, "crudata", "cmip5_rec_co2_rcp4p5_1901_2100.txt")
LPJsettings$file.cru <- file.path(mainDir, "crudata", "Cruncep_1901_2015.bin")
LPJsettings$file.cru.misc <- file.path(mainDir, "crudata", "Cruncep_1901_2015misc.bin")
LPJsettings$file.ndep <- file.path(mainDir, "crudata", "GlobalNitrogenDeposition.bin")
LPJsettings$file.temp <- file.path(mainDir, "temp.nc")
LPJsettings$file.prec <- file.path(mainDir, "prec.nc")
LPJsettings$file.insol <- file.path(mainDir, "rad.nc")

# Get the parameter lists with ranges from literature

LPJparameters <- read.csv("./LPJParameters_calibrate.csv", sep = ";")
LPJparameters <- LPJparameters[LPJparameters$type == "calibrate", ]

# include only F. sylvatica in simulations

parameterStandard <- as.matrix(LPJparameters$Standard_Values_test)
rownames(parameterStandard) <- LPJparameters$name
spp <- getParameterList(scaleLPJ, list = F)
spp <- as.matrix(spp[rownames(spp)[grep("_include", rownames(spp))],
                     ])
speciesRows <- c()
for (i in species) speciesRows <- c(speciesRows, grep(i, rownames(spp)))
spp[-speciesRows, 1] <- 0
parameterStandard <- rbind(parameterStandard, spp)

# add two more parameters, which care for non-degeneracy

for (i in calibList) assign(x = paste("errpar", i, sep = "_"),
                            value = sd(referenceData[[i]], na.rm = T)/20)
for (i in ls()[grep("errpar", ls())]) {
  x <- get(i)
  names(x) <- i
  assign(i, x)
}

# add error parameters to table

for (i in ls()[grep("errpar", ls())]) {
  errp <- get(i)
  err <- as.data.frame(matrix(NA, ncol = ncol(LPJparameters)))
  colnames(err) <- colnames(LPJparameters)
  err$name <- names(errp)
  err$Standard_Values_test <- as.numeric(errp)
  err$Range_min_test <- 0.01
  err$Range_max_test <- 0.03
  LPJparameters <- rbind(LPJparameters, err)
}


# generate noisy data

referenceDataNoise <- referenceData
for (i in calibList) {
  referenceDataNoise[[i]] <- referenceData[[i]] + c(rnorm(length(referenceData[[i]]),
                                                          mean = 0, sd = get(ls()[grep(paste("errpar", i, sep = "_"),
                                                                                       ls())])))
  referenceDataNoise[[i]][which(referenceDataNoise[[i]] < 0)] <- 0
}

## performing the real OAT Analysis####


# LPJ -Setup with 8 cores for parallelization on my laptop
numCores <- 2

LPJsetup <- setupLPJParallel(numCores = numCores, clusterType = "SOCK",
                             mainDir = mainDir)

#define target function which computes variation in output w.r.t. reference data
# Take the C biomass and LIA

coordinate_list = list("grid-cell_c(10.45, 51.08)","grid-cell_c(87.25, 46.25)")

predict <- function(parms) {

  parDefault <- as.matrix(parameterStandard[!rownames(parameterStandard) %in%
                                              colnames(parms), ])
  parDefault <- do.call(rbind, replicate(nrow(parms), t(parDefault),
                                         simplify = FALSE))
  par <- cbind(parms, parDefault)

  x2 <- par[, !colnames(par) %in% colnames(par)[grep("errpar",
                                                     colnames(par))]]


  outputLPJ <- runLPJ(x = LPJsetup, parameterList = x2, typeList = typeList,
                      settings = LPJsettings)
  likelihoods <- c()
  nruns <- nrow(x2)
  for (nparticles in 1:nruns) {
    output <- outputLPJ[[nparticles]]@dataTypes
    output <- output[calibList]
    print(output)
    sumll <- NA
    for (grid_cells in coordinate_list){
      predicted_cell <- as.data.frame(subset(output[[grid_cells]]))
      #ToDo the same with the observed fuck
      for (out in calibList){
      predicted <- as.data.frame(subset(predicted_cell[out][,
                                                               species]))
      observed <- as.data.frame(subset(referenceDataNoise[[out]]))
      print(predicted)
      print(observe)
      llobj <- merge(observed, predicted, all = FALSE,
                     by = "row.names")[, -1]
      res <- llobj[, 2] - llobj[, 1]
      err <- as.numeric(par[nparticles, grep(out, colnames(par))])
      x <- dnorm(res, sd = err, log = T)
      sumll <- sum(sumll, x, na.rm = TRUE)
      }
    }
    likelihoods[nparticles] <- sumll
  }
  result <- vector("list", length = length(calibList) + 1)
  names(result) <- c("response", calibList)
  result[["response"]] <- likelihoods

  for (i in calibList) {
    values <- lapply(outputLPJ, function(x) {
      value <- x[i][, species]
      return(value)
    })
    values <- Reduce(cbind, values)
    values <- t(data.matrix(values))
    result[[i]] <- values
  }
  return(result)
}

# utility function for plotting the results

plotResponse <- function(response, parameters, responseName,
                         ylim, ylab, xlab, cex.main) {
  parameterName <- unique(response[, "parameter"])
  plot(response[, c(parameterName, responseName)], main = parameterName,
       cex.main = cex.main, type = "l", col = "red", xlab = xlab,
       ylab = ylab, ylim = ylim)
}

# run the sensitivity analysis

parameterNames <- LPJparameters$name
parameterData <- vector("list", length(parameterNames))
names(parameterData) <- parameterNames

steps <- numCores
for (i in 1:length(parameterNames)) {

  parameterTest <- as.matrix(seq(LPJparameters[LPJparameters$name ==
                                                parameterNames[i], "Range_min_test"], LPJparameters[LPJparameters$name ==
                                                parameterNames[i], "Range_max_test"], length.out = steps))
  colnames(parameterTest) <- parameterNames[i]

  response <- predict(parameterTest)
  # if ("try-error" %in% class(response)) {
  #   cat("Error!!")
  #   cat("\n")
  #   response <- list(response = NA)
  # }

  parameter <- as.character(parameterNames[i])
  parameterData[[i]] <- cbind(parameter, parameterTest, response$response)
  colnames(parameterData[[i]])[3] <- "response"
}


rangeResponse <- range(unlist(lapply(parameterData, function(x) {
  rangeResponse <- range(as.numeric(as.character(x[, c("response")])),
                         na.rm = T)
})))

## Plot the provided results ##

par(mfrow = c(10, 6), las = 0, mar = c(3, 2, 2, 1), oma = c(2,
                                                            4, 0, 0))

for (i in 1:length(parameterData)) {
  plotResponse(response = parameterData[[i]], parameters = LPJparameters,
               responseName = "response", ylim = rangeResponse, ylab = "",
               xlab = "", cex.main = 0.9)
}
mtext("log-likelihood", side = 2, line = 1, outer = T)
mtext(paste("Parameter values (n =", steps, ")", sep = ""), side = 1,
      line = 1, outer = T)


###### Morris Sensitivity analysis ##########

library(sensitivity)

predict <- function(parms) {
  parDefault <- as.matrix(parameterStandard[!rownames(parameterStandard) %in%
                                              colnames(parms), ])
  parDefault <- do.call(rbind, replicate(nrow(parms), t(parDefault),
                                         simplify = FALSE))
  par <- cbind(parms, parDefault)

  x2 <- par[, !colnames(par) %in% colnames(par)[grep("errpar",
                                                     colnames(par))]]

  outputLPJ <- runLPJ(x = LPJsetup, parameterList = x2, typeList = typeList,
                      settings = LPJsettings)
  likelihoods <- c()
  nruns <- nrow(x2)
  pb <- txtProgressBar(min = 0, max = nruns, style = 3)
  for (nparticles in 1:nruns) {
    output <- outputLPJ[[nparticles]]@dataTypes
    output <- output[calibList]
    sumll <- NA
    for (out in calibList) {
      predicted <- as.data.frame(subset(output[[out]][,
                                                      species]))
      observed <- as.data.frame(subset(referenceDataNoise[[out]]))
      llobj <- merge(observed, predicted, all = FALSE,
                     by = "row.names")[, -1]
      res <- llobj[, 2] - llobj[, 1]
      err <- as.numeric(par[nparticles, grep(out, colnames(par))])
      x <- dnorm(res, sd = err, log = T)
      sumll <- sum(sumll, x, na.rm = TRUE)
    }
    likelihoods[nparticles] <- sumll
  }
  likelihoods <- data.matrix(likelihoods)
  return(likelihoods)
}


levels <- 20
grid.jump <- round(levels/2)
N <- 2

factors <- as.character(LPJparameters$name)

binf <- as.numeric(LPJparameters$Range_min_test)
bsup <- as.numeric(LPJparameters$Range_max_test)

design <- list(type = "oat", levels = levels, grid.jump = grid.jump)

morrisOut <- morris(model = predict, factors = factors, r = N,
                    design = design, binf, bsup, scale = TRUE)

morrisOut

mu <- apply(morrisOut$ee, 3, function(M) {
  apply(M, 2, mean)
})
mu.star <- apply(abs(morrisOut$ee), 3, function(M) {
  apply(M, 2, mean)
})
sigma <- apply(morrisOut$ee, 3, function(M) {
  apply(M, 2, sd)
})

df <- cbind(mu.star, sigma)
colnames(df) <- c("mu.star", "sigma")
df <- df[order(df[, 1]), ]


par(mar = c(8, 4, 1, 1), cex.axis = 0.7)
M <- t(df)
pdf("./morris.pdf")
barplot(M, beside = T, col = c("black", "red"), border = c("black",
                                                           "red"), las = 2, names.arg = colnames(M))

legend("topleft", c(expression(sigma), expression(mu * "*")),
       pch = 15, col = c("red", "black"), bty = "n", cex = 1)

dev.off()
