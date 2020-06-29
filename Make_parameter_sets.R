sites <- readRDS("sites_data.rds")

gridlist = sites[,c("Longitudinal","Latitudinal")]
write.table(gridlist, file = "LPJrunTest/gridlist_sensi.txt",sep="\t",row.names=FALSE, col.names = F)

mastersheet = data.frame(matrix(nrow = 800, ncol = 6))

colnames(mastersheet) = c("site-id","lon","lat","Fag_syl","Pic_abi", "Pin_syl")

mastersheet[,c("site-id","lon","lat")] = sites[,c("Site","Longitudinal","Latitudinal")]

mastersheet[1:200,4:6] = matrix(rep(c(1,0,0),200),ncol = 3, nrow = 200, byrow = T )
mastersheet[201:400,4:6] = matrix(rep(c(0,1,0),200),ncol = 3, nrow = 200, byrow = T )
mastersheet[401:600,4:6] = matrix(rep(c(0,0,1),200),ncol = 3, nrow = 200, byrow = T )
mastersheet[601:800,4:6] = matrix(rep(c(1,1,1),200),ncol = 3, nrow = 200, byrow = T )

parameterranges = read.csv("./SyntheseParameterRanges.csv", header = T, sep = ",")

parameterranges_common = parameterranges[2:12,] 

parameterranges_pinsyl = parameterranges[15:36,]

parameterranges_picabi = parameterranges[38:59,]

parameterranges_fagsyl = parameterranges[61:82,]





## Technical implementation : 

# think about drivers ( how to get the ranges )
# think about how to implement them in lpj ( + - seems to be appropriate since 
# temperature is in celcius )
# how to draw the parameters: 
# - for each experiment (3xmono, 1x mixed) same parameters 
# - 




generate_R_scripts <- function(mastertemplate, )