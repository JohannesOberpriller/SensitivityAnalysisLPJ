library(purrr)
#function to resclae the parameters to 0 to 1
rescaling <- function(x){(x-min(x))/(max(x)-min(x))}

rescaling_abs <- function(x){abs(x)/max(abs(x))}

## getting the ordered parameter names ##

parameters = readRDS("Parameter_list.rds")
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



### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects_lin.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects_lin.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects_lin.rds")


#### Analysis for Pic abi ####


mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)

interaction_matrix_Pic_abi_cpool = effects_Pic_abi[["cpool_interactions"]][["complete"]]
interaction_matrix_Pic_abi_cpool = interaction_matrix_Pic_abi_cpool[growth_sites_Pic_abi]

interaction_matrix_Pic_abi_cpool_rel = Reduce("+", interaction_matrix_Pic_abi_cpool) / length(interaction_matrix_Pic_abi_cpool)
interaction_matrix_Pic_abi_cpool_abs = Reduce("+", map(interaction_matrix_Pic_abi_cpool, abs)) / length(interaction_matrix_Pic_abi_cpool)




interaction_matrix_Pic_abi_cflux = effects_Pic_abi[["cflux_interactions"]][["complete"]]
interaction_matrix_Pic_abi_cflux = interaction_matrix_Pic_abi_cflux[growth_sites_Pic_abi]

interaction_matrix_Pic_abi_cflux_rel = Reduce("+", interaction_matrix_Pic_abi_cflux) / length(interaction_matrix_Pic_abi_cflux)
interaction_matrix_Pic_abi_cflux_abs = Reduce("+", map(interaction_matrix_Pic_abi_cflux, abs)) / length(interaction_matrix_Pic_abi_cflux)


interaction_matrix_Pic_abi_agpp = effects_Pic_abi[["agpp_interactions"]][["complete"]]
interaction_matrix_Pic_abi_agpp = interaction_matrix_Pic_abi_agpp[growth_sites_Pic_abi]

interaction_matrix_Pic_abi_agpp_rel = Reduce("+", interaction_matrix_Pic_abi_agpp) / length(interaction_matrix_Pic_abi_agpp)
interaction_matrix_Pic_abi_agpp_abs = Reduce("+", map(interaction_matrix_Pic_abi_agpp, abs)) / length(interaction_matrix_Pic_abi_agpp)


#### Analysis for Fag syl ####


mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
growth_sites_Fag_syl = which(mean_effects_Fag_syl_cpool[,39]/200 > 2)

interaction_matrix_Fag_syl_cpool = effects_Fag_syl[["cpool_interactions"]][["complete"]]
interaction_matrix_Fag_syl_cpool = interaction_matrix_Fag_syl_cpool[growth_sites_Fag_syl]

interaction_matrix_Fag_syl_cpool_rel = Reduce("+", interaction_matrix_Fag_syl_cpool) / length(interaction_matrix_Fag_syl_cpool)
interaction_matrix_Fag_syl_cpool_abs = Reduce("+", map(interaction_matrix_Fag_syl_cpool, abs)) / length(interaction_matrix_Fag_syl_cpool)




interaction_matrix_Fag_syl_cflux = effects_Fag_syl[["cflux_interactions"]][["complete"]]
interaction_matrix_Fag_syl_cflux = interaction_matrix_Fag_syl_cflux[growth_sites_Fag_syl]

interaction_matrix_Fag_syl_cflux_rel = Reduce("+", interaction_matrix_Fag_syl_cflux) / length(interaction_matrix_Fag_syl_cflux)
interaction_matrix_Fag_syl_cflux_abs = Reduce("+", map(interaction_matrix_Fag_syl_cflux, abs)) / length(interaction_matrix_Fag_syl_cflux)


interaction_matrix_Fag_syl_agpp = effects_Fag_syl[["agpp_interactions"]][["complete"]]
interaction_matrix_Fag_syl_agpp = interaction_matrix_Fag_syl_agpp[growth_sites_Fag_syl]

interaction_matrix_Fag_syl_agpp_rel = Reduce("+", interaction_matrix_Fag_syl_agpp) / length(interaction_matrix_Fag_syl_agpp)
interaction_matrix_Fag_syl_agpp_abs = Reduce("+", map(interaction_matrix_Fag_syl_agpp, abs)) / length(interaction_matrix_Fag_syl_agpp)


#### Analysis for Pin syl ####


mean_effects_Pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
growth_sites_Pin_syl = which(mean_effects_Pin_syl_cpool[,39]/200 > 2)

interaction_matrix_Pin_syl_cpool = effects_Pin_syl[["cpool_interactions"]][["complete"]]
interaction_matrix_Pin_syl_cpool = interaction_matrix_Pin_syl_cpool[growth_sites_Pin_syl]

interaction_matrix_Pin_syl_cpool_rel = Reduce("+", interaction_matrix_Pin_syl_cpool) / length(interaction_matrix_Pin_syl_cpool)
interaction_matrix_Pin_syl_cpool_abs = Reduce("+", map(interaction_matrix_Pin_syl_cpool, abs)) / length(interaction_matrix_Pin_syl_cpool)




interaction_matrix_Pin_syl_cflux = effects_Pin_syl[["cflux_interactions"]][["complete"]]
interaction_matrix_Pin_syl_cflux = interaction_matrix_Pin_syl_cflux[growth_sites_Pin_syl]

interaction_matrix_Pin_syl_cflux_rel = Reduce("+", interaction_matrix_Pin_syl_cflux) / length(interaction_matrix_Pin_syl_cflux)
interaction_matrix_Pin_syl_cflux_abs = Reduce("+", map(interaction_matrix_Pin_syl_cflux, abs)) / length(interaction_matrix_Pin_syl_cflux)


interaction_matrix_Pin_syl_agpp = effects_Pin_syl[["agpp_interactions"]][["complete"]]
interaction_matrix_Pin_syl_agpp = interaction_matrix_Pin_syl_agpp[growth_sites_Pin_syl]

interaction_matrix_Pin_syl_agpp_rel = Reduce("+", interaction_matrix_Pin_syl_agpp) / length(interaction_matrix_Pin_syl_agpp)
interaction_matrix_Pin_syl_agpp_abs = Reduce("+", map(interaction_matrix_Pin_syl_agpp, abs)) / length(interaction_matrix_Pin_syl_agpp)


parameternames_ordered_pic_abi2 = substring(parameternames_ordered_pic_abi,regexpr("_", parameternames_ordered_pic_abi) + 1)
parameternames_ordered_fag_syl2 = substring(parameternames_ordered_fag_syl,regexpr("_", parameternames_ordered_fag_syl) + 1)
parameternames_ordered_pin_syl2 = substring(parameternames_ordered_pin_syl,regexpr("_", parameternames_ordered_pin_syl) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  parameternames_ordered_pic_abi2 = gsub(noisy_things[i],"",parameternames_ordered_pic_abi2)
  parameternames_ordered_fag_syl2 = gsub(noisy_things[i],"",parameternames_ordered_fag_syl2)
  parameternames_ordered_pin_syl2 = gsub(noisy_things[i],"",parameternames_ordered_pin_syl2)
}

ordering_fag_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_fag_syl2)
ordering_pin_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_pin_syl2)

## reordering of effects such that all are the same ##

interaction_matrix_Pin_syl_agpp_rel = interaction_matrix_Pin_syl_agpp_rel[ordering_pin_syl,ordering_pin_syl]/max(abs(interaction_matrix_Pin_syl_agpp_rel), na.rm =T)
interaction_matrix_Pin_syl_cflux_rel = interaction_matrix_Pin_syl_cflux_rel[ordering_pin_syl,ordering_pin_syl]/max(abs(interaction_matrix_Pin_syl_cflux_rel), na.rm =T)
interaction_matrix_Pin_syl_cpool_rel = interaction_matrix_Pin_syl_cpool_rel[ordering_pin_syl,ordering_pin_syl]/max(abs(interaction_matrix_Pin_syl_cpool_rel), na.rm =T)

interaction_matrix_Fag_syl_agpp_rel = interaction_matrix_Fag_syl_agpp_rel[ordering_fag_syl,ordering_fag_syl]/max(abs(interaction_matrix_Fag_syl_agpp_rel), na.rm =T)
interaction_matrix_Fag_syl_cflux_rel = interaction_matrix_Fag_syl_cflux_rel[ordering_fag_syl,ordering_fag_syl]/max(abs(interaction_matrix_Fag_syl_cflux_rel), na.rm =T)
interaction_matrix_Fag_syl_cpool_rel = interaction_matrix_Fag_syl_cpool_rel[ordering_fag_syl,ordering_fag_syl]/max(abs(interaction_matrix_Fag_syl_cpool_rel), na.rm =T)

interaction_matrix_Pic_abi_agpp_rel = interaction_matrix_Pic_abi_agpp_rel/max(abs(interaction_matrix_Pic_abi_agpp_rel), na.rm =T)
interaction_matrix_Pic_abi_cflux_rel = interaction_matrix_Pic_abi_cflux_rel/max(abs(interaction_matrix_Pic_abi_cflux_rel), na.rm =T)
interaction_matrix_Pic_abi_cpool_rel = interaction_matrix_Pic_abi_cpool_rel/max(abs(interaction_matrix_Pic_abi_cpool_rel), na.rm =T)


interaction_matrix_Pin_syl_agpp_abs = interaction_matrix_Pin_syl_agpp_abs[ordering_pin_syl,ordering_pin_syl]/max(abs(interaction_matrix_Pin_syl_agpp_abs), na.rm =T)
interaction_matrix_Pin_syl_cflux_abs = interaction_matrix_Pin_syl_cflux_abs[ordering_pin_syl,ordering_pin_syl]/max(abs(interaction_matrix_Pin_syl_cflux_abs), na.rm =T)
interaction_matrix_Pin_syl_cpool_abs = interaction_matrix_Pin_syl_cpool_abs[ordering_pin_syl,ordering_pin_syl]/max(abs(interaction_matrix_Pin_syl_cpool_abs), na.rm =T)

interaction_matrix_Fag_syl_agpp_abs = interaction_matrix_Fag_syl_agpp_abs[ordering_fag_syl,ordering_fag_syl]/max(abs(interaction_matrix_Fag_syl_agpp_abs), na.rm =T)
interaction_matrix_Fag_syl_cflux_abs = interaction_matrix_Fag_syl_cflux_abs[ordering_fag_syl,ordering_fag_syl]/max(abs(interaction_matrix_Fag_syl_cflux_abs), na.rm =T)
interaction_matrix_Fag_syl_cpool_abs = interaction_matrix_Fag_syl_cpool_abs[ordering_fag_syl,ordering_fag_syl]/max(abs(interaction_matrix_Fag_syl_cpool_abs), na.rm =T)


interaction_matrix_Pic_abi_agpp_abs = interaction_matrix_Pic_abi_agpp_abs/max(abs(interaction_matrix_Pic_abi_agpp_abs), na.rm =T)
interaction_matrix_Pic_abi_cflux_abs = interaction_matrix_Pic_abi_cflux_abs/max(abs(interaction_matrix_Pic_abi_cflux_abs), na.rm =T)
interaction_matrix_Pic_abi_cpool_abs = interaction_matrix_Pic_abi_cpool_abs/max(abs(interaction_matrix_Pic_abi_cpool_abs), na.rm =T)

parameters = readRDS("Parameter_list.rds")
drivers = c("co2","ndep","insol","temp","ph","prec")

variablenames = c(as.character(parameters$Pic_abi[,"NameRLPJ"]),drivers)
variablenames2 = substring(variablenames,regexpr("_", variablenames) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  variablenames2 = gsub(noisy_things[i],"",variablenames2)
}

grouping = read.csv("./Grouping_Fag_syl.csv", header = T, sep = ";")
grouping = c(as.character(grouping[,1]), rep("Drivers",6))


order_grouping = match(parameternames_ordered_pic_abi2,variablenames2)
grouping2 = data.frame(grouping[order_grouping])


interaction_matrix_agpp_rel = (interaction_matrix_Fag_syl_agpp_rel + interaction_matrix_Pic_abi_agpp_rel +
                             interaction_matrix_Pin_syl_agpp_rel)/3



interaction_matrix_cflux_rel = (interaction_matrix_Fag_syl_cflux_rel + interaction_matrix_Pic_abi_cflux_rel +
                             interaction_matrix_Pin_syl_cflux_rel)/3



interaction_matrix_cpool_rel = (interaction_matrix_Fag_syl_cpool_rel + interaction_matrix_Pic_abi_cpool_rel +
                             interaction_matrix_Pin_syl_cpool_rel)/3



interaction_matrix_agpp_abs = (interaction_matrix_Fag_syl_agpp_abs + interaction_matrix_Pic_abi_agpp_abs +
                                 interaction_matrix_Pin_syl_agpp_abs)/3



interaction_matrix_cflux_abs = (interaction_matrix_Fag_syl_cflux_abs + interaction_matrix_Pic_abi_cflux_abs +
                                  interaction_matrix_Pin_syl_cflux_abs)/3



interaction_matrix_cpool_abs = (interaction_matrix_Fag_syl_cpool_abs + interaction_matrix_Pic_abi_cpool_abs +
                                  interaction_matrix_Pin_syl_cpool_abs)/3





interaction_matrix_agpp_rel = interaction_matrix_agpp_rel[order(grouping2),order(grouping2)]
rowsums_agpp_rel = apply(abs(interaction_matrix_agpp_rel), 2,sum, na.rm = T)
interaction_matrix_agpp_rel[upper.tri(interaction_matrix_agpp_rel)] <- NA

interaction_matrix_cpool_rel = interaction_matrix_cpool_rel[order(grouping2),order(grouping2)]
rowsums_cpool_rel = apply(abs(interaction_matrix_cpool_rel), 2,sum, na.rm = T)
interaction_matrix_cpool_rel[upper.tri(interaction_matrix_cpool_rel)] <- NA

interaction_matrix_cflux_rel = interaction_matrix_cflux_rel[order(grouping2),order(grouping2)]
rowsums_cflux_rel = apply(abs(interaction_matrix_cflux_rel), 2,sum, na.rm = T)
interaction_matrix_cflux_rel[upper.tri(interaction_matrix_cflux_rel)] <- NA

interaction_matrix_agpp_abs = interaction_matrix_agpp_abs[order(grouping2),order(grouping2)]
rowsums_agpp_abs = apply(abs(interaction_matrix_agpp_abs), 2,sum, na.rm = T)
interaction_matrix_agpp_abs[upper.tri(interaction_matrix_agpp_abs)] <- NA

interaction_matrix_cpool_abs = interaction_matrix_cpool_abs[order(grouping2),order(grouping2)]
rowsums_cpool_abs = apply(abs(interaction_matrix_cpool_abs), 2,sum, na.rm = T)
interaction_matrix_cpool_abs[upper.tri(interaction_matrix_cpool_abs)] <- NA

interaction_matrix_cflux_abs = interaction_matrix_cflux_abs[order(grouping2),order(grouping2)]
rowsums_cflux_abs = apply(abs(interaction_matrix_cflux_abs), 2,sum, na.rm = T)
interaction_matrix_cflux_abs[upper.tri(interaction_matrix_cflux_abs)] <- NA


rownames(interaction_matrix_agpp_rel) = parameternames_ordered_pic_abi2[order(grouping2)]
colnames(interaction_matrix_agpp_rel) = parameternames_ordered_pic_abi2[order(grouping2)]
rownames(interaction_matrix_cpool_rel) = parameternames_ordered_pic_abi2[order(grouping2)]
colnames(interaction_matrix_cpool_rel) = parameternames_ordered_pic_abi2[order(grouping2)]
rownames(interaction_matrix_cflux_rel) = parameternames_ordered_pic_abi2[order(grouping2)]
colnames(interaction_matrix_cflux_rel) = parameternames_ordered_pic_abi2[order(grouping2)]

rownames(interaction_matrix_agpp_abs) = parameternames_ordered_pic_abi2[order(grouping2)]
colnames(interaction_matrix_agpp_abs) = parameternames_ordered_pic_abi2[order(grouping2)]
rownames(interaction_matrix_cpool_abs) = parameternames_ordered_pic_abi2[order(grouping2)]
colnames(interaction_matrix_cpool_abs) = parameternames_ordered_pic_abi2[order(grouping2)]
rownames(interaction_matrix_cflux_abs) = parameternames_ordered_pic_abi2[order(grouping2)]
colnames(interaction_matrix_cflux_abs) = parameternames_ordered_pic_abi2[order(grouping2)]


grouping2 = data.frame(grouping2[order(grouping2),])
rownames(grouping2) = rownames(interaction_matrix_agpp_rel)
colnames(grouping2) = "Processes"

# rowsums_agpp_abs = rowsums_agpp_abs[order(grouping2)]
# rowsums_cflux_abs = rowsums_cflux_abs[order(grouping2)]
# rowsums_cpool_abs = rowsums_cpool_abs[order(grouping2)]
#
# rowsums_agpp_rel = rowsums_agpp_rel[order(grouping2)]
# rowsums_cflux_rel = rowsums_cflux_rel[order(grouping2)]
# rowsums_cpool_rel = rowsums_cpool_rel[order(grouping2)]

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggplotify)

colorings = c("gold","burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")
names(colorings) = as.character(unique(grouping2)[,1])[order(as.character(unique(grouping2)[,1]))]
colors = list("Processes" = colorings)

breaksList = cumsum(table(grouping2))

counter = 0
parameternames_plot = as.grob(~{
for(i in 1:38){
  if(i %in% (breaksList+1)){
    counter = counter +1
  }
  grid.text(rownames(interaction_matrix_agpp_rel)[i],
              y= 0.782-i*0.0144 - counter*0.006 ,
              x = 1.,gp=gpar(fontsize=7), just = "right")
}
})
p1_bar <- as.grob(~barplot(-rowsums_agpp_rel[38:1],
        horiz = T, las = 1, xaxt = 'n', yaxt = "n",
        space = c(rep(0.1,table(grouping2)[7]),0.4,
                  rep(0.1,table(grouping2)[6]-1),0.4,
                  rep(0.1,table(grouping2)[5]-1),0.4,
                  rep(0.1,table(grouping2)[4]-1),0.4,
                  rep(0.1,table(grouping2)[3]-1),0.4,
                  rep(0.1,table(grouping2)[2]-1),0.4,
                  rep(0.1,table(grouping2)[1]-1))
        ))

p1_rel <- pheatmap(interaction_matrix_agpp_rel, annotation_row = grouping2, scale ="none",
                   main = "", annotation_colors = colors,
                   angle_col= 45, annotation_legend = T,cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                      "RdYlBu")))(1000),
                   breaks = seq(-1,1,2/1000),
                   legend_labels = c("-0.5","0","0.5","Interaction strength"), legend.position = c(1,3),
                   gaps_row = cumsum(table(grouping2)), gaps_col = cumsum(table(grouping2)),cellheight=10.9,cellwidth=10.9,
                   fontsize = 6.9, na_col = "white",show_rownames = F, border_color = 'white', legend = T)

P1_rel = grid.arrange(arrangeGrob(ggplot()+ theme_void(), p1_bar,ggplot()+ theme_void(), nrow = 3, heights= c(0.8,5.4,0.8)),
             arrangeGrob(arrangeGrob(ggplot() +theme_void(), parameternames_plot,ggplot() +theme_void(), nrow = 3, heights = c(0.2,6,0.8)),
             p1_rel[[4]],ncol =2, widths = c(1,4)), ncol = 2 , widths = c(1,4))

p2_bar <- as.grob(~barplot(-rowsums_cpool_rel[38:1],
                           horiz = T, las = 1, xaxt = 'n', yaxt = "n",
                           space = c(rep(0.1,table(grouping2)[7]),0.4,
                                     rep(0.1,table(grouping2)[6]-1),0.4,
                                     rep(0.1,table(grouping2)[5]-1),0.4,
                                     rep(0.1,table(grouping2)[4]-1),0.4,
                                     rep(0.1,table(grouping2)[3]-1),0.4,
                                     rep(0.1,table(grouping2)[2]-1),0.4,
                                     rep(0.1,table(grouping2)[1]-1))
))



p2_rel <- pheatmap(interaction_matrix_cpool_rel, annotation_row = grouping2, scale ="none",
                   main = "", annotation_colors = colors,
                   angle_col= 45, annotation_legend = T,cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                      "RdYlBu")))(1000),
                   breaks = seq(-1,1,2/1000),
                   legend_labels = c("-0.5","0","0.5","Interaction strength"), legend.position = c(1,3),
                   gaps_row = cumsum(table(grouping2)), gaps_col = cumsum(table(grouping2)),cellheight=10.9,cellwidth=10.9,
                   fontsize = 6.9, na_col = "white",show_rownames = F, border_color = 'white', legend = T)

P2_rel = grid.arrange(arrangeGrob(ggplot()+ theme_void(), p2_bar,ggplot()+ theme_void(), nrow = 3, heights= c(0.8,5.4,0.8)),
                      arrangeGrob(arrangeGrob(ggplot() +theme_void(), parameternames_plot,ggplot() +theme_void(), nrow = 3, heights = c(0.2,6,0.8)),
                                  p2_rel[[4]],ncol =2, widths = c(1,4)), ncol = 2 , widths = c(1,4))


p3_bar <- as.grob(~barplot(-rowsums_cflux_rel[38:1],
                           horiz = T, las = 1, xaxt = 'n', yaxt = "n",
                           space = c(rep(0.1,table(grouping2)[7]),0.4,
                                     rep(0.1,table(grouping2)[6]-1),0.4,
                                     rep(0.1,table(grouping2)[5]-1),0.4,
                                     rep(0.1,table(grouping2)[4]-1),0.4,
                                     rep(0.1,table(grouping2)[3]-1),0.4,
                                     rep(0.1,table(grouping2)[2]-1),0.4,
                                     rep(0.1,table(grouping2)[1]-1))
))



p3_rel <- pheatmap(interaction_matrix_cflux_rel, annotation_row = grouping2, scale ="none",
                   main = "", annotation_colors = colors,
                   angle_col= 45, annotation_legend = T,cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                      "RdYlBu")))(1000),
                   breaks = seq(-1,1,2/1000),
                   legend_labels = c("-0.5","0","0.5","Interaction strength"), legend.position = c(1,3),
                   gaps_row = cumsum(table(grouping2)), gaps_col = cumsum(table(grouping2)),cellheight=10.9,cellwidth=10.9,
                   fontsize = 6.9, na_col = "white",show_rownames = F, border_color = 'white', legend = T)

P3_rel = grid.arrange(arrangeGrob(ggplot()+ theme_void(), p3_bar,ggplot()+ theme_void(), nrow = 3, heights= c(0.8,5.4,0.8)),
                      arrangeGrob(arrangeGrob(ggplot() +theme_void(), parameternames_plot,ggplot() +theme_void(), nrow = 3, heights = c(0.2,6,0.8)),
                                  p3_rel[[4]],ncol =2, widths = c(1,4)), ncol = 2 , widths = c(1,4))



p1_bar <- as.grob(~barplot(-rowsums_agpp_abs[38:1],
                           horiz = T, las = 1, xaxt = 'n', yaxt = "n",
                           space = c(rep(0.1,table(grouping2)[7]),0.4,
                                     rep(0.1,table(grouping2)[6]-1),0.4,
                                     rep(0.1,table(grouping2)[5]-1),0.4,
                                     rep(0.1,table(grouping2)[4]-1),0.4,
                                     rep(0.1,table(grouping2)[3]-1),0.4,
                                     rep(0.1,table(grouping2)[2]-1),0.4,
                                     rep(0.1,table(grouping2)[1]-1))
))

p1_abs <- pheatmap(interaction_matrix_agpp_abs, annotation_row = grouping2, scale ="none",
                   main = "", annotation_colors = colors,
                   angle_col= 45, annotation_legend = T,cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                      "RdYlBu")))(1000),
                   breaks = seq(-1,1,2/1000),
                   legend_labels = c("-0.5","0","0.5","Interaction strength"), legend.position = c(1,3),
                   gaps_row = cumsum(table(grouping2)), gaps_col = cumsum(table(grouping2)),cellheight=10.9,cellwidth=10.9,
                   fontsize = 6.9, na_col = "white",show_rownames = F, border_color = 'white', legend = T)

P1_abs = grid.arrange(arrangeGrob(ggplot()+ theme_void(), p1_bar,ggplot()+ theme_void(), nrow = 3, heights= c(0.8,5.4,0.8)),
                      arrangeGrob(arrangeGrob(ggplot() +theme_void(), parameternames_plot,ggplot() +theme_void(), nrow = 3, heights = c(0.2,6,0.8)),
                                  p1_abs[[4]],ncol =2, widths = c(1,4)), ncol = 2 , widths = c(1,4))

p2_bar <- as.grob(~barplot(-rowsums_cpool_abs[38:1],
                           horiz = T, las = 1, xaxt = 'n', yaxt = "n",
                           space = c(rep(0.1,table(grouping2)[7]),0.4,
                                     rep(0.1,table(grouping2)[6]-1),0.4,
                                     rep(0.1,table(grouping2)[5]-1),0.4,
                                     rep(0.1,table(grouping2)[4]-1),0.4,
                                     rep(0.1,table(grouping2)[3]-1),0.4,
                                     rep(0.1,table(grouping2)[2]-1),0.4,
                                     rep(0.1,table(grouping2)[1]-1))
))



p2_abs <- pheatmap(interaction_matrix_cpool_abs, annotation_row = grouping2, scale ="none",
                   main = "", annotation_colors = colors,
                   angle_col= 45, annotation_legend = T,cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                      "RdYlBu")))(1000),
                   breaks = seq(-1,1,2/1000),
                   legend_labels = c("-0.5","0","0.5","Interaction strength"), legend.position = c(1,3),
                   gaps_row = cumsum(table(grouping2)), gaps_col = cumsum(table(grouping2)),cellheight=10.9,cellwidth=10.9,
                   fontsize = 6.9, na_col = "white",show_rownames = F, border_color = 'white', legend = T)

P2_abs = grid.arrange(arrangeGrob(ggplot()+ theme_void(), p2_bar,ggplot()+ theme_void(), nrow = 3, heights= c(0.8,5.4,0.8)),
                      arrangeGrob(arrangeGrob(ggplot() +theme_void(), parameternames_plot,ggplot() +theme_void(), nrow = 3, heights = c(0.2,6,0.8)),
                                  p2_abs[[4]],ncol =2, widths = c(1,4)), ncol = 2 , widths = c(1,4))


p3_bar <- as.grob(~barplot(-rowsums_cflux_abs[38:1],
                           horiz = T, las = 1, xaxt = 'n', yaxt = "n",
                           space = c(rep(0.1,table(grouping2)[7]),0.4,
                                     rep(0.1,table(grouping2)[6]-1),0.4,
                                     rep(0.1,table(grouping2)[5]-1),0.4,
                                     rep(0.1,table(grouping2)[4]-1),0.4,
                                     rep(0.1,table(grouping2)[3]-1),0.4,
                                     rep(0.1,table(grouping2)[2]-1),0.4,
                                     rep(0.1,table(grouping2)[1]-1))
))



p3_abs <- pheatmap(interaction_matrix_cflux_abs, annotation_row = grouping2, scale ="none",
                   main = "", annotation_colors = colors,
                   angle_col= 45, annotation_legend = T,cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                                                      "RdYlBu")))(1000),
                   breaks = seq(-1,1,2/1000),
                   legend_labels = c("-0.5","0","0.5","Interaction strength"), legend.position = c(1,3),
                   gaps_row = cumsum(table(grouping2)), gaps_col = cumsum(table(grouping2)),cellheight=10.9,cellwidth=10.9,
                   fontsize = 6.9, na_col = "white",show_rownames = F, border_color = 'white', legend = T)

P3_abs = grid.arrange(arrangeGrob(ggplot()+ theme_void(), p3_bar,ggplot()+ theme_void(), nrow = 3, heights= c(0.8,5.4,0.8)),
                      arrangeGrob(arrangeGrob(ggplot() +theme_void(), parameternames_plot,ggplot() +theme_void(), nrow = 3, heights = c(0.2,6,0.8)),
                                  p3_abs[[4]],ncol =2, widths = c(1,4)), ncol = 2 , widths = c(1,4))


