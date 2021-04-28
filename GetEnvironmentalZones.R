library(sp)
library(sf)
library(tidyverse)
library(ggplot2)


aoi_boundary_HARV <- st_read("./EnSv8/ens_v8.shp")

sites <- readRDS("sites_data.rds")

pnts = sites[,3:4]

pnts$region <- apply(pnts, 1, function(row) {
  # transformation to palnar is required, since sf library assumes planar projection
  tt1_pl <- st_transform(aoi_boundary_HARV, 2163)
  coords <- as.data.frame(matrix(row, nrow = 1,
                                 dimnames = list("", c("x", "y"))))
  pnt_sf <- st_transform(st_sfc(st_point(row),crs = 4326), 2163)
  # st_intersects with sparse = FALSE returns a logical matrix
  # with rows corresponds to argument 1 (points) and
  # columns to argument 2 (polygons)
  tt1_pl[which(st_intersects(pnt_sf, tt1_pl, sparse = F)), ]$EnZ_name
})

out_of_scope = which(is.na((pnts$region == F)))

sites_enz = sites[!(1:200 %in% out_of_scope),]
sites_enz$VegetationZone = unlist(pnts$region)

saveRDS(sites_enz,"sites_cleaned_with_vegzone.Rds")

sites_enz <- st_as_sf(sites_enz, coords = c("Longitudinal", "Latitudinal"),
                  crs = 4326, agr = "constant")

ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.3, color = "black", aes(fill = EnZ_name)) +
  coord_sf()+
  geom_sf(data = sites_enz, colour = "black", size = 2.5)
