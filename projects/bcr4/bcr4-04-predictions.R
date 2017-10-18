TODO:
- use existing pred grid (intersect with AK grid)
- resample pred grid at AK grid (use slp as is from AK)
- write code for:
  - raster to table
  - table to pred
  - pred to raster

library(sp)
library(raster)

slp <- raster("e:/peter/bam/bcr4/alfresco_data/iem_prism_slope_1km.tiff")
veg <- raster("e:/peter/bam/bcr4/alfresco_data/iem_vegetation_model_input_v0_4.tiff")
lcc <- raster("e:/peter/bam/bcr4/alfresco_data/LandCover_iem_ALFRESCO_2005.tiff")
lcc[lcc==255] <- NA

veg1 <- raster("e:/peter/bam/bcr4/alfresco_data/Veg_171_2050.tif")
age1 <- raster("e:/peter/bam/bcr4/alfresco_data/Age_171_2050.tif")
age1[age1 < 0] <- 100



V4.0
0 | Not Modeled
1 | Black Spruce Forest
2 | White Spruce Forest
3 | Deciduous Forest
4 | Shrub Tundra
5 | Graminoid Tundra
6 | Wetland Tundra
7 | Barren lichen-moss
8 | Heath
9 | Maritime Upland Forest
10 | Maritime Forested Wetland
11 | Maritime Fen
12 | Maritime Alder Shrubland

v2.0
Final Legend:
0 - Not Modeled
1 - Black Spruce
2 - White Spruce
3 - Deciduous
4 - Shrub Tundra
5 - Graminoid  Tundra
6 - Wetland Tundra
7 - Barren lichen-moss
8 - Temperate Rainforest
255 - out of bounds
