library(maps)
library(mapdata)
library(tess3r)
library(wesanderson)

# set working directory
setwd('/Users/Evan/Dropbox/Code/')

# load local bugfixed TESS3 code
insertSource("TESS3/R/spatialPlots.R", package = "tess3r", functions = 'PlotInterpotationMax')

# load the Sulawesi RData
load('sulawesi/Sulawesi.RData')

render_map <- function(species.data, species.colours, species.clades) {
  
  # drop any rows with missing components or locations
  species.data <- na.omit(species.data)
  
  # drop Zoo samples and those of unknown origin
  species.data <- species.data[-which(species.data$Loc_Abbrev_LF %in% c("UN", "ZO", "SU_BU", "S")),]
  
  # get the STRUCTURE components ready for tess3
  species.struct <- as.matrix(species.data[,1:(ncol(species.data)-3)])
  class(species.struct) <- c("tess3Q", "matrix")
  
  # get the coordinates of the remaining samples
  species.coord <- species.data[,c('Longitude', 'Latitude')]  
  
  # make the colour palette
  my.palette <- CreatePalette(species.colours, palette.length = 50)
  
  # plot the interpolation surface
  plot(x=species.struct, coord=species.coord,
       method = "map.max",
       interpol = FieldsTpsModel(),
       # interpol = FieldsKrigModel(10),
       map.polygon = map.polygon,
       main = NA, xlab = NA, ylab = NA, xaxt='n', yaxt='n', bty = 'n',
       resolution = resolution, cex = .4,
       window = c(xlim, ylim),
       xlim =  xlim, ylim = ylim,
       col.palette = my.palette)
  
  # render the actual map
  # plot(map.polygon, xlim = c(118.5, 127.4), ylim = c(-5.9, 1.8), add=TRUE)
  
  # add the regions, labels and mask some of the peripheral islands
  annotate_map(Sula_Lab, Sula_Seg, species.colours, species.clades)
}

annotate_map <- function(segment.labels, segment.lines, species.colours, species.clades) {
  
  # mask the edges of the adjacent islands using white polygons
  polygon(c(par()$usr[1], par()$usr[1], 119.25, 119.25), c(0, 2, 2, 0), border = "white", col = "white")
  polygon(c(120, 120, 121, 121), c(par()$usr[3], -5.7, -5.7, par()$usr[3]), border = "white", col = "white")
  polygon(c(127, 127, par()$usr[2], par()$usr[2]), c(-1.8, 2, 2, -1.8), border = "white", col = "white")
  
  # segment the island into regions (using the coordinates in segment.lines)
  with(segment.lines, segments(Long_1, Lat_1, Long_2, Lat_2, lwd = 3))
  
  # remove the UN and ZO rows from the matrix of label coordinates
  segment.labels <- segment.labels[-which(rownames(segment.labels) %in% c("UN", "ZO", "SU_BU", "S")),]
  
  # label the regions
  with(segment.labels, text(Long, Lat, labels = dimnames(segment.labels)[[1]], cex = 1.2))
  
  # add arrows for ambiguous label positions
  arrows(122.674413, -0.411239, 122.35, -0.394760, length = 0.1, code = 2, lwd = 2)
  arrows(123.649449, -1.328482, 124.001012, -1.314753, code = 1, lwd = 2, length = 0.1)
  
  # add the legend
  legend(125.75, 1, legend = species.clades, fill = alpha(species.colours, 0.75), cex = 1.2, bty = "n", title = "Population")
  # legend(125.75, 1.8, legend = species.clades, fill = species.colours, cex = 0.8, bty = "n", title = "Clade")
}

# get the custom colour palette
WesAndersonCol <- c(wes_palette("GrandBudapest")[c(1:2, 4)],
 			wes_palette("GrandBudapest2")[c(1:2, 4)],
 			wes_palette("Royal1")[c(1, 2)],
 			wes_palette("Moonrise3")[c(1, 3, 2, 5)],
 			wes_palette("Moonrise2")[c(4, 2)],
 			wes_palette("Darjeeling")[c(1)],
 			wes_palette("Royal2")[c(5)],
 			wes_palette("Chevalier")[c(3)])

# display the colour list
# pie(rep(1,length(WesAndersonCol)), col=WesAndersonCol)

# set the lat/long boundaries for the map
xlim <- c(118.5, 127.4)
ylim <- c(-5.9, 2)

# use the worldHires map
map.hires <- map("worldHires", xlim = xlim, ylim = ylim, plot=FALSE, fill=TRUE)
map.polygon <- map2SpatialPolygons(map.hires, map.hires$names)
# map.polygon <- rworldmap::getMap(resolution = "high")

# set the resolution of the interpolation surface
# resolution <- c(100,100)
# resolution <- c(300,300)
resolution <- c(1000,1000)

locat.columns <- c('Loc_Abbrev_LF', 'Longitude', 'Latitude')

# extract the STRUCTURE components and the sample locations
# anoa.data <- subset(Anoa_All_Genetics, select = c(A4_NE_WC:A1_BT, Loc_Abbrev_LF, Longitude, Latitude))
# baby.data <- subset(Baby_All_Genetics, select = c(B3_SU_BU:B2_SE, Loc_Abbrev_LF, Longitude, Latitude))
# susc.data <- subset(Sus_cel_All_Genetics, select = c(S5_WC_SW:S3, Loc_Abbrev_LF, Longitude, Latitude))

anoa.data <- Anoa_All_Genetics[,c('A1_BT', 'A2_NW', 'A3_SE', 'A4_NE_WC', 'A5_NC_EC', locat.columns)]
baby.data <- Baby_All_Genetics[,c('B1_WC_NW', 'B2_SE', 'B3_SU_BU', 'B4_NE', 'B5_TO', locat.columns)]
susc.data <- Sus_cel_All_Genetics[,c('S1_NW', 'S2_PE', 'S3', 'S4_SE_BT', 'S5_WC_SW', 'S7_EC', locat.columns)]

# setup the colour lists
anoa.colours <- WesAndersonCol[c(3, 6, 8, 2, 5)]
baby.colours <- WesAndersonCol[c(6, 1, 2, 10, 9)]
susc.colours <- WesAndersonCol[c(1, 2, 3, 4, 5, 6)]

# set the clade names
anoa.clades <- c('A1','A2','A3','A4','A5')
baby.clades <- c('B1','B2','B3','B4','B5','B6')
susc.clades <- c('S1','S2','S3','S4','S5')

pdf(file = "sulawesi/pdf/Sula-Maps.pdf", width = (xlim[2]-xlim[1])*3*.8, height = (ylim[2]-ylim[1])*2*.8)

# display the maps in a 2x3 grid
par(mfrow = c(2,3), mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1), bty = 'n')

# STRUCTURE maps
render_map(anoa.data, anoa.colours, anoa.clades)
render_map(baby.data, baby.colours, baby.clades)
render_map(susc.data, susc.colours, susc.clades)

# TODO add the clade maps

dev.off()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# library(maps)
# 
# par(mfrow = c(1,1), mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1), bty = 'n')
# 
# # print the map with the admix components
# plot(species.struct, species.coord,
#      method = "map.max", 
#      # background = FALSE,
#      # interpol = FieldsKrigModel(10),
#      interpol = FieldsTpsModel(), 
#      main = NA, xlab = NA, ylab = NA, xaxt='n', yaxt='n', bty = 'n',
#      resolution = resolution, cex = .4,
#      window = c(118.5, 127.4, -5.9, 1.8),
#      xlim = xlim, ylim = ylim,
#      col.palette = my.palette)

# -----------------------------------------------------------------------------

# asc.raster <- tempfile()
# download.file("http://membres-timc.imag.fr/Olivier.Francois/RasterMaps/Europe.asc", asc.raster)
# 
# plot(species.struct, species.coord,
#      method = "map.max",
#      interpol = FieldsKrigModel(10),
#      raster.filename = asc.raster,
#      main = NA, xlab = NA, ylab = NA, xaxt='n', yaxt='n', bty = 'n',
#      resolution = resolution, cex = .4,
#      window = c(118.5, 126, -5.9, 1.8),
#      xlim = xlim, ylim = ylim,
#      col.palette = my.palette)


# -----------------------------------------------------------------------------

# plot(species.struct, species.coord,
#      method = "map.max", 
#      interpol = FieldsTpsModel(),
#      main = NA, xlab = NA, ylab = NA, xaxt='n', yaxt='n', bty = 'n',
#      resolution = resolution, cex = .4,
#      xlim = xlim, ylim = ylim,
#      col.palette = my.palette)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

## ggplot2 version...

# library(ggplot2)
# library(maps)
# library(mapdata)
# library(maptools)
# 
# # get the worldHires map and convert to SpatialPolygons
# map.hires <- map("worldHires", xlim = xlim, ylim = ylim, plot=FALSE, fill=TRUE)
# map.polygon <- map2SpatialPolygons(map.hires, map.hires$names)
# 
# pl <- ggtess3Q(species.struct, species.coord,
#          map.polygon = map.polygon, 
#          col.palette = my.palette,
#          window = c(118.5, 127.4, -5.9, 1.8),
#          interpolation.model = FieldsKrigModel(1),
#          resolution = resolution)
# 
# pl + geom_path(data = map.polygon, 
#             aes(x = long, y = lat, group = group), size = 0.4) +
#   xlim(118.5, 127.4) + ylim(-5.9, 1.8) +
#   coord_equal() +
#   geom_point(data = species.coord, aes(x = Longitude, y = Latitude), size = 0.4) +
#   theme_classic() +
#   theme(
#     # turn off the axis
#     axis.line=element_blank(),
#     axis.text.x=element_blank(),
#     axis.text.y=element_blank(),
#     axis.ticks=element_blank(),
#     axis.title.x=element_blank(),
#     axis.title.y=element_blank()
#   )

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 
# require('dplyr')
# 
# # make the matrices that we need for the maps
# 
# # extract the data frame
# Anoa_HapRegion <- as.data.frame.matrix(with(Anoa_All_Genetics, table(Loc_Abbrev_LF, LF_Anoa_Clust)))
# 
# # add a row sum column AND add a Loc_Abbrev_LF column (so we can join the lat/long data)
# Anoa_HapRegion <- Anoa_HapRegion %>%
#   dplyr::mutate(Sum = rowSums(Anoa_HapRegion)) %>%
#   dplyr::mutate(Loc_Abbrev_LF = dimnames(Anoa_HapRegion)[[1]])
# 
# # join the lat/long data using Loc_Abbrev_LF as the key
# Anoa_HaploLF_Region <- dplyr::left_join(Anoa_HapRegion, Sula_Reg, by = "Loc_Abbrev_LF")
# 
# # add a row sum column AND add a Loc_Abbrev_LF column (so we can join the lat/long data)
# Anoa_HapRegion <- Anoa_HapRegion %>%
# 	dplyr::mutate(Sum = rowSums(Anoa_HapRegion)) %>%
# 	dplyr::mutate(Loc_Abbrev_LF = dimnames(Anoa_HapRegion)[[1]])
# 
# # join the lat/long data using Loc_Abbrev_LF as the key
# Anoa_HaploLF_Region <- dplyr::left_join(Anoa_HapRegion, Sula_Reg, by = "Loc_Abbrev_LF")
# 
# # add % columns for each type (remove rows with NaN in the Prop_A4a_WC column)
# Anoa_HaploLF_Region <- Anoa_HaploLF_Region %>%
# 	dplyr::mutate(Prop_A4a_WC = round(A4a_WC / Sum, 3)) %>%
# 	dplyr::mutate(Prop_A1_NE = round(A1_NE / Sum, 3)) %>%
# 	dplyr::mutate(Prop_A5b_SE = round(A5b_SE / Sum, 3)) %>%
# 	dplyr::mutate(Prop_A2_NW = round(A2_NW / Sum, 3)) %>%
# 	dplyr::mutate(Prop_A3_EC = round(A3_EC / Sum, 3)) %>%
# 	dplyr::mutate(Prop_A5a_BT = round(A5a_BT / Sum, 3)) %>%
# 	dplyr::mutate(Prop_A4b = round(A4b / Sum, 3)) %>%
# 	dplyr::filter(!is.na(Prop_A4a_WC))
# 
# # remove the Zoo row and the A4b columns
# Anoa_HaploLF_Region_NoZO <- Anoa_HaploLF_Region %>%
# 	dplyr::filter(Loc_Abbrev_LF != "ZO") %>%
# 	dplyr::select(A4a_WC:A5a_BT, Sum:Prop_A5a_BT)
# 
# require('mapdata')
# require('maps')
# require('scales')
# require('plotrix')
# 
# ## mtDNA map
# 
# # pdf(file = "Anoa\\Anoa_Map_mtDNA_haplogroups_NoZO.pdf", width = 11.69, height = 8.27)
# 
# # par(mfrow = c(1,1), mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1))
# 
# # map of Sulawesi (using lat/long coordinates)
# # map("worldHires", xlim = c(118.5, 127.4), ylim = c(-5.9, 1.8))
# 
# 
# for(i in 1:nrow(Anoa_HaploLF_Region_NoZO)) {
#   # get all the non-zero % components
# 	xxx <-with(Anoa_HaploLF_Region_NoZO,
# 		c(Prop_A4a_WC[i], Prop_A1_NE[i], Prop_A5b_SE[i], Prop_A2_NW[i], Prop_A3_EC[i], Prop_A5a_BT[i]))[which(Anoa_HaploLF_Region_NoZO[i, 1:6] > 0)]
# 
# 	# add a pie chart of the % components (preserving colour order for each component) using a 0.75 transparency
# 	with(Anoa_HaploLF_Region_NoZO, floating.pie(Long_Reg_LF[i], Lat_Reg_LF[i], x = xxx, col =
# 		alpha(WesAndersonCol[c(2, 3, 5, 6, 8, 7)][which(Anoa_HaploLF_Region_NoZO[i, 1:6] > 0)], 0.99), radius = 0.4, border = NA))
# }
# 
# # add the count of samples ontop of each pie chart
# with(Anoa_HaploLF_Region_NoZO, text(Long_Reg_LF, Lat_Reg_LF, Sum))
# 
# legend(125.75, 1, legend = c("A1", "A2", "A3", "A4", "A5a", "A5b"),
# 	fill = alpha(WesAndersonCol[c(3, 6, 8, 2, 7, 5)], 0.75), cex = 1.2, bty = "n", title = "Clade")
# 
# dev.off()



