## need to set your own working directory
setwd('/Users/Evan/Dropbox/Documents/Oxford/DPhil/Sulawesi/')

# load the Sulawesi RData
load('Sulawesi.RData')

# get the custom colour palette
require('wesanderson')
WesAndersonCol <- c(wes_palette("GrandBudapest")[c(1:2, 4)],
 			wes_palette("GrandBudapest2")[c(1:2, 4)],
 			wes_palette("Royal1")[c(1, 2)],
 			wes_palette("Moonrise3")[c(1, 3, 2, 5)],
 			wes_palette("Moonrise2")[c(4, 2)],
 			wes_palette("Darjeeling")[c(1)],
 			wes_palette("Royal2")[c(5)],
 			wes_palette("Chevalier")[c(3)])

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# TODO which STRUCTURE fields to use? Az4_NE_WC:Az3_BT/ A4_NE_WC:A1_BT

# extract the STRUCTURE components and the sample locations
Anoa.data <- subset(Anoa_All_Genetics, select = c(Az4_NE_WC:Az3_BT, Loc_Abbrev_LF, Longitude, Latitude))

# drop any rows with missing components or locations
Anoa.data <- na.omit(Anoa.data)

# drop Zoo samples
Anoa.data <- Anoa.data[Anoa.data$Loc_Abbrev_LF != "ZO",]

# get the STRUCTURE components ready for tess3
Anoa.struct <- as.matrix(Anoa.data[,1:5])
class(Anoa.struct) <- c("tess3Q", "matrix")

# get the coordinates of the remaining samples
Anoa.coord <- Anoa.data[,c('Longitude', 'Latitude')]

library(tess3r)

# TODO fix this colour palette
my.palette <- CreatePalette(WesAndersonCol[c(2, 3, 5, 6, 8, 7)], palette.length = 20)

xlim <- c(118.5, 127.4)
ylim <- c(-5.9, 1.8)
resolution <- c(400,400)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

library(maps)

par(mfrow = c(1,1), mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1), bty = 'n')

# print the map with the admix components
plot(Anoa.struct, Anoa.coord,
     method = "map.max", 
     interpol = FieldsKrigModel(10),
     main = NA, xlab = NA, ylab = NA, xaxt='n', yaxt='n', bty = 'n',
     resolution = resolution, cex = .4,
     xlim = xlim, ylim = ylim,
     col.palette = my.palette)

# plot(Anoa.struct, Anoa.coord,
#      method = "map.max", 
#      interpol = FieldsTpsModel(),
#      main = NA, xlab = NA, ylab = NA, xaxt='n', yaxt='n', bty = 'n',
#      resolution = resolution, cex = .4,
#      xlim = xlim, ylim = ylim,
#      col.palette = my.palette)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

library(ggplot2)
library(maps)
library(mapdata)
library(maptools)

# get the worldHires map and convert to SpatialPolygons
map.hires <- map("worldHires", xlim = xlim, ylim = ylim, plot=FALSE, fill=TRUE)
map.polygon <- map2SpatialPolygons(map.hires, map.hires$names)

pl <- ggtess3Q(Anoa.struct, Anoa.coord,
         map.polygon = map.polygon, 
         col.palette = my.palette,
         interpolation.model = FieldsTpsModel(),
         resolution = resolution)

pl + geom_path(data = map.polygon, 
            aes(x = long, y = lat, group = group), size = 0.4) +
  xlim(xlim) + ylim(ylim) +
  coord_equal() +
  geom_point(data = Anoa.coord, aes(x = Longitude, y = Latitude), size = 0.5) +
  theme_classic() +
  theme(
    # turn off the axis
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  )

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
require('dplyr')

# make the matrices that we need for the maps

# extract the data frame
Anoa_HapRegion <- as.data.frame.matrix(with(Anoa_All_Genetics, table(Loc_Abbrev_LF, LF_Anoa_Clust)))

# add a row sum column AND add a Loc_Abbrev_LF column (so we can join the lat/long data)
Anoa_HapRegion <- Anoa_HapRegion %>%
  dplyr::mutate(Sum = rowSums(Anoa_HapRegion)) %>%
  dplyr::mutate(Loc_Abbrev_LF = dimnames(Anoa_HapRegion)[[1]])

# join the lat/long data using Loc_Abbrev_LF as the key
Anoa_HaploLF_Region <- dplyr::left_join(Anoa_HapRegion, Sula_Reg, by = "Loc_Abbrev_LF")

# add a row sum column AND add a Loc_Abbrev_LF column (so we can join the lat/long data)
Anoa_HapRegion <- Anoa_HapRegion %>%
	dplyr::mutate(Sum = rowSums(Anoa_HapRegion)) %>%
	dplyr::mutate(Loc_Abbrev_LF = dimnames(Anoa_HapRegion)[[1]])

# join the lat/long data using Loc_Abbrev_LF as the key
Anoa_HaploLF_Region <- dplyr::left_join(Anoa_HapRegion, Sula_Reg, by = "Loc_Abbrev_LF")

# add % columns for each type (remove rows with NaN in the Prop_A4a_WC column)
Anoa_HaploLF_Region <- Anoa_HaploLF_Region %>%
	dplyr::mutate(Prop_A4a_WC = round(A4a_WC / Sum, 3)) %>%
	dplyr::mutate(Prop_A1_NE = round(A1_NE / Sum, 3)) %>%
	dplyr::mutate(Prop_A5b_SE = round(A5b_SE / Sum, 3)) %>%
	dplyr::mutate(Prop_A2_NW = round(A2_NW / Sum, 3)) %>%
	dplyr::mutate(Prop_A3_EC = round(A3_EC / Sum, 3)) %>%
	dplyr::mutate(Prop_A5a_BT = round(A5a_BT / Sum, 3)) %>%
	dplyr::mutate(Prop_A4b = round(A4b / Sum, 3)) %>%
	dplyr::filter(!is.na(Prop_A4a_WC))

# remove the Zoo row and the A4b columns
Anoa_HaploLF_Region_NoZO <- Anoa_HaploLF_Region %>%
	dplyr::filter(Loc_Abbrev_LF != "ZO") %>%
	dplyr::select(A4a_WC:A5a_BT, Sum:Prop_A5a_BT)

# remove the UN and ZO rows from the matrix of label coordinates
Sula_Lab_NoUNZO <- Sula_Lab[-which(rownames(Sula_Lab) %in% c("UN", "ZO")),]

#Sula_Lab_NoUNZO <- Sula_Lab %>%
#  dplyr::filter(!(rownames(Sula_Lab) %in% c("UN", "ZO")))
  


require('mapdata')
require('maps')
require('scales')
require('plotrix')

## mtDNA map

pdf(file = "Anoa\\Anoa_Map_mtDNA_haplogroups_NoZO.pdf", width = 11.69, height = 8.27)

par(mfrow = c(1,1), mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1))

# map of Sulawesi (using lat/long coordinates)
map("worldHires", xlim = c(118.5, 127.4), ylim = c(-5.9, 1.8))

# mask the edges of the adjacent islands using white polygons
polygon(c(par()$usr[1], par()$usr[1], 119.25, 119.25), c(0, 2, 2, 0), 
	border = "white", col = "white")
polygon(c(120, 120, 121, 121), 
	c(par()$usr[3], -5.7, -5.7, par()$usr[3]), 
	border = "white", col = "white")
polygon(c(127, 127, par()$usr[2], par()$usr[2]), c(-1.8, 2, 2, -1.8), 
	border = "white", col = "white")

# segment the island into regions (using the coordinates in Sula_Seg)
with(Sula_Seg, segments(Long_1, Lat_1, Long_2, Lat_2, lwd = 3))

# label the regions
with(Sula_Lab_NoUNZO, text(Long, Lat, labels = dimnames(Sula_Lab_NoUNZO)[[1]], cex = 0.9))

# add arrows for ambiguous label positions
arrows(122.674413, -0.411239, 122.35, -0.394760, length = 0.1, code = 2, lwd = 2)
arrows(123.649449, -1.328482, 124.001012, -1.314753, code = 1, lwd = 2, length = 0.1)

# partially mask one of the labels (why not exclude before rendering?)
polygon(c(125, 125, 125.8, 125.8), c(-3.5, -3.1, -3.1, -3.5), 
	border = "white", col = "white")

for(i in 1:nrow(Anoa_HaploLF_Region_NoZO)) {
  # get all the non-zero % components
	xxx <-with(Anoa_HaploLF_Region_NoZO, 
		c(Prop_A4a_WC[i], Prop_A1_NE[i], Prop_A5b_SE[i], Prop_A2_NW[i], Prop_A3_EC[i], Prop_A5a_BT[i]))[which(Anoa_HaploLF_Region_NoZO[i, 1:6] > 0)]
	
	# add a pie chart of the % components (preserving colour order for each component) using a 0.75 transparency
	with(Anoa_HaploLF_Region_NoZO, floating.pie(Long_Reg_LF[i], Lat_Reg_LF[i], x = xxx, col = 
		alpha(WesAndersonCol[c(2, 3, 5, 6, 8, 7)][which(Anoa_HaploLF_Region_NoZO[i, 1:6] > 0)], 0.99), radius = 0.4, border = NA))
}

# add the count of samples ontop of each pie chart
with(Anoa_HaploLF_Region_NoZO, text(Long_Reg_LF, Lat_Reg_LF, Sum))

legend(125.75, 1, legend = c("A1", "A2", "A3", "A4", "A5a", "A5b"),
	fill = alpha(WesAndersonCol[c(3, 6, 8, 2, 7, 5)], 0.75), cex = 1.2, bty = "n", title = "Clade")

dev.off()



