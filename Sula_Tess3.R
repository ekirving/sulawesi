library(maps)
library(mapdata)
library(maptools)
library(tess3r)
library(wesanderson)

# set working directory
setwd('/Users/Evan/Dropbox/Code/')

# load local bugfixed TESS3 code
insertSource("TESS3/R/spatialPlots.R", package = "tess3r", functions = 'PlotInterpotationMax')

# load the Sulawesi RData
load('sulawesi/Sulawesi.RData')

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

# how many graduated steps to show in the colour palette (larger = smoother gradient)
palette.length = 10

# Tess3 interpolation method
interpol = FieldsKrigModel()
# interpol = FieldsTpsModel()

# Tess3 resolution for the interpolation surface
resolution <- c(200,200)
# resolution <- c(1000,1000)

# display the outline of the map, overlaid on the admixture components
map.display = TRUE

# set the lat/long boundaries for the map
xlim <- c(118.5, 127.4)
ylim <- c(-5.9, 2)

# use the worldHires map
map.hires <- map("worldHires", xlim = xlim, ylim = ylim, plot=FALSE, fill=TRUE)
map.polygon <- map2SpatialPolygons(map.hires, map.hires$names)
# map.polygon <- rworldmap::getMap(resolution = "high")

render_map <- function(species.data, species.colours, map.display = TRUE, legend.title = NA) {

  # drop any rows with missing components or locations
  species.data <- na.omit(species.data)

  # drop Zoo samples and those of unknown origin
  species.data <- species.data[-which(species.data$Loc_Abbrev_LF %in% c("UN", "ZO", "SU_BU", "S")),]

  # drop any columns with all zero values (i.e. any orphan components for Zoo samples)
  species.data <- species.data[, c(!colSums(species.data[1:(ncol(species.data)-3)])==0, TRUE, TRUE, TRUE)]

  # extract the clade/population names
  species.clades <- names(species.data)[1:(ncol(species.data)-3)]

  # get the STRUCTURE components ready for tess3
  species.struct <- as.matrix(species.data[,1:(ncol(species.data)-3)])
  class(species.struct) <- c("tess3Q", "matrix")

  # get the coordinates of the remaining samples
  species.coord <- species.data[,c('Longitude', 'Latitude')]

  # make the colour palette (and reverse the order because tess3 reads this backwards)
  my.palette <- rev(CreatePalette(species.colours[1:length(species.clades)], palette.length = palette.length))

  # plot the interpolation surface
  plot(x=species.struct, coord=species.coord,
       method = "map.max",
       map.polygon = map.polygon,
       resolution = resolution,
       interpol = interpol,
       xlab = NA, ylab = NA, xaxt='n', yaxt='n',
       main = NA, bty = 'n', cex = .4,
       window = c(xlim, ylim),
       xlim =  xlim, ylim = ylim,
       col.palette = my.palette)

  # render the actual map
  if (map.display) {
    plot(map.polygon, xlim = xlim, ylim = ylim, add=TRUE)
  }

  # add the regions, labels and mask some of the peripheral islands
  annotate_map(Sula_Lab, Sula_Seg, species.colours, species.clades, legend.title)
}

annotate_map <- function(segment.labels, segment.lines, species.colours, species.clades, legend.title = NA) {

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
  if (!is.na(legend.title)) {
    # trim any trailing suffix
    species.clades <- sub("_.+", "", species.clades)
    legend(125.75, 1, legend = species.clades, fill = species.colours, cex = 1.2, bty = "n", title = legend.title)
  }
}

clade2q <- function(species.data) {

  clust.column <- names(species.data)[1]

  # extract the clades
  clades <- sort(levels(species.data[[clust.column]]))

  # convert clades into a Q matrix with 100% identity
  for (i in 1:length(clades)) {
    species.data[[clades[i]]] <- with(species.data, as.integer(clades[i] == species.data[clust.column]))
  }

  # put the columns in the right order
  species.data[,c(5:ncol(species.data), 2:4)]
}

# start preparing the data
location.columns <- c('Loc_Abbrev_LF', 'Longitude', 'Latitude')

# extract the STRUCTURE components and the sample locations
anoa.data <- Anoa_All_Genetics[,c('A1_BT', 'A2_NW', 'A3_SE', 'A4_NE_WC', 'A5_NC_EC', location.columns)]
baby.data <- Baby_All_Genetics[,c('B1_WC_NW', 'B2_SE', 'B3_SU_BU', 'B4_NE', 'B5_TO', location.columns)]
susc.data <- Sus_cel_All_Genetics[,c('S1_NW', 'S2_PE', 'S3', 'S4_SE_BT', 'S5_WC_SW', 'S6_BU', 'S7_EC', location.columns)]

# anoa.data <- subset(Anoa_All_Genetics, select = c(A4_NE_WC:A1_BT, Loc_Abbrev_LF, Longitude, Latitude))
# baby.data <- subset(Baby_All_Genetics, select = c(B3_SU_BU:B2_SE, Loc_Abbrev_LF, Longitude, Latitude))
# susc.data <- subset(Sus_cel_All_Genetics, select = c(S5_WC_SW:S3, Loc_Abbrev_LF, Longitude, Latitude))

# TODO setup the colour lists
anoa.colours <- WesAndersonCol[1:(ncol(anoa.data)-3)]
baby.colours <- WesAndersonCol[1:(ncol(baby.data)-3)]
susc.colours <- WesAndersonCol[1:(ncol(susc.data)-3)]

# get the cluster data
anoa.data2 <- clade2q(Anoa_All_Genetics[,c('LF_Anoa_Clust', location.columns)])
baby.data2 <- clade2q(Baby_All_Genetics[,c('LF_Baby_Clust', location.columns)])
susc.data2 <- clade2q(Sus_cel_All_Genetics[,c('LF_Sus_cel_Clust', location.columns)])

# TODO setup the colour lists
anoa.colours2 <- WesAndersonCol[1:(ncol(anoa.data2)-3)]
baby.colours2 <- WesAndersonCol[1:(ncol(baby.data2)-3)]
susc.colours2 <- WesAndersonCol[1:(ncol(susc.data2)-3)]

pdf(file = "sulawesi/pdf/Sula-Maps.pdf", width = (xlim[2]-xlim[1])*3*.8, height = (ylim[2]-ylim[1])*2*.8)

# display the maps in a 2x3 grid
par(mfrow = c(2,3), mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1), bty = 'n')

# Clade maps
render_map(anoa.data2, anoa.colours2, map.display, 'Clades')
render_map(baby.data2, baby.colours2, map.display, 'Clades')
render_map(susc.data2, susc.colours2, map.display, 'Clades')

# STRUCTURE maps
render_map(anoa.data, anoa.colours, map.display, 'Populations')
render_map(baby.data, baby.colours, map.display, 'Populations')
render_map(susc.data, susc.colours, map.display, 'Populations')

# dev.off()
