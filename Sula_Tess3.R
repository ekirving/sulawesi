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

# text sizing
text_cex = 1.8
point_cex = 1

# how many graduated steps to show in the colour palette (larger = smoother gradient)
palette.length = 10

# Tess3 interpolation method
interpol = FieldsKrigModel()
# interpol = FieldsTpsModel()

# Tess3 resolution for the interpolation surface
# resolution <- c(100,100)
resolution <- c(1000,1000)

# display the outline of the map, overlaid on the admixture components
map.display = TRUE

# set the lat/long boundaries for the map
# xlim <- c(118.5, 127.4)
xlim <- c(118.6, 127.3)
# xlim <- c(119.1, 127.2)
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
       main = NA, bty = 'n', cex = point_cex,
       window = c(xlim, ylim),
       xlim =  xlim, ylim = ylim,
       col.palette = my.palette)

  # render the actual map
  if (map.display) {
    plot(map.polygon, xlim = xlim, ylim = ylim, add=TRUE)
  }

  # add the count of samples for each zone
  segment.labels <- cbind(name=rownames(Sula_Lab),Sula_Lab)
  counts <- as.data.frame(table(species.data$Loc_Abbrev_LF))
  segment.labels <- merge(x=segment.labels, y=counts, by.x="name", by.y="Var1", fill=0, all.x = TRUE)
  segment.labels[is.na(segment.labels)] <- 0

  # add the regions, labels and mask some of the peripheral islands
  annotate_map(segment.labels, Sula_Seg, species.colours, species.clades, legend.title)
}

annotate_map <- function(segment.labels, segment.lines, species.colours, species.clades, legend.title = NA) {

  # mask the edges of the adjacent islands using white polygons
  polygon(c(par()$usr[1], par()$usr[1], 119.25, 119.25), c(0, 2, 2, 0), border = "white", col = "white")
  polygon(c(120, 120, 121, 121), c(par()$usr[3], -5.7, -5.7, par()$usr[3]), border = "white", col = "white")
  polygon(c(127, 127, par()$usr[2], par()$usr[2]), c(-1.8, 2, 2, -1.8), border = "white", col = "white")

  # segment the island into regions (using the coordinates in segment.lines)
  with(segment.lines, segments(Long_1, Lat_1, Long_2, Lat_2, lwd = 3))

  # remove the UN and ZO rows from the matrix of label coordinates
  segment.labels <- segment.labels[-which(segment.labels$name %in% c("UN", "ZO", "SU_BU", "S")),]

  # add the count of n for each label
  segment.labels$label <- paste(segment.labels$name, " (n=", segment.labels$Freq, ")", sep="")

  # label the regions
  with(segment.labels, text(Long, Lat, labels = segment.labels$label, cex = text_cex, pos = 4))

  # add arrows for ambiguous label positions
  arrows(122.674413, -0.411239, 122.35, -0.394760, length = 0.1, code = 2, lwd = 2)
  arrows(123.649449, -1.328482, 124.001012, -1.314753, code = 1, lwd = 2, length = 0.1)

  if (!is.na(legend.title)) {
    # trim any trailing suffix
    species.clades <- sub("_.+", "", species.clades)

    # handle special case of A4a
    species.clades <- sub("A4a", "A4", species.clades)

    # add the legend
    # legend(125.75, 1.7, legend = species.clades, fill = species.colours, cex = text_cex, bty = "n", title = legend.title)
    legend(125.5, 2.1, legend = species.clades, fill = species.colours, cex = text_cex, bty = "n", title = legend.title)
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

# adjust the label positions (pretty cludgy but cleaner than changing the Rdata)
Sula_Lab$Long <- Sula_Lab$Long-0.3
Sula_Lab['WC','Lat'] <- Sula_Lab['WC','Lat'] + 2.2 # up
Sula_Lab['WC','Long'] <- Sula_Lab['WC','Long'] - 0.3 # left
Sula_Lab['NE','Lat'] <- Sula_Lab['NE','Lat'] + 0.35 # up
Sula_Lab['NE','Long'] <- Sula_Lab['NE','Long'] - 0.6 # left
Sula_Lab['NC','Lat'] <- Sula_Lab['NC','Lat'] + 0.10 # up
Sula_Lab['NC','Long'] <- Sula_Lab['NC','Long'] - 0.5 # left
Sula_Lab['NW','Lat'] <- Sula_Lab['NW','Lat'] + 0.30 # up
Sula_Lab['NW','Long'] <- Sula_Lab['NW','Long'] - 0.75 # left
Sula_Lab['SW','Lat'] <- Sula_Lab['SW','Lat'] + 0.20 # up
Sula_Lab['SW','Long'] <- Sula_Lab['SW','Long'] + 1.05 # right
Sula_Lab['SU','Long'] <- Sula_Lab['SU','Long'] - 0.2 # left
Sula_Lab['BU','Long'] <- Sula_Lab['BU','Long'] - 0.5 # left
Sula_Lab['EC','Lat'] <- Sula_Lab['EC','Lat'] - 0.5 # down
Sula_Lab['SE','Lat'] <- Sula_Lab['SE','Lat'] + 0.2 # up
Sula_Lab['SE','Long'] <- Sula_Lab['SE','Long'] + 0.3 # right
Sula_Lab['BT','Long'] <- Sula_Lab['BT','Long'] + 0.2 # right

# start preparing the data
location.columns <- c('Loc_Abbrev_LF', 'Longitude', 'Latitude')

# extract the STRUCTURE components and the sample locations
anoa.data <- Anoa_All_Genetics[,c('A1_BT', 'A2_NW', 'A3_SE', 'A4_NE_WC', 'A5_NC_EC', location.columns)]
baby.data <- Baby_All_Genetics[,c('B1_WC_NW', 'B2_SE', 'B3_SU_BU', 'B4_NE', 'B5_TO', location.columns)]
susc.data <- Sus_cel_All_Genetics[,c('S1_NW', 'S2_PE', 'S3', 'S4_SE_BT', 'S5_WC_SW', 'S6_BU', 'S7_EC', location.columns)]

# setup the STRUCTURE colour lists
anoa.colours <- c('#a6b6bb', '#95afdf', '#d4daf9', '#fe8b8d', '#d7664d')
baby.colours <- c('#ecb8d3', '#d4daf9', '#f5cc9c', '#e09568', '#fe8b8d')
susc.colours <- c('#95afdf', '#a4dfea', '#f7c8ce', '#d4daf9', '#fe8b8d', '#f5cc9c', '#d7664d')

# get the mtDNA clades and transpose them into a Q matrix
anoa.data2 <- clade2q(Anoa_All_Genetics[,c('LF_Anoa_Clust', location.columns)])
baby.data2 <- clade2q(Baby_All_Genetics[,c('LF_Baby_Clust', location.columns)])
susc.data2 <- clade2q(Sus_cel_All_Genetics[,c('LF_Sus_cel_Clust', location.columns)])

# setup the mtDNA colour lists
anoa.colours2 <- c('#e09568', '#95afdf', '#d7664d', '#fe8b8d', '#a6b6bb', '#d4daf9')
baby.colours2 <- c('#f5cc9c', '#ecb8d3', '#e09568', '#fe8b8d', '#d4daf9', '#95afdf')
susc.colours2 <- c('#95afdf', '#f5cc9c', '#a4dfea', '#b5b077', '#fe8b8d')

pdf(file = "sulawesi/pdf/Sula-Maps.pdf", width = (xlim[2]-xlim[1])*3*.8, height = (ylim[2]-ylim[1])*2*.8)
# png(file = "sulawesi/pdf/Sula-Maps.png", width = (xlim[2]-xlim[1])*3*.8, height = (ylim[2]-ylim[1])*2*.8, units = 'in', res=300)

# display the maps in a 2x3 grid
par(mfrow = c(2,3), mar = c(6, 5, 6, 5), bty = 'n')

# Clade maps
render_map(anoa.data2, anoa.colours2, map.display, 'Clades')
render_map(baby.data2, baby.colours2, map.display, 'Clades')
render_map(susc.data2, susc.colours2, map.display, 'Clades')

# STRUCTURE maps
render_map(anoa.data, anoa.colours, map.display, 'Ancestral (K)')
render_map(baby.data, baby.colours, map.display, 'Ancestral (K)')
render_map(susc.data, susc.colours, map.display, 'Ancestral (K)')

dev.off()
