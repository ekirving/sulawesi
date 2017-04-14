
library(mapdata)
library(maps)
library(scales)
library(plotrix)

SulawesiMapRegionTmp <- map("worldHires",  xlim = c(122.05, 122.15), ylim = c(-0.425, -0.25), col = "white")
SulawesiMapRegionTmp$x[SulawesiMapRegionTmp$x > 122.15] <- NA
SulawesiMapRegionTmp$x[SulawesiMapRegionTmp$x < 122.05] <- NA

SulawesiMapRegion <- map("worldHires", c("Indonesia:Sulawesi", "Indonesia:Buton", "Indonesia:Muna", 
	"Indonesia:Wowoni", "Indonesia:Kabaena", "Indonesia:Manui", "Indonesia:Togin", "Indonesia:Waleabahi",
	"Indonesia:Batudaka", "Indonesia:Unauna", "Indonesia:Peleng", "Indonesia:Bangkulu", "Indonesia:Melilis",
	"Indonesia:Mansaleang", "Indonesia:Buru", "Indonesia:Sanana", "Indonesia:Mangole", "Indonesia:Taliabu",
	"Indonesia:Selayar", "Indonesia:Klabat"), fill = T, col = "white")
# missing togian island
with(SulawesiMapRegionTmp, polygon(x, y, col = "blue"))

SulawesiMapRegion_Mainland <- map("worldHires", c("Indonesia:Sulawesi"), add = T, 
	fill = T, col = "white")
SulawesiMapRegion_Buton <- map("worldHires", c("Indonesia:Buton",  "Indonesia:Muna", "Indonesia:Wowoni",
	"Indonesia:Kabaena", "Indonesia:Manui"), add = T, fill = T, col = "white")
SulawesiMapRegion_Togian <- map("worldHires", c("Indonesia:Togin", "Indonesia:Waleabahi",
	"Indonesia:Batudaka", "Indonesia:Unauna"), add = T, 
	fill = T, col = "white")
with(SulawesiMapRegionTmp, polygon(x, y, col = "white"))
SulawesiMapRegion_Banggai <- map("worldHires", c("Indonesia:Peleng", "Indonesia:Bangkulu", "Indonesia:Melilis",
	"Indonesia:Mansaleang"), 
	add = T, fill = T, col = "white")
SulawesiMapRegion_Buru <- map("worldHires", c("Indonesia:Buru"), add = T, 
	fill = T, col = "white")
SulawesiMapRegion_Sula <- map("worldHires", c("Indonesia:Sanana", "Indonesia:Mangole", "Indonesia:Taliabu"), add = T, 
	fill = T, col = "white")
SulawesiMapRegion_Selayer <- map("worldHires", c("Indonesia:Selayar", "Indonesia:Klabat"), add = T, 
	fill = T, col = "white")