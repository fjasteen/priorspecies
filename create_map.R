# 4. GRIIS Intersect met terreinen in beheer van DVW ####
## Libraries ####
# Dit kan je ook allemaal met sf maar ik ben old skool en gebruik sp 
library(rgeos)
library(raster)
library(sp)
library(rgdal)

crs_wgs <- "EPSG:4326"

## Shape DVW inlezen ####

#Should be done from drive... https://drive.google.com/drive/folders/1rs1YXrU_wBifxH_1Z4DfMkj7000cAEcU?usp=sharing
sectoren_dvw <- readOGR(dsn="./Data/priorspecies/Sectoren", layer="DVW_SECTOREN_BASIS")

## Fgdb dVW inlezen ####

# The input file geodatabase
# List all feature classes in a file geodatabase
# Read the feature class
#fgdb <- "./Data/priorspecies/Percelen/PercelenAzsCentraalGIS.gdb"
#subset(ogrDrivers(), grepl("GDB", name))
#fc_list <- ogrListLayers(fgdb)
#print(fc_list)
#perc_eig_dvw <- readOGR(dsn=fgdb,layer="AZS_PercelenInEigendom")

# de shape van DVW transformeren
crs_wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs") #WGS84
sectoren_dvw_wgs <- spTransform(sectoren_dvw, crs_wgs)
#perc_eig_dvw_wgs <- spTransform(perc_eig_dvw, crs_wgs)


# gbif data omzetten naar spatial
latlon <- rawdata_set1_imported %>%
  dplyr::select(decimalLongitude, decimalLatitude,speciesKey) 


Gbif_sp <- SpatialPointsDataFrame(latlon,
                                  data = rawdata_set1_imported,
                                  proj4string = crs_wgs)


# Behoud enkel de datapunten die binnen de sectoren liggen (in KERN type=KERN)
gbif_dvw <- raster::intersect(Gbif_sp, sectoren_dvw_wgs)
gbif_dvw <- gbif_dvw[gbif_dvw@data$TYPE == "KERN",] 

# Histogram van # waarnemingen per taxonKey -> om te evalueren welke soorten we weerhouden?
ggplot(gbif_dvw@data %>% count(taxonKey) %>% arrange(desc(n)), aes(x=n)) +
  geom_histogram(binwidth=10) 

## gbif data in beheer spatial joinen met km - hokken ####
# als test neem be10grid ipv be1grid
download.file("https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2/gis-files/belgium-spatialite/at_download/file", destfile = file.path(tempdir(), "Belgium_spatialite.zip"), mode = "wb")
unzip(zipfile = file.path(tempdir(), "Belgium_spatialite.zip"), exdir = tempdir())  

be1grid <- readOGR(dsn = file.path(tempdir(), "Belgium.sqlite"), 
                   layer = "be_1km")
be1grid <- spTransform(be1grid, crs_wgs)
sectoren_grid <-raster::intersect(sectoren_dvw_wgs[sectoren_dvw_wgs$TYPE == "KERN",], be1grid)

gbif_over<-over(gbif_dvw,be1grid)
gbif_be1grid <-bind_cols(gbif_over, gbif_dvw@data)

list_cells_kern <-sectoren_grid$cellcode
kern_grid<-subset(be1grid, cellcode %in% list_cells_kern)
freq_occ <- count(gbif_be1grid, cellcode)
kern_grid_occ <- merge(kern_grid,freq_occ,by="cellcode")



### Plot occurences ###

library(bcmaps)
library(ows4R)
library(httr)
library(leaflet)
library(leaflet.extras)
library(mapview)
library(tidyverse)
library(RColorBrewer)

wfs_regions <- "https://eservices.minfin.fgov.be/arcgis/services/R2C/Regions/MapServer/WFSServer"
regions_client <- WFSClient$new(wfs_regions, 
                                serviceVersion = "2.0.0")
regions_client$getFeatureTypes(pretty = TRUE)
url <- parse_url(wfs_regions)
url$query <- list(service = "wfs",
                  version = "2.0.0", # optional
                  request = "GetFeature",
                  typename = "regions",
                  srsName = "EPSG:3857"
)
request <- build_url(url)
bel_regions <- readOGR(request, stringsAsFactors = FALSE) #Lambert2008
flanders <- subset(bel_regions, bel_regions$NameDUT == "Vlaams Gewest")
flanders <- spTransform(flanders, crs_wgs)
flanders <- fix_geo_problems(flanders, tries = 5)
bbox <- as.data.frame(flanders@bbox)

mybins <- c(0,10,20,50,100,500,Inf)
mypalette <- colorBin( palette="YlOrRd", domain=kern_grid_occ@data$n, na.color="transparent", bins=mybins)

map <- leaflet(gbif_dvw) %>% 
  addPolylines(data = flanders,
               weight = 1,
               color = "black") %>%
  addPolygons(data = kern_grid_occ,
              stroke = TRUE,
              fillOpacity = 0.5,
              color="black",
              opacity = 1,
              weight = 0.1,
              fillColor = ~mypalette(n)) %>%
  addCircleMarkers(data = gbif_dvw,
                   radius = 0.1,
                   fillColor = "red",
                   color = "red",
                   opacity = 1,
                   fillOpacity = 1,
                   weight = 0.5) %>% 
  setMaxBounds(lng1 = bbox$min[1],
               lat1 = bbox$min[2],
               lng2 = bbox$max[1],
               lat2 = bbox$max[2]) %>%
  setMapWidgetStyle(list(background = "white")) %>%
  addLegend( pal=mypalette, values=~n, opacity=0.9, title = "# occurences", position = "bottomleft" )


# export map
mapshot(map, 
        file = "./Data/Kaarten/Verspreiding.png",
        remove_controls = c("zoomControl"))