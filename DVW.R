# Download occurrence data from potentially invasive species from GBIF #
# Libraries #
library(rgbif)
library(tidyverse)
library(ggplot2)
library(INBOtheme)
library(googledrive)
theme_set(theme_inbo())
library(readxl)

# Credentials #
# in .Renviron save: # gbif_user = "gbif username", gbif_pwd = "gbif password", gbif_email = "gbif email"

gbif_user <- Sys.getenv("gbif_user")
gbif_pwd <- Sys.getenv("gbif_pwd")
gbif_email <- Sys.getenv("gbif_email")


# 1. Prepare list of potentially invasive species (present - not yet present) #
# 1a. List of taxonIDs from species in GRIIS checklist Belgium (https://www.gbif.org/dataset/6d9e952f-948c-4483-9807-575348147c7e) ###

df <-read_delim("./priorspecies/dwca-unified-checklist-v1.12/taxon.txt")
df <- df[df$kingdom =="Plantae",] %>% 
      mutate_at("taxonID", str_replace, "https://www.gbif.org/species/", "")

# 1b. List of taxonIDs from horizonscanningtool CABI #

cabi <- read_delim("./priorspecies/Horizonscan/Horizon Scanning_20220125154122643.csv")
cabi <- cabi[-c(1:5),]
colnames(cabi) <- cabi[1,]
cabi <- cabi[which(cabi$Phylum=='Spermatophyta'|
                  cabi$Phylum =='Pteridophyta'),]
#select with ones with only genus & species name
cabi<-cabi[grep("[ ]", cabi$`Preferred scientific name`), ]
                
# Get GBIF backbone taxon ID from taxonomic names, based on species rank, accepted status, matchtype exact #

library(taxize)
cabi_id<-as.data.frame(get_gbifid(
  na.omit(cabi$`Preferred scientific name`),
  ask = TRUE,
  messages = TRUE,
  rows = 1,
  rank = "species",
  method = "lookup",
  sciname = NULL,
  ))

# problem: not found: (Senecio squalidus subsp. rupestris) How to solve? Log output and grep not Found, lookup seperately?

# 1c Make a unified list of all species (horizon + present)
#list_taxonID <- c(cabi_id$ids,df$taxonID) %>% unique() %>% na.omit()

#Compare cabi-GRISS list with list frome Scheers et al.
Scheers<-read_excel("./priorspecies/Lijst_Handel_ScheersK.xlsx", range = "A4:D290")
Scheers<-Scheers[Scheers$`Taxon category`=="species",]
Scheers<-Scheers[grep("[ ]{1}", Scheers$'Horticultural name'), ]

library(taxize)
Scheers_id<-as.data.frame(get_gbifid(
  na.omit(Scheers$`Horticultural name`),
  ask = TRUE,
  messages = TRUE,
  rows = 1,
  rank = "species",
  method = "lookup",
  sciname = NULL,
))

# 2. Download occurence data from GBIF
# 2a. for CABI-species

cabi_taxonkeys <-cabi_id$ids %>% na.omit
occurence_counts <- vector()
for (i in 1:length(cabi_taxonkeys)){
  occurence_counts[i]<-occ_count(cabi_taxonkeys[i])
}


num_occ<-data.frame(cabi_taxonkeys,occurence_counts)
num_occ_large<-num_occ[num_occ$occurence_counts >= 100000,]
num_occ_small<-num_occ[num_occ$occurence_counts <= 100000,]

taxonkey_set1 <- pred_in("taxonKey", cabi_taxonkeys)
taxonkey_set1 <- pred_in("taxonKey", num_occ_small$cabi_taxonkeys)
taxonkey_set2 <- pred_in("taxonKey", num_occ_large$cabi_taxonkeys)

# initiate download

set1 <- occ_download(taxonkey_set1,
                     pred("hasCoordinate", TRUE), # Georefereerde gegevens
                     pred_gte("year", 1900), # alle waarnemingen vanaf 1900
                     pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "UNKNOWN")),
                     user = gbif_user, 
                     pwd = gbif_pwd, 
                     email = gbif_email)

set2 <- occ_download(taxonkey_set2,
                     pred("hasCoordinate", TRUE), # Georefereerde gegevens
                     pred_gte("year", 1900), # alle waarnemingen vanaf 1900
                     pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "UNKNOWN")),
                     user = gbif_user, 
                     pwd = gbif_pwd, 
                     email = gbif_email)

# download data from GBIF

repeat{
  Sys.sleep(time = 5*length(acceptedKeys1))
  test_set1 <- occ_download_meta(set1)
  if(test_set1$status == "SUCCEEDED"){
    rawdata_set1_imported <- occ_download_get(set1, 
                                              path = "./priorspecies",
                                              overwrite = TRUE) %>% 
      occ_download_import()
    break
  }
  print(test_set1$status)
}

# 2. Download occurence data from GBIF
# 2b. for GRIIS-species

acceptedKeys <-df$taxonID
taxonkey_set3 <- pred_in("taxonKey", acceptedKeys)

# intiate download

set3 <- occ_download(taxonkey_set3,
                     pred("hasCoordinate", TRUE), # Georefereerde gegevens
                     #pred_gte("year", 2000), # alle waarnemingen vanaf 2000
                     user = gbif_user, 
                     pwd = gbif_pwd, 
                     email = gbif_email)

# download data from GBIF

repeat{
  Sys.sleep(time = 5*length(acceptedKeys))
  test_set2 <- occ_download_meta(set1)
  if(test_set2$status == "SUCCEEDED"){
    rawdata_set2_imported <- occ_download_get(set2, 
                                              path = "./priorspecies",
                                              overwrite = TRUE) %>% 
      occ_download_import()
    break
  }
  print(test_set2$status)
}

#2b. Compare list CABI with list Kevin Scheers 280 additional plants


# 3. Climate matching CABI species
# Testcase 1 species from CABI: C:\Users\frederique_steen\Documents\GitHub\priorspecies\0120469-210914110416597.zip
library(readr)
library(googlesheets4)
library(dplyr)
library(rgbif)
library(sp)
library(devtools)
install_github("trias-project/trias", ref = "73_climate_matching_function")
library(trias)

check <- function(x){
  tryCatch(if(class(x) == 'logical') 1 else 1, error=function(e) 0)
} 

get_cred <- function(x){
  library(svDialogs)
  
  cred <- Sys.getenv(x)
  
  if(cred == ""){
    input <- dlgInput(paste0("What is your ", x, "?"))
    cred <- input$res
    Sys.setenv(x = cred)
  }
  return(cred)
}

#vector with taxonkeys
#taxonkeys <- num_occ_small$cabi_taxonkeys#as.integer(unique(cabi_id$ids %>% na.omit))
taxonkeys <- num_occ_large$cabi_taxonkeys
#climate matching function

zipfile <- "C:\Users\frederique_steen\Documents\GitHub\priorspecies\large_occ\0144617-210914110416597_large.zip"
#target <- 75 
target <- 1
cuts <- ceiling(length(taxonkeys)/target)
gc(reset = TRUE)
if(memory.limit() <= 3000 * length(taxonkeys)){
  memory.limit(size = 3000 * length(taxonkeys))
}
for(i in 1:cuts){
  gc(reset = TRUE)
  if(i == 1){
    start <- 1
    end <- target
  }else{
    start <- (i-1)*target
    end <- i*target
    if(end > length(taxonkeys)){
      end <- length(taxonkeys)
    }
  }
  taxonkeys_sub <- taxonkeys[start:end]
  output <- climate_match(region = "Belgium",
                          taxonkey = taxonkeys_sub,
                          BasisOfRecord = c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "UNKNOWN"),
                          cm_limit = 0.2,
                          n_limit = 90,
                          maps = FALSE)
  
  assign(paste0("output_", i), output)
}

#merge all the resulting data in one file (& write as csv)
final_unfiltered <- data.frame()
final_filtered <- data.frame()
final_cm <- data.frame()
raw_data_final <- data.frame()
for(i in 1:cuts){
  data_list <- get(paste0("output_", i))
  data_filtered <- data_list$filtered
  data_cm <- data_list$cm
  data_unfiltered <- data_list$unfiltered
  data_raw <- data_list$spatial@data
  
  if(nrow(final_filtered)==0){
    final_filtered <- data_filtered
  }else{
    final_filtered <- rbind(final_filtered, data_filtered)
  }
  
  if(nrow(final_cm)==0){
    final_cm <- data_cm
  }else{
    final_cm <- rbind(final_cm, data_cm)
  }
  
  if(nrow(final_unfiltered)==0){
    final_unfiltered <- data_unfiltered
  }else{
    final_unfiltered <- rbind(final_unfiltered, data_unfiltered)
  }
  
  if(nrow(raw_data_final)==0){
    raw_data_final <- data_raw
  }else{
    raw_data_final <- rbind(raw_data_final, data_raw)
  }
}
final_filtered <- final_filtered %>% 
  distinct(taxonKey, Classification, scenario, .keep_all = TRUE)
write_csv(final_filtered, "./plants_data_overlay_future_filtered.csv")
final_cm <- final_cm %>% 
  distinct(taxonKey, Classification, scenario, .keep_all = TRUE) 
write_csv(final_cm, "./plants_data_overlay_future.csv")
final_unfiltered <- final_unfiltered %>% 
  distinct(taxonKey, Classification, .keep_all = TRUE) 
write_csv(final_unfiltered, "./plants_data_overlay_present.csv")
write_csv(raw_data_final, gsub(".zip", ".csv", zipfile))

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
