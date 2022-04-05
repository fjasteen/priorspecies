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

# 1. List of taxonIDs from horizonscanningtool CABI #
cabi <- read_delim("./Horizonscan/Horizon Scanning_20220125154122643.csv")
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

cabi_taxonID<-cabi_id$ids

# problem: not found: (Senecio squalidus subsp. rupestris) How to solve? Log output and grep not Found, lookup seperately?

#Compare list_taxonID with list from Scheers et al.
Scheers<-read_delim("./List_Scheers/Lijst_Handel_ScheersK.csv")
Scheers<-Scheers[Scheers$`Taxon category`=="species",]
Scheers<-Scheers[grep("[ ]{1,}", Scheers$'Horticultural name'), ]
Scheers_list <- Scheers$`Horticultural name`

# Get GBIF backbone taxon ID from taxonomic names, based on species rank, accepted status, matchtype exact #
library(taxize)
Scheers_id<-as.data.frame(get_gbifid(
  Scheers_list,
  ask = TRUE,
  messages = TRUE,
  rows = 1,
  rank = "species",
  method = "lookup",
  sciname = NULL,
))

Scheers_taxonID<-Scheers_id$ids

# Get GBIF backbone taxon ID from RIPARIAS taxonomic names, based on species rank, accepted status, matchtype exact #
Riparias <-read_delim("./lists_horizonscan/Lijst_Handel_ScheersK.csv")
Riparias <-Riparias[Riparias$is_plant==TRUE,]
Riparias_taxonID <-Riparias$gbif_key

# All horizonscan species Scheers et al. + cabi
horizon_taxonID <- unique(c(Scheers_taxonID,cabi_taxonID, Riparias_taxonID))

# Duplicate Species in GRIIS list and in horizonscan list?
horizon_taxonKeys <- setdiff(horizon_taxonID, GRIIS_taxonID) %>% na.omit

# 2. Download occurence data from GBIF 
# count occurence data as indication (occ_count does not take multiple values for basisOfRecord)
occurence_counts_HO <- vector()
occurence_counts_PR <- vector()
occurence_counts_UK <- vector()
for (i in 1:length(horizon_taxonKeys)){
  occurence_counts_HO[i]<-occ_count(horizon_taxonKeys[i], 
                                    georeferenced =TRUE, 
                                    basisOfRecord= "HUMAN_OBSERVATION", 
                                    from= 1900, to=2022) #didn't add this one the first time
}
for (i in 1:length(horizon_taxonKeys)){
  occurence_counts_PR[i]<-occ_count(horizon_taxonKeys[i], 
                                    georeferenced =TRUE, 
                                    basisOfRecord= "PRESERVED_SPECIMEN", 
                                    from= 1900, to=2022)
}
for (i in 1:length(horizon_taxonKeys)){
  occurence_counts_UK[i]<-occ_count(horizon_taxonKeys[i], 
                                    georeferenced =TRUE, 
                                    basisOfRecord= "UNKNOWN", 
                                    from= 1900, to=2022)
}

occurence_counts<-occurence_counts_HO + occurence_counts_PR + occurence_counts_UK
zipfile <- "C:/Users/frederique_steen/Documents/Data/priorspecies/occurence_data/occurence_data_horizon_cm/0147366-210914110416597.zip"


#There is a slight difference in the way records are counted here vs. results from occ_search(). For equivalent outcomes, in the occ_search() function use hasCoordinate=TRUE, and hasGeospatialIssue=FALSE to have the same outcome for this function using georeferenced=TRUE.!!!!
num_occ<-data.frame(horizon_taxonKeys,occurence_counts)
num_occ<- num_occ[order(num_occ$occurence_counts, decreasing = T),]
View(num_occ)
#num_occ_large<-num_occ[num_occ$occurence_counts >= 100000,]
#num_occ_small<-num_occ[num_occ$occurence_counts <= 100000,]

#taxonkey_set1 <- pred_in("taxonKey", cabi_taxonkeys)
taxonkey_set1 <- pred_in("taxonKey", num_occ$horizon_taxonKeys)

# initiate download

set1 <- occ_download(taxonkey_set1,
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
taxonkey_set2 <- pred_in("taxonKey", acceptedKeys)

# intiate download

set2 <- occ_download(taxonkey_set2,
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


# 3. Climate matching horizonscan species
# Testcase 1 species from CABI: C:/Users/frederique_steen/Documents/GitHub/priorspecies/0120469-210914110416597.zip
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
taxonkeys <- num_occ$horizon_taxonKeys
#climate matching function

zipfile <- "C:/Users/frederique_steen/Documents/Data/priorspecies/occurence_data/occurence_data_horizon_cm/0147366-210914110416597.zip"
#target <- 75 
target <- 2
cuts <- ceiling(length(taxonkeys)/target)
gc(reset = TRUE)
if(memory.limit() <= 300 * length(taxonkeys)){
  memory.limit(size = 300 * length(taxonkeys))
}
for(i in 25:cuts){
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

#running for 25:cuts:
# Download file size: 0 MB
# On disk at C:\Users\frederique_steen\Documents\GitHub\priorspecies\0148301-210914110416597.zip
# Error in climate_match(region = "Belgium", taxonkey = taxonkeys_sub, BasisOfRecord = c("HUMAN_OBSERVATION",  : 
# No useable data for 7445953,6063723,7329389 left after filters. Try omiting or changing the filter setup.
# Stopped at cut 240 for the above taxonKeys  (480)                                                                                   In addition: There were 50 or more warnings (use warnings() to see the first 50)

#merge all the resulting data in one file (& write as csv)
final_unfiltered <- data.frame()
final_filtered <- data.frame()
final_cm <- data.frame()
raw_data_final <- data.frame()
for(i in 25:cuts){
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
final_filtered2 <- final_filtered %>% 
  distinct(taxonKey, Classification, scenario, .keep_all = TRUE)
write_csv(final_filtered, "./plants_data_overlay_future_filtered2.csv")
final_cm2 <- final_cm %>% 
  distinct(taxonKey, Classification, scenario, .keep_all = TRUE) 
write_csv(final_cm, "./plants_data_overlay_future2.csv")
final_unfiltered2 <- final_unfiltered %>% 
  distinct(taxonKey, Classification, .keep_all = TRUE) 
write_csv(final_unfiltered, "./plants_data_overlay_present2.csv")
write_csv(raw_data_final, gsub(".zip", ".csv", zipfile))


