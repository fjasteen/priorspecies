gbifcleaner <- function(zipfile){
  
  # Setup ####
  library(tidyverse)
  
  # Checks ####
  ## Zipfile exists ?? ####
  if(!file.exists(zipfile)){
    stop(paste(zipfile, " does not exist"))
  }
  
  ## Free some RAM ####
  fileinfo <- file.info(zipfile)
  
  memory.limit(size = fileinfo$size/1000)
  gc(reset = TRUE)
  
  ## Read data ####
  data <- read_tsv(unz(zipfile, "occurrence.txt"), 
                   col_types = c(decimalLatitude = col_number(),
                                 decimalLongitude = col_number(),
                                 establishmentMeans = col_character()))
  
  ## Select required_columns ####
  required_columns <- c("gbifID",
                        "basisOfRecord",
                        "taxonRank",
                        "taxonomicStatus",
                        "acceptedTaxonKey",
                        "genus",
                        "specificEpithet",
                        "eventDate",
                        "year",
                        "month",
                        "day",
                        "decimalLatitude",
                        "decimalLongitude",
                        "coordinateUncertaintyInMeters",
                        "occurrenceStatus",
                        "countryCode")
  
  data_redux <- data %>% 
    select(required_columns)
  
  ## rezip file ####
  ### Export data ####
  write_tsv(data_redux, "./data/data_redux.txt")
  ### zipfile ####
  zip(zipfile = zipfile,
      files = "./data/data_redux.txt",
      zip = Sys.getenv("R_ZIPCMD", "zip"))
  ### Remove temp txt ####
  file.remove("./data/data_redux.txt")
}