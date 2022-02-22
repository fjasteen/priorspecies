gbifcleaner <- function(zipfile){
  
  # Setup ####
  library(tidyverse)
  
  # Checks ####
  ## Zipfile exists ?? ####
  if(!file.exists(zipfile)){
    stop(paste(zipfile, " does not exist"))
  }
  
  # Free some RAM ####
  fileinfo <- file.info(zipfile)
  
  memory.limit(size = fileinfo$size/10000)
  gc(reset = TRUE)
  
  data <- read_tsv(unz(zipfile, "occurrence.txt"), 
                   col_types = c(decimalLatitude = col_number(),
                                 decimalLongitude = col_number(),
                                 establishmentMeans = col_character()))
  
}