library(rgbif)
library(tidyverse)
library(ggplot2)
library(INBOtheme)
library(googledrive)
theme_set(theme_inbo())

## Credentials ####
# in .Renviron in sla je op:
# gbif_user = "je gbif username"
# gbif_pwd = "je gbif password"
# gbif_email = "je gbif email"

gbif_user <- Sys.getenv("gbif_user")
gbif_pwd <- Sys.getenv("gbif_pwd")
gbif_email <- Sys.getenv("gbif_email")



## Lijst van taxonkeys obv gedownloade GRIIS checklist Belgium (https://www.gbif.org/dataset/6d9e952f-948c-4483-9807-575348147c7e) ####
df <-read_delim("./priorspecies/dwca-unified-checklist-v1.12/taxon.txt")
df <- df[df$kingdom =="Plantae",] %>%
  mutate_at("taxonID", str_replace, "https://www.gbif.org/species/", "")


cabi <- read_delim("./priorspecies/Horizonscan/Horizon Scanning_20220125154122643.csv")
cabi <- cabi[-c(1:5),]
colnames(cabi) <- cabi[1,]
cabi <- cabi[cabi$`Taxonomic group` == 'Plants', ] 


#Get the GBIF backbone taxon ID from taxonomic names, selecting species rank & accepted status + matchtype exact (first row is accepted match)
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

# Make a list of all species (horizon + present)
list_taxonID <- c(cabi_id$ids,df$taxonID) %>% unique() %>% na.omit()


## Taxonkeys of gbifkey #### 
# Download in R (meest up to date) versus alternatief, download via gbif (sneller & eenvoudiger te ref)

### Initieer de download in R ####
# unieke taxonkey op gbif:: 
# Ludwigia longifolia (DC.) Hara_ is de species pagina: https://www.gbif.org/species/5421003 -> taxonkey is 5421003
# acceptedKeys <- <integer lijst van de taxonkeys van de soorten waar je interesse in hebt>
# gebruik df$taxonID[1:30] als test, n in df is > 2000

acceptedKeys <-list_taxonID
taxonkey_set1 <- pred_in("taxonKey", acceptedKeys)

# Vervolgens initiëren we de download
set1 <- occ_download(taxonkey_set1,
                     pred("hasCoordinate", TRUE), # Georefereerde gegevens
                     #pred_gte("year", 2000), # alle waarnemingen vanaf 2000
                     user = gbif_user, 
                     pwd = gbif_pwd, 
                     email = gbif_email)

# Wanneer de download klaar is downloaden we de data 
repeat{
  Sys.sleep(time = 5*length(acceptedKeys))
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