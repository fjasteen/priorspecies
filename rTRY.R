library(rtry)
library(dplyr)
library(tidyr)
library(jsonlite)
library(curl)
library(data.table)
library(taxize)
library(tidyverse)    # To do datascience
library(rgbif)        # To lookup names in the GBIF backbone taxonomy
library(inborutils)   # To wrap GBIF API data
library(knitr)
library(traitdataform)

library(devtools)
install_github("ropensci")

### GOAL: make list of AccessionID from TRY database from horizon_species & GRIIS_list ###
packageVersion("rtry")
trait <- read.delim("./Traits/TRY_traits.txt", sep = "\t")
species <- read.delim("./Traits/TRY_species.txt", sep = "\t")

# Which plants to consider from horizonscan/climatematch?
final_filtered <-
  read.csv("./plants_data_overlay_future_filtered.csv")
final_filtered_crit <-
  subset(final_filtered, n_totaal >= 90 & perc_climate >= 0.2)
horizonscan_species <-
  unique(final_filtered_crit$acceptedScientificName)
horizonscan_species_taxonKey <- unique(final_filtered_crit$taxonKey)

# Which plants to consider from GRIIS
df <- read_delim("./lists/dwca-unified-checklist-v1.12/taxon.txt")
df <- df[df$kingdom == "Plantae", ] %>%
  mutate_at("taxonID", str_replace, "https://www.gbif.org/species/", "")
GRIIS_taxonKeys <- df$taxonID
df$speciesname <- str_c(df$genus, " ", df$specificEpithet)
GRIIS_list <- df[, c("taxonID", "scientificName", "speciesname")]

# match species name with accepted species name on Gbif/ usage key
TRY_gbif <-
  gbif_species_name_match(df = species, name = "AccSpeciesName")
write_csv(TRY_gbif, "./Traits/TRY_gbif.csv")
TRY_gbif <- TRY_gbif[TRY_gbif$kingdom == "Plantae", ]
TRY_gbif <-
  TRY_gbif[TRY_gbif$rank %in% c("SPECIES", "SUBSPECIES", "VARIETY"), ]

# For all synonymous entries in TRY_gbif$status == TRUE lookup accepted taxonkey & scientificname
synonyms <- TRY_gbif[TRY_gbif$synonym == TRUE, ]
syn_Acc <- data.frame(matrix(0, ncol = 24, nrow = 0))
colnames(syn_Acc) <-
  colnames(as.data.frame(name_backbone(name = synonyms$scientificName[1])))
for (i in 1:nrow(synonyms)) {
  syn_Acc[i, ] <-
    as.data.frame(name_backbone(name = synonyms$scientificName[i]))
}
#replace the usageKey in TRY_gbif (test) by the usagekey of the accepted name
for (i in syn_Acc$usageKey) {
  TRY_gbif["usageKey"][TRY_gbif["usageKey"] == i] <-
    as.integer(syn_Acc$acceptedUsageKey[match(i, syn_Acc$usageKey)])
}

# For all varieties/subspecies entries in TRY_gbif$rank lookup accepted parent specieskey
subsp_var<- TRY_gbif[TRY_gbif$rank == "SUBSPECIES" | TRY_gbif$rank == "VARIETY", ]
subsp_var_spKey <-vector()
for (i in subsp_var$usageKey){
  subsp_var_spKey[i]<-name_usage(key = i, data = 'all', rank = 'SPECIES')$data$speciesKey
}


#replace the usageKey in TRY_gbif by the speciesKey instead of key of subspecies & variety
for (i in syn_Acc$usageKey) {
  TRY_gbif["usageKey"][TRY_gbif["usageKey"] == i] <-
    as.integer(syn_Acc$acceptedUsageKey[match(i, syn_Acc$usageKey)])
}
name_usage(key = XXXXX[i], data = 'all', rank = 'SPECIES')$data$speciesKey


# which species of the horizonscan are within TRY if matched with usageKey?
horizon_TRY_tK <-
  TRY_gbif[TRY_gbif$usageKey %in% unique(final_filtered_crit$taxonKey), ] # 333 found out of 374

# Which taxa from horizonscan are within TRY if matched with species? 46 species are lacking
horizon_TRY_abs_tK <-
  setdiff(horizonscan_species_taxonKey, horizon_TRY_tK$usageKey)
horizon_TRY_abs_sp <-
  unique(final_filtered_crit[final_filtered_crit$taxonKey %in% horizon_TRY_abs_tK, ]$acceptedScientificName)
horizon_TRY_sp <-
  TRY_gbif[TRY_gbif$AccSpeciesName %in% horizon_TRY_abs_sp, ]

#36 species are found, so still 5 species lacking: which?


#add records matched through species names to horizon_TRY 333 from taxonKey, and 36 from species
horizon_TRY <- rbind(horizon_TRY, horizon_TRY_sp)

#Which of the records are not found from the horizon list
horizon_absent <-
  setdiff(horizon_TRY_abs_sp, horizon_TRY_sp$AccSpeciesName)

#NEED to perform the same action with the GRIIS species...(get GRIISlist grom DVW_horizon.R)
GRIIS_TRY <-
  TRY_gbif[TRY_gbif$usageKey %in% unique(GRIIS_list$taxonID), ]

GRIIS_TRY_abs_tK <-
  setdiff(GRIIS_list$taxonID, unique(GRIIS_TRY$usageKey))
GRIIS_TRY_abs_sp <-
  unique(GRIIS_list[GRIIS_list$taxonID %in% GRIIS_TRY_abs_tK, ]$speciesname)
GRIIS_TRY_sp <-
  TRY_gbif[TRY_gbif$AccSpeciesName %in% GRIIS_TRY_abs_sp, ]

GRIIS_absent <-
  setdiff(GRIIS_TRY_abs_sp, GRIIS_TRY_sp$AccSpeciesName)

# Lookup synonyms from GBIF
#tck <- taxonsearch(scientificname = GRIIS_absen, dataresourcekey=1)


# Add origin of species to both horizon_TRY & GRIIS_TRY & MERGE
rbind(horizon_TRY, GRIIS_TRY)



# GOAL: Obtain traits from selected Acc_numbers from TRY database