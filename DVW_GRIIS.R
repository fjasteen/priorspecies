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


# 1. List of taxonIDs from species in GRIIS checklist Belgium (https://www.gbif.org/dataset/6d9e952f-948c-4483-9807-575348147c7e) ###

df <-read_delim("./dwca-unified-checklist-v1.12/taxon.txt")
df <- df[df$kingdom =="Plantae",] %>% 
  mutate_at("taxonID", str_replace, "https://www.gbif.org/species/", "")
GRIIS_taxonID <- df$taxonID