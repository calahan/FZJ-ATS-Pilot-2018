#
# Candidate Experiments for an Algal Floway Research and Development Program at Forschungszentrum-Jülich, 2019-2021"
# CEfaAFRaDPaFZJ"
#
source("Code/R/Settings.R")

library(tidyverse)
library(readxl)

stations <- read_excel(GEMss)

deu <- stations %>% filter(Country_code == "DEU")

type <- unique(deu$Station_characteristic)
name <- unique(deu$Catchment_name)
