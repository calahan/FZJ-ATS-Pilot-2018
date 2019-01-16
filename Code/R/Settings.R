# Libraries
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(knitr)
library(PeriodicTable)
library(ggplot2)

# Directories
data_dir <- "Data/"
csv_dir <- paste0(data_dir, "csv/")
ss_dir <- paste0(data_dir, "Spreadsheets/")

# Conversion factors
yeardays <- 365     # days per year
yearweeks <- 52     # weeks per year
C2K <- 273.15       # convert °C to °K

# INPUT FILENAMES
# csv files exported from HOBO downloads
HOBO_fn <- paste0(csv_dir, c("ATS_1.11.09.18.csv",
                            "ATS_2.17.09.18.csv",
                            "ATS_1.24.09.18.csv",
                            "ATS_2.01.10.18.csv",
                            "ATS_1.08.10.18.csv",
                            "ATS_2.15.10.18.csv",
                            "ATS_1.22.10.18.csv"))

#fig_fn <- c("Visual Elements/Figures/Figure 1.png", "blah")
#fig_fn <- c("Figures/Figure1.png", "Visual Elements/Figures/Figure2.png")

# Spreadsheet of water chemistry and biomass
sswc_fn <- paste0(ss_dir, "ATS Treatment.xlsx")

# Spreadsheet of biomass composition
ssbc_fn <- paste0(ss_dir, "Analysen ZEA-3.xlsx")

# Output Files
pre <- "FZJ WWTP ATS Pilot "
ext <- ".csv"

# Cleaned temperature and illumination
ti_fn <- paste0(csv_dir, pre, "Temperature and Illumination", ext)

# Cleaned water chemistry
wc_fn <- paste0(csv_dir, pre, "Chemistry and Biomass", ext)

# Cleaned biomass productivity
bp_fn <- paste0(csv_dir, pre, "Biomass Productivity", ext)

# Cleaned biomass composition
bc_fn <- paste0(csv_dir, pre, "Pilot Biomass Composition", ext)