library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(knitr)
library(PeriodicTable)
library(ggplot2)

source("Code/R/Processing.R")

# data directories
orig_data_dir <- file.path("Research", "Data")
data_dir <- file.path("Data")

# atomic and molecular weights
awS <- mass("S")
awP <- mass("P")
awK <- mass("K")
awCa <- mass("Ca")
awMg <- mass("Mg")
awMn <- mass("Mn")
awC <- mass("C")
awN <- mass("N")

# # filenames
# # i: csv files exported from HOBO downloads
# HOBO_fn <- paste0(csv_dir, c("ATS_1.11.09.18.csv",
#                              "ATS_2.17.09.18.csv",
#                              "ATS_1.24.09.18.csv",
#                              "ATS_2.01.10.18.csv",
#                              "ATS_1.08.10.18.csv",
#                              "ATS_2.15.10.18.csv",
#                              "ATS_1.22.10.18.csv"))
#
# # i: spreadsheets
# sswc_fn <- paste0(xlsx_dir, "ATS Treatment.xlsx")  # water chemistry and biomass
# ssbc_fn <- paste0(xlsx_dir, "Analysen ZEA-3.xlsx") # of biomass composition
#
# # o: csv files
# pre <- "FZJAP2018 "
# csv <- ".csv"
# dataset <- list("ti" = "Temperature and Irradiance",
#                 "wc" = "Water Chemistry",
#                 "bp" = "Biomass Productivity",
#                 "bc" = "Biomass Composition")
#
# ti_fn <- paste0(csv_dir, pre, dataset[["ti"]], csv)     # temperature and irradiance
# wc_fn <- paste0(csv_dir, pre, dataset[["wc"]], csv)     # water chemistry
# bp_fn <- paste0(csv_dir, pre, dataset[["bp"]], csv)     # biomass productivity
# bc_fn <- paste0(csv_dir, pre, dataset[["bc"]], csv)     # biomass composition
#
