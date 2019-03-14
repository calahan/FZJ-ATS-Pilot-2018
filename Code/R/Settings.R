# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Copyright 2018-2019 by Forschungszentrum-Jülich (FZ-J)
#
# This file is part of the publication "FZJ ATS Pilot 2018" (FZJAP2018), a data
# set documenting the results of an algal turf scrubbing pilot project performed
# 15.07-15.09 2018, treating secondary sewage at FZ-J's campus wastewater treatment
# plant.
#
# The software component of FZJAP2018 is open access: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# The software component of FZJAP2018 is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with the
# software component of FZJAP2018. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# libraries
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(knitr)
library(PeriodicTable)
library(ggplot2)

# directories
data_dir <- "Data/"
csv_dir <- paste0(data_dir, "csv/")
ss_dir <- paste0(data_dir, "xlsx/")

# conversion factors
#yr2d <- 365     # days per year
#yearweeks <- 52     # weeks per year
C2K <- 273.15       # convert °C to °K

# filenames
# i: csv files exported from HOBO downloads
HOBO_fn <- paste0(csv_dir, c("ATS_1.11.09.18.csv",
                             "ATS_2.17.09.18.csv",
                             "ATS_1.24.09.18.csv",
                             "ATS_2.01.10.18.csv",
                             "ATS_1.08.10.18.csv",
                             "ATS_2.15.10.18.csv",
                             "ATS_1.22.10.18.csv"))

# i: spreadsheets
sswc_fn <- paste0(ss_dir, "ATS Treatment.xlsx")  # water chemistry and biomass
ssbc_fn <- paste0(ss_dir, "Analysen ZEA-3.xlsx") # of biomass composition

# o: csv files
pre <- "FZJAP2018 "
csv <- ".csv"
dataset <- list("ti" = "Temperature and Irradiance",
                "wc" = "Water Chemistry",
                "bp" = "Biomass Productivity",
                "bc" = "Biomass Composition")

ti_fn <- paste0(csv_dir, pre, dataset[["ti"]], csv)     # temperature and irradiance
wc_fn <- paste0(csv_dir, pre, dataset[["wc"]], csv)     # water chemistry
bp_fn <- paste0(csv_dir, pre, dataset[["bp"]], csv)     # biomass productivity
bc_fn <- paste0(csv_dir, pre, dataset[["bc"]], csv)     # biomass composition

# function call parameters
bcd <- list("fn" = ssbc_fn,
            "sheet" = "Tabelle1",
            "skip" = 18,
            "col_names" = c("id",
                            "date",
                            "Smw",
                            "Ssd",
                            "Pmw",
                            "Psd",
                            "Kmw",
                            "Ksd",
                            "Camw",
                            "Casd",
                            "Mgmw",
                            "Mgsd",
                            "Mnmw",
                            "Mnsd",
                            "Cmw",
                            "Csd",
                            "Nmw",
                            "Nsd"))