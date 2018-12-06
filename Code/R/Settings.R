# Directories
datadir <- "Data/"
csvdir <- paste0(datadir, "csv/")
ssdir <- paste0(datadir, "Spreadsheets/")

# Conversion factors
yeardays <- 365     # days per year
yearweeks <- 52     # weeks per year
C2K <- 273.15       # convert °C to °K

# Input filenames
# csv files exported from HOBO downloads
HOBO_fn <- paste0(csvdir, c("ATS_1.11.09.18.csv",
                            "ATS_2.17.09.18.csv",
                            "ATS_1.24.09.18.csv",
                            "ATS_2.01.10.18.csv",
                            "ATS_1.08.10.18.csv",
                            "ATS_2.15.10.18.csv",
                            "ATS_1.22.10.18.csv"))

# Spreadsheet of water chemistry and biomass observations
ss_fn <- paste0(ssdir, "ATS Treatment.xlsx")

# Output filenames
# Cleaned temperature and illumination observations
ti_fn <- paste0(datadir, "FZJ WWTP ATS Pilot Temperature and Illumination")

# Cleaned water chemistry and biomass observations
cb_fn <- paste0(datadir, "FZJ WWTP ATS Pilot Chemistry and Biomass")