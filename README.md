# FZ-J-WWTP-ATS-Pilot  

Author: Dean Calahan

This repository contains code and data for creating cleaned data sets from the ATS
Pilot operated by Forschungzentrum-Jülich's (FZ-J) IBG-2/Alternative Biomass department
at FZ-J's campus wastewater treatment plant from 8/18 to 10/18.

The raw data comprises two spreadsheets, and seven csv files. The spreadsheet `ATS Treatment.xlsx`
was compiled manually by IBG-2 personnel. The spreadsheet `Analysen ZEA-3.xlsx`
was provided by the central analytical laboratory. The csv files were downloaded
from HOBO data loggers.

The code comprises `FZJWWTPP.R`, a collection of functions for creating the data
sets from the raw data, `Settings.R`, a collection of variable definitions used
by the code, and `FZ-J WWTP ATS Pilot Data Cleaning.Rmd`, the file that is knit
to create the data cleaning summary.