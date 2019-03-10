# FZ-J-WWTP-ATS-Pilot  

Author: Dean Calahan

This repository contains code and data for creating cleaned data sets from the algal
turf scrubbing (ATS) pilot operated by Forschungzentrum-Jülich's (FZ-J) IBG-2/Alternative
Biomass department at FZ-J's campus wastewater treatment plant from August to October
2018.

The input data comprises two spreadsheets and seven csv files. The spreadsheet `ATS Treatment.xlsx`
was compiled manually by IBG-2 personnel. The spreadsheet `Analysen ZEA-3.xlsx`
was provided by FZ-J's central analytical laboratory. The csv files were downloaded
from HOBO data loggers using proprietary software.

The code comprises `FZJAP2018.R`, a collection of functions for creating the data
sets, `Settings.R`, a collection of variable definitions used by the code, and `FZJAP2018Data.Rmd`,
the file that is knit to clean the data and create the data documentation.