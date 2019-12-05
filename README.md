# FZ-J-WWTP-ATS-Pilot  

Author: Dean Calahan

This repository contains code and data for creating cleaned data sets from the algal turf scrubbing (ATS) pilot operated by the IBG-2/Alternative Biomass department of Forschungszentrum-Jülich (FZ-J) at the campus wastewater treatment plant from August to October 2018 and from July-November 2018, and in the Novagreen Greenhouse from November-June 2019.

The input data comprises spreadsheets and csv files.

The code comprises `FZJAP2018.R`, a collection of functions for creating the data sets, `Settings.R`, a collection of variable definitions used by the code, and `FZJAP2018Data.Rmd`.

Knit `FZJWWTPATSPilot.Rmd` to clean the data and create the documentation.

# File Descriptions
The data files are located in Data/2018/ and Data/2019/. The file descriptions follow:

Filename        |Description
----------------|--------------------------------------
ATS_*.xls       |Data downloaded from HOBO data loggers
Composition*    |Elemental analysis (2018: combined values; 2019: ICPOES and EA)
Treatment       |Water analysis results from kits
Weights         |Fresh and dry weights

# Notes
Due to a possibly neurotic desire for potentially unnecessary consistency in file naming, the names of the files from the original zip archive were changed in the following ways: for the files from the HOBO data loggers (ATS*.xls) the year 2019 was truncated to 19; the files containing the TN and P assay data were renamed to "Treatment"; the files containing weight data were renamed to "Weight"; the files containing the composition data were changed to "Composition", and for 2019 were appended with (EA) or (ICPOES). The original files, including PDFs with redundant information and information on methods, are contained in a zip archive in Research/.