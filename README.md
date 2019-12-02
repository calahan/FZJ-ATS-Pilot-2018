# FZ-J-WWTP-ATS-Pilot  

Author: Dean Calahan

This repository contains code and data for creating cleaned data sets from the algal
turf scrubbing (ATS) pilot operated by Forschungszentrum-Jülich's (FZ-J) IBG-2/Alternative
Biomass department at FZ-J's campus wastewater treatment plant from August to October
2018 and from July-November 2018, and in the Novagreen Greenhouse from November-June 2019.

The input data comprises spreadsheets and csv files. These were copied from their locations on the FZ-J server K:\\Projekte\ATS\Data\WWTP Pilot\ to Research/. See the table for the file descriptions.

The code comprises `FZJAP2018.R`, a collection of functions for creating the data
sets, `Settings.R`, a collection of variable definitions used by the code, and `FZJAP2018Data.Rmd`.

Knit `FZJAP2018Data.Rmd` to clean the data and create the documentation.

# File Descriptions
Path                |Filename                                    |Description
--------------------|--------------------------------------------|--------------------------
2018/Analysen ZEA-3/|Analyzen ZEA-3.xlsx                         |Elemental analysis (ICPOES)
&nbsp;              |Elemental Analysis for C, H, N, S and O.docx|Elemental analysis (#)
2018/               |ATS Treatment_last version.xlsx             |Water analysis
2018/HOBOware       |*.xslx                                      |Temperature & lux
&nbsp;              |Elemental Analysis for C, H, N, S and O.docx|
2019/               |Treatment data.xlsx                         |
2019/               |Weight data.xlsx                            |Biomass fresh & dry weight
2019/HOBOware       |*.xslx                                      |Temperature & lux
2019/ZEA-3/191206   |191206_Ergebnisse_ICPOES.xlsx               |Elemental analysis (ICPOES)
2019/ZEA-3/191206   |auf191206_EA.xlsx                           |Elemental analysis (#)
