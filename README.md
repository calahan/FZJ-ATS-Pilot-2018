# FZJ ATS Pilot  

Author: Dean Calahan

This repository contains code and data for creating a tidy data set from the algal turf scrubbing (ATS) pilot operated by the IBG-2/Alternative Biomass & Biorefineries department of Forschungszentrum-Jülich (FZ-J). This pilot was performed at the campus wastewater treatment plant from August to October 2018 and from July-November 2019, in the Novagreen greenhouse from November-June 2019 and from November 2019 to January 2020.

The original non-tidy data is in the directory *Research/Data/*. Files there are duplicated under a regular naming convention in the directory *Data/*. To rebuild the data set from scratch, delete the directory *Data/* and knit *BuildDataSet.Rmd*. The directory *Data/* will be created and the original files will be copied there and the data set will be built anew.

The code comprises `FZJAP2018.R`, a collection of functions for creating the data sets, `Settings.R`, a collection of variable definitions used by the code, and `BuildDataSet.Rmd`, knitting which creates the data set

## Data Dictionary


# Notes
