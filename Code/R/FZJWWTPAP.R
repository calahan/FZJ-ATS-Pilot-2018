# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Copyright 2018-2019 by Forschungszentrum-Jülich (FZ-J)
#
# This file is part of the internal publication "FZ-J WWTP Pilot" (FZJWWTPP), a
# literate program that documents data cleaning procedures for the creation of several
# data sets arising from a project, running from July 1 through August 31 2018,
# that demonstrated operation of an algal turf scrubber at FZ-J's Jülich campus
# wastewater treatment plant.
#
# The software component of FZJWWTPP is open access: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.

# The software component of FZJWWTPP is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with the
# software component of FZJWWTPP If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Functions to support analysis of data from FZ-J WWTP ATS Pilot

HOBOData <- function(fn, zone="Europe/Berlin", trim=c(NA, NA), window=c(NA, NA), view=TRUE, text=FALSE, ret=FALSE, week=NA) {
# Load and optionally trim or plot data from a HOBO data logger.
#
# i:
# fn        Name of the HOBO csv file to load
# zone      Time zone in which the HOBO readings were recorded
# trim      Number of observations to remove from beginning and end of data
# window    Beginning and end indices of data to plot (if view==TRUE)
# view      Whether to plot the data
# text      Whether to display text with count of observations in the csv file
# ret       Whether to return the trimmed dataset
# week      Which week is represented by this file
#
# o:        If requested, a tibble containing the HOBO data from the specified file, with first
#           and last data elided.
#
# s:        A plot is generated if requested, displaying the data within the specified window
#
    tb <- read_csv(fn, skip=2, col_names=c(
        "id",
        "date",
        "time",
        "tmp",
        "lx")) %>% transmute(
            datetime = make_datetime(year(mdy(date)),
                                     month(mdy(date)),
                                     day(mdy(date)),
                                     hour(time),
                                     minute(time),
                                     second(time),
                                     tz=zone),
            temp=tmp,
            lux=lx
        )

    row_ct <- nrow(tb)

    # Adjust trim[2] to reflect the index, not the number of observations to remove
    if(is.na(trim[1])) {
        trim <- c(0, row_ct)
    } else {
    trim <- c(trim[1], row_ct - trim[2])
    }

    # If window was not supplied, set it to include the entire data set
    if(is.na(window[1])) {
        window <- c(1, row_ct)
    }

    # Output the text if desired
    if(text) {
        cat(paste0(fn, ": ", nrow(tb), " ", "observations\n"))
    }

    # Plot the data if desired
    if(view) {
        par(mfrow=c(1,2))
        par(ps = 12, cex = 1, cex.main = 1)

        # Plot illumination data
        plot(x=tb[window[1]:window[2],]$datetime,
             y=tb[window[1]:window[2],]$lux,
             type="l",
             main=paste0("lux ", fn),
             xlab="Date",
             ylab="Lux")

        # Add vertical lines to show the window
        abline(v=c(tb[trim[1],]$datetime,
                   tb[trim[2],]$datetime),
               col=c("blue", "blue"),
               lty=c(1, 1),
               lwd=c(1, 1))

        # Plot temperature data
        plot(x=tb[window[1]:window[2],]$datetime,
             y=tb[window[1]:window[2],]$temp,
             type="l",
             main=paste0("temp ", fn),
             xlab="Date",
             ylab="Temperature")

        # Add vertical lines to show the window
        abline(v=c(tb[trim[1],]$datetime,
                   tb[trim[2],]$datetime),
               col=c("blue", "blue"),
               lty=c(1, 1),
               lwd=c(1, 1))

    }

    # If week number was supplied, create a column with it
    if(!is.na(week)) {
        tb$week <- week
    }

    # Return the data within the window.
    if(ret) {
        return(tb[trim[1]:trim[2],])
    }
}

CleanHOBOData <- function(fn, l, r, w) {
#
# input:
# fn    name of the csv file to clean
# l     index of the first observaton to keep
# r     index of the last observation to keep
# w     2-element vector of the window to display
#
# side effects:
# Three pairs of plots are produced: one pair showing all data with vertical lines
# at the cutoff observations, one pair zoomed in to the early cutoff zone, and one
# pair zoomed in to the late cutoff zone
#
# value:
# A tibble of the cleaned data
#
    tb <- HOBOData(fn, view=FALSE, ret=TRUE, week=w)
    row_ct <- nrow(tb)
    HOBOData(fn, trim=c(l, r))

    # Uncomment the next two lines to print higher detail
    # HOBOData(fn, trim=c(l, r), window=c(l-50, l+50))
    # HOBOData(fn, trim=c(l, r), window=c(row_ct-r-50, row_ct-r+50))
    return(tb[l:(row_ct-r),])
}

IntegrateObs <- function(d, df1, df2, diag=FALSE) {
#
# i:
# d   the observation date
# diag  whether to print diagnostic information
#
# v:
# The number of degree minutes between the first observation time of day and the
# first observation time of the previous day.
#
    # Day 1 is an error because there was no previous sample
    if(d == 1) {
        stop("Day 1 illegal")
    }

    # If d is larger than number of observations, it is an error
    if(d > length(nutdates)) {
        stop("Day too large")
    }

    # What are the bracketing datetimes
    t2 <- min(df1[which(year(df1$datetime) == year(nutdates[d]) &
                              month(df1$datetime) == month(nutdates[d]) &
                              day(df1$datetime) == day(nutdates[d])),]$datetime)

    t1 <- min(df1[which(year(df1$datetime) == year(nutdates[d-1]) &
                              month(df1$datetime) == month(nutdates[d-1]) &
                              day(df1$datetime) == day(nutdates[d-1])),]$datetime)

    # What are the temperature and lux measurements between these brackets?
    temps <- df2[df2$datetime >= t1 & df2$datetime < t2,]$temp + C2K
    luxs <-  df2[df2$datetime >= t1 & df2$datetime < t2,]$lux

    if(diag) {
        cat(paste0("Number of observations: ", length(temps), "\n"))
    }

    return(c(sum(temps), sum(luxs)))
}

# Create a cleaned data frame from the spreadsheet "ATS Treatment"
WaterChemistryBiomass <- function() {
# i: NA
# v: A cleaned data frame containing results from water chemistry analsysis
# s: Creation of the csv file "FZJ WWTP ATS Pilot Chemistry and Biomass.csv"

    # Read the spreadsheet
    t_df <- read_excel(sswc_fn,
                       sheet="Tabelle1",
                       skip=3,
                       col_names=c("date",
                                   "time",
                                   "beforeafter",
                                   "PO4dil",
                                   "PO4P",
                                   "PO4Pxdil",
                                   "PO4Pred",
                                   "dilTN",
                                   "TNxdil",
                                   "TN",
                                   "TNred",
                                   "comment",
                                   "remark",
                                   "assaydate",
                                   "bmwet",
                                   "bmdry",
                                   "bmdrypct",
                                   "pH"))

    # Impute missing times. We take the earliest time, and simply add 24 hours to it
    # [todo] Is this even used? There is something wrong here at any rate
    no_times <- which(is.na(t_df$time))
    min_time <- min(ti_df[ti_df$week == 1,]$datetime)
    t_df[which(is.na(t_df$time)),]$time <- make_datetime(2018,
                                                         12,
                                                         31,
                                                         8,
                                                         0,
                                                         0)

    # Create a new data frame with date and time combined into datetime for each observation
    wqb_df <- data.frame(datetime=make_datetime(year(t_df$date),
                                                month(t_df$date),
                                                day(t_df$date),
                                                hour(t_df$time),
                                                minute(t_df$time),
                                                second(t_df$time),
                                                tz="Europe/Berlin"))

    # Make before/after a logical
    wqb_df$before <- TRUE
    wqb_df[which(t_df$beforeafter != "before"),]$before <- FALSE

    # Nutrient concentrations & pH
    wqb_df$PO4P <- t_df$PO4Pxdil
    wqb_df$TN <- t_df$TN
    wqb_df$pH <- t_df$pH

    # Make frozen/not frozen a logical
    wqb_df$frozen <- TRUE
    wqb_df[which(t_df$comment != "frozen"),]$frozen <- FALSE

    # Date of assay
    wqb_df$assaydate <- make_date(year(t_df$assaydate),
                                  month(t_df$assaydate),
                                  day(t_df$assaydate))

    # Biomass measurements
    wqb_df$wet_biomass <- t_df$bmwet
    wqb_df$dry_biomass <- t_df$bmdry
    wqb_df$solids <- t_df$bmdrypct

    # Need observation date separate from datetime to make life easier
    wqb_df$obsdate <- make_date(year(wqb_df$datetime), month(wqb_df$datetime), day(wqb_df$datetime))

    return(wqb_df)
}

# Create a data frame with nutrient removal metrics for samples that were assayed
# both frozen and fresh.
WaterChemistryBiomassDupes <- function() {
    dupe_df <- wcb_df[c(3:6, 13:16),] %>% select(datetime, before, PO4P, TN, frozen)
    dupedeltas_df <- data.frame(before = c(FALSE, TRUE, FALSE, TRUE),
                                PO4Pdelta = c(dupe_df[2,]$PO4P - dupe_df[1,]$PO4P,
                                              dupe_df[4,]$PO4P - dupe_df[3,]$PO4P,
                                              dupe_df[6,]$PO4P - dupe_df[5,]$PO4P,
                                              dupe_df[8,]$PO4P - dupe_df[7,]$PO4P),
                                TNdelta = c(dupe_df[2,]$TN - dupe_df[1,]$TN,
                                            dupe_df[4,]$TN - dupe_df[3,]$TN,
                                            dupe_df[6,]$TN - dupe_df[5,]$TN,
                                            dupe_df[8,]$TN - dupe_df[7,]$TN))
    dupedeltas_df$PO4Ppct <- 100*dupedeltas_df$PO4Pdelta/dupe_df[c(2, 4, 6, 8),]$PO4P
    dupedeltas_df$TNpct <- 100*dupedeltas_df$TNdelta/dupe_df[c(2, 4, 6, 8),]$TN

    return(dupedeltas_df)
}

BiomassCompositionData <- function() {
    return(read_excel(ssbc_fn,
                      sheet="Tabelle1",
                      skip=18,
                      col_names=c("id",
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
                                  "Nsd")))
}

PlotBiomassCompositionData <- function(df) {
    atoms <- c("C", "Ca", "K","Mg", "Mn", "N", "P", "S")
    masses <- setNames(mass(atoms), atoms)

    val_cols <- c("Cmw", "Nmw", "Mgmw", "Pmw", "Smw", "Kmw", "Camw", "Mnmw")
    sel_cols <- c("date", "reading", "val")

    # "NAs introduced by coercion" warning is OK.
    sd <- suppressWarnings(as.numeric(with(df, c(Csd, Nsd, Mgsd, Psd, Ssd, Ksd, Casd, Mnsd))))

    plot_df <- df %>%
        gather_(key="reading", value="val", gather_cols=val_cols) %>%
        select(sel_cols)
    plot_df$sd <- sd

    # Plot weight percents of biomass composition readings with their standard deviations.
    return(ggplot(plot_df, aes(x=date, y=val, fill=reading)) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=val-sd, ymax=val+sd),
                      position=position_dodge(),
                      na.rm=TRUE))
}

MeanBiomassComposition <- function(df) {
    return(df %>% summarize(Cmean = mean(Cmw, na.rm=TRUE), Csd = sd(Cmw, na.rm=TRUE),
                            Nmean = mean(Nmw), Nsd = sd(Nmw),
                            Mgmean = mean(Mgmw), Mgsd = sd(Mgmw),
                            Pmean = mean(Pmw), Psd = sd(Pmw),
                            Smean = mean(Smw), Ssd = sd(Smw),
                            Kmean = mean(Kmw), Ksd = sd(Kmw),
                            Camean = mean(Camw), Casd = sd(Camw),
                            Mnmean = mean(Mnmw), Mnsd = sd(Mnmw)))
}

PlotMeanBiomassCompositionData <- function(df) {
    mean_df <- MeanBiomassComposition(df)

    val_cols <- c("Cmean", "Nmean", "Mgmean", "Pmean", "Smean", "Kmean", "Camean", "Mnmean")
    atom_labs <- c("C", "N", "Mg", "P", "S", "K", "Ca", "Mn")
    sd <- with(mean_df, c(Csd, Nsd, Mgsd, Psd, Ssd, Ksd, Csd, Mnsd))
    sel_cols <- c("reading", "val")

    mean_df <- mean_df %>%
        gather_(key="reading", value="val", val_cols) %>%
        select(sel_cols)

    # ggplot orders x axis alphabetically for labels, but with factors we can order at will
    # get order of readings
    ro <- order(mean_df$val, decreasing = TRUE)
    reading <- atom_labs[ro]
    mean_df$reading <- factor(atom_labs, levels=reading)

    return(ggplot(mean_df, aes(x=reading, y=val)) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=val-sd, ymax=val+sd),
                      position=position_dodge()))
}

PlotMolarRatios <- function(df) {
    # Molar ratios
    mean_df <- MeanBiomassComposition(df)
    atoms <- c("C", "N", "Mg", "P", "S", "K", "Ca", "Mn")
    masses <- setNames(mass(atoms), atoms)

    molar_v <- with(mean_df,
                     c(Cmean, Nmean, Mgmean, Pmean, Smean, Kmean, Camean, Mnmean)/
                         sort(masses))
    molar_v <- molar_v/molar_v[["P"]]

    molar_df <- data.frame(atom = names(molar_v), val = molar_v)

    # ggplot orders x axis alphabetically for labels, but with factors we can order at will
    # get order of readings
    ro <- order(molar_df$val, decreasing = TRUE)
    atoms <- atoms[ro]
    molar_df$atom <- factor(atoms, levels=atoms)
    molar_df$val <- molar_df[ro,]$val

    return(ggplot(molar_df, aes(x=atom, y=val)) +
        geom_bar(stat="identity", position=position_dodge()))

}