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

# Functions
# BiomassCompositionData : function ()
# CleanHOBOData : function (fn, l, r, w)
# HOBOData : function (fn, zone = "Europe/Berlin", trim = c(NA, NA), window = c(NA, NA), view = TRUE, text = FALSE, ret = FALSE, week = NA)
# IntegrateObs : function (d, df1, df2, diag = FALSE)
# MeanBiomassComposition : function (df)
# PlotBiomassCompositionData : function (df)
# PlotMeanBiomassCompositionData : function (df)
# PlotMolarRatios : function (df)
# WaterChemistryBiomass : function ()
# WaterChemistryBiomassDupes : function ()

# Load biomass composition data from spreadsheet
# i: args   List of function arguments maintained in Setup.R
# v:        Tidy data frame of biomass composition data
# s:        NA
BiomassCompositionData <- function(args) {
    # Load ss.
    df <- read_excel(path =      args[["fn"]],
                     sheet =     args[["sheet"]],
                     skip=       args[["skip"]],
                     col_names = args[["col_names"]])

    # Tidy ss data.
    tprop <- df %>%
        gather(key = atom, value = proportion, Smw, Pmw, Kmw, Camw, Mgmw, Mnmw, Cmw, Nmw) %>%
        select(id, date, atom, proportion)
    tprop$atom <- sub("(.+)mw$" , "\\1", tprop$atom)

    tsd <- df %>%
        gather(key = atom, value = sd, Ssd, Psd, Ksd, Casd, Mgsd, Mnsd, Csd, Nsd) %>%
        select(id, date, atom, sd)
    tsd$atom <- sub("(.+)sd$" , "\\1", tprop$atom)

    # $eplace < 0,1 with NA
    tsd[with(tsd, which(sd == "< 0,1")), ]$sd <- NA
    tsd$sd <- as.numeric(tsd$sd)

    return(plyr::join(tprop, tsd, c("id", "date", "atom")))
}

# Load and optionally trim or plot data from a HOBO data logger.
# i: fn        Name of HOBO csv file to load
#    zone      Time zone of HOBO readings
#    view      Whether to plot the data
#    window    Beginning and end indices of data to plot (if view==TRUE)
#    text      Display console text with HOBO observation count?
#    ret       Whether to return the trimmed dataset
#    trim      Number of observations to remove from beginning and end before returning data
#    week      Week index represented by fn
# v:           If ret == TRUE, a tibble containing the HOBO data from the specified file, with first
#              and last data elided according to trim.
# s:           A plot is generated if view == TRUE, text is displayed if text == TRUE
#
HOBOData <- function(fn, zone="Europe/Berlin", trim=c(NA, NA), window=c(NA, NA), view=TRUE, text=FALSE, ret=FALSE, week=NA) {
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

    # subset of data to return or plot
    if(is.na(trim[1])) {
        trim <- c(0, row_ct)
    } else {
        trim <- c(trim[1], row_ct - trim[2])
    }

    # Output diagnostic text if text == TRUE
    if(text) {
        cat(paste0(fn, ": ", nrow(tb), " ", "observations\n"))
    }

    # Plot the data if view == TRUE
    if(view) {
        if(is.na(window[1])) {
            window <- c(1, row_ct)
        }
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

    # If week index was supplied, create a column for it
    if(!is.na(week)) {
        tb$week <- week
    }

    # Return the data defined by window.
    if(ret) {
        return(tb[trim[1]:trim[2],])
    }
}

# i: fn     name of the csv file to clean
#    l      index of the first observaton to keep
#    r      index of the last observation to keep
#    w      2-element vector of the window to plot if view == TRUE
# v:        Cleaned data tibble
# s:        Three pairs of plots: one showing all data with vertical lines at clean
#           data boundaries, one pair zoomed in to the early cutoff zone, and one
#           pair zoomed in to the late one
#
CleanHOBOData <- function(fn, l, r, w) {
    tb <- HOBOData(fn, view=FALSE, ret=TRUE, week=w)
    row_ct <- nrow(tb)
    HOBOData(fn, trim=c(l, r))

    # Uncomment the next two lines to print higher detail
    # HOBOData(fn, trim=c(l, r), window=c(l-50, l+50))
    # HOBOData(fn, trim=c(l, r), window=c(row_ct-r-50, row_ct-r+50))
    return(tb[l:(row_ct-r),])
}

# Create data frames of the data dictionaries
# i: bc     biomass composition data frame
#    bp     biomass productivity data frame
#    ti     temperature and illumination data frame
#    wc     water chemistry data frame
# v:        list of data dictionary data frames
# s:        NA
DataDictionaries <- function(bc, bp, ti, wc) {

}

# Create a data frame representing integrated temperature x irradiance.
# i: d      the observation date
#    df1    data frame containing observations of one variable (e.g. temp)
#    df2    data frame containing observations of the other variable (e.g. irrad)
#    diag   whether to display diagnostic information
# v:        The sum of the integrated units (e.g. degree minutes) inclusive of the
#           first and last measurements for that day.
# s:        NA
#
IntegrateObs <- function(d, df1, df2, diag=FALSE) {
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
# [todo]
# i: NA
# v: A cleaned data frame containing results from water chemistry analsysis
# s: Creation of the csv file "FZJ WWTP ATS Pilot Chemistry and Biomass.csv"
#
WaterChemistryBiomass <- function(ti_df) {
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


MeanBiomassComposition <- function(df) {
    # return(df %>% summarize(Cmean = mean(Cmw, na.rm=TRUE), Csd = sd(Cmw, na.rm=TRUE),
    #                         Nmean = mean(Nmw), Nsd = sd(Nmw),
    #                         Mgmean = mean(Mgmw), Mgsd = sd(Mgmw),
    #                         Pmean = mean(Pmw), Psd = sd(Pmw),
    #                         Smean = mean(Smw), Ssd = sd(Smw),
    #                         Kmean = mean(Kmw), Ksd = sd(Kmw),
    #                         Camean = mean(Camw), Casd = sd(Camw),
    #                         Mnmean = mean(Mnmw), Mnsd = sd(Mnmw)))

    return(df %>%
               group_by(atom) %>%
               summarize(mean = mean(proportion),
                         sd = sd(proportion)))
}

PlotBiomassCompositionData <- function(df) {
    atoms <- sort(unique(df$atom))
    #    masses <- setNames(mass(atoms), atoms)

    plot_df <- df %>%
        gather_(key="reading", value="value", gather_cols=c("proportion")) %>%
        select(date, atom, value, sd)


    # Plot weight percents of biomass composition readings with their standard deviations.
    return(ggplot(plot_df, aes(x=date, y=value, fill=atom)) +
               geom_bar(stat="identity", position=position_dodge()) +
               geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                             position=position_dodge(),
                             na.rm=TRUE))
}

# Plot means+sd of biomass composition.
# i: df     Biomass composition data frame
# v:        The plot
# s:        NA
PlotMeanBiomassCompositionData <- function(df) {
    mbc_df <- MeanBiomassComposition(df)

    # val_cols <- c("Cmean", "Nmean", "Mgmean", "Pmean", "Smean", "Kmean", "Camean", "Mnmean")
    # sd <- with(mean_df, c(Csd, Nsd, Mgsd, Psd, Ssd, Ksd, Csd, Mnsd))
    # sel_cols <- c("reading", "val")

    # mbc_df <- mbc_df %>%
    #     gather_(key="reading", value="value", c("mean")) %>%
    #     select(sel_cols)

    # ggplot orders x axis alphabetically for labels, but with factors we can order at will
    # get order of readings
    mean_order <- order(mbc_df$mean, decreasing = TRUE)
    atoms <- mbc_df$atom[mean_order]
    mbc_df$reading <- factor(mbc_df$atom, levels=atoms)

    return(ggplot(mbc_df, aes(x=reading, y=mean)) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                      position=position_dodge()))
}

PlotMolarRatios <- function(df) {
    # Molar ratios
    mbc_df <- MeanBiomassComposition(df)
    # atoms <- c("C", "N", "Mg", "P", "S", "K", "Ca", "Mn")
    #atoms <- mbc_df$atom
    mbc_df <- mbc_df %>%
        mutate(mass = mass(mbc_df$atom)) %>%
        mutate(molar = mean/mass)

    molrat_p <- (mbc_df %>% filter(atom == "P"))$molar
    mbc_df <- mbc_df %>%
        mutate(molarP = molar / molrat_p) %>%
        mutate(molarPsd = sd * molar * molrat_p)

    # molar_v <- with(mean_df,
    #                  c(Cmean, Nmean, Mgmean, Pmean, Smean, Kmean, Camean, Mnmean)/
    #                      sort(masses))
    # molar_v <- molar_v/molar_v[["P"]]

    # molar_df <- data.frame(atom = names(molar_v), val = molar_v)

    # ggplot orders x axis alphabetically for labels, but with factors we can order at will
    # get order of readings
    reading_order <- order(mbc_df$molarP, decreasing = TRUE)
    atoms <- mbc_df$atom[reading_order]
    means <- mbc_df$molarP[reading_order]
    sds <- mbc_df$molarPsd[reading_order]

    plot_df <- data.frame(atom = factor(atoms, levels=atoms),
                          mean = mbc_df[reading_order,]$molarP,
                          sd = mbc_df[reading_order,]$molarPsd)

    return(ggplot(plot_df, aes(x=atom, y=mean)) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                      position=position_dodge()))
}