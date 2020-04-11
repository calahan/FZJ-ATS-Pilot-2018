# Functions
# BiomassCompositionData : function (fn)
# CheckData()
# CleanHOBOData : function (fn, l, r, w)
# HOBOData : function (fn, zone = "Europe/Berlin", trim = c(NA, NA), window = c(NA, NA), view = TRUE, text = FALSE, ret = FALSE, week = NA)
# IntegrateObs : function (d, df1, df2, diag = FALSE)
# LoadWCSpreadsheet : function(p)
# MeanBiomassComposition : function (df)
# PlotBiomassCompositionData : function (df)
# PlotMeanBiomassCompositionData : function (df)
# PlotMolarRatios : function (df)
# WaterChemistryBiomass : function ()
# WaterChemistryBiomassDupes : function ()

BiomassCompositionData2018 <- function(fn) {
# Load biomass composition data from spreadsheet
# i: args   List of function arguments maintained in Setup.R
# v:        Tidy data frame of biomass composition data
# s:        NA

    # Load spreadsheet.
    df <- read_excel(path = fn,
                     sheet = "Tabelle1",
                     skip = 18,
                     col_names = c("id",
                                   "date",
                                   "S",
                                   "sdS",
                                   "P",
                                   "sdP",
                                   "K",
                                   "sdK",
                                   "Ca",
                                   "sdCa",
                                   "Mg",
                                   "sdMg",
                                   "Mn",
                                   "sdMn",
                                   "C",
                                   "sdC",
                                   "N",
                                   "sdN"))

    # Tidy the data.
    tprop <- df %>%
        gather(key = element, value = wtprop, S, P, K, Ca, Mg, Mn, C, N) %>%
        select(id, date, element, wtprop)
    tprop$element <- sub("(.+)mw$" , "\\1", tprop$element)

    tsd <- df %>%
        gather(key = element, value = wtpropsd, sdS, sdP, sdK, sdCa, sdMg, sdMn, sdC, sdN) %>%
        select(id, date, element, wtpropsd)
    tsd$element <- sub("(.+)sd$" , "\\1", tprop$element)

    # Replace < 0,1 with NA, convert character to numeric
    tsd[with(tsd, which(wtpropsd == "< 0,1")), ]$wtpropsd <- NA
    tsd$wtpropsd <- as.numeric(tsd$wtpropsd)

    return(plyr::join(tprop, tsd, c("id", "date", "element")))
}

# This is data from the mixed biomass.
BiomassCompositionData2019 <- function(fn_icp, fn_ea) {
# Load biomass composition data from spreadsheet
# i: args   List of function arguments maintained in Setup.R
# v:        Tidy data frame of biomass composition data
# s:        NA

    # Load spreadsheet data
    df_icp <- read_excel(path = fn_icp,
                         sheet = "Bericht",
                         skip = 13,
                         col_names = c("element",
                                       "na",
                                       "percent",
                                       "sd"))

    ticp <- df_icp %>% select(element, percent, sd) %>%
        rename(pct = percent)

    df_ea <- read_excel(path = fn_ea,
                        sheet = "Ergebnis",
                        skip = 2,
                        n_max = 1,
                        col_names = c("bm",
                                      "C",
                                      "H",
                                      "N"))

    extract <- function(str) {
        t <- unlist(strsplit(str, "\\s"))
        tv <- as.numeric(sub(",", ".", t[1]))
        tsd <- as.numeric(sub(",", ".", t[3]))
        return(c(tv, tsd))
    }
    tl <- lapply(c(df_ea$C, df_ea$H, df_ea$N), extract)

    tdf <- data.frame(
        element = c("C", "H", "N"),
        pct = c(tl[[1]][1], tl[[2]][1], tl[[3]][1]),
        sd = c(tl[[1]][2], tl[[2]][2], tl[[3]][2])
    )
    return(
        ticp %>%
            bind_rows(tdf)
    )

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

# If the directory "Data/" doesn't exist, create it and copy the original data there.
CheckData <- function() {
    if(!file.exists(data_dir)) {
        # Create the directory structure
        dir.create(data_dir)
        dir.create(file.path(data_dir, "2018"))
        dir.create(file.path(data_dir, "2018/HOBO"))
        dir.create(file.path(data_dir, "2019"))
        dir.create(file.path(data_dir, "2019/HOBO"))

        # Copy spreadsheets
        file.copy(file.path(orig_data_dir, "2018", "ATS Treatment_last version.xlsx"),
                  file.path(data_dir, "2018", "Treatment.xlsx"), copy.date = TRUE)
        file.copy(file.path(orig_data_dir, "2018", "Biomasse.xlsx"),
                  file.path(data_dir, "2018", "Biomass.xlsx"), copy.date = TRUE)
        file.copy(file.path(orig_data_dir, "2018", "Analysen ZEA-3/Analysen ZEA-3.xlsx"),
                  file.path(data_dir, "2018", "Elemental Analysis.xlsx"), copy.date = TRUE)
        file.copy(file.path(orig_data_dir, "2019", "Treatment Data.xlsx"),
                  file.path(data_dir, "2019", "Treatment.xlsx"), copy.date = TRUE)
        file.copy(file.path(orig_data_dir, "2019", "Weight Data.xlsx"),
                  file.path(data_dir, "2019", "Biomass.xlsx"), copy.date = TRUE)
        file.copy(file.path(orig_data_dir, "2019", "ZEA-3", "191206_Ergebnisse_ICPOES.xlsx"),
                  file.path(data_dir, "2019", "ICPOES.xlsx"), copy.date = TRUE)
        file.copy(file.path(orig_data_dir, "2019", "ZEA-3", "auf191206_EA.XLS"),
                  file.path(data_dir, "2019", "EA.XLS"), copy.date = TRUE)

        # Copy HOBO csv files
        file.copy(Sys.glob(file.path("Research", "Data", "2018", "HOBOware", "*.xlsx")),
                  file.path("Data", "2018", "HOBO", "/"),
                  copy.date = TRUE
        )
        file.copy(Sys.glob(file.path("Research", "Data", "2019", "HOBOware", "*.xlsx")),
                  file.path("Data", "2019", "HOBO", "/"),
                  copy.date = TRUE
        )
    } else {
        cat(paste0(data_dir, " directory exists and is assumed correct"))
    }
}

CleanHOBOData <- function(fn, l, r, w) {
# i: fn     name of the csv file to clean
#    l      index of the first observaton to keep
#    r      index of the last observation to keep
#    w      2-element vector of the window to plot if view == TRUE
# v:        Cleaned data tibble
# s:        Three pairs of plots: one showing all data with vertical lines at clean
#           data boundaries, one pair zoomed in to the early cutoff zone, and one
#           pair zoomed in to the late one
#
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
# i: NA
# v: A cleaned data frame containing results from water chemistry analsysis
# s: Creation of the csv file "FZJ WWTP ATS Pilot Chemistry and Biomass.csv"
#
WaterChemistryBiomass2018 <- function(fn) {
    # Read the spreadsheet
    t_df <- read_excel(fn,
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
                                   "pH",
                                   "ZEA"))

    # Impute missing times. We take the earliest time, and simply add 24 hours to it
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

    # Did ZEA analyze this sample?
    wqb_df$ZEA <- !is.na(t_df$ZEA)

    # Need observation date separate from datetime to make life easier
    wqb_df$obsdate <- make_date(year(wqb_df$datetime), month(wqb_df$datetime), day(wqb_df$datetime))

    # Impute person (person column added in 2019)
    wqb_df$person <- "Isabel"

    return(wqb_df)
}

# Structural changes to the spreadsheet are easier to deal with in a separate function
WaterChemistryBiomass2019 <- function(fn) {
    # Read the spreadsheet
    t_df <- read_excel(fn,
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
                                   "person",
                                   "bmwet",
                                   "bmdry",
                                   "bmdrypct",
                                   "pH",
                                   "ZEA",
                                   "harvest"))

    # Impute missing times to 9:30
    t_df[which(is.na(t_df$time)),]$time <- (9.5/24)

    # There are some non-standard characters in t_df$time, impute it
    t_df$time[which(t_df$time == "?")] <- (9.5/24)

    # the above come through as fractions rather than numbers of seconds, requiring
    # shenanigans to eventually get something that will have a decent time value
    # date doesn't matter here, we will be extracting only the time.
    t_df$time <- as.POSIXct(as.Date("2019-01-01")) + (as.numeric(t_df$time) * 24 * 60 * 60)

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

    # Did ZEA analyze this sample?
    wqb_df$ZEA <- !is.na(t_df$ZEA)

    # Need observation date separate from datetime to make life easier
    wqb_df$obsdate <- make_date(year(wqb_df$datetime), month(wqb_df$datetime), day(wqb_df$datetime))

    # Person is a new column
    wqb_df$person <- t_df$person

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

# Compute mean biomass composition.
# i: df     Biomass composition data frame
# v:        Data frame with means and standard deviations
# s:        NA
MeanBiomassComposition <- function(df) {
    return(df %>%
               group_by(atom) %>%
               summarize(mean = mean(wtprop),
                         wtpropsd = sd(wtprop)))
}

# Plot biomass composition data by mass ratio
# i: df     Biomass composition data frame
# v:        The plot
# s:        NA
PlotBiomassCompositionData <- function(df) {
    atoms <- sort(unique(df$atom))
    #    masses <- setNames(mass(atoms), atoms)

    plot_df <- df %>%
        gather_(key="reading", value="value", gather_cols=c("wtprop")) %>%
        select(date, atom, value, wtpropsd)


    # Plot weight percents of biomass composition readings with their standard deviations.
    return(ggplot(plot_df, aes(x=date, y=value, fill=atom)) +
               geom_bar(stat="identity", position=position_dodge()) +
               geom_errorbar(aes(ymin=value-wtpropsd, ymax=value+wtpropsd),
                             position=position_dodge(),
                             na.rm=TRUE))
}

# Plot means+sd of biomass composition.
# i: df     Biomass composition data frame
# v:        The plot
# s:        NA
PlotMeanBiomassCompositionData <- function(df) {
    mbc_df <- MeanBiomassComposition(df)

    # ggplot orders x axis alphabetically, but with factors we order at will
    plot_order <- order(mbc_df$mean, decreasing = TRUE)
    atoms <- mbc_df$atom[plot_order]
    mbc_df$reading <- factor(mbc_df$atom, levels=atoms)

    return(ggplot(mbc_df, aes(x=reading, y=mean)) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=mean-wtpropsd, ymax=mean+wtpropsd),
                      position=position_dodge()))
}

# Convert atomic mass ratios to molar ratios and plot
# i: df     Biomass composition data frame
# v:        The plot
# s:        NA
PlotMolarRatios <- function(df) {
    # Calculate molar ratios and append them to df.
    df <- MeanBiomassComposition(df)
    df <- df %>%
        mutate(mass = mass(df$atom)) %>%
        mutate(molar = mean/mass)

    molrat_p <- (df %>% filter(atom == "P"))$molar
    df <- df %>%
        mutate(molarP = molar / molrat_p) %>%
        mutate(molarPsd = wtpropsd * molar * molrat_p)

    # ggplot orders x axis labels alphabetically, but with factors we order at will
    plot_order <- order(df$molarP, decreasing = TRUE)
    atoms <- df$atom[plot_order]
    means <- df$molarP[plot_order]
    sds <- df$molarPsd[plot_order]

    plot_df <- data.frame(atom = factor(atoms, levels=atoms),
                          mean = df[plot_order,]$molarP,
                          sd = df[plot_order,]$molarPsd)

    return(ggplot(plot_df, aes(x=atom, y=mean)) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                      position=position_dodge()))
}