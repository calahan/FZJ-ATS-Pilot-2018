# Functions to support analysis of data from FZ-J WWTP ATS Pilot
HOBOData <- function(fn, zone="Europe/Berlin", trim=c(NA, NA), wind=c(NA, NA), view=TRUE, text=FALSE, ret=FALSE, week=NA) {
#
# i:
# fn        Name of the HOBO csv file to load
# trim      Number of observations to remove from beginning and end of data
# wind      Beginning and end indices of data to plot (if view==TRUE)
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

    # Adjust trim[2] to reflect actual index, not number of observations to remove
    if(is.na(trim[1])) {
        trim <- c(0, row_ct)
    } else {
    trim <- c(trim[1], row_ct - trim[2])
    }

    if(is.na(wind[1])) {
        wind <- c(1, row_ct)
    }

    if(text) {
        cat(paste0(fn, ": ", nrow(tb), " ", "observations\n"))
    }

    if(view) {
        par(mfrow=c(1,2))
        par(ps = 12, cex = 1, cex.main = 1)

        plot(x=tb[wind[1]:wind[2],]$datetime,
             y=tb[wind[1]:wind[2],]$lux,
             type="l",
             main=paste0("lux ", fn),
             xlab="Date",
             ylab="Lux")

        abline(v=c(tb[trim[1],]$datetime,
                   tb[trim[2],]$datetime),
               col=c("blue", "blue"),
               lty=c(1, 1),
               lwd=c(1, 1))

        plot(x=tb[wind[1]:wind[2],]$datetime,
             y=tb[wind[1]:wind[2],]$temp,
             type="l",
             main=paste0("temp ", fn),
             xlab="Date",
             ylab="Temperature")

        abline(v=c(tb[trim[1],]$datetime,
                   tb[trim[2],]$datetime),
               col=c("blue", "blue"),
               lty=c(1, 1),
               lwd=c(1, 1))

    }

    if(!is.na(week)) {
        tb$week <- week
    }

    if(ret) {
        return(tb[trim[1]:trim[2],])
    }
}

CleanHOBOData <- function(fn, l, r, w) {
#
# i:
# fn    name of the csv file to clean
# l     index of the first observaton to keep
# r     index of the last observation to keep
# w     2-element vector of the window to display
#
# se:
# Three pairs of plots are produced: one pair showing all data with vertical lines
# at the cutoff observations, one pair zoomed in to the early cutoff zone, and one
# pair zoomed in to the late cutoff zone
#
# v:
# A tibble of the cleaned data
#
    tb <- HOBOData(fn, view=FALSE, ret=TRUE, week=w)
    row_ct <- nrow(tb)
    HOBOData(fn, trim=c(l, r))
    HOBOData(fn, trim=c(l, r), wind=c(l-50, l+50))
    HOBOData(fn, trim=c(l, r), wind=c(row_ct-r-50, row_ct-r+50))
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
