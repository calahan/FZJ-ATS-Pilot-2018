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
    tb <- read_csv(paste0(csvdir, fn), skip=2, col_names=c(
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
    tb <- HOBOData(fn, view=FALSE, ret=TRUE, week=w)
    row_ct <- nrow(tb)
    HOBOData(fn, trim=c(l, r))
    HOBOData(fn, trim=c(l, r), wind=c(l-50, l+50))
    HOBOData(fn, trim=c(l, r), wind=c(row_ct-r-50, row_ct-r+50))
    return(tb[l:(row_ct-r),])
}

# File names
csvfiles <- c("ATS_1.11.09.18.csv",
           "ATS_2.17.09.18.csv",
           "ATS_1.24.09.18.csv",
           "ATS_2.01.10.18.csv",
           "ATS_1.08.10.18.csv",
           "ATS_2.15.10.18.csv",
           "ATS_1.22.10.18.csv")

chemfile <- paste0(ssdir, "ATS Treatment.xlsx")
