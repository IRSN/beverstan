## ***************************************************************************
##' Find possible intersections between given periods of time or time
##' ranges.
##'
##' @title Find Possible Intersections among Periods of Time
##'
##' @param Periods Either a character vector with each of its elements
##' \code{\link{timeRange}}, or a \code{"Date"} matrix with one row
##' per period and wto colums corresponding to the beginning and the
##' end.
##' 
##' @param details Logical. If \code{TRUE} the results are detailed
##' and tell for all couples of periods if there is intersection or
##' not.
##' 
##' @return A logical with length \code{1} or a data frame with
##' one row per couple of periods.
##'
##' @export
##' 
##' @examples
##' per1 <- timeRange("1990_2000")
##' per2 <- timeRange("1970_1992")
##' per3 <- timeRange("1930_1960")
##' Periods <- rbind(per1, per2, per3)
##' doIntersect(Periods)
##' doIntersect(Periods, details = TRUE)
##' doIntersect(c("1990_2000", "1970_1992", "1930_1960"), details = TRUE)
##' 
doIntersect <- function(Periods, details = FALSE) {

    if (is.character(Periods)) {
        if (!is.vector(Periods)) {
            stop("'Periods'can be either a character vector",
                 " or a `Date`matrix")
        }
        Periods <- t(sapply(Periods, timeRange))
    }

    nP <- nrow(Periods)
    if (nP < 2) return(FALSE)
    ## See https://stackoverflow.com/a/3269471/8836534
    doIntersect2 <- function(x, y) {
        (x[1] < y[2]) && (y[1] < x[2])
    }
    ## Could be nicer with a suitable 'outer'...  
    Ans <- rep(FALSE, nP * (nP - 1) /2)
    k <- 0
    if (details) {
        if (is.null(nms <- rownames(Periods)))
            nms <- as.character(1:nrow(Periods))
        Period1 <- Period2 <- rep("", nP * (nP - 1) /2)
        for (i in 1:(nP-1)) {
            for (j in (i + 1):nP) {
                k <- k + 1
                Period1[k] <- nms[i]
                Period2[k] <- nms[j]
                Ans[k] <- doIntersect2(Periods[i, ], Periods[j, ])
            } 
        }
        return(data.frame(Period1, Period2, Ans))
    } else {
        for (i in 1:(nP-1)) {
            for (j in (i + 1):nP) {
                k <- k + 1
                Ans[k] <- doIntersect2(Periods[i, ], Periods[j, ])
            } 
        }
        return(any(Ans))
    }
    
}

## ***************************************************************************
##' Transform a year given as a string into a date corresponding
##' either the beginning of the year ot the begining of the following
##' year.
##'
##' @title Transform a Year into a Date
##'
##' @param str A character string.
##'
##' @param type The type of date to return. This has only an effect
##' when only the year is given, as oppposed to a full date. Then,
##' when \code{type} is \code{"start"}, the character \code{str} is
##' expected to describe the beginning of a time range and the
##' returned date will be the beginning of the year. When instead
##' \code{type} is \code{"end"}, the character \code{str} is expected
##' to describle the end of a time range and the returned date will be
##' the beginning of the next year. See \bold{Examples}.
##'
##' @return An object of class \code{"Date"} representing either the
##' beginning or the end of a time range.
##'
##' @section Caution: This function is not expected to work when
##' \code{str} has length \code{> 1}.
##'
##' @export
##' 
##' @examples
##' year2Date("2021")
##' year2Date("2021", type = "end")
##' 
year2Date <-  function(str, type = c("start", "end")) {
    type <- match.arg(type)
    isYear <- grepl("^\\d{1,4}$", str)
    if (isYear) {
        d <-  as.integer(str) 
        if (type == "start") d <- as.Date(sprintf("%d-01-01", d))
        else d <- as.Date(sprintf("%d-01-01", d + 1))
    } else {
        d <- as.Date(str)
        if (type == "end") {
            d <- seq(from = d, length.out = 2, by = "day")[2]
        }
    }
    d
}

## ***************************************************************************
##' Define a time range from a character string in suitable
##' format. The format is \eqn{start}\code{_}\eqn{end} where
##' \eqn{start} and \eqn{end} define the beginning and the end of the
##' time range.
##'
##' @title Define a Time Range
##'
##' @param str A character string in suitable format, see \bold{Examples}.
##' 
##' @return An object with class \code{"Date"} and length 2 containing
##' the beginning and the end. The time-range or "period" is formed by
##' all dates which are \code{>=} the beginning and are \code{<} the
##' end.
##'
##' @export
##' 
##' @examples
##' timeRange("2020_2222")
##' ## error: and underscore '_' is needed! 
##' try(timeRange("2020-2222"))
##' timeRange("2020-01-01_2222")
##' timeRange("2020-01-01_2222-09-21")
##' 
timeRange <- function(str) {
    if (!grepl("_", str)) {
        stop("'period' must be \"beginning_end\" where 'beg' and 'end' ",
             "are years or dates")
    }    
    end <- gsub("^.*_", "", str)
    beg <- gsub("_.*$", "", str)
    
    dEnd <- year2Date(end, type = "end")
    dBeg <- year2Date(beg)
    
    if (dEnd <= dBeg) stop("end <= beginning")

    ## if (!i) stop("'x' must be \"year1-year2\" where 'year1' and 'year2' ",
    ##              "are dates")
    c(dBeg, dEnd)
}

## ***************************************************************************
##' Define the cuts (dates) for the time blocks in a given time range.
##'
##' @title Define the Time Blocks of a Given Time Range
##'
##' @param timeRange A time range. Can either be a character string to
##' be passed to the \code{str} argument of the
##' \code{\link{timeRange}} function, or a \code{Date} object returned
##' by this function.
##'
##' @param blockDuration The duration of the blocks. Should be an integer
##' number of years as in \code{"2 years"}.
##' 
##' @return A vector with class \code{"Date"} giving the cuts
##' for the blocks. By dropping the last element we get the beginnings
##' of the blocks and by dropping the first element we get the ends of
##' the blocks.
##'
##' @export
##' 
##' @examples
##' tB <- blockCuts("1930_2020", block = "2 years")
##' 
blockCuts <- function(timeRange, blockDuration = "year") {
    if (is.character(timeRange)) tR <-  timeRange(timeRange)
    else tR <-  timeRange
    tB <- seq(from = tR[1], to = tR[2], by = blockDuration)
    ## tB <- tB[tB < tR[2]]
    attr(tB, "blockDuration") <- blockDuration
    tB
}

## ***************************************************************************
##' Define MAX data with time stamps that can be used with
##' time-varying extreme-value models.
##'
##' This data structure is inspired by the \code{MAX} data structure
##' of \emph{Renext}. The data relates to a certain variable which is
##' recorded in time but can be censored over periods. The data
##' describes a number of \emph{periods} in time (or \emph{time
##' ranges}, see \code{\link{timeRange}}), in general as so-called
##' \emph{historical period}. For each period, the largest
##' observations of the variable are given along with a year or date
##' at which the observation was made. The observations for a period
##' are given either as a \emph{named} numeric vector or as a
##' \code{named} list, the names being dates or years in both
##' cases. When a numeric vector is given, each element of the vector
##' is the observation for the year of date. When a list is given, a
##' list element can be either a vector of length \eqn{1} giving the
##' observation or a vector of length \eqn{2} giving the lower and
##' upper bounds of an interval in which the observation is known to
##' lie: this is \emph{interval censoring}. An upper bound can be
##' \code{Inf}. However a lower bound can not be \code{-Inf} in the
##' block maxima framework unless the period is a one-block period.
##'
##' @title Define MAX Data with Time Stamps
##' 
##' @param ... A \emph{named} collection of objects each being
##' a list or a numeric vector. See \bold{Examples}.
##'
##' @param blocks Logical. If \code{TRUE}, the results contain a
##' a list of \emph{block time series} corresponding to the periods.
##' More precisely, each row of the \code{MAXdata} data frame correspond
##' to a block, usually with a one-year duration.
##' 
##' @return An object with class \code{timeMAXdata}. This is mainly a
##' list with the following two elements.
##' \item{MAXinfo }{
##'     A data frame with columns \code{period}, \code{start},
##'     \code{end}, \code{duration}. Each row corrrespond to
##'     a period or time range with beginning, end and duration
##'     indicated by the columns.  
##' }
##' \item{MAXdata }{
##'     A data frame with columns \code{block}, \code{period},
##'     \code{date}, \code{y}, \code{yL}, \code{yU}. The column
##'     \code{block} is for compatibility reason, \code{period} links do
##'     \code{MAXinfo}.  The column \code{date} give the date of an event
##'     in the period or the begining of the year. The columns \code{y},
##'     \code{yL} and \code{yU} give the observed level for the variable
##'     of interest and the bounds when the observation is censored.
##' }
##' 
##' Unfortunately the name \emph{block} in \code{MAXdata} does not
##' correspond to a block in the block maxima meaning, but rather to a
##' period or time range.
##' 
##' @importFrom lubridate Date
##' @export
##'
##' @section Caution: Within a given period, we find observations (I)
##' for which some information is given explicitly, and observations
##' (N) for which no information explicitely given. The observations
##' (N) are assumed to be smaller than all the observations (I), hence
##' to be smaller than the smallest lower bound provided for the
##' observations (I), say \code{minyL}. As a result, each observation
##' (N) is considered as surely lying in the interval \code{(-Inf,
##' minyL)}. So as a rule, an explicit information \code{(yL, yU)}
##' should correspond to a quite large value of \code{yL} because this
##' applies to one of the (few) largest observation in the period.
##'
##' @seealso \code{\link{autoplot.timeMAXdata}} and
##' \code{\link{autolayer.timeMAXdata}}.
##' 
##' @examples
##' ## block maxima style. There should be only one observation by block (year)
##' timeMAXdata("1961_1980" = list("1961" = 7, "1973" = 40, "1979" = c(3, 30)),
##'             "1930_1960" = list("1940" = c(14, Inf), "1943" = 22),
##'             "1900_1929" = list("1910" = c(10, Inf), "1929" = c(20, 30))) 
##' 
##' ## marked process (POT) style. The observations are attached to a day,
##' ## so there can 
##' timeMAXdata("1961_1980" = list("1961-12-02" = 7, "1973-12-20" = 40,
##'                                "1976-01-21" = c(3, 30), "1979-02-19" = 27),
##'             "1930_1960" = list("1940" = c(14, Inf), "1943" = 22),
##'             "1900_1929" = list("1910" = c(10, Inf), "1929" = c(12, 30))) 
##' 
timeMAXdata <-  function(..., blocks = FALSE) {

    L <- list(...)
    nms <- names(L)

    trs <- character(0)
    d <- Date(0)
    y <- yL <- yU <- numeric(0)
    Periods <- character(0)
    BegsU <- EndsU <-  Date(0)
    nmsU <- character(0)
        
    for (i in seq_along(L)) {
        Period <- timeRange(nms[i])
        item <-  as.list(L[[i]])
        nms2 <- names(L[[i]])
        nMax <-  0
        for (j in seq_along(item)) {
            nMax <- nMax + 1
            dj <- year2Date(nms2[[j]])
            if (dj < Period[1] || dj >= Period[2]) {
                stop("bad date/block indication in time range ",
                     nms[i])
            }
            d <-  c(d, dj)
            if (length(item[[j]]) == 1) {
                y <- c(y, item[[j]][1])
                yL <- c(yL, item[[j]][1])
                yU <- c(yU, item[[j]][1])
            } else if (length(item[[j]]) == 2) {
                if (item[[j]][1] >= item[[j]][2]) {
                    stop("Please give the data in order: min and max")
                }
                y <- c(y, NA)
                yL <- c(yL, item[[j]][1])
                yU <- c(yU, item[[j]][2])
            } else stop("Bad length for the data of block ", i)
        }
        Periods <- c(Periods, rep(nms[i], nMax))
        BegsU <- c(BegsU, Period[1L])
        EndsU <- c(EndsU, Period[2L])
        nmsU <-  c(nmsU, nms[i])
    }
    
    PeriodsU <- data.frame(start = BegsU, end = EndsU)
    
    ## rownames(PeriodsU) <- nmsU
    ## print(PeriodsU)
    
    if (doIntersect(PeriodsU)) {
        print(doIntersect(PeriodsU, details = TRUE))
        stop("The given periods intersect")
    }
    
    MAXinfo <-
        data.frame(period = nmsU, PeriodsU,
                   duration = round(as.numeric(EndsU - BegsU) / 365.25,
                       digits = 2))

    MAXinfo <- MAXinfo[order(MAXinfo$start), ]
    rownames(MAXinfo) <- NULL
    
    Periods <-  factor(Periods)
    MAXdata <- data.frame(period = Periods,
                          date = d, y = y, yL = yL, yU = yU)

    MAXdata <- MAXdata[order(MAXdata$period, MAXdata$date), ]
    MAXdata <- data.frame(block = as.integer(MAXdata$period), MAXdata)
    rownames(MAXdata) <- NULL
    
    res <- list(MAXinfo = MAXinfo, MAXdata = MAXdata)
    class(res) <- "timeMAXdata"
    res
 
}

## ***************************************************************************
##' Extract information for block data.
##'
##' @title Extract Information for Block Data
##' 
##' @param object An object with class \code{timeMAXdata}. This is
##' mainly a list with two data frame elements called \code{MAXinfo}
##' and \code{MAXdata}.
##'
##' @param blockDuration The duration of the blocks. This is most
##' generally the default value of one year, but we could aslos use
##' \code{"2 years"} doe instance.
##'
##' @return A data frame with one row by block and with the variables
##' \code{y}, \code{yL} and \code{yU} giving the value of the variable
##' (if known), a lower bound and an upper bound for the value of the
##' variable.
##' 
##' @importFrom dplyr bind_rows
##' @export
##' 
##' @examples
##' ## block maxima style. There should be only one observation by block (year)
##' M <- timeMAXdata("1961_1980" = list("1961" = 7, "1973" = 40, "1979" = c(3, 30)),
##'                  "1930_1960" = list("1940" = c(10, Inf), "1943" = 22),
##'                  "1900_1929" = list("1910" = c(14, Inf), "1929" = c(20, 30))) 
##' byBlock <- byBlock(M)
##' 
byBlock <- function(object, blockDuration = "year") {

    MAXinfo <- object$MAXinfo
    MAXdata <- object$MAXdata
    nPeriods <- nrow(MAXinfo)

    for (i in 1:nPeriods) {
        p <- MAXinfo[i, "period"]
        b <- blockCuts(c(MAXinfo[i, "start"],
                         MAXinfo[i, "end"]),
                       blockDuration = blockDuration)
        nBlocks <-  length(b) - 1
        df <- subset(MAXdata, period == p)
        ## avoid warning for the all missing case
        ind <- !is.na(df$y)
        if (any(ind)) m1 <-  min(df$y[ind])
        else m1 <- Inf
        ind <- !is.na(df$yL)
        if (any(ind)) m2 <-  min(df$yL[ind])
        else m2 <- Inf
        m <-  min(m1, m2)
        BI <- data.frame(period = rep(p, nBlocks),
                         date = b[1:nBlocks],
                         y = NA, yL = NA, yU = rep(m, nBlocks))
        
        fi <- findInterval(x = df$date, b)
        if (any(duplicated(fi))) {
            stop("'Several MAXdata obs. in the same block.",
                 " Not allowed yet.")
        }
        BI[fi, "y"] <- df[ , "y"]
        BI[fi, "yL"] <- df[ , "yL"]
        BI[fi, "yU"] <- df[ , "yU"]
        if (i == 1) BlockInfo <- BI
        else BlockInfo <- dplyr::bind_rows(BlockInfo, BI)
        
    }
    
    BlockInfo <- within(BlockInfo, period <- factor(period))
    BlockInfo
    
}

##  **************************************************************************
##' Autoplot method for \code{timeMAXdata} objects representing censored
##' observations of a same variable in time. 
##' 
##' @title Autoplot Method for \code{timeMAXdata} Objects 
##'
##' @param object An object with class \code{"timeMAXdata"} representing
##' censored observations of a same variable.
##'
##' @param ... Not used yet.
##'
##' @import ggplot2
##' @export
##' @method autoplot timeMAXdata
##' 
##' @return A graphical object inheriting from \code{"ggplot"}.
##'
##' 
##' @examples
##' ## block maxima style. There should be only one observation by block (year)
##' M <- timeMAXdata("1961_1980" = list("1961" = 7, "1973" = 40, "1979" = c(3, 30)),
##'                  "1930_1960" = list("1940" = c(10, Inf), "1943" = 22),
##'                  "1900_1929" = list("1910" = c(14, Inf), "1929" = c(20, 30))) 
##' g <- autoplot(M)
##'
autoplot.timeMAXdata <- function(object, ...) {

    y <- NULL
    bB <- byBlock(object)
    
    ind <- is.na(bB$y) & is.na(bB$yL)
    bB[ind, "yL"] <- -Inf    
    
    g <- ggplot() + 
        geom_point(data = subset(bB, !is.na(y)),
                   mapping = aes_string(x = "date", y = "y", colour = "period")) +
            geom_segment(data = subset(bB, is.na(y)),
                         mapping = aes_string(x = "date", xend = "date",
                             y = "yL", yend = "yU",
                             colour = "period")) + xlab("")
    g
}

##  **************************************************************************
##' Autolayer method for \code{timeMAXdata} objects representing censored
##' observations of a same variable in time. 
##' 
##' @title Autolayer Method for \code{timeMAXdata} Objects 
##'
##' @param object An object with class \code{"timeMAXdata"} representing
##' censored observations of a same variables.
##'
##' @param which Character vector giving the layer to add: \code{"points"}
##' shows the non-censored observations as points and \code{"segments"} 
##' shows the censored observations as vertical segements.
##'
##' @param ... Not used yet.
##'
##' @importFrom ggplot2 autolayer
##' @export
##' @method autolayer timeMAXdata
##' 
##' @return A graphical layer to be added to an existing \code{ggplot}
##' object.
##' 
##' 
autolayer.timeMAXdata <- function(object, which = c("segments", "points"), ... ) {
    y <- NULL
    which <- match.arg(which)
    
    bB <- byBlock(object)
    ind <- is.na(bB$y) & is.na(bB$yL)
    bB[ind, "yL"] <- -Inf    

    if (which == "points") {
        geom_point(data = subset(bB, !is.na(y)),
                   mapping = aes_string(x = "date", y = "y", colour = "period"))
    } else {
        geom_segment(data = subset(bB, is.na(y)),
                     mapping = aes_string(x = "date", xend = "date", y = "yL", yend = "yU",
                         colour = "period"))
    }
        
}
