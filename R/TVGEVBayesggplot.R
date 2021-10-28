## ********************************************************************
##' Autoplot a \code{TVGEVBayes} object, trying to show the dependence
##' of the GEV marginal distribution on the time/date variable.
##'
##' @title Autoplot Method for \code{TVGEVBayes} Objects
##' 
##' @param object A \code{TVGEVBayes} object.
##'
##' @param ... Further arguments passed to
##' \code{\link{fitted.TVGEVBayes}}, such as \code{wich} or
##' \code{level}.
##'
##' @return A graphical object inheriting from \code{"ggplot"}.
##'
##' @seealso \code{\link{TVGEVBayes}}.
##'
##' @export
##' @method autoplot TVGEVBayes
##' 
autoplot.TVGEVBayes <- function(object, ...) {

    f <- fitted(object, ...)
    df <- data.frame(date = object$data$Date, f)

    dnm <- object$TVGEV$date
    
    g <- ggplot() + xlab("") + ylab(attr(f, "which"))
    g <- g + geom_ribbon(data = df,
                         mapping = aes_string(x = "date", ymin = "lower", ymax = "upper"),
                         fill = "SteelBlue2", alpha = 0.4)
    g <- g + geom_line(data = df,
                      mapping = aes_string(x = "date", y = "mean"),
                      colour = "NavyBlue")
    g <- g + geom_point(data = object$data,
                        mapping = aes_string(x = dnm, y = ".y0"),
                        col = "orangered")
    g <- g + geom_segment(data = object$data,
                          mapping = aes_string(x = dnm, xend = "Date",
                              y = ".yL", yend = ".yU"),
                          col = "orangered", size = 1.3)
    if (TRUE) {
        g <- g + geom_segment(data = object$data,
                              mapping = aes_string(x = dnm, xend = dnm,
                                  y = -Inf, yend = ".y0"),
                              col = "orangered", alpha = 0.4)
    }
    g 
}

## **************************************************************************
##' Autoplot the predictive distribution computed from a
##' \code{TVGEVBayes} object, i.e. the result returned by the
##' \code{predict} method for such an object.
##' 
##' @title Autoplot a \code{PredRLTVGEVBayes} object
##'
##' @param object A \code{predRLTVGEVBayes} object.
##'
##' @param ... Not used yet.
##'
##' @return A graphical object inheriting from \code{"ggplot"}.
##' 
##' @importFrom ggplot2 autoplot
##' @export
##' @method autoplot predRLTVGEVBayes
##' 
autoplot.predRLTVGEVBayes <- function(object, ...) {
    
    p <- NULL ## avoid warning
    
    aL <- attributes(object)
    ## points <- match.arg(points)
    
    g <- ggplot(data = object)
    g <-  g + geom_line(mapping = aes_string(x = "Prob", y = "Quant"))
    g <- g + scale_x_continuous(trans = bever:::.gumbel_trans_p,
                                breaks = bever:::.gumBreaks_p,
                                minor_breaks = bever:::.gumBreaks_p) +
        xlab("prob of exceedance")


    ## ========================================================================
    ## Plot points if required. These are stored in the attribute
    ## `yMax`of the prediction object.
    ## ========================================================================
    
    ## if ((points != "none") && is.null(aL$yMax)) {
    ##     warning("'object' does not have an element 'yMax' attached so no ",
    ##             "points will be shown")
    ##     points <- "none"
    ## }
    
    ## if (points != "none") {
    ##     if (aL$newDuration != aL$blockDuration) {
    ##         stop("'points' can be != \"none\" only when the attributes ",
    ##              "\"newDuration\" and \"blockDuration\" of 'object are equal")
    ##     }
    ##     nMax <- length(aL$yMax)
    ##     df <- data.frame(p = (1 - ppoints(nMax, a = a)),    
    ##                      x = sort(aL$yMax),
    ##                      source = "data")
    ##     g <- g +
    ##         geom_point(data = dplyr::filter(df, p < 1),
    ##                    mapping = aes_string(x = "p", y = "x"),
    ##                    alpha = 0.7)
    ## }
    
    g <- g + theme_gray()
    g
    
}
