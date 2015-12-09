#' Lattice plot of the raw, seasonal, trend, and remainder components
#'
#' Lattice plot of the raw, seasonal, trend, and remainder components.  If post-trend smoothing was done, these components will be plotted instead of the trend component.
#'
#' @param x object of class \code{"stlplus"}.
#' @param scales,type,as.table,strip,strip.left,between,layout,\ldots parameters to be passed to xyplot.
#' @return object of class \code{"trellis"}.
#'
#' @references R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @seealso \code{\link{stlplus}}
#' @export
#' @importFrom lattice xyplot panel.abline panel.xyplot panel.loess panel.segments
plot.stlplus <- function(x,
  scales = list(y = list(relation = "sliced")),
  type = "l", as.table = TRUE, strip = FALSE, strip.left = TRUE,
  between = list(y = 0.5), layout = NULL, ...) {

  if(x$pars$fc.number == 0) {
    d <- data.frame(
      time = rep(time.stlplus(x), 4),
      values = c(getraw(x), seasonal(x), trend(x), remainder(x)),
      ind = factor(rep(c(1:4), each = x$n)))
    levels(d$ind) <- c("raw", "seasonal", "trend", "remainder")
    d$which <- "isnotNA"

    if(is.null(layout)) layout <- c(1, 4)

    # if there are any NA values, plot the part of the seasonal and trend
    # components with a different color
    if(any(is.na(getraw(x)))) {
      d$values[d$ind == "seasonal" & is.na(getraw(x))] <- NA
      d$values[d$ind == "trend" & is.na(getraw(x))] <- NA
      d <- rbind(d,
        data.frame(
          time = rep(time.stlplus(x), 2),
          values = c(seasonal(x), trend(x)),
          ind = rep(c("seasonal", "trend"), each = x$n),
          which = "isxNA")
      )
      d$values[d$ind == "seasonal" & d$which == "isxNA" & !is.na(getraw(x))] <- NA
      d$values[d$ind == "trend" & d$which == "isxNA" & !is.na(getraw(x))] <- NA
    }

  } else {
    fc.number <- x$pars$fc.number
    nvar <- 3 + fc.number
    fc.name <- as.character(x$pars$fc$fc.name)
    if(fc.number == 1) {
      fcdat <- x$fc[,1]
    } else {
      fcdat <- stack(x$fc[,1:fc.number])$values
    }

    d <- data.frame(
      time = rep(time.stlplus(x), nvar),
      values = c(getraw(x), seasonal(x), fcdat, remainder(x)),
      ind = factor(rep(c(1:nvar), each = x$n)))
    levels(d$ind) <- c("raw", "seasonal", fc.name, "remainder")
    d$which <- "isnotNA"

    if(is.null(layout)) layout <- c(1, nvar)

    # if there are any NA values, plot the part of the seasonal and trend
    # components with a different color
    if(any(is.na(getraw(x)))) {
      d$values[d$ind == "seasonal" & is.na(getraw(x))] <- NA
      for(i in 1:fc.number) {
        d$values[d$ind == fc.name[i] & is.na(getraw(x))] <- NA
      }
      d <- rbind(d,
        data.frame(
          time = rep(time.stlplus(x), fc.number + 1),
          values = c(seasonal(x), fcdat),
          ind = rep(c("seasonal", fc.name), each = x$n),
          which = "isxNA")
      )
      d$values[d$ind == "seasonal" & d$which == "isxNA" & !is.na(getraw(x))] <- NA
      for(i in 1:fc.number) {
        d$values[d$ind == fc.name[i] & d$which == "isxNA" & !is.na(getraw(x))] <- NA
      }
    }
  }

  p <- lattice::xyplot(values ~ time | ind,
    data = d,
    groups = which,
    type = type,
    layout = layout,
    scales = scales,
    as.table = as.table,
    strip = strip, strip.left = strip.left,
    between = between,
    ...
  )
  p
}

.midmean <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  isInner <- x < q[2] & x > q[1]
  mean(x[isInner], na.rm = TRUE)
}

#' Cycle-Subseries Plot for an stlplus Object
#'
#' Plots the seasonal component by cycle-subseries, with lines emanating
#' from the midmean of the values within each cycle-subseries.
#'
#' @param x object of class \code{"stlplus"}.
#' @param layout,col,xlab,ylab,panel,\ldots parameters to be passed to \code{xyplot()}.
#' @return object of class \code{"trellis"}.
#' @references R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @seealso \code{\link{stlplus}}
#' @export
plot_cycle <- function(x, layout = c(x$pars$n.p, 1),
  col = "#0080ff", xlab = "Time", ylab = "Seasonal",
  panel= function(x, y, ...) {
    lattice::panel.segments(x, rep(.midmean(y), length(x)), x, y, col = col)
  }, ...) {

  seas <- seasonal(x)
  t <- time.stlplus(x)

  cycleSubIndices <- x$data$sub.labels

  if(x$pars$periodic) {
    p <- lattice::xyplot(seas ~ t | cycleSubIndices,
      layout = layout,
      type = "l",
      xlab = xlab,
      ylab = ylab,
      ...
    )
  } else {
    p <- lattice::xyplot(seas ~ t | cycleSubIndices,
      layout = layout,
      type = "l",
      panel = panel,
      xlab = xlab,
      ylab = ylab,
      ...
    )
  }

  p
}

#' Seasonal Diagnostic Plot for an stlplus Object
#'
#' Plots each cycle-subseries of the detrended data (or equivalently,
#' seasonal plus remainder), with the mean subtracted.  The fitted seasonal
#' component is plotted as a line through the points.
#'
#' @param x object of class \code{"stlplus"}.
#' @param col,lwd,xlab,ylab,\ldots parameters to be passed to \code{xyplot()}.
#'
#' @details Helps decide how much of the variation in the data other than the trend should go into the seasonal component, and how much in the remainder.
#' @return object of class \code{"trellis"}.
#' @references R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @seealso \code{\link{stlplus}}
#' @export
plot_seasonal <- function(x, col = c("darkgray", "black"),
  lwd = 2, xlab = "Time", ylab = "Centered Seasonal + Remainder", ...) {

  dat <- by(
    data.frame(
      v1 = seasonal(x) + remainder(x),
      v2 = seasonal(x), t = time.stlplus(x),
      v3 = x$data$sub.labels), list(x$data$sub.labels),
    function(dd) {
      mn <- mean(dd$v1, na.rm = TRUE)
      data.frame(a = dd$v1 - mn, b = dd$v2 - mn, t = dd$t, sub.labels = dd$v3)
    }
  )

  dat <- do.call(rbind, dat)

  p <- lattice::xyplot(a + b ~ t | sub.labels, data = dat,
    type = c("p", "l"),
    col = col,
    lwd = lwd,
    distribute.type = TRUE,
    as.table = TRUE,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  p
}

#' Plot of Remainder Component by Cycle-Subseries
#'
#' Plots the remainder component by cycle-subseries with a loess line.
#'
#' @param x object of class \code{"stlplus"}.
#' @param col,locol,lolwd,xlab,ylab,\ldots parameters to be passed to \code{xyplot()}. \code{locol} and \code{lolwd} are the line color and width for the loess line.
#'
#' @return object of class \code{"trellis"}.
#' @references R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @seealso \code{\link{stlplus}}
#' @export
plot_rembycycle <- function(x, col = "darkgray", locol = "black",
  lolwd = 2, xlab = "Time", ylab = "Remainder", ...) {

  vals2 <- data.frame(
    values = remainder(x),
    ind = x$data$sub.labels,
    t = time.stlplus(x))

  vals2 <- vals2[!is.na(vals2$values),]

  p <- lattice::xyplot(values ~ t | ind,
    data = vals2,
    type = "p",
    panel = function(x, y, ...) {
      lattice::panel.abline(h = 0, lty = 2, col = "darkgray")
      lattice::panel.xyplot(x, y, ...)
      lattice::panel.loess(x, y, col = locol, lwd = lolwd)
    },
    xlab = xlab,
    ylab = ylab,
    col = col,
    as.table = TRUE,
    ...
  )

  p
}

#' Trend Diagnostic Plot for an stlplus Object
#'
#' Plots the trend+remainder with the trend component overlaid, and the
#' remainder component in a separate panel.
#'
#' @param x object of class \code{"stlplus"}.
#' @param xlab,ylab,span,type,scales,lwd,col,layout parameters to be passed to
#'   xyplot.
#' @param between,strip,strip.left,as.table,\ldots parameters to be passed to
#' xyplot.
#'
#' @return object of class \code{"trellis"}.
#' @references R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @seealso \code{\link{stlplus}}
#' @export
plot_trend <- function(x, xlab = "Time", ylab = "Trend",
  span = 0.3, type = c("p", "l"),
  scales = list(y = list(relation = "free")),
  lwd = c(1, 1), col = c("darkgray", "black", "darkgray"),
  layout = c(1, 2), between = list(y = 0.5),
  strip = FALSE, strip.left = TRUE, as.table = TRUE, ...) {

  dat <- rbind(
    data.frame(
      x = time.stlplus(x),
      y = remainder(x) + trend(x),
      type = "p",
      pan = "Trend"),
    data.frame(
      x = time.stlplus(x),
      y = trend(x),
      type = "l",
      pan = "Trend"),
    data.frame(
      x = time.stlplus(x),
      y = remainder(x),
      type = "p",
      pan = "Remainder"),
    data.frame(
      x = time.stlplus(x),
      y = predict(loess(remainder(x) ~ c(1:length(remainder(x))),
        span = span, weights = x$data$weights),
        newdata = c(1:length(remainder(x)))),
      type = "l",
      pan = "Remainder")
  )

  p <- lattice::xyplot(y ~ x | pan,
    groups = type,
    data = dat,
    type = type,
    col = col,
    as.table = as.table,
    layout = layout,
    lwd = lwd,
    scales = scales,
    distribute.type = TRUE,
    between = between,
    strip = strip, strip.left = strip.left,
    xlab = xlab, ylab = ylab,
    ...
  )

  p
}

