#' Seasonal Decomposition of Time Series by Loess
#'
#' Decompose a time series into seasonal, trend and irregular components using \code{loess}, acronym STL.  A new implementation of STL.  Allows for NA values, local quadratic smoothing,  post-trend smoothing, and endpoint blending.  The usage is very similar to that of R's built-in \code{stl()}.
#'
#' @param x vector of time series values, in order of time.  If \code{x} is a time series object, then \code{t} and \code{n.p} do not need to be specified, although they still can be.
#' @param t times at which the time series values were observed.  Not required.
#' @param n.p periodicity of the seasonal component.  In R's \code{stl} function, this is the frequency of the time series.
#' @param s.window either the character string \code{"periodic"} or the span (in lags) of the loess window for seasonal extraction, which should be odd.  This has no default.
#' @param s.degree degree of locally-fitted polynomial in seasonal extraction.  Should be 0, 1, or 2.
#' @param t.window the span (in lags) of the loess window for trend extraction, which should be odd.  If \code{NULL}, the default, \code{nextodd(ceiling((1.5*period) / (1-(1.5/s.window))))}, is taken.
#' @param t.degree degree of locally-fitted polynomial in trend extraction.  Should be 0, 1, or 2.
#' @param l.window the span (in lags) of the loess window of the low-pass filter used for each subseries.  Defaults to the smallest odd integer greater than or equal to \code{n.p} which is recommended since it prevents competition between the trend and seasonal components.  If not an odd integer its given value is increased to the next odd one.
#' @param l.degree degree of locally-fitted polynomial for the subseries low-pass filter.  Should be 0, 1, or 2.
#' @param critfreq the critical frequency to use for automatic calculation of smoothing windows for the trend and high-pass filter.
#' @param fc.window vector of lengths of windows for loess smoothings for other trend frequency components after the original STL decomposition has been obtained.  The smoothing is applied to the data with the STL seasonal component removed.  A frequency component is computed by a loess fit with the window length equal to the first element of fc.window, the component is removed, another component is computed with the window length equal to the second element of fc.window, and so forth. In most cases, the values of the argument should be decreasing, that is, the frequency bands of the fitted components should increase. The robustness weights from original STL are used as weights in the loess fitting if specified.
#' @param fc.degree vector of degrees of locally-fitted polynomial in the loess smoothings for the frequency components specified in fc.window. Values of 0, 1 and 2 are allowed. If the length of fc.degree is less than that of fc.window, the former is expanded to the length of the latter using rep; thus, giving the value 1 specifies a degree of 1 for all components.
#' @param fc.name vector of names of the post-trend smoothing operations specified by \code{fc.window} and \code{fc.degree} (optional).
#' @param inner integer; the number of \sQuote{inner} (backfitting) iterations; usually very few (2) iterations suffice.
#' @param outer integer; the number of \sQuote{outer} robustness iterations.  Default is 0, but Recommended if outliers are present.
#' @param sub.labels optional vector of length n.p that contains the labels of the subseries in their natural order (such as month name, day of week, etc.), used for strip labels when plotting.  All entries must be unique.
#' @param sub.start which element of sub.labels does the series begin with.  See details.
#' @param zero.weight value to use as zero for zero weighting
#' @param details if \code{TRUE}, returns a list of the results of all the intermediate iterations.
#' @param s.jump,t.jump,l.jump,fc.jump integers at least one to increase speed of the respective smoother.  Linear interpolation happens between every \code{*.jump}th value.
#' @param s.blend,t.blend,l.blend,fc.blend vectors of proportion of blending to degree 0 polynomials at the endpoints of the series.
#' @param \ldots additional parameters
#' @details The seasonal component is found by \emph{loess} smoothing the seasonal sub-series (the series of all January values, \ldots); if \code{s.window = "periodic"} smoothing is effectively replaced by taking the mean. The seasonal values are removed, and the remainder smoothed to find the trend. The overall level is removed from the seasonal component and added to the trend component. This process is iterated a few times.  The \code{remainder} component is the residuals from the seasonal plus trend fit.
#'
#' Cycle-subseries labels are useful for plotting and can be specified through the sub.labels argument.  Here is an example for how the sub.labels and sub.start parameters might be set for one situation.  Suppose we have a daily series with n.p=7 (fitting a day-of-week component).  Here, sub.labels could be set to c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat").  Now, if the series starts with a Wednesday value, then one would specify sub.labels=4, since Wednesday is the fourth element of sub.labels.  This ensures that the labels in the plots to start the plotting with Sunday cycle-subseries instead of Wednesday.
#' @return
#' returns an object of class \code{"stlplus"}, containing
#' \item{data}{data frame containing all of the components: \code{raw}, \code{seasonal}, \code{trend}, \code{remainder}, \code{weights}.}
#' \item{pars}{list of parameters used in the procedure.}
#' \item{fc.number}{number of post-trend frequency components fitted.}
#' \item{fc}{data frame of the post-trend frequency components.}
#' \item{time}{vector of time values corresponding to the raw values, if specified.}
#' \item{n}{the number of observations.}
#' \item{sub.labels}{the cycle-subseries labels.}
#' @references R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @author Ryan Hafen
#' @note This is a complete re-implementation of the STL algorithm, with the loess part in C and the rest in R.  Moving a lot of the code to R makes it easier to experiment with the method at a very minimal speed cost.  Recoding in C instead of using R's built-in loess results in better performance, especially for larger series.
#' @seealso \code{\link{plot.stlplus}} for plotting the components, \code{\link{getraw}}, \code{\link{seasonal}}, \code{\link{trend}}, \code{\link{remainder}} for accessing the components.
#' @example man-roxygen/ex-stlplus.R
#' @importFrom stats frequency loess median predict quantile weighted.mean time
#' @importFrom utils head stack tail
#' @export
#' @rdname stlplus
stlplus <- function(x, t = NULL, n.p, s.window, s.degree = 1,
  t.window = NULL, t.degree = 1,
  fc.window = NULL, fc.degree = NULL, fc.name = NULL,
  l.window = NULL, l.degree = t.degree,
  s.jump = ceiling(s.window / 10),
  t.jump = ceiling(t.window / 10),
  l.jump = ceiling(l.window / 10),
  fc.jump = NULL, critfreq = 0.05,
  s.blend = 0, t.blend = 0, l.blend = t.blend,
  fc.blend = NULL, inner = 2, outer = 1,
  sub.labels = NULL, sub.start = 1, zero.weight = 1e-6,
  details = FALSE, ...) {

  UseMethod("stlplus")
}

#' @export
#' @rdname stlplus
stlplus.ts <- function(x, t = as.numeric(stats::time(x)), n.p = frequency(x),
  s.window, s.degree = 1,
  t.window = NULL, t.degree = 1,
  fc.window = NULL, fc.degree = NULL, fc.name = NULL,
  l.window = NULL, l.degree = t.degree,
  s.jump = ceiling(s.window / 10),
  t.jump = ceiling(t.window / 10),
  l.jump = ceiling(l.window / 10),
  fc.jump = NULL, critfreq = 0.05,
  s.blend = 0, t.blend = 0, l.blend = t.blend,
  fc.blend = NULL, inner = 2, outer = 1,
  sub.labels = NULL, sub.start = 1, zero.weight = 1e-6,
  details = FALSE, ...) {

  if (is.matrix(x))
     stop("only univariate series are allowed")

  if(missing(n.p)) n.p <- frequency(x)
  if(!is.null(t)) {
    if(length(t) != length(x))
      stop("t must be same length as time series")
  } else {
    t <- as.vector(stats::time(x))
  }

  stlplus.default(x, t = t, n.p = n.p,
    s.window = s.window, s.degree = s.degree,
    t.window = t.window, t.degree = t.degree,
    fc.window = fc.window, fc.degree = fc.degree, fc.name = fc.name,
    l.window = l.window, l.degree = l.degree,
    s.jump = s.jump, t.jump = t.jump, l.jump = l.jump,
    fc.jump = fc.jump, critfreq = 0.05,
    s.blend = s.blend, t.blend = t.blend, l.blend = l.blend,
    fc.blend = NULL, inner = inner, outer = outer,
    sub.labels = sub.labels, sub.start = sub.start,
    details = details, ...)
}

#' @export
#' @rdname stlplus
stlplus.zoo <- function(...) {
  stlplus.ts(...)
}

#' @export
stlplus.default <- function(x, t = NULL, n.p, s.window, s.degree = 1,
  t.window = NULL, t.degree = 1,
  fc.window = NULL, fc.degree = NULL, fc.name = NULL,
  l.window = NULL, l.degree = t.degree,
  s.jump = ceiling(s.window / 10),
  t.jump = ceiling(t.window / 10),
  l.jump = ceiling(l.window / 10),
  fc.jump = NULL, critfreq = 0.05,
  s.blend = 0, t.blend = 0, l.blend = t.blend,
  fc.blend = NULL, inner = 2, outer = 1,
  sub.labels = NULL, sub.start = 1, zero.weight = 1e-6,
  details = FALSE, ...) {

  if(missing(n.p)) stop("must specify periodicity of seasonal (either explicitly or through a time series object)")
  n.p <- as.integer(n.p)

  if(n.p < 4)
    stop(paste("Parameter n.p was set to ", n.p, ".  Must be at least 4.", sep = ""))

  # family <- ifelse(robust, "symmetric", "gaussian")

  Y <- as.vector(x)
  n <- length(Y)
  nextodd <- function(x) {
    x <- round(x)
    x2 <- ifelse(x %% 2 == 0, x + 1, x)
    # if(any(x != x2))
    #   warning("A smoothing span was not odd, was rounded to nearest odd.  Check final object parameters to see which spans were used.")
    as.integer(x2)
  }


  wincheck <- function(x) {
    x <- nextodd(x)
    if(any(x <= 0)) stop("Window lengths must be positive.")
    x
  }

  degcheck <- function(x) {
    if(! all(x == 0 | x == 1 | x == 2)) stop("Smoothing degree must be 0, 1, or 2")
  }

  get.t.window <- function(t.dg, s.dg, n.s, n.p, omega) {
    if(t.dg == 0) t.dg <- 1
    if(s.dg == 0) s.dg <- 1

    coefs_a <- data.frame(
      a = c(0.000103350651767650, 3.81086166990428e-6),
      b = c(-0.000216653946625270, 0.000708495976681902))
    coefs_b <- data.frame(
      a = c(1.42686036792937, 2.24089552678906),
      b = c(-3.1503819836694, -3.30435316073732),
      c = c(5.07481807116087, 5.08099438760489))
    coefs_c <- data.frame(
      a = c(1.66534145060448, 2.33114333880815),
      b = c(-3.87719398039131, -1.8314816166323),
      c = c(6.46952900183769, 1.85431548427732))

    # estimate critical frequency for seasonal
    betac0 <- coefs_a$a[s.dg] + coefs_a$b[s.dg] * omega
    betac1 <- coefs_b$a[s.dg] + coefs_b$b[s.dg] * omega + coefs_b$c[s.dg] * omega^2
    betac2 <- coefs_c$a[s.dg] + coefs_c$b[s.dg] * omega + coefs_c$c[s.dg] * omega^2
    f_c <- (1 - (betac0 + betac1 / n.s + betac2 / n.s^2)) / n.p

    # choose
    betat0 <- coefs_a$a[t.dg] + coefs_a$b[t.dg] * omega
    betat1 <- coefs_b$a[t.dg] + coefs_b$b[t.dg] * omega + coefs_b$c[t.dg] * omega^2
    betat2 <- coefs_c$a[t.dg] + coefs_c$b[t.dg] * omega + coefs_c$c[t.dg] * omega^2

    betat00 <- betat0 - f_c

    n.t <- nextodd((-betat1 - sqrt(betat1^2 - 4 * betat00 * betat2)) / (2 * betat00))

    n.t
  }

  y_idx <- !is.na(Y)
  noNA <- all(y_idx)

  csSmooth <- function(d) {
    nn <- length(d$y)

    # check to make sure there is at least one valid value
    if(all(is.na(d$y))) {
      a <- rep(NA, 1:(nn + 2))
    } else {
      cs.ev <- seq(1, length(d$y), by = s.jump)
      if(tail(cs.ev, 1) != nn) cs.ev <- c(cs.ev, nn)
      cs.ev <- c(0, cs.ev, nn + 1)
      # print(m)
      a <- .loess_stlplus(y = d$y, span = s.window, degree = s.degree,
        weights = d$w, m = cs.ev, noNA = noNA, jump = s.jump,
        at = c(1:(nn + 2)))
    }
    c(a, rep(NA, csLength - nn - 2))
  }

  if(is.null(l.window)) {
    l.window <- nextodd(n.p)
  } else {
    l.window <- wincheck(l.window)
  }

  if(is.null(sub.labels)) {
    sub.labels <- paste("subseries", 1:n.p)
  } else {
    if(length(sub.labels) != n.p) stop("sub.labels must be of length n.p")
    if(length(unique(sub.labels)) != n.p) stop("sub.labels must be unique")
  }
  tmp <- c(sub.start:n.p, rep(1:n.p, ceiling(n / n.p)))[1:n]
  sub.labels <- factor(tmp, labels = sub.labels)

  periodic <- FALSE
  if (is.character(s.window)) {
    if (is.na(pmatch(s.window, "periodic")))
      stop("unknown string value for s.window")
    else {
      periodic <- TRUE
      s.window <- 10 * n + 1
      s.degree <- 0
      s.jump <- ceiling(s.window / 10)
    }
  } else {
    s.window <- wincheck(s.window)
  }

  degcheck(s.degree)
  degcheck(t.degree)
  degcheck(l.degree)

  if (is.null(t.window)) {
    # t.window <- nextodd(ceiling(1.5 * n.p/(1 - 1.5 / s.window)))
    t.window <- get.t.window(t.degree, s.degree, s.window, n.p, critfreq)
  } else {
    t.window <- wincheck(t.window)
  }

  if(is.null(s.jump) || length(s.jump)==0) s.jump <- ceiling(s.window / 10)
  if(is.null(t.jump) || length(t.jump)==0) t.jump <- ceiling(t.window / 10)
  if(is.null(l.jump) || length(l.jump)==0) l.jump <- ceiling(l.window / 10)

  # cat(s.degree, " ", t.degree, " ", l.degree, "\n")
  # cat(s.window, " ", t.window, " ", l.window, "\n")

  # Trend vector - initialize to 0 or NA, depending on what's in Y
  trend <- 0

  # start and end indices for after adding in extra n.p before and after
  st <- n.p + 1
  nd <- n + n.p

  # cycleSubIndices will keep track of what part of the
  # seasonal each observation belongs to
  cycleSubIndices <- rep(c(1:n.p), ceiling(n / n.p))[1:n]

  if(any(by(Y, list(cycleSubIndices), function(x) all(is.na(x)))))
    stop("There is at least one subseries for which all values are missing.")

  C <- rep(NA, n + 2 * n.p)

  dtls <- NULL

  w <- rep(1, n)

  for(o_iter in 1:outer) {

    for(iter in 1:inner) {

      # step 1: detrending...
      Y.detrended <- Y - trend

      csLength <- ceiling(n / n.p) + 2

      # step 2: smoothing of cycle-subseries
      for(i in 1:n.p) {
        cycleSub <- Y.detrended[cycleSubIndices == i]
        subWeights <- w[cycleSubIndices == i]
        cycleSub.length <- length(cycleSub)

        cs1 <- head(cycleSubIndices, n.p)
        cs2 <- tail(cycleSubIndices, n.p)

        notEnoughData <- length(cycleSub[!is.na(cycleSub)]) < s.window / 2
        # if(notEnoughData || periodic) {
        if(periodic) {
          C[c(cs1, cycleSubIndices, cs2) == i] <- rep(weighted.mean(cycleSub,
            w = w[cycleSubIndices == i], na.rm = TRUE), cycleSub.length + 2)
        } else {
          cs.ev <- seq(1, cycleSub.length, by = s.jump)
          if(tail(cs.ev, 1) != cycleSub.length) cs.ev <- c(cs.ev, cycleSub.length)
          cs.ev <- c(0, cs.ev, cycleSub.length + 1)
          tmps <- .loess_stlplus(y = cycleSub, span = s.window, degree = s.degree,
            m = cs.ev, weights = w[cycleSubIndices == i], blend = s.blend,
            jump = s.jump,  at = c(0:(cycleSub.length + 1)))
          C[c(cs1, cycleSubIndices, cs2) == i] <- tmps
          # approx(x = cs.ev, y = tmps, xout = c(0:(cycleSub.length + 1)))$y
        }
      }

      # Step 3: Low-pass filtering of collection of all the cycle-subseries
      # moving averages
      ma3 <- c_ma(C, n.p)

      l.ev <- seq(1, n, by = l.jump)
      if(tail(l.ev, 1) != n) l.ev <- c(l.ev, n)
      L <- .loess_stlplus(y = ma3, span = l.window, degree = l.degree,
        m = l.ev, weights = w, y_idx = y_idx, noNA = noNA,
        blend = l.blend, jump = l.jump, at = c(1:n))

      # L <- predict(loess(ma3 ~ c(1:n), degree = l.degree,
      #   span = l.window / n, family = family), newdata = c(1:n))

      # Step 4: Detrend smoothed cycle-subseries
      seasonal <- C[st:nd] - L

      # Step 5: Deseasonalize
      D <- Y - seasonal

      # Step 6: Trend Smoothing

      t.ev <- seq(1, n, by = t.jump)
      if(tail(t.ev, 1) != n) t.ev <- c(t.ev, n)
      trend <- .loess_stlplus(y = D, span = t.window, degree = t.degree,
        m = t.ev, weights = w, y_idx = y_idx, noNA = noNA,
        blend = t.blend, jump = t.jump, at = c(1:n))
      # if(blend$t) {
      #   # TODO: validate blend parameters
      #   trend0 <- .loess_stlplus(y = D, span = t.window, degree = t.degree,
      #     m = t.ev, weights = w, y_idx = y_idx, noNA = noNA)
      # }

      # trend <- predict(loess(D ~ c(1:length(D)), degree = t.degree,
      #   span = t.window / length(D), family = family), newdata = c(1:length(D)))

      if(details)
        dtls <- c(dtls, list(list(trend = trend, seasonal = seasonal,
          C = C, D = D, L = L, Y.detrended = Y.detrended, weights = w)))
    }

    if(outer > 1) {
      mid1 <- floor(n / 2 + 1)
      mid2 <- n - mid1 + 1
      R <- Y - seasonal - trend
      R.abs <- abs(R)
      h <- 3 * sum(sort(R.abs)[mid1:mid2])
      h9 <- 0.999 * h
      h1 <- 0.001 * h
      w <- (1 - (R.abs / h)^2)^2
      w[R.abs <= h1] <- 1
      w[R.abs >= h9] <- 0
      w[w == 0] <- zero.weight
      w[is.na(w)] <- 1
    }
  }

  ## post-seasonal smoothing, if specified
  fc <- NULL
  fc.number <- 0
  fc.res <- NULL
  if(!is.null(fc.window)) {
    fc.number <- length(fc.window)
    fc.window <- wincheck(fc.window)
    if(is.null(fc.degree))
      fc.degree <- 1
    if(length(fc.degree) < fc.number)
      fc.degree <- c(fc.degree,
        rep(fc.degree[length(fc.degree)], fc.number - length(fc.degree)))
    fc.cumulative <- rep(0, n) # keep sum of all previous fc smoothings

    degcheck(fc.degree)

    if(is.null(fc.name))
      fc.name <- paste("fc.", fc.window, sep = "")

    if(is.null(fc.jump))
      fc.jump <- ceiling(fc.window / 10)

    if(is.null(fc.blend))
      fc.blend <- rep(t.blend, fc.number)

    if(length(fc.blend) < fc.number)
      fc.blend <- c(fc.blend, rep(0, fc.number - length(fc.blend)))

    if(length(fc.jump) < fc.number)
      fc.jump <- c(fc.jump,
        rep(fc.jump[length(fc.jump)], fc.number - length(fc.jump)))

    fc.res <- data.frame(matrix(nrow = n, ncol = fc.number, data = 0))

    for(ii in 1:fc.number) {
      fc.ev <- seq(1, n, by = fc.jump[ii])
      if(tail(fc.ev, 1) != n) fc.ev <- c(fc.ev, n)
      tmp <- .loess_stlplus(y = Y - seasonal - fc.cumulative,
        span = fc.window[ii], degree = fc.degree[ii],
        m = fc.ev, weights = w, y_idx = y_idx, noNA = noNA,
        blend = fc.blend[ii], jump = fc.jump[ii], at = c(1:n))

      fc.cumulative <- fc.cumulative + tmp
      fc.res[,ii] <- tmp
    }

    if(any(is.null(fc.name)) || length(fc.name) < fc.number)
      fc.name <- paste("trend", c(1:fc.number))

    names(fc.res) <- fc.name
    fc.res$remainder <- Y - seasonal - apply(fc.res, 1, sum)

    fc <- data.frame(
      fc.window = fc.window, fc.degree = fc.degree,
      fc.name = fc.name, fc.jump = fc.jump, fc.blend = fc.blend)
  }

  # compute remainder
  R <- Y - seasonal - trend

  dat <- data.frame(raw = Y, seasonal = seasonal, trend = trend,
    remainder = R, weights = w, sub.labels = sub.labels)

  pars <- list(
    deg = data.frame(s.degree = s.degree, t.degree = t.degree, l.degree = l.degree),
    win = data.frame(s.window = s.window, t.window = t.window, l.window = l.window),
    blend = data.frame(s.blend = s.blend, t.blend = t.blend, l.blend = l.blend),
    jump = data.frame(s.jump = s.jump, t.jump = t.jump, l.jump = l.jump),
    fc.number = fc.number, fc = fc, n.p = n.p,
    inner = inner, outer = outer, periodic = periodic)

  res <- list(data = dat, fc = fc.res, pars = pars, time = t, n = nrow(dat),
    details = dtls)

  class(res) <- "stlplus"

  res
}

