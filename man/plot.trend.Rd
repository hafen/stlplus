\name{plot.trend}
\alias{plot.trend}
\title{Trend Diagnostic Plot for an stl2 Object}
\description{
Plots the trend+remainder with the trend component overlaid, and the 
remainder component in a separate panel.
}
\usage{
plot.trend(x, xlab = "Time", ylab = "Trend", span = 0.3, 
   type = c("p", "l"), scales = list(y = list(relation = "free")), 
   lwd = c(1, 1), col = c("darkgray", "black", "darkgray"), 
   layout = c(1, 2), between = list(y = 0.5), strip = FALSE, 
   strip.left = TRUE, as.table = TRUE, ...)
}
\arguments{
  \item{x}{object of class \code{"stl2"}.}
  \item{xlab, ylab, span, type, scales, lwd, col, layout}{parameters to be passed to 
  xyplot.}
  \item{between, strip, strip.left, as.table, ...}{parameters to be passed to 
xyplot.}
}
% \details{
% }
\value{
object of class \code{"trellis"}.
}
\references{
R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
}
\author{Ryan Hafen}
% \note{ ~~further notes~~ 
%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{stl2}}}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--   or do  help(data=index)  for the standard data sets.
% 
% }
