\name{plot.seasonal}
\alias{plot.seasonal}
\title{Seasonal Diagnostic Plot for an stl2 Object}
\description{
Plots each cycle-subseries of the detrended data (or equivalently, 
seasonal plus remainder), with the mean subtracted.  The fitted seasonal 
component is plotted as a line through the points.
}
\usage{
plot.seasonal(x, col = c("black", "red"), lwd = 2, xlab = "Time", 
    ylab = "Centered Seasonal + Remainder", ...)
}
\arguments{
  \item{x}{object of class \code{"stl2"}.}
  \item{col, lwd, xlab, ylab, ...}{parameters to be passed to \code{xyplot()}.}
}
\details{
Helps decide how much of the variation in the data other than the trend should go into the seasonal component, and how much in the remainder.
}
\value{
object of class \code{"trellis"}.
}
\references{
R. B. Cleveland, W. S. Cleveland, J.E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
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
