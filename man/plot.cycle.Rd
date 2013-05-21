\name{plot.cycle}
\alias{plot.cycle}
\title{Cycle-Subseries Plot for an stl2 Object}
\description{
Plots the seasonal component by cycle-subseries, with lines emanating 
from the midmean of the values within each cycle-subseries.
}
\usage{
plot.cycle(x, layout=c(x$pars$n.p, 1), col="#0080ff", 
   xlab="Time", ylab="Seasonal", 
   panel = function(x, y, ...) {
      panel.segments(x, rep(.midmean(y), length(x)), x, y, col=col)
   }, ...)
}
\arguments{
   \item{x}{object of class \code{"stl2"}.}
   \item{layout, col, xlab, ylab, panel, ...}{parameters to be passed to \code{xyplot()}.}
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
