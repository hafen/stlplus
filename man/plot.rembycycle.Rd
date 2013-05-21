\name{plot.rembycycle}
\alias{plot.rembycycle}
\title{Plot of Remainder Component by Cycle-Subseries}
\description{
Plots the remainder component by cycle-subseries with a loess line.
}
\usage{
plot.rembycycle(x, col = "black", locol = "red", lolwd = 2, 
   xlab = "Time", ylab = "Remainder", ...)
}
\arguments{
   \item{x}{object of class \code{"stl2"}.}
   \item{col, locol, lolwd, xlab, ylab, ...}{parameters to be passed to \code{xyplot()}. \code{locol} and \code{lolwd} are the line color and width for the loess line.}
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
