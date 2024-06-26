% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IDETECT.R
\name{sol.idetect}
\alias{sol.idetect}
\title{Solution path generation via the Isolate-Detect method}
\usage{
sol.idetect(
  x,
  type = "const",
  thr_ic_cons = 0.9,
  thr_ic_lin = 1.25,
  points = 3
)
}
\arguments{
\item{x}{A numeric vector containing the data to be processed.}

\item{type}{The model type considered. \code{type = "const"}, \code{type = "lin.cont"}, \code{type = "lin.discont"} mean, respectively, that the signal (mean of \code{x}) is piecewise constant, 
piecewise linear and continuous, and piecewise linear but not necessarily continuous. If not given, the default is \code{type = "const"}}

\item{thr_ic_cons}{A positive real number with default value equal to 0.9. It is used to create the solution path for the piecewise-constant model. The lower the value, the longer the solution path.}

\item{thr_ic_lin}{A positive real number with default value 1.25. Used to create the solution path if \code{type == "lin.cont"} or \code{type == "lin.discont"}}

\item{points}{A positive integer with default value equal to 3. It defines the distance between two consecutive end- or start-points of the right- or
left-expanding intervals, as described in the Isolate-Detect methodology.}
}
\value{
An S3 object of class \code{cptpath}, which contains the following fields: 
\item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
\item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
\item{solution.set}{Empty list}
\item{x}{Input vector \code{x}}
\item{type}{The input parameter \code{type}}
\item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
\item{method}{The method used, which has value "idetect" here}
}
\description{
This function arranges all possible change-points in the mean of the input vector, or in its linear trend, in the order of importance, via the Isolate-Detect (ID) method.
It is developed to be used with the sdll and information criterion (ic) model selection rules.
}
\details{
The Isolate-Detect method and its algorithm is described in 
"Detecting multiple generalized change-points by isolating single ones", A. Anastasiou & P. Fryzlewicz (2022), Metrika, https://doi.org/10.1007/s00184-021-00821-6.
}
\examples{
r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
sol.idetect(r3)
}
\references{
A. Anastasiou & P. Fryzlewicz (2022). Detecting multiple generalized change-points by isolating single ones. \emph{Metrika}, https://doi.org/10.1007/s00184-021-00821-6.
}
\seealso{
\code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{sol.tguh}}
}
