% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IDETECT_seq.R
\name{sol.idetect_seq}
\alias{sol.idetect_seq}
\title{Solution path generation using the sequential approach of the Isolate-Detect method}
\usage{
sol.idetect_seq(x, points = 3)
}
\arguments{
\item{x}{A numeric vector containing the data to be processed}

\item{points}{A positive integer with default value equal to 3. It defines the distance between two consecutive end- or start-points of the right- or
left-expanding intervals, as described in the Isolate-Detect methodology.}
}
\value{
An S3 object of class \code{cptpath}, which contains the following fields: 
\item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
\item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
\item{solution.set}{Empty list}
\item{x}{Input vector \code{x}}
\item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
\item{method}{The method used, which has value "idetect_seq" here}
}
\description{
This function uses the Isolate-Detect method in its original sequential way in order to create
the solution path. It is developed to be used with the thresholding model selection rule.
}
\details{
The Isolate-Detect method and its algorithm is described in 
"Detecting multiple generalized change-points by isolating single ones", A. Anastasiou & P. Fryzlewicz (2021), Metrika, https://doi.org/10.1007/s00184-021-00821-6.
}
\examples{
r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
sol.idetect_seq(r3)
}
\references{
A. Anastasiou & P. Fryzlewicz (2021). Detecting multiple generalized change-points by isolating single ones. \emph{Metrika}, https://doi.org/10.1007/s00184-021-00821-6.
}
\seealso{
\code{\link{sol.idetect}}, \code{\link{sol.not}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{sol.tguh}},
}
