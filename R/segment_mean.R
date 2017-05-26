#' Multiple change-point detection in the mean of a vector
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, using a user-specified method. 
#' It also estimates the constant means between each pair of neighbouring change-points.
#' It works best when the noise in the input vector is independent and identically distributed Gaussian.
#'
#' If \code{method=="tguh"}, then the change-point detection algorithm used is the 
#' Tail-Greedy Unbalanced Haar method as described in "Tail-greedy bottom-up data
#' decompositions and fast multiple change-point detection", P. Fryzlewicz (2017),
#' preprint. This paper describes two optional post-processing steps; neither of
#' them is implemented in this package.
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param method Specifies the method to be used for change-point detection. In this
#' version of the package, the only option (and the default) is \code{"tguh"}: the Tail-Greedy Unbalanced
#' Haar method, see Details for the relevant reference.
#' @param sigma The estimate or estimator of the standard deviation of the noise in \code{x}; 
#' the default is the Median Absolute Deviation of \code{x} computed under the assumption that
#' the noise is independent and identically distributed Gaussian.
#' @param th.const Tuning parameter. If \code{method=="tguh"}, then change-points are
#' estimated by connected thresholding (of the Tail-Greedy Unbalanced Haar decomposition of \code{x})
#' in which the threshold has magnitude \code{sigma * sqrt(2 * (1 + 0.01) * log(n)) * th.const},
#' where \code{n} is the length of \code{x}. The default value of \code{th.const} is 1.
#' @param p Only relevant if \code{method=="tguh"}. Specifies the number of region pairs merged 
#' in each pass through the data, as the proportion of all remaining region pairs. The default is
#' 0.01.
#' @param minseglen The minimum permitted length of each segment of constancy in the estimated mean of \code{x}; the
#' default is 1.
#' @param bal Only relevant if \code{method=="tguh"}. Specifies the minimum ratio of the length
#' of the shorter wing of each Unbalanced Haar wavelet whose coefficient survives the thresholding, to the length
#' of its support. The default is 0.05.
#' @param num.zero Numerical zero; the default is 0.00001.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}.}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{tguh.cpt}}
#' @examples
#' stairs <- rep(1:50, each=10)
#' stairs.noisy <- stairs + rnorm(500)/5
#' stairs.cleaned <- segment.mean(stairs.noisy)
#' ts.plot(stairs.cleaned$est)
#' stairs.cleaned$no.of.cpt
#' stairs.cleaned$cpt
#' @importFrom stats mad
#' @export


segment.mean <- function(x, method="tguh", sigma = stats::mad(diff(x)/sqrt(2)), th.const = 1, p = .01, minseglen = 1, bal = 1/20, num.zero = 10^(-5)) {

	res <- list(est = x, no.of.cpt = 0, cpt = integer(0))

	if (method == "tguh") res <- tguh.cpt(x, sigma, th.const, p, minseglen, bal, num.zero)

	list(est = res$est, no.of.cpt = res$no.of.cpt, cpt = res$cpt)

}
