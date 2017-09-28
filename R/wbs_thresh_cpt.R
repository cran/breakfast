#' Multiple change-point detection in the mean of a vector using the (Adaptive) WBS method, with the number of change-points chosen by thresholding
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, using the (Adaptive) Wild Binary Segmentation
#' method (see Details for the relevant literature references). The number of change-points
#' is chosen via a thresholding-type criterion. The constant means between each pair 
#' of neighbouring change-points are also estimated. The method works best when the noise in the 
#' input vector is independent and identically distributed Gaussian.
#'
#' The change-point detection algorithms used in \code{wbs.thresh.cpt} are: standard
#' Wild Binary Segmentation [see "Wild Binary Segmentation for multiple 
#' change-point detection", P. Fryzlewicz (2014), Annals of Statistics, 42, 2243-2281]
#' and Adaptive Wild Binary Segmentation [see "Data-adaptive Wild Binary Segmentation",
#' P. Fryzlewicz (2017), in preparation as of September 28th, 2017].
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param sigma The estimate or estimator of the standard deviation of the noise in \code{x}; 
#' the default is the Median Absolute Deviation of \code{x} computed under the assumption that
#' the noise is independent and identically distributed Gaussian.
#' @param universal If \code{TRUE}, then \code{M} and \code{th.const} (see below) are chosen automatically
#' in such a way that if the mean of \code{x} is constant (i.e. if there are no change-points),
#' the probability of no detection (i.e. \code{est} being constant) is approximately
#' \code{lambda}. When \code{universal} is \code{TRUE}, then \code{M=1000} for longer signals
#' and \code{M<1000} for shorter signals to avoid \code{th.const} being larger than \code{1.3},
#' which empirically appears to be too high a value. If \code{universal} is \code{FALSE}, then
#' both \code{M} and \code{th.const} must be specified.
#' @param M The number of randomly selected sub-segments of the data on which to build
#' the CUSUM statistics in the (Adaptive) Wild Binary Segmentation algorithm.
#' If you are using Adaptive Wild
#' Binary Segmentation (\code{adapt=TRUE}) and do not wish to set \code{universal}
#' to \code{TRUE} (and therefore have \code{M} chosen for you), try \code{M=1000}. If you are
#' using standard Wild Binary Segmentation (\code{adapt=TRUE}), try \code{M=20000} or higher. 
#' @param th.const Tuning parameter. Change-points are
#' estimated by thresholding [of the (Adaptive) WBS CUSUMs of \code{x}]
#' in which the threshold has magnitude \code{th.const * sqrt(2 * log(n)) * sigma},
#' where \code{n} is the length of \code{x}. There is an extra twist if \code{adapt=TRUE}, see
#' \code{th.const.min.mult} below.
#' @param th.const.min.mult If \code{adapt=TRUE}, then the threshold gradually decreases in each
#' recursive pass through the data, but in such a way that in never goes below
#' \code{th.const.min.mult * th.const * sqrt(2 * log(n)) * sigma}.
#' @param adapt If \code{TRUE} (respectively, \code{FALSE}), then Adaptive (respectively, standard)
#' Wild Binary Segmentation is used.
#' @param lambda See the description for the \code{universal} parameter above. Currently, the only
#' permitted values are \code{0.9} and \code{0.95}.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}.}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{segment.mean}}, \code{\link{wbs.bic.cpt}},
#' \code{\link{wbs.cpt}}, \code{\link{tguh.cpt}}, \code{\link{hybrid.cpt}}, \code{\link{wbs.K.cpt}}
#' @examples
#' teeth <- rep(rep(0:1, each=5), 20)
#' teeth.noisy <- teeth + rnorm(200)/5
#' teeth.cleaned <- wbs.thresh.cpt(teeth.noisy)
#' ts.plot(teeth.cleaned$est)
#' teeth.cleaned$no.of.cpt
#' teeth.cleaned$cpt
#' @importFrom stats mad
#' @export


wbs.thresh.cpt <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.825, adapt = TRUE, lambda = 0.9) {

	n <- length(x)

	if (universal) {
		u <- universal.M.th(n, lambda)
		th.const <- u$th.const
		M <- u$M
	}
	else if (is.null(M) || is.null(th.const))
		stop("If universal is FALSE, then M and th.const must be specified.")
	
	th.const.min <- th.const * th.const.min.mult
	th <- th.const * sqrt(2 * log(n)) * sigma
	th.min <- th.const.min * sqrt(2 * log(n)) * sigma

	est <- wbs.thresh.int(x, M, th, th.min, adapt)
	
	cpt <- which(abs(diff(est)) > 0)
	no.of.cpt <- length(cpt)

	list(est=est, no.of.cpt=no.of.cpt, cpt=cpt)
}
