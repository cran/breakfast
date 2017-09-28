#' Detecting exactly \code{K} change-points in the mean of a vector using the Adaptive WBS method
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, using the Adaptive Wild Binary Segmentation
#' method (see Details for the relevant literature reference). The number of change-points
#' is exactly \code{K}. The constant means between each pair 
#' of neighbouring change-points are also estimated. The method works best when the noise in the 
#' input vector is independent and identically distributed Gaussian. As a by-product, the function
#' also computes the entire solution path, i.e. all estimated \code{n-1} change-point locations 
#' (where \code{n} is the length of the input data) sorted from the most to the least important.
#'
#' This function should only be used if (a) you know exactly how many change-points you wish
#' to detect, or (b) you wish to order all possible change-points from the most to the least
#' important. If you need a function to estimate the number of change-points for you, 
#' try \code{\link{segment.mean}} (for a
#' default recommended estimation technique), \code{\link{wbs.thresh.cpt}}, \code{\link{wbs.bic.cpt}},
#' \code{\link{wbs.cpt}} (if you require an (Adaptive) WBS-based technique), \code{\link{tguh.cpt}}
#' (if you require a TGUH-based technique), or \code{\link{hybrid.cpt}}
#' (to use a hybrid between TGUH and Adaptive WBS). If you are unsure where to start, try
#' \code{\link{segment.mean}}.
#'
#' The change-point detection algorithm used in \code{wbs.K.cpt} is the 
#' Adaptive Wild Binary Segmentaton method as described in 
#' "Data-adaptive Wild Binary Segmentation",
#' P. Fryzlewicz (2017), in preparation as of September 28th, 2017.
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param K The number of change-points you wish to detect.
#' @param M The number of randomly selected sub-segments of the data on which to build
#' the CUSUM statistics on each recursively identified interval in the Adaptive Wild Binary Segmentation algorithm.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}; the minumum of \code{K} and \code{n-1}, where \code{n} is the length of \code{x}}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' \item{cpt.sorted}{The list of all possible change-point locations, sorted from the most to the least likely}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{segment.mean}}, \code{\link{wbs.thresh.cpt}},
#' \code{\link{wbs.cpt}}, \code{\link{tguh.cpt}}, \code{\link{hybrid.cpt}}, \code{\link{wbs.bic.cpt}}
#' @examples
#' teeth <- rep(rep(0:1, each=5), 20)
#' teeth.noisy <- teeth + rnorm(200)/5
#' teeth.cleaned <- wbs.K.cpt(teeth.noisy, 39)
#' teeth.cleaned$cpt
#' teeth.cleaned <- wbs.K.cpt(teeth.noisy, 78)
#' teeth.cleaned$cpt
#' teeth.cleaned$cpt.sorted
#' @export

wbs.K.cpt <- function(x, K, M = 1000) {
	
	n <- length(x)
	
	cpt <- integer(0)

	if (n == 1) {

		est <- x
		no.of.cpt <- 0
		cpt.sorted <- integer(0)

	} else {
	
		rc <- t(wbs.K.int(x, M))
		i <- order(abs(rc[,4]), decreasing=T)
		cpt.sorted <- rc[i,3]
		no.of.cpt <- min(K, n-1)
		if (K) cpt <- sort(cpt.sorted[1:no.of.cpt])
		est <- mean.from.cpts(x, cpt)
	}
	
	list(est=est, no.of.cpt=no.of.cpt, cpt=cpt, cpt.sorted = cpt.sorted)

}
