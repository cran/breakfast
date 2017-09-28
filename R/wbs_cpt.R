#' Multiple change-point detection in the mean of a vector using the (Adaptive) WBS method.
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, using the (Adaptive) Wild Binary Segmentation
#' method (see Details for the relevant literature references).
#' The constant means between each pair 
#' of neighbouring change-points are also estimated. The method works best when the noise in the 
#' input vector is independent and identically distributed Gaussian.
#'
#' This is a hybrid method, which returns the result of \code{\link{wbs.thresh.cpt}} or
#' \code{\link{wbs.bic.cpt}}, whichever of the two detect the larger number of change-points.
#' If there is a tie, \code{\link{wbs.bic.cpt}} is returned.
#'
#' The change-point detection algorithms used in \code{wbs.thresh.cpt} are: standard
#' Wild Binary Segmentation [see "Wild Binary Segmentation for multiple 
#' change-point detection", P. Fryzlewicz (2014), Annals of Statistics, 42, 2243-2281]
#' and Adaptive Wild Binary Segmentation [see "Data-adaptive Wild Binary Segmentation",
#' P. Fryzlewicz (2017), in preparation as of September 28th, 2017].
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param sigma Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @param M.bic Only relevant to the \code{\link{wbs.bic.cpt}} part (see Details below); the same as the \code{M}
#' parameter in \code{\link{wbs.bic.cpt}}.
#' @param Kmax Only relevant to the \code{\link{wbs.bic.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.bic.cpt}}.
#' @param universal Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @param M.thresh Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the \code{M}
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @param th.const Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @param th.const.min.mult Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @param adapt Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @param lambda Only relevant to the \code{\link{wbs.thresh.cpt}} part (see Details below); the same as the corresponding
#' parameter in \code{\link{wbs.thresh.cpt}}.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}.}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{segment.mean}}, \code{\link{wbs.bic.cpt}},
#' \code{\link{wbs.thresh.cpt}}, \code{\link{tguh.cpt}}, \code{\link{hybrid.cpt}}, \code{\link{wbs.K.cpt}}
#' @examples
#' teeth <- rep(rep(0:1, each=5), 20)
#' teeth.noisy <- teeth + rnorm(200)/5
#' teeth.cleaned <- wbs.cpt(teeth.noisy)
#' ts.plot(teeth.cleaned$est)
#' @importFrom stats mad
#' @export


wbs.cpt <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), M.bic = 20000, Kmax = ceiling(length(x)/5), universal = TRUE, M.thresh = NULL, th.const = NULL, th.const.min.mult = 0.825, adapt = TRUE, lambda = 0.9) {

	wbs.bic <- wbs.bic.cpt(x, M.bic, Kmax)
	wbs.thresh <- wbs.thresh.cpt(x, sigma, universal, M.thresh, th.const, th.const.min.mult, adapt, lambda)
	
	if (wbs.thresh$no.of.cpt > wbs.bic$no.of.cpt) return(wbs.thresh) else return(wbs.bic)
	
}
