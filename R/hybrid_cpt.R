#' Multiple change-point detection in the mean of a vector using a hybrid between the TGUH and Adaptive WBS methods.
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, combining the Tail-Greedy Unbalanced Haar and Adaptive Wild Binary Segmentation
#' methods (see Details for the relevant literature references).
#' The constant means between each pair 
#' of neighbouring change-points are also estimated. The method works best when the noise in the 
#' input vector is independent and identically distributed Gaussian.
#'
#' This is a hybrid method, which first estimates the number of change-points using 
#' \code{\link{tguh.cpt}} and then estimates their locations using \code{\link{wbs.K.cpt}}.
#'
#' The change-point detection algorithms used in \code{tguh.cpt} are: the 
#' Tail-Greedy Unbalanced Haar method as described in "Tail-greedy bottom-up data
#' decompositions and fast multiple change-point detection", P. Fryzlewicz (2017),
#' preprint, and Adaptive Wild Binary Segmentation as described in "Data-adaptive Wild Binary Segmentation",
#' P. Fryzlewicz (2017), in preparation as of September 28th, 2017.
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param M The same as the corresponding parameter in \code{\link{wbs.K.cpt}}.
#' @param sigma The same as the corresponding parameter in \code{\link{tguh.cpt}}.
#' @param th.const The same as the corresponding parameter in \code{\link{tguh.cpt}}.
#' @param p The same as the corresponding parameter in \code{\link{tguh.cpt}}.
#' @param minseglen The same as the corresponding parameter in \code{\link{tguh.cpt}}.
#' @param bal The same as the corresponding parameter in \code{\link{tguh.cpt}}.
#' @param num.zero The same as the corresponding parameter in \code{\link{tguh.cpt}}.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}.}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{segment.mean}}, \code{\link{wbs.bic.cpt}},
#' \code{\link{wbs.thresh.cpt}}, \code{\link{wbs.cpt}}, \code{\link{tguh.cpt}}, \code{\link{wbs.K.cpt}}
#' @examples
#' teeth <- rep(rep(0:1, each=5), 20)
#' teeth.noisy <- teeth + rnorm(200)/5
#' teeth.cleaned <- hybrid.cpt(teeth.noisy)
#' ts.plot(teeth.cleaned$est)
#' @importFrom stats mad
#' @export


hybrid.cpt <- function(x, M = 1000, sigma = stats::mad(diff(x)/sqrt(2)), th.const = 1, p = 0.01, minseglen = 1, bal = 1/20, num.zero = 10^(-5)) {
	
	x.tguh <- tguh.cpt(x, sigma, th.const, p, minseglen, bal, num.zero)
	wbs.K.cpt(x, x.tguh$no.of.cpt, M)
	
}
