#' Multiple change-point detection in the mean of a vector
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, using a method that puts more emphasis
#' either on "speed" (i.e. is faster but possibly less accurate) or on "accuracy" (i.e. is
#' possibly more accurate but slower).
#' It also estimates the constant means between each pair of neighbouring change-points.
#' It works best when the noise in the input vector is independent and identically distributed Gaussian.
#'
#' In the current version of the package, \code{attribute="speed"} triggers the function
#' \code{\link{tguh.cpt}} and \code{attribute="accuracy"} triggers the function
#' \code{\link{hybrid.cpt}}. \strong{Warning:} this can change in future versions of the
#' package. Note that \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}} return the same number
#' of change-points and the only difference lies in their estimated locations.
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param attribute As described in the Details section of this help file.
#' @param M The same as the corresponding parameter in \code{\link{hybrid.cpt}}.
#' @param sigma The same as the corresponding parameter in \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}}.
#' @param th.const The same as the corresponding parameter in \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}}.
#' @param p The same as the corresponding parameter in \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}}.
#' @param minseglen The same as the corresponding parameter in \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}}.
#' @param bal The same as the corresponding parameter in \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}}.
#' @param num.zero The same as the corresponding parameter in \code{\link{tguh.cpt}} and \code{\link{hybrid.cpt}}.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}.}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{tguh.cpt}}, \code{\link{hybrid.cpt}}, \code{\link{wbs.cpt}}
#' @examples
#' stairs <- rep(1:50, each=10)
#' stairs.noisy <- stairs + rnorm(500)/5
#' stairs.cleaned <- segment.mean(stairs.noisy)
#' ts.plot(stairs.cleaned$est)
#' stairs.cleaned$no.of.cpt
#' stairs.cleaned$cpt
#' @importFrom stats mad
#' @export

segment.mean <- function(x, attribute="speed", M = 1000, sigma = stats::mad(diff(x)/sqrt(2)), th.const = 1, p = 0.01, minseglen = 1, bal = 1/20, num.zero = 10^(-5)) {

	if (attribute == "speed") return(tguh.cpt(x, sigma, th.const, p, minseglen, bal, num.zero))
	else if (attribute == "accuracy") return(hybrid.cpt(x, M, sigma, th.const, p, minseglen, bal, num.zero))
	else return(tguh.cpt(x, sigma, th.const, p, minseglen, bal, num.zero))

}
