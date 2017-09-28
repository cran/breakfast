#' Multiple change-point detection in the mean of a vector using the WBS method, with the number of change-points chosen by BIC
#'
#' This function estimates the number and locations of change-points in the 
#' piecewise-constant mean of the noisy input vector, using the Wild Binary Segmentation
#' method (see Details for the relevant literature reference). The number of change-points
#' is chosen via the Bayesian Information Criterion. The constant means between each pair 
#' of neighbouring change-points are also estimated. The method works best when the noise in the 
#' input vector is independent and identically distributed Gaussian, and when the number
#' change-points is small.
#'
#' The BIC penalty is unsuitable as a model selection tool in long signals
#' with frequent change-points; if you need a more versatile function that works well 
#' regardless of the number of change-points, try \code{\link{segment.mean}} (for a
#' default recommended estimation technique), \code{\link{wbs.thresh.cpt}},
#' \code{\link{wbs.cpt}} (if you require an (Adaptive) WBS-based technique), \code{\link{tguh.cpt}}
#' (if you require a TGUH-based technique), or \code{\link{hybrid.cpt}}
#' (to use a hybrid between TGUH and Adaptive WBS). If you are unsure where to start, try
#' \code{\link{segment.mean}}. (If you know how many change-points you wish to detect,
#' try \code{\link{wbs.K.cpt}}.)
#'
#' The change-point detection algorithm used in \code{wbs.bic.cpt} is the 
#' Wild Binary Segmentaton method as described in "Wild Binary Segmentation for multiple 
#' change-point detection", P. Fryzlewicz (2014), Annals of Statistics, 42, 2243-2281.
#'
#' @param x A vector containing the data in which you wish to find change-points.
#' @param M The number of randomly selected sub-segments of the data on which to build
#' the CUSUM statistics in the Wild Binary Segmentation algorithm; generally, the larger
#' the value of M, the more accurate but slower the algorithm - but see the remarks below
#' about the BIC penalty.
#' @param Kmax The maximum number of change-points that can be detected.
#' @return A list with the following components:
#' \item{est}{The estimated piecewise-constant mean of \code{x}.}
#' \item{no.of.cpt}{The estimated number of change-points in the piecewise-constant mean of \code{x}.}
#' \item{cpt}{The estimated locations of change-points in the piecewise-contant mean of \code{x} (these
#' are the final indices \emph{before} the location of each change-point).}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{segment.mean}}, \code{\link{wbs.thresh.cpt}},
#' \code{\link{wbs.cpt}}, \code{\link{tguh.cpt}}, \code{\link{hybrid.cpt}}, \code{\link{wbs.K.cpt}}
#' @examples
#' teeth <- rep(rep(0:1, each=5), 20)
#' teeth.noisy <- teeth + rnorm(200)/5
#' teeth.cleaned <- wbs.bic.cpt(teeth.noisy)
#' ts.plot(teeth.cleaned$est)
#' teeth.cleaned$no.of.cpt
#' teeth.cleaned$cpt
#' @importFrom stats var
#' @export


wbs.bic.cpt <- function(x, M=20000, Kmax = ceiling(length(x)/5)) {

	n <- length(x)
	
	if (n == 1) {

		est <- x
		no.of.cpt <- 0
		cpt <- integer(0)

	} else {

		rc <- random.cusums(x, M)
		w <- wbs.sort(rc)

		K <- dim(w)[1]
	
		if (is.null(Kmax)) K.eff <- K else K.eff <- min(Kmax, K)

		bic.criterion <- min.log.lik.K <- rep(0, K.eff+1)

		min.log.lik.K[1] <- bic.criterion[1] <- n/2 * log(stats::var(x))

		if (K.eff) for (i in 1:K.eff) {	
			tmp <- mean.from.cpts(x, w[1:i, 3])	
			min.log.lik.K[i+1] <- n/2 * log(stats::var(x - tmp))
			bic.criterion[i+1] <- min.log.lik.K[i+1] + i * log(n)
		}

		K <- which.min(bic.criterion)[1]-1

		if (K == 0) est <- rep(mean(x), n) else 
			est <- mean.from.cpts(x, w[1:K, 3])
			
		cpt <- which(abs(diff(est)) > 0)
		no.of.cpt <- length(cpt)	
			
	}

	list(est=est, no.of.cpt=no.of.cpt, cpt=cpt)
}
