#' The inverse Tail-Greedy Unbalanced Haar transformation
#'
#' This function performs the inverse Tail-Greedy Unbalanced Haar transformation,
#' also referred to as reconstruction.
#' 
#' The Tail-Greedy Unbalanced Haar decomposition and reconstruction algorithms are described in 
#' "Tail-greedy bottom-up data decompositions and fast multiple change-point 
#' detection", P. Fryzlewicz (2017), preprint.
#'
#' @param tguh.decomp.obj A variable returned by \code{tguh.decomp} or \code{tguh.denoise}.
#' @return A vector being the result of the inverse Tail-Greedy Unbalanced Haar transformation of \code{tghu.decomp.obj}.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{tguh.cpt}}, \code{\link{tguh.decomp}}, \code{\link{tguh.denoise}}
#' @examples
#' rnoise <- rnorm(10)
#' rnoise.tguh <- tguh.decomp(rnoise)
#' print(rnoise.tguh)
#' rnoise.denoise <- tguh.denoise(rnoise.tguh, 3)
#' rnoise.clean <- tguh.reconstr(rnoise.denoise)
#' print(rnoise.clean)
#' @export




tguh.reconstr <- function(tguh.decomp.obj) {

	n <- tguh.decomp.obj$n

	for (i in (n-1):1) {

		ind <- tguh.decomp.obj$decomp.hist[1,,i]

		tguh.decomp.obj$decomp.hist[3,2,i] <- tguh.decomp.obj$tguh.coeffs[min(ind)]

		tguh.decomp.obj$tguh.coeffs[ind[2]] <- sum(tguh.decomp.obj$decomp.hist[2,2:1,i] *
			tguh.decomp.obj$decomp.hist[3,,i])
		tguh.decomp.obj$tguh.coeffs[ind[1]] <- sum(tguh.decomp.obj$decomp.hist[2,,i] *
			c(1,-1) * tguh.decomp.obj$decomp.hist[3,,i])


	}

	tguh.decomp.obj$tguh.coeffs

}
