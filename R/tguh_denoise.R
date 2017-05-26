#' Noise removal from Tail-Greedy Unbalanced Haar coefficients via connected thresholding
#'
#' This function performs the connected thresholding of the Tail-Greedy Unbalanced Haar 
#' coefficients.
#' 
#' Typically, the first parameter of \code{tguh.denoise} will be an object returned by
#' \code{tguh.decomp}. The function \code{tguh.denoise} performs the "connected thresholding"
#' of this object, in the sense that if a Tail-Greedy Unbalanced Haar detail coefficient does
#' not have any surviving children coefficients, then it gets set to zero if it falls under
#' the threshold, or if the corresponding Unbalanced Haar wavelet is too unbalanced or has too
#' short a wing. See "Tail-greedy bottom-up data decompositions and fast multiple change-point 
#' detection", P. Fryzlewicz (2017), preprint, for details.
#'
#' @param tguh.decomp.obj A variable returned by \code{tguh.decomp} or \code{tguh.denoise}.
#' @param lambda The threshold value.
#' @param minseglen The minimum permitted length of either wing of any Unbalanced Haar wavelet whose
#' corresponding coefficient survives the thresholding.
#' @param bal The minimum permitted ratio of the length of either wing to the sum of the lengths of both wings of
#' any Unbalanced Haar wavelet whose corresponding coefficient survives the thresholding.
#' @return Modified object \code{tguh.decomp.obj}; the modification is that the detail coefficients
#' in the \code{decomp.hist} field that do not survive the thresholding get set to zero.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{tguh.cpt}}, \code{\link{tguh.decomp}}, \code{\link{tguh.reconstr}}
#' @examples
#' rnoise <- rnorm(10)
#' rnoise.tguh <- tguh.decomp(rnoise)
#' print(rnoise.tguh)
#' rnoise.denoise <- tguh.denoise(rnoise.tguh, 3)
#' rnoise.clean <- tguh.reconstr(rnoise.denoise)
#' print(rnoise.clean)
#' @export



tguh.denoise <- function(tguh.decomp.obj, lambda, minseglen = 1, bal = 1/20) {

	n <- tguh.decomp.obj$n

	protected <- rep(0, n)

	for (i in 1:(n-1)) {

		if (!protected[tguh.decomp.obj$decomp.hist[1,1,i]] & !protected[tguh.decomp.obj$decomp.hist[1,2,i]])
			tguh.decomp.obj$decomp.hist[3,1,i] <- tguh.decomp.obj$decomp.hist[3,1,i] * ((abs(tguh.decomp.obj$decomp.hist[3,1,i]) > lambda) & 
				(round(tguh.decomp.obj$decomp.hist[4,1,i]) >= minseglen) & (round(tguh.decomp.obj$decomp.hist[4,2,i]) >= minseglen) & 
				(tguh.decomp.obj$decomp.hist[2,1,i]^2 > bal) & (tguh.decomp.obj$decomp.hist[2,2,i]^2 > bal))

		if (abs(tguh.decomp.obj$decomp.hist[3,1,i]) > 0) protected[tguh.decomp.obj$decomp.hist[1,1,i]] <- 1

	}

	tguh.decomp.obj

}
