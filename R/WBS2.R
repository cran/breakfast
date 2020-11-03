#' @title Solution path generation via the Wild Binary Segmentation 2 method
#' @description This function arranges all possible change-points in the mean of the input vector in the order of importance, via the Wild Binary Segmentation 2 method.
#' @details 
#' The Wild Binary Segmentation 2 algorithm is described in 
#' "Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection", P. Fryzlewicz (2020), Journal of the Korean Statistical Society, to appear.
#'
#' @param x A numeric vector containing the data to be processed.
#' @param M The maximum number of data sub-samples drawn at each recursive stage of the algorithm. The default is
#' \code{M = 100}. Setting \code{M = 0} executes the standard binary segmentation.
#' @param cusums Whether data sub-intervals for CUSUM computation are drawn systematically (start- and end-points taken from an approximately equispaced grid) or randomly (obtained uniformly with replacement). Use "systematic" (default) or "random".
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{M}{Input parameter \code{M}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "wbs2" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}
#' @references P. Fryzlewicz (2020). Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{Journal of the Korean Statistical Society (to appear)}.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.wbs2(r3)
#' @export
sol.wbs2 <- function(x, M=100, cusums = "systematic") {
	
	# M = 0 means standard Binary Segmentation; M >= 1 means WBS2
	
	solutions.nested <- TRUE
	
	solution.set <- list()
	
	n <- length(x)
	
	sorted.cusums <- matrix(NA, 0, 4)
	
	if (n <= 1) solution.path <- integer()
	
	else {

		if (cusums == "random") cusum.sampling <- random.cusums else if (cusums == "systematic") cusum.sampling <- systematic.cusums

		rc <- t(wbs.K.int(x, M, cusum.sampling))
		
#		rc <- t(wbs.K.int(x, M))

		ord <- order(abs(rc[,4]), decreasing=T)

		sorted.cusums <- abs(rc[ord,, drop=F])

		solution.path <- sorted.cusums[,3]
				
	}	
	
	ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, M = M, cands = sorted.cusums, method = "wbs2")

	class(ret) <- "cptpath"
	
	ret
	
}


wbs.K.int <- function(x, M, cusum.sampling) {

	

	n <- length(x)

	if (n == 1) return(matrix(NA, 4, 0))

	else {

		cpt <- t(cusum.sampling(x, M)$max.val)

		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M, cusum.sampling), wbs.K.int(x[(cpt[3]+1):n], M, cusum.sampling) + c(rep(cpt[3], 3), 0)            ))

	}

	

}





#wbs.K.int <- function(x, M) {

	# M = 0 means standard Binary Segmentation; M >= 1 means WBS2	

#	n <- length(x)

#	if (n == 1) return(matrix(NA, 4, 0))

#	else {

#		cpt <- t(random.cusums(x, M)$max.val)

#		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M), wbs.K.int(x[(cpt[3]+1):n], M) + c(rep(cpt[3], 3), 0)            ))

#	}

	

#}




mean.from.cpt <- function(x, cpt) {



	n <- length(x)

	len.cpt <- length(cpt)

	if (len.cpt) cpt <- sort(cpt)

	beg <- endd <- rep(0, len.cpt+1)

	beg[1] <- 1

	endd[len.cpt+1] <- n

	if (len.cpt) {

		beg[2:(len.cpt+1)] <- cpt+1

		endd[1:len.cpt] <- cpt

	}

	means <- rep(0, len.cpt+1)

	for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])

	rep(means, endd-beg+1)

}


all.intervals.flat <- function(n) {

	

	if (n == 2) ind <- matrix(1:2, 2, 1) else {

		M <- (n-1)*n/2	

		ind <- matrix(0, 2, M)

		ind[1,] <- rep(1:(n-1), (n-1):1)

		ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)

	}

	ind



}



grid.intervals <- function(n, M) {

		

	if (n==2) ind <- matrix(c(1, 2), 2, 1)

	

	else if (M >= (n-1)*n/2) ind <- all.intervals.flat(n)

	

	else {

		k <- 1

		while (k*(k-1)/2 < M) k <- k+1

		ind2 <- all.intervals.flat(k)

		ind2.mx <- max(ind2)

		ind <- round((ind2 - 1) * ((n-1) / (ind2.mx-1)) + 1)

	}	

	

	ind	

}
