#' @title Solution path generation via the Wild Binary Segmentation 2 method
#' @description This function arranges all possible change-points in the mean of the input vector in the order of importance, via the Wild Binary Segmentation 2 method.
#' @details 
#' The Wild Binary Segmentation 2 algorithm is described in 
#' "Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection", P. Fryzlewicz (2020), Journal of the Korean Statistical Society, 49, 1027-1070.
#' @param x A numeric vector containing the data to be processed.
#' @param type The model type considered. \code{type = "const"} means piecewise-constant; this is the only type currently supported in \code{sol.wbs2}
#' @param M The maximum number of data sub-samples drawn at each recursive stage of the algorithm. The default is
#' \code{M = 1000}. Setting \code{M = 0} executes the standard binary segmentation.
#' @param systematic.intervals Whether data sub-intervals for CUSUM computation are drawn systematically (TRUE; start- and end-points taken from an approximately equispaced grid) or randomly (FALSE; obtained uniformly with replacement). The default is TRUE.
#' @param seed If a random scheme is used, a random seed can be provided so that every time the same sets of random sub-intervals would be drawn. The default is \code{seed = NULL}, which means that this option is not set
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{type}{Input parameter \code{type}}
#' \item{M}{Input parameter \code{M}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "wbs2" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}
#' @references P. Fryzlewicz (2020). Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{Journal of the Korean Statistical Society}, 49, 1027-1070.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.wbs2(r3)
#' @export
sol.wbs2 <- function(x, type = "const", M = 1000, systematic.intervals = TRUE, seed = NULL) {
	
	# M = 0 means standard Binary Segmentation; M >= 1 means WBS2
	x <- as.numeric(x)
	storage.mode(x) <- "double"
	n <- length(x)
	check.input(x)
	
	stopifnot(type %in% c('const'))
  
  #veryfing the input parameters - M
	M <- as.integer(M)
	if(any(is.na(M))) stop("M cannot be NA")
	if(length(M)> 1)  stop("M should be a single integer.")
	if(M<0)  stop("M should be an integer > 0.")

	solutions.nested <- TRUE
	
	solution.set <- list()
		
	sorted.cusums <- matrix(NA, 0, 4)
	
	if (n <= 1) solution.path <- integer()
	
	else {
	  
	  set.seed(seed)

		if (!systematic.intervals) cusum.sampling <- random_cusums else cusum.sampling <- systematic_cusums

		rc <- t(wbs_K_int(x, M, cusum.sampling))
		
		ord <- order(abs(rc[,4]), decreasing=T)

		sorted.cusums <- abs(rc[ord,, drop=F])

		solution.path <- sorted.cusums[,3]
				
	}
	ret <- list()	  
	ret$solutions.nested = TRUE
	ret$solution.path <- solution.path
	ret$solution.set <- solution.set
	ret$x = x
	ret$type = type
	ret$M = M
	ret$cands = as.matrix(sorted.cusums)
	ret$method = "wbs2"
	class(ret) <- "cptpath"
	
	return(ret)
	
}
