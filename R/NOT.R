#' @title Solution path generation via the Narrowest-Over-Threshold method
#' @description This function arranges all possible change-points in the mean of the input vector in the order of importance, via the Narrowest-Over-Threshold (NOT) method.
#' @details 
#' The Narrowest-Over-Threshold method and its algorithm is described in 
#' "Narrowest-over-threshold detection of multiple change points and change-point-like features", R. Baranowski, Y. Chen and P. Fryzlewicz (2019), Journal of Royal Statistical Society: Series B, 81(3), 649--672.
#' @param x A numeric vector containing the data to be processed
#' @param M The maximum number of all data sub-samples at the beginning of the algorithm. The default is
#' \code{M = 10000}
#' @param rand.intervals When drawing the sub-intervals, whether to use a random or a fixed scheme. The default is \code{rand.intervals = TRUE}
#' @param seed If a random scheme is used, a random seed can be provided so that every time the same sets of random sub-intervals would be drawn. The default is \code{seed = NULL}, which means that this option is not taken
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{FALSE}, i.e., the change-point outputs are not nested}
#' \item{solution.path}{Empty list}
#' \item{solution.set}{Locations of possible change-points in the mean of \code{x} for each threshold level (in the decreasing order), arranged in the form of a list of lists}
#' \item{solution.set.th}{A list that contains threshold levels corresponding to the detections in \code{solution.set}}
#' \item{x}{Input vector \code{x}}
#' \item{M}{Input parameter \code{M}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column resulted from applying NOT to all threshold levels. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows reflect the strength of each detection in decreasing order. To avoid repetition, each possible location would appear at most once in the matrix (with the sub-interval that carries its highest possible strength)}
#' \item{method}{The method used, which has value "not" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}},  \code{\link{sol.wbs2}}
#' @references R. Baranowski, Y. Chen & P. Fryzlewicz (2019). Narrowest-over-threshold detection of multiple change points and change-point-like features. \emph{Journal of the Royal Statistical Society: Series B}, 81(3), 649--672.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.not(r3)
#' @export
sol.not <- function(x, M=10000, rand.intervals = TRUE, seed = NULL){


  #veryfing the input parameters - x
  x <- as.numeric(x)
  storage.mode(x) <- "double"
  n <- length(x)
  check.input(x)
  
  #veryfing the input parameters - M
  M <- as.integer(M)
  if(any(is.na(M))) stop("M cannot be NA")
  if(length(M)> 1)  stop("M should be a single integer.")
  if(M<0)  stop("M should be an integer > 0.")

  #veryfing the input parameters - rand.intervals
  rand.intervals <- as.logical(rand.intervals)
  
  #drawing the intervals over which the contrast function will be computed
  if(rand.intervals) intervals <-  matrix(random.intervals(n,M,seed),ncol=2)
  else {
    intervals <- matrix(fixed.intervals(n,M),ncol=2)
  }
  #making sure that not_r_wrapper gets the right type
  storage.mode(intervals) <- "integer"
  
  #picking the correct type of the contrast function - pcwsConstMean for now
  contrast <- as.integer(0)
  
  # method = 0 means not (instead of the "wbs" type approach - when method = 1)
  method <- as.integer(0)
  
  # parallel option and augmented option suppressed for now
  parallel <- as.integer(0)
  augmented <- as.integer(0)
  
  
  out <- call_not_r_wrapper (x, intervals, method, contrast, parallel, augmented)
  
  #calling the C code - old way without Rcpp
  #out <-.Call("not_r_wrapper", x, intervals, method, contrast, parallel, augmented, PACKAGE="breakfast")
  index <- unique(unlist(out$solution.path$index,use.names=F))
  cands <- out$contrasts[index,, drop = FALSE]
  cands <- cands[order(cands$max.contrast, decreasing = TRUE),c(1,2,4,5), drop = FALSE]
  rownames(cands) <- NULL
  colnames(cands) <- NULL
  cands <- as.matrix(cands)
  cands <- cands [!duplicated(cands[,3]),, drop = FALSE]

  ret <- list()	  
  ret$solutions.nested = FALSE
  ret$solution.path = integer(0)
  ret$solution.set = out$solution.path$cpt
  ret$solution.set.th = out$solution.path$th
  ret$x = x
  ret$M = M
  ret$cands = cands  
  ret$method = "not"
  class(ret) <- "cptpath"

  return(ret)
}


