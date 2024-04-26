#' @title Solution path generation via the Wild Binary Segmentation method
#' @description This function arranges all possible change-in-mean features of the input vector in the order of importance, via the Wild Binary Segmentation (WBS) method.
#' @details 
#' The Wild Binary Segmentation algorithm is described in 
#' "Wild binary segmentation for multiple change-point detection", P. Fryzlewicz (2014), The Annals of Statistics, 42: 2243--2281.
#' @param x A numeric vector containing the data to be processed
#' @param type The model type considered. Currently \code{type = "const"} is the only accepted value. This assumes that the mean of the input vector is piecewise-constant.
#' @param M The maximum number of all data sub-samples at the beginning of the algorithm. The default is
#' \code{M = 10000}
#' @param systematic.intervals When drawing the sub-intervals, whether to use a systematic (and fixed) or random scheme. The default is \code{systematic.intervals = TRUE}
#' @param seed If a random scheme is used, a random seed can be provided so that every time the same sets of random sub-intervals would be drawn. The default is \code{seed = NULL}, which means that this option is not set
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{type}{The input parameter \code{type}}
#' \item{M}{Input parameter \code{M}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "wbs" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs2}}
#' @references P. Fryzlewicz (2014). Wild binary segmentation for multiple change-point detection. \emph{The Annals of Statistics}, 42(6), 2243--2281.
#' @references R. Baranowski, Y. Chen & P. Fryzlewicz (2019). Narrowest-over-threshold detection of multiple change points and change-point-like features. \emph{Journal of the Royal Statistical Society: Series B}, 81(3), 649--672.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.wbs(r3)
#' @export
sol.wbs <- function(x, type="const", M=10000, systematic.intervals = TRUE, seed = NULL){
  
  #verifying the input parameters - x
  x <- as.numeric(x)
  storage.mode(x) <- "double"
  n <- length(x)
  check.input(x)
  
  #verifying the input parameters - type
  #picking the correct type of the contrast function:
  # 0 for const, 1 for lin.cont and 2 for lin.discont for now
#  stopifnot(type %in% c('const', 'lin.cont', 'lin.discont'))
  stopifnot(type %in% c('const'))
    if(type=="lin.cont")   contrast <- as.integer(1)
  else {
    if(type=="lin.discont") contrast <- as.integer(2)
    else contrast <- as.integer(0)
  }
  
  #verifying the input parameters - M
  M <- as.integer(M)
  if(any(is.na(M))) stop("M cannot be NA")
  if(length(M)> 1)  stop("M should be a single integer.")
  if(M<0)  stop("M should be an integer > 0.")

  #veriyfing the input parameters - systematic.intervals 
  systematic.intervals  <- as.logical(systematic.intervals)
  
  #drawing the intervals over which the contrast function will be computed
  if(!(systematic.intervals)) intervals <-  matrix(random_intervals(n,M,seed),ncol=2)
  else {
    intervals <- matrix(fixed_intervals(n,M),ncol=2)
  }
  #making sure that not_r_wrapper gets the right type
  storage.mode(intervals) <- "integer"
  
  # method = 1 means wbs (instead of the "not" type approach - when method = 0)
  method <- as.integer(1)
  
  # parallel option and augmented option suppressed for now
  parallel <- as.integer(0)
  augmented <- as.integer(0)
  
  out <- call_not_r_wrapper (x, intervals, method, contrast, parallel, augmented)
  
  #calling the C code - old way without Rcpp - suppressed for now
  #out <-  .Call("not_r_wrapper", x, intervals, method, contrast, parallel, augmented, PACKAGE="breakfast")

  if (length(out$solution.path$cpt)==0){
  	index = NULL 
  }
  else{
    index <- out$solution.path$index[[length(out$solution.path$cpt)]]
  }
  cands <- out$contrasts[index,, drop = FALSE]
  cands <- cands[order(cands$max.contrast, decreasing = TRUE),c(1,2,4,5), drop = FALSE]
  path <-  cands[,3]
  rownames(cands) <- NULL
  colnames(cands) <- NULL
  
 
  
  ret <- list()	  
  ret$solutions.nested = TRUE
  ret$solution.path <- path
  ret$solution.set <- list()
  ret$x = x
  ret$type = type
  ret$M = M
  ret$cands = as.matrix(cands)
  ret$method = "wbs"
  class(ret) <- "cptpath"

  return(ret)
}
