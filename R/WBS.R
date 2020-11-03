#' @title Solution path generation via the Wild Binary Segmentation method
#' @description This function arranges all possible change-points in the mean of the input vector in the order of importance, via the Wild Binary Segmentation (WBS) method.
#' @details 
#' The Wild Binary Segmentation algorithm is described in 
#' "Wild binary segmentation for multiple change-point detection", P. Fryzlewicz (2014), The Annals of Statistics, 42: 2243--2281.
#'
#' @param x A numeric vector containing the data to be processed
#' @param M The maximum number of all data sub-samples at the beginning of the algorithm. The default is
#' \code{M = 10000}
#' @param rand.intervals When drawing the sub-intervals, whether to use a random or a fixed scheme. The default is \code{rand.intervals = TRUE}
#' @param seed If a random scheme is used, a random seed can be provided so that every time the same sets of random sub-intervals would be drawn. The default is \code{seed = NULL}, which means that this option is not taken
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{M}{Input parameter \code{M}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "wbs" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs2}}
#' @references P. Fryzlewicz (2014). Wild binary segmentation for multiple change-point detection. \emph{The Annals of Statistics}, 42(6), 2243--2281.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.wbs(r3)
#' @export
sol.wbs <- function(x, M=10000, rand.intervals = TRUE, seed = NULL){
  
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
  
  # method = 1 means wbs (instead of the "not" type approach)
  method <- as.integer(1)
  
  # parallel option and augmented option suppressed for now
  parallel <- as.integer(0)
  augmented <- as.integer(0)
  
  #calling the C code
  out <-  .Call("not_r_wrapper", x, intervals, method, contrast, parallel, augmented, PACKAGE="breakfast")
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
  ret$M = M
  ret$method = "wbs"
  ret$cands = as.matrix(cands)
  class(ret) <- "cptpath"

  return(ret)
}
