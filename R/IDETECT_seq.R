#' @title  Solution path generation using the sequential approach of the Isolate-Detect method
#' @description This function uses the Isolate-Detect method in its original sequential way in order to create the solution path.
#' It is developed to be used with the thresholding model selection rule.
#' @details 
#' The Isolate-Detect method and its algorithm is described in 
#' "Detecting multiple generalized change-points by isolating single ones", A. Anastasiou & P. Fryzlewicz (2019), arXiv preprint arXiv:1901.10852.
#' 
#' @param x A numeric vector containing the data to be processed
#' @param points A positive integer with default value equal to 3. It defines the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, as described in the Isolate-Detect methodology.
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "idetect_seq" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.not}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{sol.tguh}}, 
#' @references A. Anastasiou & P. Fryzlewicz (2019). Detecting multiple generalized change-points by isolating single ones. \emph{arXiv preprint arXiv:1901.10852}.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.idetect_seq(r3)
#' @export
sol.idetect_seq <- function(x, points = 3) {
  solutions.nested <- TRUE
  solution.set <- list()
  cands <- matrix(NA, 0, 4)
  lx <- length(x)
  if (lx < points) {solution.path <- integer()
  ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, cands = cands, method = "idetect_seq")
  
  class(ret) <- "cptpath"
  
  ret
  }
  else{
    points <- as.integer(points)
    step1 <- window.idetect.th(x, thr_con = 1.15, w_points = points)
    s2 <- step1$changepoints
    if(length(s2) == 0){
      m <- matrix(0, 1, 4)
      right_points <- seq(points, lx, points)
      l_r <- length(right_points)
      k_r <- 1
      while (k_r <= l_r) {
        m <- rbind(m, c(1,right_points[k_r],max.cusum(ind = c(1,right_points[k_r]), y = step1$y)))
        k_r <- k_r + 1
      }
      m <- m[which(m[,3] != 0),,drop = FALSE]
      ord <- order(m[,4], decreasing = T)
      cands <- m[ord, ,drop=FALSE]
      cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
      solution.path <- cands[,3]
      
      ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, cands = cands, method = "idetect_seq")
      
      class(ret) <- "cptpath"
      
      ret
    }
    else{Res <- step1$full_information
    ord <- order(Res[,4,drop = FALSE], decreasing = T)
    cands <- Res[ord, ,drop=FALSE]
    cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
    solution.path <- cands[,3]
  ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, cands = cands, method = "idetect_seq")
  
  class(ret) <- "cptpath"
  
  ret
    }
  }
}
