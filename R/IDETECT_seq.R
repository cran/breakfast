#' @title  Solution path generation using the sequential approach of the Isolate-Detect method
#' @description This function arranges all possible change-points in the mean of the input vector, or in its linear trend, in the order of importance, via the Isolate-Detect (ID) method. 
#' It is developed to be used with the thresholding model selection rule.
#' @details
#' The Isolate-Detect method and its algorithm is described in 
#' "Detecting multiple generalized change-points by isolating single ones", A. Anastasiou & P. Fryzlewicz (2022), Metrika, https://doi.org/10.1007/s00184-021-00821-6.
#' @param x A numeric vector containing the data to be processed
#' @param type The model type considered. \code{type = "const"}, \code{type = "lin.cont"}, \code{type = "lin.discont"} mean, respectively, that the signal (mean of \code{x}) is piecewise constant, 
#'   piecewise linear and continuous, and piecewise linear but not necessarily continuous. If not given, the default is \code{type = "const"}
#' @param points A positive integer with default value equal to 4. It defines the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, as described in the Isolate-Detect methodology.
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{type}{The input parameter \code{type}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "idetect_seq" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.not}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{sol.tguh}} 
#' @references A. Anastasiou & P. Fryzlewicz (2022). Detecting multiple generalized change-points by isolating single ones. \emph{Metrika}, https://doi.org/10.1007/s00184-021-00821-6.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.idetect_seq(r3)
#' @export
sol.idetect_seq <- function(x, type = "const", points = 4) {
  
  #veryfing the input parameters - x
  x <- as.numeric(x)
  storage.mode(x) <- "double"
  lx <- length(x)
  check.input(x)
  
  stopifnot(type %in% c('const', 'lin.cont', 'lin.discont'))
  solutions.nested <- TRUE
  solution.set <- list()
  cands <- matrix(NA, 0, 4)
  if (lx < points) {solution.path <- integer()
  ret <- list()
  ret$solutions.nested = TRUE
  ret$solution.path <- solution.path
  ret$solution.set <- solution.set
  ret$x = x
  ret$type = type
  ret$cands = as.matrix(cands)
  ret$method = "idetect_seq"
  class(ret) <- "cptpath"
  
  return(ret)
  }
  else{
    points <- as.integer(points)
  if(type=="lin.cont")  {step1 <- window.idetect.th_lin(x, thr_con = 1.4, w_points = points)}
  if(type=="lin.discont")  {step1 <- window.idetect.th_discont_lin(x, thr_con = 1.4, w_points = points)}
  if(type == "const"){step1 <- window.idetect.th(x, thr_con = 1.15, w_points = points)}
    s2 <- step1$changepoints
    if(length(s2) == 0){
      m <- matrix(0, 1, 4)
      right_points <- seq(points, lx, points)
      pos_r <- numeric()
      CUSUM_r <- numeric()
      l_r <- length(right_points)
      k_r <- 1
      while (k_r <= l_r) {
        if(type=="lin.cont") {x_temp_r <- x[1:right_points[k_r]]
        ipcr <- IDetect_cumsum_lin(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr))
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r]])
        m <- rbind(m, c(1,right_points[k_r],pos_r[k_r], CUSUM_r[k_r]))}
        if(type=="lin.discont") {x_temp_r <- x[1:right_points[k_r]]
        l_temp_r <- right_points[k_r]
        ipcr <- IDetect_linear_discontr_one(x_temp_r, s = rep(1, l_temp_r - 2), b = 1:(l_temp_r - 2), e = rep(l_temp_r, l_temp_r - 2))
        pos_r[k_r] <- which.max(abs(ipcr))
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r]])
        m <- rbind(m, c(1,right_points[k_r],pos_r[k_r], CUSUM_r[k_r]))}
        else{m <- rbind(m, c(1,right_points[k_r],max_cusum(ind = c(1,right_points[k_r]), y = step1$y)))
        }
        k_r <- k_r + 1
      }
      m <- m[which(m[,3] != 0),,drop = FALSE]
      ord <- order(m[,4], decreasing = T)
      cands <- m[ord, ,drop=FALSE]
      cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
      solution.path <- cands[,3]
      ret <- list()
      ret$solutions.nested = TRUE
      ret$solution.path <- solution.path
      ret$solution.set <- solution.set
      ret$x = x
      ret$type = type
      ret$cands = as.matrix(cands)
      ret$method = "idetect_seq"
      class(ret) <- "cptpath"
      
      return(ret)}
    else{Res <- step1$full_information
    ord <- order(Res[,4,drop = FALSE], decreasing = T)
    cands <- Res[ord, ,drop=FALSE]
    cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
    solution.path <- cands[,3]
    ret <- list()
    ret$solutions.nested = TRUE
    ret$solution.path <- solution.path
    ret$solution.set <- solution.set
    ret$x = x
    ret$type = type
    ret$cands = as.matrix(cands)
    ret$method = "idetect_seq"
    class(ret) <- "cptpath"
    return(ret)
    }
  }
}
