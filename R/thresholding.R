#' @title Estimating change-points in the piecewise-constant mean of a noisy data sequence via thresholding
#' @description This function estimates the number and locations of change-points in the piecewise-constant mean of a noisy data sequence via thresholding.
#' 
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. Note that the field \code{sols.object$x} contains the input data sequence.
#' @param sigma An estimate of the standard deviation of the noise in the data \code{cptpath.object$x}. Can be a functional of \code{cptpath.object$x} or a specific value if known. The default is the Median Absolute Deviation of the vector \code{diff(cptpath.object$x)/sqrt(2)}, tuned to the Gaussian distribution. Note that \code{model.thresh} works particularly well when the noise is i.i.d. Gaussian.
#' @param th_const A positive real number with default value equal to 1. It is used to define the threshold for the detection process.
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{model}{The model selection method used to return the final change-point estimators object, here its value is \code{"thresh"}}
#' \item{no.of.cpt}{The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#' \item{cpts}{The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#' \item{est}{An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#' @seealso \code{\link{sol.idetect_seq}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{breakfast}}
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f))
#' model.thresh(sol.idetect_seq(x))
#' @export
model.thresh <- function(cptpath.object,  sigma = stats::mad(diff(cptpath.object$x)/ sqrt(2)), th_const = 1.15){
if(class(cptpath.object) != 'cptpath') stop('cptpath.object must be in cptpath class')
x <- cptpath.object$x
lx <- length(x)
if (lx <= 1) {
  est <- x
  no.of.cpt <- 0
  cpts <- integer(0)
}
else {
  if (sigma == 0) {
    s0 <- all.shifts.are.cpts(x)
    est <- s0$est
    no.of.cpt <- s0$no.of.cpt
    cpts <- s0$cpts
    }
else{
if(cptpath.object$method == 'idetect') warning('To get the Isolate-Detect method results when the model selection is "thresh" we recommend that you use sol.idetect_seq to create the cptpath.object instead of sol.idetect')
if(cptpath.object$method == 'tguh'){th_const <- 1}
th_fin <- sigma * th_const * sqrt(2*log(lx))
 if (cptpath.object$solutions.nested == FALSE){
   critical <- length(which(cptpath.object$solution.set.th > th_fin))
   if (critical == 0){cpts <- 0
   no.of.cpt <- 0
   }
  else{cpts <- sort(cptpath.object$solution.set[[critical]])
  no.of.cpt <- length(cpts)
  }
 }
else{
  positions <- which(cptpath.object$cands[,4] > th_fin)
  critical <- length(positions)
  if (critical == 0){cpts <- 0
  no.of.cpt <- 0
  }
  else{cpts <- sort(cptpath.object$cands[positions,3])
  no.of.cpt <- length(cpts)
  }
}
est <- mean.from.cpt(x, cpts)
}
}

ret <- list(solution.path=cptpath.object$method, model="thresh", no.of.cpt=no.of.cpt, cpts=cpts, est=est)

class(ret) <- "cptmodel"

return(ret)
}
