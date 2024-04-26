#' @title Estimating change-points in the piecewise-constant or piecewise-linear mean of a noisy data sequence via thresholding
#' @description This function estimates the number and locations of change-points in the piecewise-constant or piecewise-linear mean of a noisy data sequence via thresholding.
#' 
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. The \code{cptpath.object$type} variable decides the model type: piecewise-constant (\code{type == "const"}),
#' piecewise-linear and continuous (\code{type == "lin.cont"}) or piecewise-linear and discontinuous (\code{type == "lin.discont"}). In the piecewise-linear model (whether continuous or not), the output of \code{sol.idetect_seq} or \code{sol.not} should be supplied as
#' \code{cptpath.object}. Note that the field \code{cptpath.object$x} contains the input data sequence.
#' @param sigma An estimate of the standard deviation of the noise in the data \code{cptpath.object$x}. Can be a functional of \code{cptpath.object$x} or a specific value if known.
#' The default in the piecewise-constant model is the Median Absolute Deviation of the vector \code{diff(cptpath.object$x)/sqrt(2)}, tuned to the Gaussian distribution.
#' In the piecewise-linear models, \code{diff(cptpath.object$x, differences = 2)/sqrt(6)} is used by default.
#' Note that \code{model.thresh} works particularly well when the noise is i.i.d. Gaussian.
#' @param th.const A positive real number used to define the threshold for the detection process. The default used in the piecewise-constant model is 1.15, while in the piecewise-linear model, the value is taken equal to 1.4.
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{type}{The model type used, inherited from the given \code{cptpath.object}}
#' \item{model.selection}{The model selection method used to return the final change-point estimators object, here its value is \code{"thresh"}}
#' \item{no.of.cpt}{The number of estimated change-points}
#' \item{cpts}{The locations of estimated change-points}
#' \item{est}{An estimate of the mean of the vector \code{cptpath.object$x}}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{breakfast}}
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f))
#' model.thresh(sol.idetect_seq(x))
#' @export
model.thresh <- function(cptpath.object,  sigma = NULL, th.const = NULL){
  
  if(!("cptpath" %in%  class(cptpath.object))) stop("A cptmodel class object has to be supplied in the first argument.")
  
  x <- cptpath.object$x
  type <- cptpath.object$type
  lx <- length(x)
  if (lx <= 1) {
    est <- x
    no.of.cpt <- 0
    cpts <- integer(0)
  }
  else {
    if (type == "const"){sigma <- stats::mad(diff(cptpath.object$x)/ sqrt(2))}
    else{
      sigma <- stats::mad(diff(diff(cptpath.object$x))/ sqrt(6))
    }
    if (sigma == 0) {
      if(type=="const"){s0 <- all_shifts_are_cpts(x)}
      else{
        s0 <- all_slopechanges_are_cpts(x)
}
      est <- s0$est
      no.of.cpt <- s0$no.of.cpt
      cpts <- s0$cpts
    }
    else{
      th.const <- ifelse(type == "const", 1.15,1.4)
      if(cptpath.object$method == 'idetect') warning('To get the Isolate-Detect method results when the model selection is "thresh" we recommend that you use sol.idetect_seq to create the cptpath.object instead of sol.idetect')
      if((cptpath.object$method == 'tguh') & (type=="const")){th.const <- 1}
      if ((cptpath.object$method != "idetect_seq") & (type != "const") & (cptpath.object$method != "not")) warning("In the piecewise-linear change-point model, please consider supplying a cptpath.object produced with sol.idetect_seq or with sol.not.")
      th_fin <- sigma * th.const * sqrt(2*log(lx))
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
      if(type == "const"){est <- mean_from_cpt(x, cpts)}
      if(type == "lin.cont"){est <- slope_from_cpt(x, cpts)}
      if(type == "lin.discont"){est <- prediction(cptpath.object, cpts)}
    }
  }
  ret <- list(solution.path=cptpath.object$method, type = type, model.selection="thresh", no.of.cpt=no.of.cpt, cpts=cpts, est=est)
  
  class(ret) <- "cptmodel"
  
  return(ret)
}
