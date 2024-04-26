#' @title Estimate the location of change-points when the number of them is fixed
#' @description Return a solution with the given number of change-points or change-point-type features from the solution path 
#' @details 
#' The model selection method which returns results with a given number of change-points or change-point-type features. If there are 
#' multiple such elements on the solution path, the one with the smaller residual sum of squares will be returned. On the other hand,
#' if no such element exists, an empty set (i.e. with no change-points) will be returned.
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. Note that the field \code{cptpath.object$x} contains the input data sequence. 
#' @param fixednum The number of change-points or change-point-type features 
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{type}{The model type used, inherited from the given \code{cptpath.object}}
#' \item{model.selection}{The model selection method used to return the final change-point or change-point-type feature estimators object, here its value is \code{"ic"}}
#' \item{no.of.cpt}{The number of estimated features in the mean of the vector \code{cptpath.object$x} based on the given \code{type} of the model}
#' \item{cpts}{The locations of estimated features in the mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean or constant-slope intervals}
#' \item{est}{An estimate of the mean of the vector \code{cptpath.object$x}; for piecewise-constant signals, the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points; 
#' for piecewise-linear but discontinuous signals, the values are the estimated linear trend (replicated a suitable number of times) between each pair of consecutive detected change of slopes; 
#' for piecewise-linear and continuous signals, it is similar to the previous case but with the continuity constraint enforced, which envolves solving a global least squares problem.}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{sol.wcm}}, \code{\link{breakfast}}
#' @examples
#' x <- c(rep(0, 100), rep(1, 100), rep(0, 100)) + rnorm(300)
#' model.fixednum(sol.wbs(x),2)
#' model.fixednum(sol.not(x),2)
model.fixednum <- function(cptpath.object, fixednum=NULL) {
  
    ret <- list()
    class(ret) <- "cptmodel"
    
    if(!("cptpath" %in%  class(cptpath.object)))  stop("A cptmodel class object has to be supplied in the first argument.")
    if(cptpath.object$method == 'idetect_seq') warning('To get the Isolate-Detect method results when the model selection is "ic"
                                                        we recommend that you use sol.idetect to create the cptpath.object instead of
                                                        sol.idetect_seq')
    ret$solution.path = cptpath.object$method
    ret$model.selection <- "fixednum"
    
     
    if(is.null(fixednum)) {
      fixednum = 0
      warning("The number of change-points not supplied, zero is used by default")
    }
    fixednum <- as.integer(fixednum)
    if (!(is.integer(fixednum)) || fixednum < 0)
      stop("The number of change-points must be non-negative")
     
    if (length(cptpath.object$x)==0){
      ret$cpts<-integer(0)
      ret$no.of.cpt<-integer(0)
      ret$est<-numeric(0)
    }

    if (cptpath.object$solutions.nested == FALSE){
      id <- which(sapply(cptpath.object$solution.set,length) == fixednum)
      if(length(id)==0){
        if(fixednum != 0) warning("Cannot find the required element on the solution path; a solution with no change-point shall be returned.")
        fixednum = 0
      }
    } else {
      if(length(cptpath.object$solution.path)<fixednum){
        fixednum = 0
        warning("Cannot find the required element on the solution path;  a solution with no change-point shall be returned.")
      }
    }
      
      
    if(fixednum == 0){
      ret$cpts <- 0
      ret$no.of.cpt <- 0
      ret$est <- prediction(cptpath.object, 0)
    }
    else{
      if(cptpath.object$solutions.nested == FALSE){
        lh <- rep(0, length(id))
        for (j in 1:length(id)) {
          cpt.cand <- sort(cptpath.object$solution.set[[id[j]]])
          lh[j] <- logLik(cptpath.object, cpts = cpt.cand)
        }
        K.min.rss <- which.max(lh)
        ret$cpts <- sort(cptpath.object$solution.set[[id[K.min.rss]]])  
      }
      else{
        ret$cpts <-  sort(cptpath.object$solution.path[1:fixednum])
      }
      ret$no.of.cpt <- length(ret$cpts)
      ret$est <- prediction(cptpath.object, ret$cpts)
    }
      
    return(ret)
}
