#' @title Estimating change-points in the piecewise-constant mean of a noisy data sequence via the strengthened Schwarz information criterion
#' @description This function estimates the number and locations of change-points in the piecewise-constant mean of a noisy data sequence via the sSIC (strengthened Schwarz information criterion) method.
#' @details 
#' The model selection method for algorithms that produce nested solution path is described in 
#' "Wild binary segmentation for multiple change-point detection", P. Fryzlewicz (2014), The Annals of Statitics, 42: 2243--2281.
#' The corresponding description for those that produce non-nested solution set can be found in 
#' "Narrowest-over-threshold detection of multiple change points and change-point-like features", R. Baranowski, Y. Chen and P. Fryzlewicz (2019), Journal of Royal Statistical Society: Series B, 81(3), 649--672.
#' 
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. Note that the field \code{cptpath.object$x} contains the input data sequence. 
#' @param alpha The parameter associated with the sSIC. The default value is 1.01. Note that the SIC is recovered when alpha = 1.
#' @param q.max The maximum number of change-points allowed. If nothing or \code{NULL} is provided, the default value of around \code{min(25, n/log(n))} will be used.
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{model}{The model selection method used to return the final change-point estimators object, here its value is \code{"ic"}}
#' \item{no.of.cpt}{The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#' \item{cpts}{The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#' \item{est}{An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{breakfast}}
#' @references P. Fryzlewicz (2014). Wild binary segmentation for multiple change-point detection. \emph{The Annals of Statistics}, 42(6), 2243--2281.
#' @references R. Baranowski, Y. Chen & P. Fryzlewicz (2019). Narrowest-over-threshold detection of multiple change points and change-point-like features. \emph{Journal of the Royal Statistical Society: Series B}, 81(3), 649--672.
#' @examples
#' x <- c(rep(0, 100), rep(1, 100), rep(0, 100)) + rnorm(300)
#' model.ic(sol.wbs(x))
#' model.ic(sol.not(x))
model.ic <- function(cptpath.object, alpha=1.01, q.max = NULL) {
  
    ret <- list()
    class(ret) <- "cptmodel"
    
    if(class(cptpath.object)!="cptpath") stop("A cptmodel class object has to be supplied in the first argument.")
    if(cptpath.object$method == 'idetect_seq') warning('To get the Isolate-Detect method results when the model selection is "ic"
                                                        we recommend that you use sol.idetect to create the cptpath.object instead of
                                                        sol.idetect_seq')
    ret$solution.path = cptpath.object$method
    ret$model <- "ic"
    
    if(is.null(q.max)) q.max.null=TRUE
    else q.max.null=FALSE
    
    if(is.null(q.max)) q.max = min(floor(length(cptpath.object$x)/max(1, log(length(cptpath.object$x)))),100)
    q.max <- as.integer(q.max)
    if (length(q.max) != 1)
      stop("q.max must be a single number.")
     
    if (length(cptpath.object$x)==0){
      ret$cpts<-integer(0)
      ret$no.of.cpt<-integer(0)
      ret$est<-numeric(0)
    }
    
    if(q.max > length(cptpath.object$x)-1) q.max = length(cptpath.object$x)-1
    
    if( floor(length(cptpath.object$x)/max(1, log(length(cptpath.object$x))))<q.max){
      warning("The information criterion approach might not work well for short signals. Please consider reducing the value of q.max.")
    }
    
    penalty.fun <- sic.penalty 
    
    if (cptpath.object$solutions.nested == FALSE){
      id <- which(sapply(cptpath.object$solution.set,length) <= q.max)
    } else id <- 1:q.max
    
    resic <- rep(0, length(id) + 1)
    
    lh <- logLik(cptpath.object, cpts = c())
    n.param <- attr(lh, "df")
    
    resic[1] <-
      -2 * lh + penalty.fun(n = length(cptpath.object$x),
                            n.param = n.param, alpha=alpha)
    
    #calculate ic = log.lik + pen
    if(length(resic)>=2){
      
    for (j in 2:length(resic)) {
      if (cptpath.object$solutions.nested == FALSE){
        cpt.cand <- sort(cptpath.object$solution.set[[id[j - 1]]])
      } else cpt.cand <- sort(cptpath.object$solution.path[1:(j-1)])
      
      if(length(cpt.cand) == 0){ret$cpts <- 0
      ret$no.of.cpt <- 0
      ret$est <- prediction(cptpath.object, 0)
      return(ret)}
      
      lh <- logLik(cptpath.object, cpts = cpt.cand)
      n.param <- attr(lh, "df")
      resic[j] <-
        -2 * lh + penalty.fun(n = length(cptpath.object$x),
                              n.param = n.param, alpha=alpha)
    }
    }
    K.min.ic <- which.min(resic) - 1
    ret$cpts <- list()
    
    if (K.min.ic == 0) {
      ret$cpts <- 0
      ret$no.of.cpt <- 0
      ret$est <- prediction(cptpath.object, 0)
    }
    else {
      if (cptpath.object$solutions.nested == FALSE) {
        ret$cpts <- sort(cptpath.object$solution.set[[id[K.min.ic]]])
      } else 
      ret$cpts <- sort(cptpath.object$solution.path[1:id[K.min.ic]])
      
      ret$no.of.cpt <- length(ret$cpts)
      ret$est <- prediction(cptpath.object, ret$cpts)
    }
    
    if (q.max == ret$no.of.cpt & q.max.null) warning('IC estimated the maximum default permitted number of change-points. It is possible that there are more. Please use other methods or specifically supply an upper bound of the number of change-points in the options of model.ic.')
    if (q.max == ret$no.of.cpt & !q.max.null) warning('IC estimated the maximum permitted number of change-points. It is possible that there are more. Please use other methods or increase the upper bound of the number of change-points in model.ic.')
    
    return(ret)
    
}

sic.penalty <- function(n, n.param, alpha=1.0){
  alpha <- as.numeric(alpha)
  
  pen <- log(n)^alpha
  
  return(n.param*pen)
}



logLik <- function(cptpath.object, cpts){
  
  if(missing(cpts)) cpts <- model.ic(cptpath.object)$cpts
  if(!is.null(cpts)) if(any(is.na(cpts))) cpts <- cpts[!is.na(cpts)]
  
  
  n.cpt <- length(cpts)
  
  residuals <- cptpath.object$x -  prediction(cptpath.object, cpts = cpts)
  
  ret <- - length(cptpath.object$x)/2 * log(mean(residuals^2))
    
  attr(ret, "df") <- 2*n.cpt + 2
  
  return(ret)
  
}


prediction <- function(cptpath.object, cpts) {
  if (length(cptpath.object$x)==0) return(numeric(0))
  if (missing(cpts))
    cpts <- model.ic(cptpath.object)$cpts
  if (!is.null(cpts))
    if (any(is.na(cpts)))
      cpts <- cpts[!is.na(cpts)]
    
    cpts <- sort(unique(c(cpts, 0, length(cptpath.object$x))))
    
    fit <- rep(0, length(cptpath.object$x))
      
    for (i in 1:(length(cpts) - 1)) {
        fit[(cpts[i] + 1):cpts[i + 1]] <- mean(cptpath.object$x[(cpts[i] + 1):cpts[i + 1]])
        
    }
    
    return(fit)
    
}


