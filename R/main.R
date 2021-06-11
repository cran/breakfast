#' @title Methods for fast multiple change-point detection and estimation
#' @description This function estimates the number and locations of change-points in a data sequence,
#' which is modelled as a piecewise-constant function plus i.i.d. Gaussian noise. 
#' This is carried out via a two-stage procedure combining solution path generation and model selection methodologies.
#' @details Please also take a look at the vignette for tips/suggestions/examples of using the breakfast package.
#' @param x A numeric vector containing the data to be processed
#' @param solution.path A string or a vector of strings containing the name(s) of solution path generating method(s);
#' if individual methods are accessed via this option, default tuning parameters are used. 
#' Alternatively, you can directly access each solution path generating method via \code{sol.[method]}, see below.
#' If both \code{solution.path} and \code{model.selection} are unspecified, we return the output from the suggested combinations based on their performance, which are: \code{("idetect", "ic")}, \code{("idetect_seq", "thresh")}, \code{("not", "ic")}, \code{("tguh", "lp")}, \code{("wbs", "ic")} and \code{("wbs2", "sdll")}.
#' If \code{solution.path} is specified but \code{model.selection} is not, we return the output from the specified \code{solution.path} methods combined with the suggested model selection methods (respectively) as above.
#' \itemize{
#' \item{"idetect"}{ IDetect, see \link[breakfast]{sol.idetect}} 
#' \item{"idetect_seq"}{ Sequential IDetect, see \link[breakfast]{sol.idetect_seq}} 
#' \item{"not"}{ Narrowest-Over-Threshold, see \link[breakfast]{sol.not}}
#' \item{"tguh"}{ Tail-Greedy Unbalanced Haar, see \link[breakfast]{sol.tguh}}
#' \item{"wbs"}{ Wild Binary Segmentation, see \link[breakfast]{sol.wbs}}
#' \item{"wbs2"}{ Wild Binary Segmentation 2, see \link[breakfast]{sol.wbs2}}
#' \item{"all"}{ All of the above }
#' }
#' @param model.selection A string or a vector of strings containing the name(s) of model selection method(s);
#' if individual methods are accessed via this option, default tuning parameters are used. 
#' Alternatively, you can directly access each model selection method via \code{model.[method]}, see below.
#' If both \code{solution.path} and \code{model.selection} are unspecified, we return the output from the suggested combinations based on their performance, which are: \code{("idetect", "ic")}, \code{("idetect_seq", "thresh")}, \code{("not", "ic")}, \code{("tguh", "lp")}, \code{("wbs", "ic")} and \code{("wbs2", "sdll")}.
#' If \code{model.selection} is specified but \code{solution.path} is not, we return the output from the specified \code{model.selection} methods combined with the suggested solution path methods (respectively) as above.
#' \itemize{
#' \item{"ic"}{ Strengthened Schwarz information criterion, see \link[breakfast]{model.ic}}
#' \item{"lp"}{ Localised pruning, see \link[breakfast]{model.lp}}
#' \item{"sdll"}{ Steepest Drop to Low Levels method, see \link[breakfast]{model.sdll}}
#' \item{"thresh"}{ Thresholding, see \link[breakfast]{model.thresh}}
#' \item{"all"}{ All of the above }
#' }
#' @return An S3 object of class \code{breakfast.cpts}, which contains the following fields:
#' \itemize{
#' \item{x}{ Input vector \code{x}}
#' \item{cptmodel.list}{ A list containing S3 objects of class \code{cptmodel}; each contains the following fields:}
#' \itemize{
#'   \item{solution.path}{ The solution path method used}
#'   \item{model.selection}{ The model selection method used to return the final change-point estimators object}
#'   \item{no.of.cpt}{ The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#'   \item{cpts}{ The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#'   \item{est}{ An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#'   }
#' }
#' @references A. Anastasiou & P. Fryzlewicz (2019). Detecting multiple generalized change-points by isolating single ones. \emph{arXiv preprint arXiv:1901.10852}.
#' @references R. Baranowski, Y. Chen & P. Fryzlewicz (2019). Narrowest-over-threshold detection of multiple change points and change-point-like features. \emph{Journal of the Royal Statistical Society: Series B}, 81(3), 649--672.
#' @references H. Cho & C. Kirch (2021) Two-stage data segmentation permitting multiscale change points, heavy tails and dependence. \emph{arXiv preprint arXiv:1910.12486}.
#' @references P. Fryzlewicz (2014). Wild binary segmentation for multiple change-point detection. \emph{The Annals of Statistics}, 42(6), 2243--2281.
#' @references P. Fryzlewicz (2020). Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{To appear in Journal of the Korean Statistical Society}.
#' @references P. Fryzlewicz (2018). Tail-greedy bottom-up data decompositions and fast multiple change-point detection. \emph{The Annals of Statistics}, 46(6B), 3390--3421.
#' @examples 
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f)) * .5
#' breakfast(x)
#' @export
breakfast <- function(x, solution.path = NULL, model.selection = NULL) {
  
  #veryfing the input parameters - x
  x <- as.numeric(x)
  storage.mode(x) <- "double"
  n <- length(x)
  check.input(x)

  # length
  # if piecewise const length >= 3
  if(n <= 2) stop("x should contain at least three elements")
  # if piecewise linear length
  # if(n <= 5) stop("x should contain at least five elements")

  # calculate variance
  # check for zero variance
  # if(!is.null(sd.fn)){
  #   # sigma contains the sqrt of variance to be passed on
  #   if(is.numeric(sd.fn) && sd.fn > 0) sigma <- sd.fn else sigma <- sd.fn(x) 
  # } else{
  #   sigma <- mad(diff(x))/sqrt(2)
  # }
  # for model selection methods requiring sigma (thresh, sdll), we pass on this value (I have crudely implemented this)
  
  # input argument checks
  # stopifnot(is.integer(seed))
  
  option <- 'user'
  
  if(is.null(solution.path) && is.null(model.selection)) option <- 'suggested'
  
  if(option == 'user'){
    if(!is.null(solution.path)) stopifnot(solution.path %in% c('idetect', 'idetect_seq', 'wbs', 'wbs2', 'not', 'tguh', 'all'))
    if(!is.null(model.selection)) stopifnot(model.selection %in% c('ic', 'thresh', 'sdll', 'lp', 'all'))
  }
  
  if(option %in% c('suggested', 'all')) solution.path <- model.selection <- 'all'
  
  if(!is.null(solution.path) && is.null(model.selection)){
    option <- 'suggested'; model.selection <- 'all'
  }
  if(is.null(solution.path) && !is.null(model.selection)){
    option <- 'suggested'; solution.path <- 'all'
  }
  
  if('all' %in% solution.path) solution.path <- c('idetect', 'idetect_seq', 'wbs', 'wbs2', 'not', 'tguh')
  if('all' %in% model.selection) model.selection <- c('ic', 'thresh', 'sdll', 'lp')

  no.sol.path <- length(solution.path)
  no.mod.select <- length(model.selection)
  
  # generate all solution paths
  all.sol.paths <- list()
  for(l in 1:no.sol.path){
    sp <- solution.path[l]
    args <- list()
    args$x <- x
    out <- list(do.call(paste("sol.", sp, sep = ""), args))
    all.sol.paths <- c(all.sol.paths, out) # each solution path to contain name
  }
    
  # apply all model selection methods to each of the solution paths
  cptmodel.list <- list() 
  no.lp.warning <- 0
  for(l in 1:no.sol.path){
    sp <- all.sol.paths[[l]]
    args <- list()
    args$cptpath.object <- sp
    
    for(m in 1:no.mod.select){
      ms <- model.selection[m]
      
      # exclude non-suggested combo 
      if(option == 'suggested') if(ms == 'ic' & sp$method %in% c('tguh', 'wbs2')) next
      if(option == 'suggested') if(ms == 'lp' & sp$method %in% c('not', 'wbs', 'wbs2')) next
      if(option == 'suggested') if(ms == 'sdll' & sp$method %in% c('idetect', 'tguh', 'wbs')) next
      if(option == 'suggested') if(ms == 'thresh' & sp$method %in% c('not', 'wbs', 'wbs2', 'tguh')) next
      
      # exclude bad combo for suggested/all; if they are run for option = 'user' appropriate warnings generated
      if(option %in% c('suggested', 'all')){
        if(ms == 'ic' & sp$method %in% c('idetect_seq')) next
        if(ms == 'lp' & sp$method %in% c('idetect', 'idetect_seq')) next
        if(ms == 'sdll' & sp$method %in% c('idetect_seq', 'not')) next
        if(ms == 'thresh' & sp$method %in% c('idetect')) next
      }
      
      if(ms == 'lp' & n < 11){
        if(no.lp.warning == 0) warning(paste0('x contains too few elements for model.lp, skipping the model selection method'))
        no.lp.warning <- no.lp.warning + 1
        next
      }
      
      out <- list(do.call(paste("model.", ms, sep = ""), args))
      cptmodel.list <- c(cptmodel.list, out)
    }
  }
  
  # output
  ret <- structure(
    list(x = x,
        # cptpath.list = solution.path,
         cptmodel.list = cptmodel.list),
        class = 'breakfast.cpts')
  return(ret)
}

# recommend <- function(){
#   # make recommendation
# }

#' #' Summary of change-points estimated by BREAKFAST MAIN ROUTINE
#' #' 
#' #' Summary method for objects of class \code{breakfast.cpts}
#' #' @method summary breakfast.cpts
#' #' @param object a \code{breakfast.cpts} object
#' #' @param ... not in use
#' #' @details Provide information about each estimated change-point, 
#' #' including associated CUSUM statistic and its detection interval.
#' #' @examples 
#' summary.breakfast.cpts <- function(object, ...) { 
#'   
#' }
#' 
