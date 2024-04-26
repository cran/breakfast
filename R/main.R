#' @title Methods for fast detection of multiple change-points
#' @description This function estimates the number and locations of change-points in a univariate data sequence,
#' which is modelled as (i) a piecewise-constant function plus i.i.d. Gaussian noise, 
#' (ii) a piecewise-constant function plus autoregressive time series,  
#' (iii) a piecewise-linear and continuous function plus i.i.d. Gaussian noise, or
#' (iv) a piecewise-linear and discontinuous function plus i.i.d. Gaussian noise.
#' This is carried out via a two-stage procedure combining solution path generation and model selection methodologies.
#' @details Please also take a look at the vignette for tips/suggestions/examples of using the breakfast package.
#' @param x A numeric vector containing the data to be processed
#' @param type The type of change-point models fitted to the data; currently supported models are: piecewise constant signals (\code{type = "const"}, chosen by default), piecewise linear and continuous signals (\code{type = "lin.cont"}) and piecewise linear and discontinuous signals (\code{type = "lin.discont"}).
#' @param solution.path A string or a vector of strings containing the name(s) of solution path generating method(s);
#' if individual methods are accessed via this option, default tuning parameters are used. 
#' Alternatively, you can directly access each solution path generating method via \code{sol.[method]}.
#' If both \code{solution.path} and \code{model.selection} are unspecified, we return the output from the suggested combinations based on their performance, which depends on \code{type} as below:
#' 
#' When \code{type = "const"}: \code{("idetect", "ic")}, \code{("idetect_seq", "thresh")}, 
#' \code{("not", "ic")}, \code{("tguh", "lp")}, \code{("wbs", "ic")}, \code{("wbs2", "sdll")} and \code{("wcm", "gsa")}.
#' 
#' When \code{type = "lin.cont"} or \code{type = "lin.discont"}: \code{("idetect_seq", "thresh")}, 
#' \code{("not", "ic")} and \code{("idetect", "sdll")}.
#' 
#' If \code{solution.path} is specified but \code{model.selection} is not, we return the output from the specified \code{solution.path} methods combined with the suggested model selection methods (respectively) as above.
#' \describe{
#' \item{"idetect"}{ IDetect, supporting \code{type = "const", "lin.cont", "lin.discont"}, see \link[breakfast]{sol.idetect}} 
#' \item{"idetect_seq"}{ Sequential IDetect, supporting \code{type = "const", "lin.cont", "lin.discont"}, see \link[breakfast]{sol.idetect_seq}} 
#' \item{"not"}{ Narrowest-Over-Threshold, supporting \code{type = "const", "lin.cont", "lin.discont"}, see \link[breakfast]{sol.not}}
#' \item{"tguh"}{ Tail-Greedy Unbalanced Haar, supporting \code{type = "const"}, see \link[breakfast]{sol.tguh}}
#' \item{"wbs"}{ Wild Binary Segmentation, supporting \code{type = "const"}, see \link[breakfast]{sol.wbs}}
#' \item{"wbs2"}{ Wild Binary Segmentation 2, supporting \code{type = "const"}, see \link[breakfast]{sol.wbs2}}
#' \item{"wcm"}{ Wild Contrast Maximisation, supporting \code{type = "const"} in combination with \link[breakfast]{model.gsa} handling model (ii), see \link[breakfast]{sol.wcm}}
#' \item{"all"}{ All of the above that support the \code{type} }
#' }
#' 
#' @param model.selection A string or a vector of strings containing the name(s) of model selection method(s);
#' if individual methods are accessed via this option, default tuning parameters are used. 
#' Alternatively, you can directly access each model selection method via \code{model.[method]}.
#' If both \code{solution.path} and \code{model.selection} are unspecified, we return the output from the suggested combinations based on their performance, see \code{solution.path}.
#' If \code{model.selection} is specified but \code{solution.path} is not, we return the output from the specified \code{model.selection} methods combined with the suggested solution path methods (respectively).
#' Not all \code{solution.path} methods are supported by all \code{model.selection} methods; check the individual functions for more information.
#' \describe{
#' \item{"ic"}{ Strengthened Schwarz information criterion, supporting \code{type = "const", "lin.cont", "lin.discont"}, see \link[breakfast]{model.ic}}
#' \item{"lp"}{ Localised pruning, supporting \code{type = "const"}, see \link[breakfast]{model.lp}}
#' \item{"sdll"}{ Steepest Drop to Low Levels method, supporting \code{type = "const", "lin.cont", "lin.discont"}, see \link[breakfast]{model.sdll}}
#' \item{"thresh"}{ Thresholding, supporting \code{type = "const", "lin.cont", "lin.discont"}, see \link[breakfast]{model.thresh}}
#' \item{"gsa"}{ gappy Schwarz algorithm, supporting \code{type = "const"} in combination with \link[breakfast]{sol.wcm} handling model (ii), see \link[breakfast]{model.gsa}}
#' \item{"all"}{ All of the above that support the given \code{type}}
#' }
#' 
#' @return An S3 object of class \code{breakfast.cpts}, which contains the following fields:
#' \describe{
#' \item{x}{ Input vector \code{x}}
#' \item{cptmodel.list}{ A list containing S3 objects of class \code{cptmodel}; each contains the following fields:}
#' \describe{
#'   \item{solution.path}{ The solution path method used}
#'   \item{model.selection}{ The model selection method used to return the final change-point estimators object}
#'   \item{no.of.cpt}{ The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#'   \item{cpts}{ The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#'   \item{est}{ An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#'   }
#' }
#' @references A. Anastasiou & P. Fryzlewicz (2022) Detecting multiple generalized change-points by isolating single ones. \emph{Metrika}, 85(2), 141--174.
#' @references R. Baranowski, Y. Chen & P. Fryzlewicz (2019) Narrowest-over-threshold detection of multiple change points and change-point-like features. \emph{Journal of the Royal Statistical Society: Series B}, 81(3), 649--672.
#' @references H. Cho & C. Kirch (2022) Two-stage data segmentation permitting multiscale change points, heavy tails and dependence. \emph{Annals of the Institute of Statistical Mathematics}, 74(4), 653--684.
#' @references H. Cho & P. Fryzlewicz (2024) Multiple change point detection under serial dependence: Wild contrast maximisation and gappy Schwarz algorithm. \emph{Journal of Time Series Analysis}, 45(3): 479--494.
#' @references P. Fryzlewicz (2014) Wild binary segmentation for multiple change-point detection. \emph{The Annals of Statistics}, 42(6), 2243--2281.
#' @references P. Fryzlewicz (2018) Tail-greedy bottom-up data decompositions and fast multiple change-point detection. \emph{The Annals of Statistics}, 46(6B), 3390--3421.
#' @references P. Fryzlewicz (2020) Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{Journal of the Korean Statistical Society}, 49(4), 1027--1070.
#' @examples 
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f)) * .5
#' breakfast(x)
#' @export
breakfast <- function(x, type = c('const', 'lin.cont', 'lin.discont'), solution.path = NULL, model.selection = NULL) {
  
  # verifyng the input parameters - x
  x <- as.numeric(x)
  storage.mode(x) <- "double"
  n <- length(x)
  check.input(x)
  
  type <- match.arg(type)
  
  stopifnot(type == 'const' || type == "lin.cont" || type == "lin.discont")
  
  # length
  # if piecewise const length >= 3
  if(type == "const" & n <= 2) stop("x should contain at least three elements")
  # if piecewise linear length
  if((type == "lin.cont" || type == "lin.discont") & n <= 5) stop("x should contain at least five elements")
  
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
  
  if(is.null(solution.path) || is.null(model.selection)) {
    option <- 'suggested'
    if(is.null(solution.path)) solution.path <- 'all'
    if(is.null(model.selection)) model.selection <- 'all'
  } else if(solution.path == 'all' || model.selection == 'all'){
    option <- 'all'
  } else option <- 'user'
  
  if(type == 'const'){
    stopifnot(solution.path %in% c('idetect', 'idetect_seq', 'wbs', 'wbs2', 'not', 'tguh', 'wcm', 'all'))
    stopifnot(model.selection %in% c('ic', 'thresh', 'sdll', 'lp', 'gsa', 'all'))
  } else if(type == 'lin.cont' || type == 'lin.discont'){
    stopifnot(solution.path %in% c('idetect', 'idetect_seq', 'not', 'all'))
    stopifnot(model.selection %in% c('ic', 'thresh', 'sdll', 'all'))
  }
  
  if('all' %in% solution.path){
    if(type == 'const'){
      solution.path <- c('idetect', 'idetect_seq', 'wbs', 'wbs2', 'not', 'tguh', 'wcm')
      # if(option == 'suggested') solution.path <- setdiff(solution.path, 'wcm')
    }
    if(type == 'lin.cont' || type == 'lin.discont') solution.path <- c('idetect_seq', 'idetect', 'not')
  } 
  
  if('all' %in% model.selection){
    if(type == 'const') model.selection <- c('ic', 'thresh', 'sdll', 'lp', 'gsa')
    if(type == 'lin.cont' || type == 'lin.discont') model.selection <- c('ic', 'thresh', 'sdll')
  }
  
  no.sol.path <- length(solution.path)
  no.mod.select <- length(model.selection)
  
  # generate all solution paths
  all.sol.paths <- list()
  no.wcm.warning <- 0
  for(l in 1:no.sol.path){
    sp <- solution.path[l]
    
    if(sp == 'wcm' && n < 2 * max(20, 10 + ceiling(log(n)^1.1)) + 1){
      if(no.wcm.warning == 0) warning(paste0('x is too short for sol.wcm, skipping the solution path method'))
      no.wcm.warning <- no.wcm.warning + 1
      next
    }
    
    args <- list()
    args$x <- x
    args$type <- type
    out <- list(do.call(paste("sol.", sp, sep = ""), args))
    all.sol.paths <- c(all.sol.paths, out) # each solution path to contain name
  }
  
  # apply all model selection methods to each of the solution paths
  cptmodel.list <- list() 
  no.lp.warning <- 0
  if(length(all.sol.paths) >= 1){
    for(l in 1:length(all.sol.paths)){
      sp <- all.sol.paths[[l]]
      # args <- list()
      # args$cptpath.object <- sp
      # args$type <- NULL
      args <- sp
      
      for(m in 1:no.mod.select){
        ms <- model.selection[m]
        
        # exclude non-suggested combo 
        if(option == 'suggested'){
          if(type == 'const'){
            if(ms == 'ic' & !(sp$method %in% c('idetect', 'not', 'wbs'))) next
            if(ms == 'lp' & !(sp$method %in% c('tguh'))) next
            if(ms == 'sdll' & !(sp$method %in% c('wbs2'))) next
            if(ms == 'thresh' & !(sp$method %in% c('idetect_seq'))) next
            if(ms == 'gsa' & !(sp$method %in% c('wcm'))) next
          } else if(type %in% c('lin.cont', 'lin.discont')){
            if(ms == 'ic' & !(sp$method %in% c('not'))) next
            if(ms == 'lp') next
            if(ms == 'sdll' & !(sp$method %in% c('idetect'))) next
            if(ms == 'thresh' & !(sp$method %in% c('idetect_seq'))) next
            if(ms == 'gsa') next
          }
        }
        
        # exclude bad combo for all; if they are run for option = 'user' appropriate warnings generated
        if(option == 'all'){
          if(type == 'const'){
            if(ms == 'ic' & sp$method %in% c('idetect_seq', 'wcm')) next
            if(ms == 'lp' & sp$method %in% c('idetect', 'idetect_seq', 'wcm')) next
            if(ms == 'sdll' & sp$method %in% c('idetect_seq', 'not', 'wcm')) next
            if(ms == 'thresh' & sp$method %in% c('idetect', 'wcm')) next
            if(ms == 'gsa' & !(sp$method %in% c('wcm'))) next
          }
          
          if(type %in% c('lin.cont', 'lin.discont')){
            if(ms == 'ic' & !(sp$method %in% c('not', 'idetect', 'idetect_seq'))) next
            if(ms == 'lp') next
            if(ms == 'sdll' & !(sp$method %in% c('idetect'))) next
            if(ms == 'thresh' & !(sp$method %in% c('not', 'idetect', 'idetect_seq'))) next
            if(ms == 'gsa') next
          }
        }
        
        if(ms == 'lp' && n < 11){
          if(no.lp.warning == 0) warning(paste0('x is too short for model.lp, skipping the model selection method'))
          no.lp.warning <- no.lp.warning + 1
          next
        }
      
        out <- list(do.call(paste("model.", ms, sep = ""), list(cptpath.object = args)))
        cptmodel.list <- c(cptmodel.list, out)
      }
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
