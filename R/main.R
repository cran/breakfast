#' @title Methods for fast multiple change-point detection and estimation
#' @description This function estimates the number and locations of change-points in a data sequence,
#' which is modelled as a piecewise-constant function plus i.i.d. Gaussian noise. 
#' This is carried out via a two-stage procedure combining solution path generation and model selection methodologies.
#' @details Please also take a look at the vignette for tips/suggestions/examples of using the breakfast package.
#' @param x A numeric vector containing the data to be processed
#' @param sd.fn Either a positive numeric constant containing the square root of the error variance, 
#' or a function for generating this quantity that takes \code{x} as an input; if \code{sd.fn = NULL}, 
#' a mean absolute deviation estimator is generated internally.
#' @param option A string specifying the option for selecting the combination of
#' solution path and model selection procedures. 
#' The default choice is \code{option = "suggested"} with combinations of solution path and model selection methods
#' suggested for their good properties.
#' \itemize{
#' \item{"suggested"}{ Suggested combinations based on their performance, which are \code{("wbs2", "sdll")}, \code{("idetect", "ic")}, \code{("idetect_seq", "thresh")} and \code{("tguh", "lp")}}
#' \item{"user"}{ User-specified combination (to be supplied by \code{solution.path} and \code{model.selection}) }
#' \item{"all"}{ All currently supported combinations; we exclude the following combinations of the solution path and model selection methods based on their performance:  \code{("not", "sdll")}, \code{("idetect", "thresh")}, \code{("idetect", "lp")}, \code{("idetect_seq", "ic")}, \code{("idetect_seq", "lp")}, \code{("idetect_seq", "sdll")} }
#' }
#' @param solution.path A list containing the name of a solution path generating method
#' and any arguments to be used with them; currently supported methods are
#' \itemize{
#' \item{"idetect"}{ IDetect, see \link[breakfast]{sol.idetect}} 
#' \item{"idetect_seq"}{ Sequential IDetect, see \link[breakfast]{sol.idetect_seq}} 
#' \item{"not"}{ Narrowest-Over-Threshold, see \link[breakfast]{sol.not}}
#' \item{"tguh"}{ Tail-Greedy Unbalanced Haar, see \link[breakfast]{sol.tguh}}
#' \item{"wbs"}{ Wild Binary Segmentation, see \link[breakfast]{sol.wbs}}
#' \item{"wbs2"}{ Wild Binary Segmentation 2, see \link[breakfast]{sol.wbs2}}
#' }
#' @param model.selection A list of containing the name of a model selection method
#' and any arguments to be used with them; currently supported methods are
#' \itemize{
#' \item{"ic"}{ Strengthened Schwarz information criterion, see \link[breakfast]{model.ic}}
#' \item{"lp"}{ Localised pruning, see \link[breakfast]{model.lp}}
#' \item{"sdll"}{ Steepest Drop to Low Levels method, see \link[breakfast]{model.sdll}}
#' \item{"thresh"}{ Thresholding, see \link[breakfast]{model.thresh}}
#' }
#' @return An S3 object of class \code{breakfast.cpts}, which contains the following fields:
#' \itemize{
#' \item{x}{ Input vector \code{x}}
#' \item{cptmodel.list}{ A list containing S3 objects of class \code{cptmodel}; each contains the following fields:}
#' \itemize{
#'   \item{solution.path}{ The solution path method used}
#'   \item{model}{ The model selection method used to return the final change-point estimators object}
#'   \item{no.of.cpt}{ The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#'   \item{cpts}{ The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#'   \item{est}{ An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#'   }
#' }
#' @references A. Anastasiou & P. Fryzlewicz (2019). Detecting multiple generalized change-points by isolating single ones. \emph{arXiv preprint arXiv:1901.10852}.
#' @references R. Baranowski, Y. Chen & P. Fryzlewicz (2019). Narrowest-over-threshold detection of multiple change points and change-point-like features. \emph{Journal of the Royal Statistical Society: Series B}, 81(3), 649--672.
#' @references H. Cho & C. Kirch (2020) Two-stage data segmentation permitting multiscale change points, heavy tails and dependence. \emph{arXiv preprint arXiv:1910.12486}.
#' @references P. Fryzlewicz (2014). Wild binary segmentation for multiple change-point detection. \emph{The Annals of Statistics}, 42(6), 2243--2281.
#' @references P. Fryzlewicz (2020). Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{To appear in Journal of the Korean Statistical Society}.
#' @references P. Fryzlewicz (2018). Tail-greedy bottom-up data decompositions and fast multiple change-point detection. \emph{The Annals of Statistics}, 46(6B), 3390--3421.
#' @examples 
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f))
#' breakfast(x, option = 'suggested')
#' @export
breakfast <- function(x, sd.fn = NULL, 
                      option = c('suggested','user', 'all')[1],
                      solution.path = list(name = 'idetect', args = list()),
                      model.selection = list(name = 'ic', args = list())) {
  
  # check x is numeric
  if(NA %in% x) stop("x vector cannot contain NA's")
  
  n <- length(x)
  # length
  # if piecewise const length >= 3
  if(n <= 2) stop("x should contain at least three elements")
  # if piecewise linear length
  # if(n <= 5) stop("x should contain at least five elements")

  # calculate variance
  # check for zero variance
  if(!is.null(sd.fn)){
    # sigma contains the sqrt of variance to be passed on
    if(is.numeric(sd.fn) && sd.fn > 0) sigma <- sd.fn else sigma <- sd.fn(x) 
  } else{
    sigma <- mad(diff(x))/sqrt(2)
  }
  # for model selection methods requiring sigma (thresh, sdll), we pass on this value (I have crudely implemented this)
  
  # input argument checks
  # stopifnot(is.integer(seed))
  stopifnot(option %in% c('suggested','user', 'all'))

  if(option == 'user'){
    stopifnot(solution.path$name %in% c('idetect', 'idetect_seq', 'wbs', 'wbs2', 'not', 'tguh'))
    stopifnot(model.selection$name %in% c('ic', 'thresh', 'sdll', 'lp'))
    solution.path <- list(solution.path)
    model.selection <- list(model.selection)
    # if(model.selection[[k]]$name %in% c('thresh', 'sdll')) model.selection[[k]]$args$sigma <- sigma
#    for(k in 1:no.sol.path) {if((solution.path[[k]] == 'idetect') & (model.selection[[k]] == 'thresh')) {solution.path[[k]] <- 'idetect2'}
#      if((solution.path[[k]] == 'idetect_seq') & (model.selection[[k]] == 'ic')) {solution.path[[k]] <- 'idetect'}
#      if((solution.path[[k]] == 'idetect_seq') & (model.selection[[k]] == 'sdll')) {solution.path[[k]] <- 'idetect'}
#      }
  }

  if(option %in% c('suggested', 'all')){
    solution.path <- list(list(name = 'wbs'), list(name = 'wbs2'), list(name = 'idetect'), list(name = 'not'), list(name = 'tguh'), list(name = 'idetect_seq'))
    model.selection <- list(list(name = 'ic'), list(name = 'thresh', args = list(sigma = sigma)), list(name = 'sdll', args = list(sigma = sigma)), list(name = 'lp'))
  }
  
  no.sol.path <- length(solution.path)
  no.mod.select <- length(model.selection)
  
  # if(option == 'fastest'){ # have we decided on a method?
  #   solution.path <- list(name = 'idetect')
  #   model.selection <- list(name = 'ic')
  # }
    
  # generate all solution paths
  all.sol.paths <- list()
  for(l in 1:no.sol.path){
    sp <- solution.path[[l]]
    args <- list()
    args$x <- x
    if(option == 'user') args <- append(args, sp$args)
    all.sol.paths <- c(all.sol.paths, list(do.call(paste("sol.", sp$name, sep = ""), args))) # each solution path to contain name
  }
    
  # apply all model selection methods to each of the solution paths
  cptmodel.list <- list() 
  no.lp.warning <- 0
  for(l in 1:no.sol.path){
    sp <- all.sol.paths[[l]]
    for(m in 1:no.mod.select){
      ms <- model.selection[[m]]
      
      # exclude non-suggested combo 
      if(option == 'suggested') if(ms$name == 'ic' & sp$method %in% c('not', 'tguh', 'wbs', 'wbs2')) next
      if(option == 'suggested') if(ms$name == 'lp' & sp$method %in% c('not', 'wbs', 'wbs2')) next
      if(option == 'suggested') if(ms$name == 'sdll' & sp$method %in% c('idetect', 'tguh', 'wbs')) next
      if(option == 'suggested') if(ms$name == 'thresh' & sp$method %in% c('not', 'wbs', 'wbs2', 'tguh')) next
      
      # exclude bad combo -- if(c(sp$name, ms$name) not supported
      if(option %in% c('suggested', 'all')) if(ms$name == 'ic' & sp$method %in% c('idetect_seq')) next
      if(option %in% c('suggested', 'all')) if(ms$name == 'lp' & sp$method %in% c('idetect', 'idetect_seq')) next
      if(option %in% c('suggested', 'all')) if(ms$name == 'sdll' & sp$method %in% c('idetect_seq', 'not')) next
      if(option %in% c('suggested', 'all')) if(ms$name == 'thresh' & sp$method %in% c('idetect')) next
      
      if(ms$name == 'thresh' || ms$name == 'sdll') args <- list(sigma = sigma) else args <- list()
      
      if(ms$name == 'lp' & n < 11){
        if(no.lp.warning == 0) warning(paste0('x contains too few elements for model.lp, skipping the model selection method'))
        no.lp.warning <- no.lp.warning + 1
        next
      }
      
      args$cptpath.object <- sp # make sure each solution path contains x, make sure each model selection inputs object
      
      if(option == 'user') args <- append(args, ms$args)
      out <- list(do.call(paste("model.", ms$name, sep = ""), args))
      cptmodel.list <- c(cptmodel.list, (out))
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
