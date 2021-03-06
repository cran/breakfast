#' @title Breakfast: Methods for Fast Multiple Change-point Detection and Estimation
#' @description  A developing software suite for multiple change-point detection/estimation (data segmentation) in data sequences. 
#' @details
#' The current version implements the Gaussian mean-shift model, 
#' in which the data are assumed to be a piecewise-constant signal observed with i.i.d. Gaussian noise. 
#' Change-point detection in breakfast is carried out in two stages: 
#' (i) computation of a solution path, and (ii) model selection along the path. 
#' A variety of solution path and model selection methods are included, which can be accessed individually,
#' or through \link[breakfast]{breakfast}.
#' Currently supported solution path methods are: \link[breakfast]{sol.idetect}, \link[breakfast]{sol.idetect_seq},
#' \link[breakfast]{sol.wbs}, \link[breakfast]{sol.wbs2}, \link[breakfast]{sol.not} and \link[breakfast]{sol.tguh}.
#' 
#' Currently supported model selection methods are: \link[breakfast]{model.ic}, \link[breakfast]{model.lp},
#' \link[breakfast]{model.sdll} \link[breakfast]{model.thresh}.
#' 
#' Check back future versions for more change-point models and further methods.
#' @useDynLib breakfast, .registration = TRUE 
#' @author \itemize{
#' \item \href{https://www.andreasanastasiou-statistics.com/}{Andreas Anastasiou}
#' \item \href{http://personal.lse.ac.uk/cheny100/}{Yining Chen}
#' \item \href{https://sites.google.com/view/haeran-cho/}{Haeran Cho}
#' \item \href{http://stats.lse.ac.uk/fryzlewicz/}{Piotr Fryzlewicz}
#' }
#' 
#' We would like to thank Shakeel Gavioli-Akilagun, Anica Kostic, Shuhan Yang and Christine Yuen for their comments and suggestions that helped improve this package.
#' 
#' @seealso \code{browseVignettes(package = "breakfast")} contains a detailed comparative simulation study of various methods 
#' implemented in \link[breakfast]{breakfast} for the Gaussian mean-shift model.
#'
#' @docType package
#' @name breakfast-package
NULL