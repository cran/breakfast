#' breakfast: Multiple change-point detection and segmentation for data sequences
#'
#' The breakfast package performs multiple change-point detection in data sequences, or sequence
#' segmentation, using computationally efficient multiscale methods. This version of the
#' package implements the "Tail-Greedy Unbalanced Haar", "Wild Binary Segmentation" and
#' "Adaptive Wild Binary Segmentation" change-point detection and segmentation
#' methodologies. To start with, see the function
#' \code{segment.mean}.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @references "Tail-greedy bottom-up data decompositions and fast multiple change-point detection", P. Fryzlewicz (2017), preprint.
#' "Wild Binary Segmentation for multiple change-point detection", P. Fryzlewicz (2014), Annals of Statistics, 42, 2243-2281.
#' "Data-adaptive Wild Binary Segmentation", P. Fryzlewicz (2017), in preparation as of September 28th, 2017.
#' @seealso \code{\link{segment.mean}}
#' @examples
#' #See Examples for segment.mean
#' @docType package
#' @name breakfast
NULL
