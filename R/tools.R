#' @title Change-points estimated by breakfast
#' @description  Print method for objects of class \code{breakfast.cpts}
#' @method print breakfast.cpts
#' @param x a \code{breakfast.cpts} object
#' @param ... not in use
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 5)
#' x <- f + rnorm(length(f))
#' print(breakfast(x, option = 'user'))
print.breakfast.cpts <- function(x, ...) {
  L <- length(x$cptmodel.list)
  brks <- rep(0, length(x$x))
  mm <- character(length(x$x))
  all.nm <- character(0)
  for(l in 1:L){
    cl <- x$cptmodel.list[[l]]
    nm <- paste(cl$solution.path, '.', cl$model, sep = '')
    if(l == 1) all.nm <- nm else all.nm <- paste(all.nm, ', ', nm, sep = '')
    if(cl$no.of.cpt > 0){
      brks[cl$cpts] <- brks[cl$cpts] + 1
      for(ii in cl$cpts) if(brks[ii] == 1) mm[ii] <- nm else mm[ii] <- paste(mm[ii], nm, sep = ', ')
    }
  }
  if(sum(brks > 0) == 0) cat(paste('No change point is found')) else{
    cat(paste('Change-point locations estimated by: ', all.nm, sep = ''))
    cat('\n\n')
    for(ii in which(brks > 0)){
      cat(paste(ii, ': ', mm[ii], sep = ''))
      cat('\n')
    }
  }
}


#' #' Change-points estimated by breakfast
#' #'
#' #' Plot method for objects of class \code{breakfast.cpts}
#' #' @method plot breakfast.cpts
#' #' @param x a \code{breakfast.cpts} object
#' #' @param ... not in use
#' #' @importFrom ggplot2
#' #' @examples
#' #' f <- rep(rep(c(0, 1), each = 50), 5)
#' #' x <- f + rnorm(length(f))
#' #' plot(breakfast(x, option = 'user'))
#' plot.breakfast.cpts <- function(x, ...) {
#'   L <- length(x$cptmodel.list)
#'   brks <- matrix(0, nrow = length(x$x), ncol = L)
#'   all.nm <- character(L)
#'   for(l in 1:L){
#'     cl <- x$cptmodel.list[[l]][[1]]
#'     nm <- paste(cl$solution.path, '.', cl$model, sep = '')
#'     all.nm[l] <- nm
#'     if(cl$no.of.cpt > 0) brks[cl$cpts, l] <- brks[cl$cpts, l] + 1
#'   }
#' 
#'   if(sum(brks > 0) == 0){
#'     df <- data.frame()
#'     ggplot(df) + geom_point() + xlim(0, length(x$x)) + ylim(0, L) + ggtitle('Estimated change-point locations') + theme_classic()
#'   } else{
#'     cpts <- which(apply(brks, 1, sum) > 0)
#'     location <- rep(cpts, each = L)
#'     method <- rep(all.nm, length(cpts))
#'     value <- c(t(brks[cpts,, drop = FALSE]))
#'     df <- data.frame(location, method, value)
#'     ggplot(df, aes(fill = method, y = value, x = location)) +
#'       geom_bar(position = "stack", stat = "identity", width = 1) +
#'       theme(axis.text.y=element_blank(), axis.title.y = element_blank()) +
#'       ggtitle('Estimated change-point locations') +
#'       theme_classic()
#'   }
#' }

#plot.cptpath <- function(x, max.cpts = 10, max.lwd = 10, cpt.col="red", ...) {
#
#	a <- list(...)
#	if (!("ylab" %in% names(a))) a$ylab <- "Data"
#	ts.plot(x$data, gpars=a)
#	m <- length(x$solution.path)
#	if (m) abline(v=x$solution.path[1:min(m, max.cpts)], col=cpt.col, lwd=ceiling(x$sorted.cusums[1:min(m, max.cpts),4]/x$sorted.cusums[1,4]*max.lwd))
#		
#}
