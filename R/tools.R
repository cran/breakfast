#' @title Change-points estimated by the "breakfast" routine
#' @description  Print method for objects of class \code{breakfast.cpts}
#' @method print breakfast.cpts
#' @param x a \code{breakfast.cpts} object
#' @param by if \code{by = 'method'}, change-point estimators are printed by method;
#' if \code{by = 'estimator'}, each change-point estimator is printed with the methods that detect it.
#' @param ... current not in use
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 5)
#' x <- f + rnorm(length(f)) * .5
#' print(breakfast(x, solution.path = 'all', model.selection = 'all'), by = 'method')
#' print(breakfast(x), by = 'estimator')
#' @export
print.breakfast.cpts <- function(x, by = c('method', 'estimator'), ...) {
  
  L <- length(x$cptmodel.list)
  if(L == 0) stop('No change point analysis is performed')
    
  by <- match.arg(by, choices = c('method', 'estimator'))
  
  if(by == 'method'){
    max.char <- 0
    for(l in 1:L){
      cl <- x$cptmodel.list[[l]]
      nm <- paste(cl$solution.path, '.', cl$model.selection, sep = '')
      max.char <- max(max.char, nchar(nm))
    }
    
    cat(paste('Change-point locations estimated by:'))
    cat('\n\n')
    for(l in 1:L){
      cl <- x$cptmodel.list[[l]]
      nm <- paste(cl$solution.path, '.', cl$model.selection, sep = '')
      buff <- character(max.char - nchar(nm) + 1)
      if(cl$no.of.cpt > 0){
        cat(paste(nm, paste(buff, collapse = ' '), ': ', paste(cl$cpts, collapse = ', '), sep = ''))
        cat('\n')
      } else{
        cat(paste(nm, paste(buff, collapse = ' '), ': none'))
        cat('\n')
      }
    }
  }
  
  if(by == 'estimator'){
    brks <- rep(0, length(x$x))
    mm <- character(length(x$x))
    all.nm <- character(0)
    for(l in 1:L){
      cl <- x$cptmodel.list[[l]]
      nm <- paste(cl$solution.path, '.', cl$model.selection, sep = '')
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
}

#' @title Change-points estimated by solution path generation + model selection methods
#' @description  Print method for objects of class \code{cptmodel}
#' @method print cptmodel
#' @param x a \code{cptmodel} object
#' @param ... current not in use
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 5)
#' x <- f + rnorm(length(f)) * .5
#' print(model.ic(sol.idetect(x)))
#' @export
print.cptmodel <- function(x, ...) {

 cat(paste('Change-point locations estimated by:'))
 cat('\n\n')
 nm <- paste(x$solution.path, '.', x$model.selection, sep = '')
 if(x$no.of.cpt > 0){
   cat(paste(nm, ': ', paste(x$cpts, collapse = ', '), sep = ''))
   cat('\n')
 } else{
   cat(paste(nm, ': none'))
   cat('\n')
 }
}


#' @title Change-points estimated by the "breakfast" routine
#' @description Plot method for objects of class \code{breakfast.cpts}
#' @method plot breakfast.cpts
#' @param x a \code{breakfast.cpts} object
#' @param display.data if \code{display.data = TRUE}, change-point estimators are plotted against the data by method.
#' If \code{display.data = FALSE}, only the estimators are plotted; this option is recommended when \code{length(x)} is large.
#' @param ... current not in use
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 5)
#' x <- f + rnorm(length(f)) * .5
#' plot(breakfast(x, solution.path = 'all', model.selection = 'all'), display.data = TRUE)
#' plot(breakfast(x), display.data = FALSE)
#' @export
plot.breakfast.cpts <- function(x, display.data = TRUE, ...) {
  L <- length(x$cptmodel.list)
  hues <- seq(15, 375, length = L + 1)
  colors <- hcl(h = hues, l = 65, c = 100)[1:L]
  points <- (1:L) - 1
  
  all.nm <- character(L)
  for(l in 1:L){
    cl <- x$cptmodel.list[[l]]
    nm <- paste(cl$solution.path, '.', cl$model.selection, sep = '')
    all.nm[l] <- nm
  }
  all.nm <- factor(all.nm, levels = all.nm)
  
  if(display.data){
    ss <- sd(x$x)
    ii <- 1:length(x$x)
    
    df <- data.frame(i = ii, x = x$x)
    g <- ggplot(df, aes_string(x = "i", y = "x")) + geom_line(color = 'grey', size = .5) + 
      ggtitle(paste('Estimated change-point locations')) +
      ylim(c(min(x$x) - ss * .22 * (L + 2), max(x$x) * 1.1)) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      theme_classic() +
      xlab("time") + ylab("") 
    
    df <- data.frame(location = rep(0, L), method = all.nm, frequency = rep(0, L))
    g <- g + geom_point(df, mapping = aes_string(fill = "method", y = "frequency", x = "location"), 
                        inherit.aes = FALSE, alpha = 0) + 
      guides(fill = guide_legend(override.aes = list(alpha = 1, color = colors, shape = points))) 
    
    for(l in 1:L){
      cl <- x$cptmodel.list[[l]]
      #if(cl$no.of.cpt > 0){
        nm <- paste(cl$solution.path, '.', cl$model.selection, sep = '')
        h <- min(x$x) - .22 * ss * (L - l + 1)
        df1 <- data.frame(x = ii[cl$cpts], y = rep(h, cl$no.of.cpt))
        g <- g + 
          geom_point(df1, mapping = aes_string(x = "x", y = "y"), inherit.aes = FALSE, 
                     color = colors[l], shape = points[l], size = 2)
      #}
    }
  } else{
    df <- data.frame(location = rep(0, L), method = all.nm, frequency = rep(0, L))
    g <- ggplot(df, aes_string(fill = "method", y = "frequency", x = "location")) + 
      geom_point(alpha = 0) + 
      xlim(0, length(x$x)) + ylim(0, L + 0.5) + 
      ggtitle('Estimated change-point locations') + 
      theme_classic() +
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank()) +
      xlab("time") + ylab("") +
      guides(fill = guide_legend(override.aes = list(alpha = 1, color = colors, shape = points))) 
   
    for(l in 1:L){
      cl <- x$cptmodel.list[[l]]
      if(cl$no.of.cpt > 0){
        h <- L - l + 1
        df1 <- data.frame(x = (1:length(x$x))[cl$cpts], y = rep(h, cl$no.of.cpt))
        g <- g + 
          geom_point(df1, mapping = aes_string(x = "x", y = "y"), 
                     inherit.aes = FALSE, alpha = 1, color = colors[l], shape = points[l], size = 2)
      }
    }
      # ggplot(df, aes_string(fill = "method", y = "frequency", x = "location")) +
      #   geom_bar(position = "stack", stat = "identity", width = 1) +
      #   theme(axis.text.y=element_blank(), axis.title.y = element_blank()) +
      #   ggtitle('Estimated change-point locations') +
      #   theme_classic()
  }
  
  return(g)

}

# #' @title Change-points estimated by solution path generation + model selection methods
# #' @description Plot method for objects of class \code{cptmodel}
# #' @method plot cptmodel
# #' @param x a \code{cptmodel} object
# #' @param data a numeric vector containing the data processed by the combined method
# #' @param display.data if \code{display.data = TRUE}, change-point estimators are plotted against the data by method.
# #' If \code{display.data = FALSE}, only the estimators are plotted; this option is recommended when \code{length(data)} is large.
# #' @param ... current not in use
# #' @importFrom ggplot2
# #' @examples
# #' f <- rep(rep(c(0, 1), each = 50), 5)
# #' x <- f + rnorm(length(f)) * .5
# #' plot(model.ic(sol.idetect(x)), x, display.data = TRUE)
# #' plot(model.lp(sol.not(x)), x, display.data = FALSE)
# #' @export
# plot.cptmodel <- function(x, data, display.data = TRUE, ...) {
#   
#   ret <- structure(
#     list(x = data,
#          cptmodel.list = list(x)),
#     class = 'breakfast.cpts')
#   
#   plot(ret, display.data = display.data)
# }
