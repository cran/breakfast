all.shifts.are.cpts <- function(x) {
	
	diff.x <- abs(diff(x))
	
	cpts <- which(diff.x > 0)
	no.of.cpt <- length(cpts)
	est <- x
	
	
	list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
	
	
}

all.slopechanges.are.cpts <- function(x) {
  
  diff.x <- abs(diff(diff(x)))
  
  cpts <- which(diff.x > 0)
  no.of.cpt <- length(cpts)
  est <- x
  
  
  list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
  
  
}




## for wbs2
# for roxygen2 to generate the help file - add ' below
# Generates a fixed number of random intervals and identifes one
# with the maximum absolute CUSUM value.
# 
# @param x A numeric vector containing the input time series
# @param M An integer indicating the maximal number of intervals to be drawn
# @param max.cusum.fun A function that calculates the max of the cusum for a chosen setting
# @param seed An integer to be passed on to \code{\link[base]{set.seed}}
# @return A list which contains the below besides the input arguments
# \itemize{
# \item{res}{A matrix where each row contains the beginning and the end (second) of a random interval,
# the maximiser of the (absolute) CUSUM statistics and the maximum CUSUM value}
# \item{max.val}{The row of \code{res} with the largest maximum CUSUM value}
# }
# @examples
# # x <- mosum::testData('teeth10')$x
# # random.cusums(x, 100)


random.cusums <- function(x, M, max.cusum.fun = max.cusum, seed = NULL) {
	
	  if(is.integer(seed)) set.seed(seed)

	# M=0 means there is a single cusum that gets taken over [1,length(x)] and therefore we are in the standard BS setting


	y <- c(0, cumsum(x))



	n <- length(x)

	

	M <- min(M, (n-1)*n/2)

		

	res <- matrix(0, max(1, M), 4)

	
	if ((n==2) || (M == 0)) ind <- matrix(c(1, n), 2, 1)
	
	else if (M == (n-1)*n/2) {

		ind <- matrix(0, 2, M)

		ind[1,] <- rep(1:(n-1), (n-1):1)

		ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)

	}

	else {

		ind <- ind2 <- matrix(floor(runif(2*M) * (n-1)), nrow=2)

		ind2[1,] <- apply(ind, 2, min)

		ind2[2,] <- apply(ind, 2, max)

		ind <- ind2 + c(1, 2)

	}



	res[,1:2] <- t(ind)

	res[,3:4] <- t(apply(ind, 2, max.cusum, y))



	max.ind <- which.max(abs(res[,4]))



	max.val <- res[max.ind,,drop=F]



	list(res=res, max.val=max.val, M.eff=max(1, M))

}




systematic.cusums <- function(x, M) {
  
  
  
  y <- c(0, cumsum(x))
  
  
  
  n <- length(x)
  
  
  
  M <- min(M, (n-1)*n/2)
  
  
  
  ind <- grid.intervals(n, M)
  
  
  
  M <- dim(ind)[2]
  
  
  
  res <- matrix(0, M, 4)
  
  
  
  res[,1:2] <- t(ind)
  
  res[,3:4] <- t(apply(ind, 2, max.cusum, y))
  
  
  
 # max.ind <- which.max(abs(res[,4]))
  
  max.ind <- max.col(matrix(abs(res[,4]), 1, length(res[,4])))
  
  max.val <- res[max.ind,,drop=F]
  
  
  
  list(res=res, max.val=max.val, M.eff=M)
  
  
  
}







max.cusum <- function(ind, y, min.d = 0) {
  
  z <- y[(ind[1]+1):(ind[2]+1)] - y[ind[1]]
  m <- ind[2] - ind[1] + 1
  if(m >= 2*(min.d) + 2){
    ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])
    mv <- max(abs(ip[(min.d + 1):(m - 1 - min.d)]))
    ip.max <- min(which(abs(ip) == mv))
  } else{
    mv <- 0
    ip.max <- 1
  }
  
  c(ip.max + ind[1] - 1, mv)

}



## for not, wbs
random.intervals <-	function(n, M, seed = NULL) {
  
  if(is.integer(seed)) set.seed(seed)
  
  n <- as.integer(n)
  M <- as.integer(M)
  M <- min(M, (n-1)*n/2)
  
  if ((n==2) || (M == 0)) intervals <- matrix(c(1, n), 2, 1)
  else if (M == (n-1)*n/2) {
    intervals <- matrix(0, nrow = M, ncol = 2)
    intervals[,1] <- rep(1:(n-1), (n-1):1)
    intervals[,2] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
  }
  else{
    intervals <- matrix(0, nrow = M, ncol = 2)
    intervals[,1] <- ceiling(runif(M)*(n - 1))
    intervals[,2] <- intervals[, 1] + ceiling(runif(M)*(n - intervals[, 1]))
    unique(intervals)
  }
  return(intervals)
}



# same two functions for interval drawing as in wbs and not
fixed.intervals <-function(n,M){
    n <- as.integer(n)
    M <- as.integer(M)
    
    m <- ceiling(0.5*(sqrt(8*M+1)+1))
    m <- min(n,m)
    M <- m*(m-1)/2
    end.points <- round(c(1,seq.int(2,n-1,length.out=(m-2)),n))
    intervals <- matrix(0,nrow=M,ncol=2)
    
    k <- 0;
    for(i in 1:(m-1)){
        tmp <- (m-i);
        intervals[(k+1):(k+tmp),1] <- rep(end.points[i],tmp)
        intervals[(k+1):(k+tmp),2] <- end.points[(i+1):m]
        
        k <- k+tmp;
    }
    
    unique(intervals)
}

# check input numermic

check.input <- function(x){
  # if(length(x) < 1) stop("Data vector x should contain at least two elements.")
  if(!is.numeric(x)) stop("Data vector x vector must be numeric")
  if(!all(is.finite(x))) stop("Data vector x vector cannot contain NA's")
}

## for Isolate-Detect
start_end_points <- function(r, l, s, e) {
  r <- sort(r)
  l <- sort(l, decreasing = TRUE)
  if (s > e){
    stop("s should be less than or equal to e")
  }
  if (!(is.numeric(c(r, l, s, e))) | (r[1] <= 0) | (l[length(l)] <= 0) | s <= 0 | e <= 0){
    stop("The input arguments must be positive integers")
  }
  if (any(abs(r - round(r)) > .Machine$double.eps ^ 0.5)){
    warning("The input for r should be a vector of positive integers. If there is at least a positive real
            number then the integer part of that number is used.")
  }
  if (any(abs(l - round(l)) > .Machine$double.eps ^ 0.5)){
    warning("The input for l should be a vector of positive integers. If there is at least a positive real
            number then the integer part of that number is used.")
  }
  if (abs(s - round(s)) > .Machine$double.eps ^ 0.5){
    warning("The input for s should be a positive integer. If it is a positive real
            number then the integer part of that number is used.")
  }
  if (abs(e - round(e)) > .Machine$double.eps ^ 0.5){
    warning("The input for e should be a positive integer. If it is a positive real
            number then the integer part of that number is used.")
  }
  r <- as.integer(r)
  l <- as.integer(l)
  e <- as.integer(e)
  s <- as.integer(s)
  e_points <- unique(c(r[which( (r > s) & (r < e))], e))
  s_points <- unique(c(l[which( (l > s) & (l < e))], s))
  return(list(e_points = e_points, s_points = s_points))
}

# The CUSUM function for the case of a piecewise-constant signal.
cusum_function <- function(x) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data
         for which the CUSUM function will be calculated.")
  }
  n <- length(x)
  y <- cumsum(x)
  res <- sqrt( ( (n - 1):1) / n / (1:(n - 1))) * y[1:(n - 1)] - sqrt( (1:(n - 1)) / n / ( (n - 1):1)) * (y[n] - y[1:(n - 1)])
  return(res)
}

# This function estimates the signal in a given data sequence x with change-points at cpt.
# The type of the signal depends on whether the change-points represent changes in a
# piecewise-constant or continuous, piecewise-linear signal.
# est_signal <- function(x, cpt, type = c("mean", "slope")) {
#   if (!(is.numeric(x))){
#     stop("The input in `x' should be a numeric vector containing the data for
#          which you want to estimate the underlying signal.")
#   }
#   x <- as.numeric(x)
#   if (NA %in% x)
#     stop("x vector cannot contain NA's")
#   n <- length(x)
#   if (type == "mean") {
#     if (missing(cpt)){
#       cpt <- ID_pcm(x)$cpt
#     }
#     if (!is.null(cpt)){
#       if (any(is.na(cpt))){
#         cpt <- cpt[!is.na(cpt)]
#       }
#     }
#     cpt <- as.integer(cpt)
#     len_cpt <- length(cpt)
#     if (len_cpt) {
#       if (min(cpt) < 0 || max(cpt) >= n)
#         stop("change-points cannot be negative or greater than and n-1")
#       cpt <- sort(cpt)
#     }
#     s <- e <- rep(0, len_cpt + 1)
#     s[1] <- 1
#     e[len_cpt + 1] <- n
#     if (len_cpt) {
#      s[2:(len_cpt + 1)] <- cpt + 1
#      e[1:len_cpt] <- cpt
#    }
#    means <- rep(0, len_cpt + 1)
#    for (i in 1:(len_cpt + 1)) {
#      means[i] <- mean(x[s[i]:e[i]])
#    }
#    fit <- rep(means, e - s + 1)
#  }
#  if (type == "slope"){
#    if (missing(cpt)){
#      cpt <- ID_cplm(x)$cpt
#    }
#    if (!is.null(cpt)){
#      if (any(is.na(cpt))){
#        cpt <- cpt[!is.na(cpt)]
#      }
#    }
#    cpt <- as.integer(cpt)
#    len_cpt <- length(cpt)
#    if (len_cpt) {
#      if (min(cpt) < 0 || max(cpt) >= n)
#        stop("change-points cannot be negative or greater than and n-1")
#      cpt <- sort(cpt)
#    }
#    cpt <- sort(unique(c(cpt, 0, n)))
#    fit <- rep(0, n)
#    cpt <- setdiff(cpt, c(0, n))
#    X <- splines::bs(1:n, knots = cpt, degree = 1, intercept = TRUE)
#    fit <- stats::lm.fit(X, x)$fitted.values
#  }
#  return(fit)
#}

# The following function finds the CUSUM value for the intervals [s[j], e[j]) at the location
# b[j], where j = 1,2,...,L, with L being the length of the vectors s,b,e.
# The input arguments are:
#  1) x which is a numeric vector containing the data.
#  2) s which is a vector of length L which has the startpoints for the CUSUM values that will be computed.
#  3) e which is a vector of length L which has the endpoints for the CUSUM values that will be computed.
#  4) b which is a vector of length L which has the location of the points where the CUSUM will be computed.
#
# The output is a numeric vector of length L which has the CUSUM values, CS[j], j=1,2,...,L calculated at the
# location b[j], when the interval [s[j], e[j]] is used.


IDetect_cusum_one <- function(x, s, e, b) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector.")
  }
  y <- cumsum(x)
  l <- numeric()
  d <- numeric()
  result <- numeric()
  if ( (length(s) != length(b)) || (length(s) != length(e)) || (length(e) != length(b))){
    stop("The vectors s, b, e, should be of the same length")
  }
  if (any(s < 1) | any(b < 1) | any(e < 1)){
    stop("The entries of the vectors s, b, e should be positive integers.")
  }
  if (any(s > b) | any(b >= e)){
    stop("The value for b should be in the interval [s,e)")
  }
  if ( (any(abs( (s - round(s))) > .Machine$double.eps ^ 0.5))
       || (any(abs( (b - round(b))) > .Machine$double.eps ^ 0.5))
       || (any(abs( (e - round(e))) > .Machine$double.eps ^ 0.5))){
    stop("The input values  for s, b, and  e should be positive integers.")
  }
  for (j in 1:length(b)) {
    l[j] <- e[j] - s[j] + 1
    d[j] <- ifelse(s[j] == 1, 0, y[s[j] - 1])
    result[j] <- abs(sqrt( (e[j] - b[j]) / (l[j] * (b[j] - s[j] + 1))) * (y[b[j]] - d[j]) - sqrt( (b[j] - s[j] + 1) / (l[j] * (e[j] - b[j]))) * (y[e[j]] - y[b[j]]))
  }
  return(result)
}

