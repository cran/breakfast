all_shifts_are_cpts <- function(x) {
	
	diff.x <- abs(diff(x))
	
	cpts <- which(diff.x > 0)
	no.of.cpt <- length(cpts)
	est <- x
	
	
	list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
}

all_slopechanges_are_cpts <- function(x) {
  
  diff.x <- abs(diff(diff(x)))
  
  cpts <- which(diff.x > 0) + 1
  no.of.cpt <- length(cpts)
  est <- x
  
  
  list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
  
  
}




## for wbs2
# for roxygen2 to generate the help file - add ' below
# Generates a fixed number of random_intervals and identifes one
# with the maximum absolute CUSUM value.
# 
# @param x A numeric vector containing the input time series
# @param M An integer indicating the maximal number of intervals to be drawn
# @param max_cusum_fun A function that calculates the max of the cusum for a chosen setting
# @param seed An integer to be passed on to \code{\link[base]{set.seed}}
# @return A list which contains the below besides the input arguments
# \itemize{
# \item{res}{A matrix where each row contains the beginning and the end (second) of a random interval,
# the maximiser of the (absolute) CUSUM statistics and the maximum CUSUM value}
# \item{max.val}{The row of \code{res} with the largest maximum CUSUM value}
# }
# @examples
# # x <- mosum::testData('teeth10')$x
# # random_cusums(x, 100)
#' @keywords internal
#' @noRd


random_cusums <- function(x, M, max_cusum_fun = max_cusum, seed = NULL, min.d = 0) {
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
	res[,3:4] <- t(apply(ind, 2, max_cusum, y, min.d))
	max.ind <- which.max(abs(res[,4]))
	max.val <- res[max.ind,,drop=FALSE]
	list(res=res, max.val=max.val, M.eff=max(1, M))
}


systematic_cusums <- function(x, M, min.d = 0) {
  y <- c(0, cumsum(x))
  n <- length(x)
  M <- min(M, (n-1)*n/2)
  ind <- grid_intervals(n, M)
  M <- dim(ind)[2]
  res <- matrix(0, M, 4)
  res[,1:2] <- t(ind)
  res[,3:4] <- t(apply(ind, 2, max_cusum, y, min.d))
 # max.ind <- which.max(abs(res[,4]))
  max.ind <- max.col(matrix(abs(res[,4]), 1, length(res[,4])))
  max.val <- res[max.ind,,drop=FALSE]
  list(res=res, max.val=max.val, M.eff=M)
}

max_cusum <- function(ind, y, min.d = 0) {
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



wbs_K_int <- function(x, M, cusum.sampling) {
  n <- length(x)
  if (n == 1) return(matrix(NA, 4, 0))
  else {
    cpt <- t(cusum.sampling(x, M)$max.val)
    return(cbind(cpt, wbs_K_int(x[1:cpt[3]], M, cusum.sampling), wbs_K_int(x[(cpt[3]+1):n], M, cusum.sampling) + c(rep(cpt[3], 3), 0)            ))
  }
}


mean_from_cpt <- function(x, cpt) {
  n <- length(x)
  len.cpt <- length(cpt)
  if (len.cpt) cpt <- sort(cpt)
  beg <- endd <- rep(0, len.cpt+1)
  beg[1] <- 1
  endd[len.cpt+1] <- n
  if (len.cpt) {
    beg[2:(len.cpt+1)] <- cpt+1
    endd[1:len.cpt] <- cpt
  }
  means <- rep(0, len.cpt+1)
  for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])
  rep(means, endd-beg+1)
}

slope_from_cpt <- function (x, cpt) {
  n <- length(x)
  if (!is.null(cpt)) {
    if (any(is.na(cpt))) {
      cpt <- cpt[!is.na(cpt)]
    }
    }
    cpt <- as.integer(cpt)
    len_cpt <- length(cpt)
    if (len_cpt) {
      if (min(cpt) < 0 || max(cpt) >= n) 
        stop("change-points cannot be negative or greater than and n-1")
      cpt <- sort(cpt)
    }
    cpt <- sort(unique(c(cpt, 0, n)))
    fit <- rep(0, n)
    cpt <- setdiff(cpt, c(0, n))
    X <- splines::bs(1:n, knots = cpt, degree = 1, intercept = TRUE)
    fit <- stats::lm.fit(X, x)$fitted.values
  return(fit)
}


all_intervals_flat <- function(n) {
  if (n == 2) ind <- matrix(1:2, 2, 1) else {
    M <- (n-1)*n/2	
    ind <- matrix(0, 2, M)
    ind[1,] <- rep(1:(n-1), (n-1):1)
    ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
  }
  return(ind)
}

grid_intervals <- function(n, M) {
  if (n==2) ind <- matrix(c(1, 2), 2, 1)
  else if (M >= (n-1)*n/2) ind <- all_intervals_flat(n)
  else {
    k <- 1
    while (k*(k-1)/2 < M) k <- k+1
    ind2 <- all_intervals_flat(k)
    ind2.mx <- max(ind2)
    ind <- round((ind2 - 1) * ((n-1) / (ind2.mx-1)) + 1)
  }	
  return(ind)	
}


## for not, wbs - mostly follow the code from WBS2 in this version
random_intervals <-	function(n, M, seed = NULL) {
  
  if(is.integer(seed)) set.seed(seed)
  n <- as.integer(n)
  M <- as.integer(M)
  M <- min(M, (n-1)*n/2)
  
  #res <- matrix(0, max(1, M), 2)
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
  res <- t(ind)
  return(res)
}

# same two functions for interval drawing as in wbs and not
fixed_intervals <-function(n,M){
    n <- as.integer(n)
    M <- as.integer(M)
    M <- min(M, (n-1)*n/2)
    ind <- grid_intervals(n, M)
    M <- dim(ind)[2]
    res <- t(ind)
    return(res)
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

IDetect_linear_contr_one <- function (x, s, e, b) {
  r <- numeric()
  for (j in 1:length(b)) {
    x1 <- x[s[j]:e[j]]
    n <- length(x1)
    if ((b[j] - s[j] + 1) == 1) {
      r[j] <- 0
    }
    else {
      y1 <- cumsum(x1 * (1:n))
      y <- cumsum(x1)
      a <- sqrt(6/((n - 1) * n * (n + 1) * (2 - 2 * (b[j] -  s[j] + 1)^2 + 2 * (b[j] - s[j] + 1) * n - 1 + 2 * (b[j] - s[j] + 1) - n)))
      be <- sqrt(((n - (b[j] - s[j] + 1) + 1) * (n - (b[j] - s[j] + 1)))/((b[j] - s[j]) * (b[j] - s[j] + 1)))
      r[j] <- a * be * ((2 * (b[j] - s[j] + 1) + n - 1) * y1[b[j] - s[j] + 1] - (n + 1) * (b[j] - s[j] + 1) * y[b[j] - s[j] + 1]) - (a/be) * ((3 * n - 2 * (b[j] - s[j] + 1) + 1) * (y1[n] - y1[b[j] - s[j] + 1]) - (n + 1) * (2 * n - (b[j] - s[j] + 1)) * (y[n] - y[b[j] - s[j] + 1]))
    }
  }
  return(abs(r))
}

IDetect_cumsum_lin <- function(x) {
  if (!(is.numeric(x))) {
    stop("The input in `x' should be a numeric vector.")
  }
  res <- numeric()
  n <- length(x)
  if (n <= 2) {
    res <- 0
  }
  else {
    b <- 2:(n - 1)
    y1 <- cumsum(x * (1:n))
    y <- cumsum(x)
    a <- sqrt(6/((n - 1) * n * (n + 1) * (2 - 2 * b^2 + 2 * b * n - 1 + 2 * b - n)))
    be <- sqrt(((n - b + 1) * (n - b))/((b - 1) * b))
    res[1] <- 0
    res[b] <- a * be * ((2 * b + n - 1) * y1[b] - (n + 1) * b * y[b]) - (a/be) * ((3 * n - 2 * b + 1) * (y1[n] - y1[b]) - (n + 1) * (2 * n - b) * (y[n] - y[b]))
  }
  return(res)
}

phi <- function(x, s, e, b){
  res <- rep(0, length(x))
  if (s<b && b<e){
    alpha <- sqrt(6/((e - s)*(e - s + 1)*(e - s + 2)*(2-2*b^2+e+2*b*e-s+2*b*s-2*e*s)))
    beta <- sqrt(((e - b + 1) * (e - b))/((b - s)*(b - s + 1)))
    for(i in 1:length(x)){
      if((i >= s) && (i <=b)){
        res[i] = alpha * beta * (i *(e + 2*b - 3*s + 2) - (b*e + b*s + 2*s - 2*s^2))
      }else if((b < i) && (i <= e)){
        res[i] = - (alpha / beta) * (i *(3*e - 2*b - s + 2)- (2*e - b*e + 2*e^2 - b*s))
      }
    }
  }
  res
}


# contrast basis for piecewise constant
# x is the vector with observations (but we only use its length in this function)
# s is the start index, e is the end index, b is the location of the feature
psi <- function(x, s, e, b){
  l <- length(x)
  y <- 1:l
  res <- rep(0, l)
  res[(y>=s) & (y<=b)] <- rep(sqrt((e-b)/((e-s+1) * (b-s+1))), b-s+1)
  res[(y>=b+1) & (y<=e)] <- -rep(sqrt((b-s+1)/((e-s+1) * (e-b))), e-b)
  res
}

# contrast basis for linear
# x is the vector with observations (but we only use its length in this function)
# s is the start index, e is the end index, b is the location of the feature
gamma <- function(x, s, e){
  vec <- 1: (length(x))
  res <- rep(0,length(x))
  if(s<e){
    vec <- ((e-s+1)*(e*e-2*e*s+2*e + s*s -2*s)/12)^(-0.5)*(vec-(e+s)/2) 
    res[s:e] <- vec[s:e]
  }
  res
}


IDetect_linear_discontr_one <- function (x, s, e, b) {
  r <- rep(0,length(b))
  n <- length(x)
  for (j in 1:length(b)) {
    x1<-rep(0,n)
    x1[s[j]:e[j]] <- x[s[j]:e[j]]
    if ((e[j] - s[j] <= 3) || (b[j]<s[j]) || (b[j]+1 > e[j])) {
      r[j] <- 0
    }
    else{
      r[j] <- sqrt(abs((sum(x * psi(x,s[j],e[j],b[j])))^2+
                         (sum(x * gamma(x,s[j],b[j])))^2+
                         (sum(x * gamma(x,b[j]+1,e[j])))^2-
                         (sum(x * gamma(x,s[j],e[j])))^2))
    }
  }
  return(r)
}