#' @title Estimating change-points in the piecewise-constant or piecewise-linear mean of a noisy data sequence via the Steepest Drop to Low Levels method
#' @description This function estimates the number and locations of change-points in the piecewise-constant or piecewise-linear mean of a noisy data sequence via the Steepest Drop to Low Levels method.
#' @details 
#' The Steepest Drop to Low Levels method is described in 
#' "Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection", P. Fryzlewicz (2020), Journal of the Korean Statistical Society, 49, 1027--1070.
#' 
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. The \code{cptpath.object$type} variable decides the model type: piecewise-constant (\code{type == "const"}),
#' piecewise-linear and continuous (\code{type == "lin.cont"}) or piecewise-linear and discontinuous (\code{type == "lin.discont"}). In the piecewise-constant model, SDLL model selection should work well 
#' when \code{cptpath.object} is an object returned by the \code{sol.wbs2} routine. In the piecewise-linear model (whether continuous or not), the output of \code{sol.idetect} should be supplied as
#' \code{cptpath.object}. Note that the field \code{cptpath.object$x} contains the input data sequence. 
#' @param sigma An estimate of the standard deviation of the noise in the data \code{cptpath.object$x}. Can be a functional of \code{cptpath.object$x} or a specific value if known. 
#' The default in the piecewise-constant model is the Median Absolute Deviation of the vector \code{diff(cptpath.object$x)/sqrt(2)}, tuned to the Gaussian distribution.
#' In the piecewise-linear models, \code{diff(cptpath.object$x, differences = 2)/sqrt(6)} is used by default.
#' Note that \code{model.sdll} works particularly well when the noise is i.i.d. Gaussian.
#' @param universal If \code{TRUE}, then the threshold that decides if there are any change-points is chosen automatically, so that the probability of type-I error (i.e. indicating change-points if there are none) 
#' is approximately \code{1 - alpha}. If \code{FALSE}, then \code{th.const} must be specified.
#' @param th.const Only relevant if \code{universal == FALSE}; in that case a numerical value must be provided. Used to create the threshold (applicable to the contrast magnitudes stored in \code{cptpath.object}) 
#' that decides if there are any change-points in the mean vector; that threshold is then \code{th.const * sqrt(2 * log(n)) * sigma}, where \code{n} is the length of the data vector \code{cptpath.object$x}.
#' @param th.const.min.mult A fractional multiple of the threshold, used to decide the lowest magnitude of contrasts from \code{cptpath.object} still considered by the SDLL model selection criterion as potentially change-point-carrying.
#' @param lambda Only relevant if \code{universal == TRUE}; can be set to 0.9 or 0.95. The approximate probability of not detecting any change-points if the truth does not contain any.
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{type}{The model type used, inherited from the given \code{cptpath.object}}
#' \item{model.selection}{The model selection method used to return the final change-point estimators object, here its value is \code{"sdll"}}
#' \item{no.of.cpt}{The number of estimated change-points}
#' \item{cpts}{The locations of estimated change-points}
#' \item{est}{An estimate of the mean of the vector \code{cptpath.object$x}}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{breakfast}}
#' @references P. Fryzlewicz (2020). Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{Journal of the Korean Statistical Society}, 49, 1027--1070.
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f))
#' model.sdll(sol.wbs2(x))
model.sdll <- function(cptpath.object, sigma = stats::mad(diff(cptpath.object$x)/sqrt(2)), universal = TRUE, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9) {

  
  if(!("cptpath" %in%  class(cptpath.object)))  stop("A cptmodel class object has to be supplied in the first argument.")
  
	x <- cptpath.object$x
	type <- cptpath.object$type

	n <- length(x)

    if (n <= 1) {

        est <- x

        no.of.cpt <- 0

        cpts <- integer(0)

    }

    else {

		if (!(type == "const")) sigma <- stats::mad(diff(cptpath.object$x, differences = 2)/sqrt(6))

		if (sigma == 0) {
			
			if (type == "const") s0 <- all_shifts_are_cpts(x) else s0 <- all_slopechanges_are_cpts(x)
			est <- s0$est
			no.of.cpt <- s0$no.of.cpt
			cpts <- s0$cpts
			
		} else {

		if (universal) {

        	if (type == "const") u <- universal.M.th.v3(n, lambda)
        	else if (type == "lin.cont") u <- universal.th.lin(n, lambda)
        	else if (type == "lin.discont") u <- universal.th.lin.discont(n, lambda)

        	th.const <- u$th.const

    	}

    	else if (is.null(th.const)) stop("If universal is FALSE, then th.const must be specified.")
    	    	
		if ((type == "const") & (cptpath.object$method %in% c("not", "idetect_seq"))) warning("model.sdll won't work well on cptpath.object produced with sol.not or sol.idetect_seq; consider using other model. functions, or produce your cptpath.object with a different sol. function")

		if (!(type == "const") & !(cptpath.object$method == "idetect")) warning("In the piecewise-linear change-point model, for the SDLL
		model selection criterion to work well, please consider supplying a cptpath.object produced with sol.idetect.")
 
    	th.const.min <- th.const * th.const.min.mult

    	th <- th.const * sqrt(2 * log(n)) * sigma

    	th.min <- th.const.min * sqrt(2 * log(n)) * sigma



 	#	rc <- t(wbs_K_int(x, M))

 		if (dim(cptpath.object$cands)[1] == 0) {
 			
    	    no.of.cpt <- 0

        	cpts <- integer(0)
 			 			
 			
 		} else
 		 		
 		if (cptpath.object$cands[1,4] < th) {

 #			est <- mean(x)

    	    no.of.cpt <- 0

        	cpts <- integer(0)



 		}

		else {

			indices <- which(cptpath.object$cands[,4] > th.min)

			if (length(indices) == 1) {

				cpts <- cptpath.object$cands[indices, 3]

				no.of.cpt <- 1

	#			est <- mean_from_cpt(x, cpts)								

			}

			else {

				rc.sel <- cptpath.object$cands[indices,,drop=F]

	#			ord <- order(abs(rc.sel[,4]), decreasing=T)

				z <- cptpath.object$cands[indices,4]




				z.l <- length(z)

				dif <- -diff(log(z))

				dif.ord <- order(dif, decreasing=T)

				j <- 1

				while ((j < z.l) & (z[dif.ord[j]+1] > th)) j <- j+1

				if (j < z.l) no.of.cpt <- dif.ord[j] else no.of.cpt <- z.l

				cpts <- sort(cptpath.object$cands[1:no.of.cpt,3])			

	#			est <- mean_from_cpt(x, cpts)						

			}

		} 

	est <- prediction(cptpath.object, cpts)

#	if (type == "const") est <- mean_from_cpt(x, cpts) else est <- slope_from_cpt(x, cpts)

    }

	}


	ret <- list(solution.path=cptpath.object$method, type = type, model.selection="sdll", no.of.cpt=no.of.cpt, cpts=cpts, est=est)
	
	class(ret) <- "cptmodel"
	
	ret

}

universal.M.th.v3 <- function(n, lambda = 0.9) {

		

	mat.90 <- matrix(0, 24, 3)

	mat.90[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.90[,2] <- c(2.091965, 1.482996, 1.410729, 1.35482, 1.335922, 1.304641, 1.283584, 1.265133, 1.260038, 1.244437, 1.254992, 1.248876, 1.228693, 1.225614, 1.214514, 1.209304, 1.197611, 1.190821, 1.187124, 1.187707, 1.182522, 1.180259, 1.175435, 1.16677)

#	mat.90[,2] <- c(1.420, 1.310, 1.280, 1.270, 1.250, 1.220, 1.205, 1.205, 1.200, 1.200, 1.200, 1.185, 1.185, 1.170, 1.170, 1.160, 1.150, 1.150, 1.150, 1.150, 1.145, 1.145, 1.135, 1.135)

	mat.90[,3] <- rep(1000, 24)

	

	mat.95 <- matrix(0, 24, 3)

	mat.95[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.95[,2] <- c(2.641451, 1.630137, 1.542994, 1.4237, 1.433023, 1.369914, 1.352619, 1.315107, 1.303483, 1.299074, 1.304927, 1.30621, 1.275396, 1.269159, 1.254369, 1.251187, 1.236001, 1.243129, 1.220056, 1.226702, 1.22321, 1.222939, 1.211347, 1.207717)


#	mat.95[,2] <- c(1.550, 1.370, 1.340, 1.320, 1.300, 1.290, 1.265, 1.265, 1.247, 1.247, 1.247, 1.225, 1.225, 1.220, 1.210, 1.190, 1.190, 1.190, 1.190, 1.190, 1.190, 1.180, 1.170, 1.170)

	mat.95[,3] <- rep(1000, 24)



	if (lambda == 0.9) A <- mat.90 else A <- mat.95



	d <- dim(A)

	if (n < A[1,1]) {

		th <- A[1,2]

		M <- A[1,3]

	}

	else if (n > A[d[1],1]) {

		th <- A[d[1],2]

		M <- A[d[1],3]

	}

	else {

		ind <- order(abs(n - A[,1]))[1:2]

		s <- min(ind)

		e <- max(ind)

		th <- A[s,2] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,2] * (n - A[s,1])/(A[e,1] - A[s,1])

		M <- A[s,3] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,3] * (n - A[s,1])/(A[e,1] - A[s,1])

	}



	list(th.const=th, M=M)

}



universal.th.lin <- function(n, lambda = 0.9) {
	
	mat.90 <- matrix(0, 13, 2)
	mat.90[,1] <- c(30, 50, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000)
	mat.90[,2] <- c(1.32, 1.24, 1.19, 1.15, 1.15, 1.14, 1.13, 1.13, 1.13, 1.12, 1.12, 1.12, 1.12)

	mat.95 <- matrix(0, 13, 2)
	mat.95[,1] <- c(30, 50, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000)
	mat.95[,2] <- c(1.5, 1.36, 1.27, 1.21, 1.19, 1.18, 1.18, 1.18, 1.17, 1.17, 1.17, 1.17, 1.16)
	
	if (lambda == 0.9) A <- mat.90 else A <- mat.95

	d <- dim(A)

	if (n < A[1,1]) th <- A[1,2]
		else if (n > A[d[1],1]) th <- A[d[1],2]
			else th <- approx(A[,1], A[,2], xout = n)$y

	list(th.const=th)
	
}



universal.th.lin.discont <- function(n, lambda = 0.9) {
  
  mat.90 <- matrix(0, 13, 2)
  mat.90[,1] <- c(30, 50, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000)
  mat.90[,2] <- c(1.70, 1.61, 1.49, 1.42, 1.37, 1.34, 1.34, 1.31, 1.28, 1.26, 1.25, 1.25, 1.23)
  
  mat.95 <- matrix(0, 13, 2)
  mat.95[,1] <- c(30, 50, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000)
  mat.95[,2] <- c(1.90, 1.73, 1.61, 1.49, 1.42, 1.40, 1.40, 1.35, 1.33, 1.30, 1.30, 1.29, 1.27)
  
  if (lambda == 0.9) A <- mat.90 else A <- mat.95
  
  d <- dim(A)
  
  if (n < A[1,1]) th <- A[1,2]
  else if (n > A[d[1],1]) th <- A[d[1],2]
  else th <- approx(A[,1], A[,2], xout = n)$y
  
  list(th.const=th)
  
}
