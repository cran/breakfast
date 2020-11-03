#' @title Estimating change-points in the piecewise-constant mean of a noisy data sequence via the Steepest Drop to Low Levels method
#' @description This function estimates the number and locations of change-points in the piecewise-constant mean of a noisy data sequence via the Steepest Drop to Low Levels method.
#' @details 
#' The Steepest Drop to Low Levels method is described in 
#' "Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection", P. Fryzlewicz (2020), Journal of the Korean Statistical Society, to appear.
#' 
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. In particular, SDLL model selection should work well when \code{cptpath.object} is an object returned by the \code{sol.wbs2} routine. Note that the field \code{cptpath.object$x} contains the input data sequence. 
#' @param sigma An estimate of the standard deviation of the noise in the data \code{cptpath.object$x}. Can be a functional of \code{cptpath.object$x} or a specific value if known. The default is the Median Absolute Deviation of the vector \code{diff(cptpath.object$x)/sqrt(2)}, tuned to the Gaussian distribution. Note that \code{model.sdll} works particularly well when the noise is i.i.d. Gaussian.
#' @param universal If \code{TRUE}, then the threshold that decides if there are any change-points is chosen automatically, so that the probability of type-I error (i.e. indicating change-points if there are none) is approximately \code{1 - alpha} when the number \code{M} of intervals drawn in the \code{sol.wbs2} solution path routine is 100. If \code{FALSE}, then \code{th.const} must be specified.
#' @param th.const Only relevant if \code{universal == FALSE}; in that case a numerical value must be provided. Used to create the threshold (applicable to the CUSUM magnitudes stored in \code{cptpath.object}) that decides if there are any change-points in the mean vector; that threshold is then \code{th.const * sqrt(2 * log(n)) * sigma}, where \code{n} is the length of the data vector \code{cptpath.object$x}.
#' @param th.const.min.mult A fractional multiple of the threshold, used to decide the lowest magnitude of CUSUMs from \code{cptpath.object} still considered by the SDLL model selection criterion as potentially change-point-carrying.
#' @param lambda Only relevant if \code{universal == TRUE}; can be set to 0.9 or 0.95. The approximate probability of not detecting any change-points if the truth does not contain any, when the number \code{M} of intervals drawn in the \code{sol.wbs2} solution path routine is 100.
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{model}{The model selection method used to return the final change-point estimators object, here its value is \code{"sdll"}}
#' \item{no.of.cpt}{The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#' \item{cpts}{The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#' \item{est}{An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{breakfast}}
#' @references P. Fryzlewicz (2020). Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection. \emph{Journal of the Korean Statistical Society (to appear)}.
#' @examples
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f))
#' model.sdll(sol.wbs2(x))
model.sdll <- function(cptpath.object, sigma = stats::mad(diff(cptpath.object$x)/sqrt(2)), universal = TRUE, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9) {

	x <- cptpath.object$x

	n <- length(x)

    if (n <= 1) {

        est <- x

        no.of.cpt <- 0

        cpts <- integer(0)

    }

    else {

		if (sigma == 0) {
			
			s0 <- all.shifts.are.cpts(x)
			est <- s0$est
			no.of.cpt <- s0$no.of.cpt
			cpts <- s0$cpts
			
		} else {

		if (universal) {

        	u <- universal.M.th.v3(n, lambda)

        	th.const <- u$th.const

#        	M <- u$M

    	}

    	else if (is.null(th.const)) stop("If universal is FALSE, then th.const must be specified.")
    	    	
		if (cptpath.object$method %in% c("not", "idetect_seq")) warning("model.sdll won't work well on cptpath.object produced with sol.not or sol.idetect_seq; consider using other model. functions, or produce your cptpath.object with a different sol. function")
 
    	th.const.min <- th.const * th.const.min.mult

    	th <- th.const * sqrt(2 * log(n)) * sigma

    	th.min <- th.const.min * sqrt(2 * log(n)) * sigma



 	#	rc <- t(wbs.K.int(x, M))

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

	#			est <- mean.from.cpt(x, cpts)								

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

	#			est <- mean.from.cpt(x, cpts)						

			}

		} 

	est <- mean.from.cpt(x, cpts)

    }

	}


	ret <- list(solution.path=cptpath.object$method, model="sdll", no.of.cpt=no.of.cpt, cpts=cpts, est=est)
	
	class(ret) <- "cptmodel"
	
	ret

}

universal.M.th.v3 <- function(n, lambda = 0.9) {

		

	mat.90 <- matrix(0, 24, 3)

	mat.90[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.90[,2] <- c(1.420, 1.310, 1.280, 1.270, 1.250, 1.220, 1.205, 1.205, 1.200, 1.200, 1.200, 1.185, 1.185, 1.170, 1.170, 1.160, 1.150, 1.150, 1.150, 1.150, 1.145, 1.145, 1.135, 1.135)

	mat.90[,3] <- rep(100, 24)

	

	mat.95 <- matrix(0, 24, 3)

	mat.95[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.95[,2] <- c(1.550, 1.370, 1.340, 1.320, 1.300, 1.290, 1.265, 1.265, 1.247, 1.247, 1.247, 1.225, 1.225, 1.220, 1.210, 1.190, 1.190, 1.190, 1.190, 1.190, 1.190, 1.180, 1.170, 1.170)

	mat.95[,3] <- rep(100, 24)



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



