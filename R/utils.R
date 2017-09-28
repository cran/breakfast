#' @importFrom stats runif


cpts.into.intervals <- function(cpt, n) {
	
	len.cpt <- length(cpt)
	if (len.cpt) cpt <- sort(cpt)
	iv <- matrix(0, 2, len.cpt+1)
	iv[1,1] <- 1
	iv[2,len.cpt+1] <- n
	if (len.cpt) {
		iv[1, 2:(len.cpt+1)] <- cpt+1
		iv[2, 1:len.cpt] <- cpt
	}
	iv
	
}


mean.from.cpts <- function(x, cpt) {

	iv <- cpts.into.intervals(cpt, length(x))
	
	len.cpt <- length(cpt)
	
	means <- rep(0, len.cpt+1)
	for (i in 1:(len.cpt+1)) means[i] <- mean(x[iv[1,i]:iv[2,i]])
	rep(means, iv[2,]-iv[1,]+1)
}


max.cusum <- function(ind, y) {
	
	z <- y[(ind[1]+1):(ind[2]+1)] - y[ind[1]]
	m <- ind[2]-ind[1]+1
	ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])
	ip.max <- which.max(abs(ip))
		
	c(ip.max + ind[1] - 1, ip[ip.max])

}


wbs.K.int <- function(x, M) {
	
	n <- length(x)
	if (n == 1) return(matrix(NA, 4, 0))
	else {
		cpt <- t(random.cusums(x, M)$max.val)
		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M), wbs.K.int(x[(cpt[3]+1):n], M) + c(rep(cpt[3], 3), 0)            ))
	}
	
}


random.cusums <- function(x, M) {

	y <- c(0, cumsum(x))

	n <- length(x)
	
	M <- min(M, (n-1)*n/2)
		
	res <- matrix(0, M, 4)
	
	if (n==2) ind <- matrix(1:2, 2, 1)
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

	list(res=res, max.val=max.val, M.eff=M)

}


wbs.sort <- function(rc) {
		
	sorted <- rc$res[order(abs(rc$res[,4]), decreasing=T), , drop=F]

	res <- matrix(0, rc$M.eff, 4)

	j <- 1

	while (dim(sorted)[1]) {
		res[j,] <- sorted[1,]
		chp.cand <- sorted[1,3]
		sorted <- matrix(sorted[!(sorted[,1] <= chp.cand & sorted[,2] > chp.cand),], ncol=4)
		j <- j + 1
	}

	cpt <- res[1:(j-1),, drop=F]

	cpt
	
}


wbs.thresh.int <- function(x, M, th, th.min, adapt = TRUE) {

	n <- length(x)
	if (n == 1) return(x)
	else {
		rc <- random.cusums(x, M)
		if (abs(rc$max.val[4]) > th) {
			cpt <- wbs.sort(rc)

			cpt <- cpt[(abs(cpt[,4]) > th), 3]

			est <- rep(0, n)

			len.cpt <- length(cpt)
			
			iv <- cpts.into.intervals(cpt, n)
			
			
			if (adapt) for (i in 1:(len.cpt+1)) {
				th.current <- max(th.min, th*sqrt(log(max(n/2, iv[2,i]-iv[1,i]+1))/log(n)))
				est[iv[1,i]:iv[2,i]] <- 
				wbs.thresh.int(x[iv[1,i]:iv[2,i]], M, th.current, th.min, adapt)

			}
			else for (i in 1:(len.cpt+1))
				est[iv[1,i]:iv[2,i]] <- mean(x[iv[1,i]:iv[2,i]]) 

			return(est)
		}
		else return(rep(mean(x), n))
	}

}


universal.M.th <- function(n, lambda = 0.9) {
		
	mat.90 <- matrix(0, 22, 3)
	mat.90[,1] <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
	mat.90[,2] <- c(rep(1.3, 6), 1.25, 1.25, 1.25, 1.25, 1.20, 1.20, 1.20, 1.20, 1.20, 1.15, 1.15, 1.15, 1.15, 1.15, 1.10, 1.10)
	mat.90[,3] <- c(20, 50, 500, rep(1000, 19))

	mat.95 <- matrix(0, 22, 3)
	mat.95[,1] <- mat.90[,1]
	mat.95[,2] <- c(rep(1.3, 10), 1.25, 1.25, 1.25, 1.25, 1.25, 1.20, 1.20, 1.20, 1.20, 1.15, 1.15, 1.15)
	mat.95[,3] <- c(5, 10, 150, 300, 500, 500, rep(1000, 16))

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
