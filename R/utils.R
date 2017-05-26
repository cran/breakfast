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

mean.from.cpt <- function(x, cpt) {

	iv <- cpts.into.intervals(cpt, length(x))
	
	len.cpt <- length(cpt)
	
	means <- rep(0, len.cpt+1)
	for (i in 1:(len.cpt+1)) means[i] <- mean(x[iv[1,i]:iv[2,i]])
	rep(means, iv[2,]-iv[1,]+1)
}
