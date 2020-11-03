#' @title Solution path generation via the Tail-Greedy Unbalanced Haar method
#' @description This function arranges all possible change-points in the mean of the input vector in the order of importance, via the Tail-Greedy Unbalanced Haar method.
#' @details 
#' The Tail-Greedy Unbalanced Haar decomposition algorithm is described in 
#' "Tail-greedy bottom-up data decompositions and fast multiple change-point 
#' detection", P. Fryzlewicz (2018), The Annals of Statistics, 46: 3390--3421.
#'
#' @param x A numeric vector containing the data to be processed
#' @param p Specifies the number of region pairs merged 
#' in each pass through the data, as the proportion of all remaining region pairs. The default is
#' \code{p = 0.01}
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{M}{Input parameter \code{M}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "tguh" here}
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}
#' @references P. Fryzlewicz (2018). Tail-greedy bottom-up data decompositions and fast multiple change-point detection. \emph{The Annals of Statistics}, 46(6B), 3390--3421.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.tguh(r3)
#' @export
sol.tguh <- function(x, p = .01) {

	solutions.nested <- TRUE
	
	solution.set <- list()
	
	n <- length(x)
	
	sorted.cusums <- matrix(NA, 0, 4)
	
	if (n <= 1) solution.path <- integer()
	
	else {

		x.tguh <- tguh.decomp(x, p)
	
		sorted.cusums <- matrix(0, n-1, 4)
	
		sorted.cusums[,1] <- x.tguh$decomp.hist[1,1,]
		sorted.cusums[,2] <- x.tguh$decomp.hist[1,2,] + x.tguh$decomp.hist[4,2,] - 1
		sorted.cusums[,3] <- x.tguh$decomp.hist[1,2,] - 1
	
		sorted.cusums[,4] <- abs(x.tguh$decomp.hist[3,1,])
	
		ord <- order(sorted.cusums[,4], decreasing = T)
	
		sorted.cusums <- sorted.cusums[ord,, drop=F]
		
		solution.path <- sorted.cusums[,3]
		
	}	
	
	ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, p = p, cands = sorted.cusums, method = "tguh")
	
	class(ret) <- "cptpath"
	
	ret
	
}




tguh.decomp <- function(x, p = .01) {

# Note tguh.decomp requires plyr::mapvalues

	n <- length(x)

	noe <- n-1

	weights <- x - x + 1

	edges <- matrix(0, noe, 2)

	edges[,1] <- 1:(n-1)
	edges[,2] <- 2:n

	decomp.hist <- array(0, dim=c(4, 2, n-1))

	tguh.coeffs <- as.numeric(x)
	vec.weights <- as.numeric(weights)

	steps.left <- n-1

	current.step <- 0

	while (dim(edges)[1]) {

		max.current.steps <- ceiling(p * steps.left)

		removable.nodes <- rep(1, max(edges))

		a <- vec.weights[edges[,1]]
		b <- vec.weights[edges[,2]]

		h1 <- 1/sqrt(1 + (a/b)^2)
		h2 <- -1/sqrt(1 + (b/a)^2)
		l1 <- -h2
		l2 <- h1

		details <- h1 * tguh.coeffs[edges[,1]] + h2 * tguh.coeffs[edges[,2]]
		
		ord.det <- order(abs(details))

		edge.indices.2b.removed <- 1
		traverse.edges.index <- 1
		removable.nodes[edges[ord.det[1],1]] <- removable.nodes[edges[ord.det[1],2]] <- 0

		while  ( (length(edge.indices.2b.removed) < max.current.steps) & (traverse.edges.index < noe) ) {
			traverse.edges.index <- traverse.edges.index + 1
			if (removable.nodes[edges[ord.det[traverse.edges.index],1]] & removable.nodes[edges[ord.det[traverse.edges.index],2]]) {
				edge.indices.2b.removed <- c(edge.indices.2b.removed, traverse.edges.index)
				removable.nodes[edges[ord.det[traverse.edges.index],1]] <- removable.nodes[edges[ord.det[traverse.edges.index],2]] <- 0
			}
		}

		details.min.ind <- ord.det[edge.indices.2b.removed]

		no.of.current.steps <- length(edge.indices.2b.removed)

		smooth.at.min <- l1[details.min.ind] * tguh.coeffs[edges[details.min.ind,1]] + 
				l2[details.min.ind] * tguh.coeffs[edges[details.min.ind,2]]

		det.weight.at.min <- h1[details.min.ind] * vec.weights[edges[details.min.ind,1]] + 
				h2[details.min.ind] * vec.weights[edges[details.min.ind,2]]
		sm.weight.at.min <- l1[details.min.ind] * vec.weights[edges[details.min.ind,1]] + 
				l2[details.min.ind] * vec.weights[edges[details.min.ind,2]]

		decomp.hist[1,,(current.step+1):(current.step+no.of.current.steps)] <- t(edges[details.min.ind,])
		decomp.hist[2,1,(current.step+1):(current.step+no.of.current.steps)] <- h1[details.min.ind]
		decomp.hist[2,2,(current.step+1):(current.step+no.of.current.steps)] <- h2[details.min.ind]
		decomp.hist[3,1,(current.step+1):(current.step+no.of.current.steps)] <- details[details.min.ind]
		decomp.hist[3,2,(current.step+1):(current.step+no.of.current.steps)] <- smooth.at.min
		decomp.hist[4,1,(current.step+1):(current.step+no.of.current.steps)] <- vec.weights[edges[details.min.ind,1]]^2
		decomp.hist[4,2,(current.step+1):(current.step+no.of.current.steps)] <- vec.weights[edges[details.min.ind,2]]^2


		eating.up <- apply(matrix(edges[details.min.ind,], no.of.current.steps, 2), 1, min)
		eaten.up <- apply(matrix(edges[details.min.ind,], no.of.current.steps, 2), 1, max)

		tguh.coeffs[eating.up] <- smooth.at.min
		tguh.coeffs[eaten.up] <- details[details.min.ind]

		vec.weights[eating.up] <- sm.weight.at.min
		vec.weights[eaten.up] <- det.weight.at.min

		edges <- plyr::mapvalues(edges, eaten.up, eating.up)

		edges <- edges[edges[,1] != edges[,2],]

		if (length(edges) == 2) edges <- matrix(edges, 1, 2)

		noe <- dim(edges)[1]
		steps.left <- steps.left - no.of.current.steps
		current.step <- current.step + no.of.current.steps

	}

	list(n = n, decomp.hist=decomp.hist, tguh.coeffs=tguh.coeffs)

}
