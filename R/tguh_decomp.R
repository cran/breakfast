#' The Tail-Greedy Unbalanced Haar decomposition of a vector
#'
#' This function performs the Tail-Greedy Unbalanced Haar decomposition of the input
#' vector.
#' 
#' The Tail-Greedy Unbalanced Haar decomposition algorithm is described in 
#' "Tail-greedy bottom-up data decompositions and fast multiple change-point 
#' detection", P. Fryzlewicz (2017), preprint.
#'
#' @param x A vector you wish to decompose.
#' @param p Specifies the number of region pairs merged 
#' in each pass through the data, as the proportion of all remaining region pairs. The default is
#' 0.01.
#' @return A list with the following components:
#' \item{n}{The length of \code{x}.}
#' \item{decomp.hist}{The decomposition history: the complete record of the \code{n}-1 steps taken to decompose \code{x}.
#' This is an array of dimensions 4 by 2 by \code{n}-1. Each of the \code{n}-1 matrices of dimensions 4 by 2
#' contains the following: first row - the indices of the regions merged, in increasing order (note: the indexing changes
#' through the transform); second row - the values of the Unbalanced Haar filter coefficients used to produce the
#' corresponding detail coefficient; third row - the (detail coefficient, smooth coefficient) of the decomposition;
#' fourth row - the lengths of (left wing, right wing) of the corresponding Unbalanced Haar wavelet.}
#' \item{tguh.coeffs}{The coefficients of the Tail-Greedy Unbalanced Haar transform of \code{x}.}	
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{tguh.cpt}}, \code{\link{tguh.denoise}}, \code{\link{tguh.reconstr}}
#' @examples
#' rnoise <- rnorm(10)
#' tguh.decomp(rnoise)
#' @export




tguh.decomp <- function(x, p = .01) {

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
