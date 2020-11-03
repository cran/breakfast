#' @title Estimating change-points in the piecewise-constant mean of a noisy data sequence via the localised pruning
#' @description This function estimates the number and locations of change-points in the piecewise-constant mean of a noisy data sequence via the localised pruning method, which performs a Schwarz criterion-based model selection on the given candidate set in a localised way.
#' @details 
#'  Further information can be found in Cho and Kirch (2020), arXiv preprint arXiv:1910.12486.
#' 
#' @param cptpath.object A solution-path object, returned by a \code{sol.[name]} routine. Note that the field \code{cptpath.object$x} contains the input data sequence. 
#' @param min.d A number specifying the minimal spacing between change points; \code{min.d = 5} by default
#' @param penalty A string specifying the type of penalty term to be used in Schwarz criterion; possible values are:
#' \itemize{
#'    \item{\code{"log"}}{ Use \code{penalty = log(length(x))^pen.exp}}
#'    \item{\code{"polynomial"}}{ Use \code{penalty = length(x)^pen.exp}}
#' }
#' @param pen.exp Exponent for the penalty term (see \code{penalty})
#' @param do.thr If \code{do.thr = TRUE}, mild threshoding on the CUSUM test statistics is performed after internal standardisation step in order to "pre-prune down" the candidates
#' @param th.const A constant multiplied to \code{sqrt(2*log(length(x)))} to form a mild threshold; if not supplied, a default value (\code{0.5*} a value suggested in Fryzlewicz (2020) is used,
#' see \code{th.const} in \code{\link{model.sdll}}
#' @return An S3 object of class \code{cptmodel}, which contains the following fields: 
#' \item{solution.path}{The solution path method used to obtain \code{cptpath.object}}
#' \item{model}{The model selection method used to return the final change-point estimators object, here its value is \code{"ip"}}
#' \item{no.of.cpt}{The number of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}}
#' \item{cpts}{The locations of estimated change-points in the piecewise-constant mean of the vector \code{cptpath.object$x}. These are the end-points of the corresponding constant-mean intervals}
#' \item{est}{An estimate of the piecewise-constant mean of the vector \code{cptpath.object$x}; the values are the sample means of the data (replicated a suitable number of times) between each pair of consecutive detected change-points}
#' @references H. Cho & C. Kirch (2020) Two-stage data segmentation permitting multiscale change points, heavy tails and dependence. \emph{arXiv preprint arXiv:1910.12486}.
#' @seealso \code{\link{sol.idetect}}, \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.tguh}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{breakfast}}
#' @examples 
#' f <- rep(rep(c(0, 1), each = 50), 10)
#' x <- f + rnorm(length(f))
#' model.lp(sol.not(x))
#' @export
model.lp <- function(cptpath.object, min.d = 5, 
                   penalty = c('log', 'polynomial')[1], pen.exp = 1.01, 
                   do.thr = TRUE, th.const = .5){ 
  
  if(class(cptpath.object) != 'cptpath') stop('object must be in cptpath class')
  if(cptpath.object$method %in% c('idetect', 'idetect_seq')) warning("model.lp won't work well on cptpath.object produced with sol.idetect or sol.idetect_seq; consider using other model. functions, or produce your cptpath.object with a different sol. function")
  
  x <- cptpath.object$x
  cands <- as.matrix(cptpath.object$cands)[, 1:4, drop = FALSE]
  cands <- dup.merge(cands)
  
  n <- length(x)
  
  if(min.d < 1) stop("min.d should be greater than or equal to one")
  if(n < 2 * min.d + 1) stop("For model.lp, x should contain at least max(2*min.d + 1, 3) elements")
  
  if (penalty == 'log') {
    log.penalty <- TRUE
  } else {
    if (penalty != 'polynomial') {
      stop('penalty has to set to log or polynomial')
    }
    log.penalty <- FALSE
  }
  
  cands <- cands[abs(cands[, 4]) > 0,, drop = FALSE]
  
  if(nrow(cands) == 0){
    ret <- structure(list(solution.path = cptpath.object$method,
                          model = 'lp',
                          no.of.cpt = 0,
                          cpts = integer(0), 
                          est = mean.from.cpt(x, integer(0))),
                     class = 'cptmodel')
    return(ret)
  }
  
  cs1 <- cumsum(x)
  cs2 <- cumsum(x^2)
  
  tmp <- cbind(cands[, 3], cands[, 3] - cands[, 1] + 1, cands[, 2] - cands[, 3], cands[, 2] - cands[, 1] + 1, cands[, 4])
  ind <- tmp[, 2] > min.d & tmp[, 3] > min.d
  tmp <- tmp[ind,, drop = FALSE]
  cands <- cands[ind, , drop = FALSE]
  
  if(nrow(tmp) == 0){
    ret <- structure(list(solution.path = cptpath.object$method,
                          model = 'lp',
                          no.of.cpt = 0,
                          cpts = integer(0), 
                          est = mean.from.cpt(x, integer(0))),
                     class = 'cptmodel')
    return(ret)
  }
  
  tmp <- cbind(tmp, t(apply(cands, 1, function(z){
    if(z[1] > 1) denom <- cs2[z[3]] - cs2[z[1] - 1] - (cs1[z[3]] - cs1[z[1] - 1])^2/(z[3] - z[1] + 1) else denom <- cs2[z[3]] - cs1[z[3]]^2/z[3]
    denom <- denom + cs2[z[2]] - cs2[z[3]] - (cs1[z[2]] - cs1[z[3]])^2/(z[2] - z[3])
    denom <- denom/(z[2] - z[1] + 1)
    c(abs(z[4]), sqrt(denom * (denom > 0)))
    })))
  tmp <- tmp[tmp[, 6] > 0, , drop = FALSE]

  if(nrow(tmp) == 0){
    ret <- structure(list(solution.path = cptpath.object$method,
                          model = 'lp',
                          no.of.cpt = 0,
                          cpts = integer(0), 
                          est = mean.from.cpt(x, integer(0))),
                     class = 'cptmodel')
    return(ret)
  }
  
  if(do.thr){
    if(is.null(th.const) || th.const < 0) th.const <- .5 * universal.M.th.v3(n, .9)$th.const
    th <- th.const * sqrt(2*log(n))
    tmp <- tmp[tmp[, 6] > tmp[, 7] * th, , drop = FALSE]
  } 
  
  if(nrow(tmp) == 0){
    ret <- structure(list(solution.path = cptpath.object$method,
                          model = 'lp',
                          no.of.cpt = 0,
                          cpts = integer(0), 
                          est = mean.from.cpt(x, integer(0))),
                     class = 'cptmodel')
    return(ret)
  }
  
  all.cpts <- cbind(as.integer(tmp[, 1]), tmp[, 2:4, drop= FALSE], NA, tmp[, 6, drop = FALSE]/tmp[, 7, drop = FALSE])
  all.cpts <- all.cpts[order(all.cpts[, 1]), , drop = FALSE]
  ac <- nrow(all.cpts)

  if(ac > 0){
    lp <- mosum.local.prune(x, all.cpts, rule = 'jump', log.penalty, pen.exp)
    est.cpts <- lp$est.cpts; est.cpts.ind <- lp$est.cpts.ind; # min.cost <- lp$min.cost
  } else{
    est.cpts.ind <- est.cpts <- integer(0)
    # min.cost <- cs2[n] - cs1[n]^2/n
  }
  # 
  # est.cpts.info <- data.frame(cpts = est.cpts, 
  #                             G.left =  all.cpts[est.cpts.ind, 2], 
  #                             G.right =  all.cpts[est.cpts.ind, 3],
  #                             jump = all.cpts[est.cpts.ind, 6])
  # if (log.penalty) {
  #   penalty_term <- length(est.cpts)*log(n)^pen.exp
  # } else {
  #   penalty_term <- length(est.cpts)*n^pen.exp
  # }
  # final.bic <- n/2*log(min.cost/n) + penalty_term
  
  ret <- structure(list(solution.path = cptpath.object$method,
                        model = 'lp',
                        no.of.cpt = length(est.cpts),
                        cpts = est.cpts, 
                        est = mean.from.cpt(x, est.cpts)), 
                   class = 'cptmodel')
  return(ret)
}

#' Removes duplicates from the candidate set
#' @keywords internal
dup.merge <- function(cands){
  
  all.unique.cpts <- unique(cands[, 3])
  out <- matrix(NA, nrow = 0, ncol = 4)
  for (k in all.unique.cpts) {
    ind <- which(cands[, 3] == k)
    ind.min <- ind[cands[ind, 4] == max(cands[ind, 4])]
    if (length(ind.min) > 1) ind.min <- ind.min[which.min(cands[ind.min, 4])]
    out <- rbind(out, cands[ind.min, ])
  }
  out

}

#' Localised pruning algorithm
#' 
#' @keywords internal
#' @noRd
mosum.local.prune <- function(x, all.cpts, rule, log.penalty, pen.exp){
  
  THRESH_MANUAL_MERGING <- 24
  
  n <- length(x)
  ac <- dim(all.cpts)[1]
  
  all.cpts <- cbind(all.cpts, 1:ac)
  cand_used <- rep(FALSE, ac)
  all.unique.cpts <- c(0, all.cpts[, 1], n)
  auc <- length(all.unique.cpts) - 2
  sums <- matrix(0, nrow = auc+1, ncol = 4) # calculated for efficient computation of rss
  for(j in 1:(auc + 1)){
    sums[j, 1:2] <- c(all.unique.cpts[j] + 1, all.unique.cpts[j + 1])
    sums[j, 3] <- sum(x[sums[j, 1]:sums[j, 2]])
    sums[j, 4] <- sum(x[sums[j, 1]:sums[j, 2]]^2)
  }
  min.cost <- sum(sums[, 4] - sums[, 3]^2/(sums[, 2] - sums[, 1]+1)) # min rss with all the candidates
  
  if(rule == 'pval'){
    u <- all.cpts[order(all.cpts[, 5], all.cpts[, 4], all.cpts[, 2], all.cpts[, 3]),, drop = FALSE]
    rule.seq <- u[, 7]; rm(u)
  }
  if(rule == 'jump'){
    u <- all.cpts[order(-all.cpts[, 6], all.cpts[, 4], all.cpts[, 2], all.cpts[, 3]),, drop = FALSE]
    rule.seq <- u[, 7]; rm(u)
  }
  if(rule == 'lr') rule.seq <- pool
  if(rule == 'rl') rule.seq <- rev(pool)
  
  current <- pool <- seq_len(ac); est.cpts.ind <- est.cpts <- integer(0)
  # current = C, pool = P, est.cpts.ind = B (index)
  while(length(pool)>0){
    # step 1
    j <- rule.seq[1]; adj <- 0
    
    # step 2
    le <- local.env(j, est.cpts.ind, all.cpts, current, ac)
    li <- le$li; li_final <- le$li_final
    ri <- le$ri; ri_final <- le$ri_final
    
    #step 3
    left <- li + 1
    right <- ri - 1
    cand.ind <- (left:right)[is.element(left:right, pool)]
    cand <- all.cpts[cand.ind, 1] # = D
    
    # left <- max(est.cpts.ind[est.cpts.ind < j]+1, j - all.cpts[j, 7])
    # right <- min(est.cpts.ind[est.cpts.ind > j]-1, j + all.cpts[j, 8])
    # if(sum(current < left) > 0) li <- max(current[current < left]) else li <- 0
    # if(sum(current > right) > 0) ri <- min(current[current > right]) else ri <- ac+1  
    
    ind_middl_tmp <- sums[(li + 1):(ri - 1), 2]
    ind_middl_tmp <- ind_middl_tmp[which(!cand_used[(li + 1):(ri - 1)])]
    ind_tmp <- c(sums[li + 1, 1] - 1, ind_middl_tmp, sums[ri, 2])
    sub.sums <- extract_sub(ind_tmp, x)
    
    doExhaustiveSearch <- TRUE
    # Too many candidates to do exhaustive search?
    if (length(cand) > THRESH_MANUAL_MERGING) {
      # Count neighbourhood size of neighbours
      # warn_msg <- paste0('Warning: ', length(cand), ' conflicting candidates')
      # warning(warn_msg)
      
      cand.rule.seq <- rule.seq[is.element(rule.seq, cand.ind)]
      cand_size <- rep(NA, length(cand))
      cand_size[1] <- length(cand)
      for (i_tmp in 2:length(cand)) {
        jj <- cand.rule.seq[i_tmp]
        le_jj <- local.env(jj, est.cpts.ind, all.cpts, current, ac)
        left_jj <- le_jj$li + 1
        right_jj <- le_jj$ri - 1
        cand.ind_jj <- (left_jj:right_jj)[is.element(left_jj:right_jj, pool)]
        #cand_jj <- all.cpts[cand.ind_jj, 1] # = D
        cand_size[i_tmp] <- length(cand.ind_jj)
      }
      
      if (any(cand_size <= THRESH_MANUAL_MERGING)) {
        # Proceed with next candidate, for which exhaustive search IS possible
        rule_tmp <- cand.rule.seq[min(which(cand_size <= THRESH_MANUAL_MERGING))]
        ind_star <- which(rule.seq == rule_tmp)
        rule.seq[ind_star] <- rule.seq[1]; rule.seq[1] <- rule_tmp
        doExhaustiveSearch <- FALSE
      } else {
        # Count neighbourhood size of remaining candidates
        cand_size <- rep(NA, length(rule.seq))
        cand_size[1] <- length(cand)
        for (i_tmp in seq(from = 2, length.out = length(rule.seq) - 1)) {
          jj <- rule.seq[i_tmp]
          le_jj <- local.env(jj, est.cpts.ind, all.cpts, current, ac)
          left_jj <- le_jj$li + 1
          right_jj <- le_jj$ri - 1
          cand.ind_jj <- (left_jj:right_jj)[is.element(left_jj:right_jj, pool)]
          #cand_jj <- all.cpts[cand.ind_jj, 1] # = D
          cand_size[i_tmp] <- length(cand.ind_jj)
        }
        
        if (any(cand_size <= THRESH_MANUAL_MERGING)) {
          # Proceed with next candidate, for which exhaustive search IS possible
          ind_star <- min(which(cand_size <= THRESH_MANUAL_MERGING))
          rule_tmp <- rule.seq[ind_star]; rule.seq[ind_star] <- rule.seq[1]; rule.seq[1] <- rule_tmp
          doExhaustiveSearch <- FALSE
        } else {
          # No more exhaustive search possible at all
          # --> Do manual merging, until exhaustive search becomes possible
          while(length(cand) > THRESH_MANUAL_MERGING) {
            warn_msg <- paste0('Warning: ', length(cand), ' conflicting candidates, thinning manually')
            warning(warn_msg)
            k <- cand[which.min(diff(cand))]
            l <- which(sub.sums[, 2] == k)
            a <- sub.sums[l, ]; b <- sub.sums[l + 1, ]
            # as change-points are merged, the minimum rss in the local environment needs to be updated
            adj <- adj + (a[2] - a[1] + 1)*(b[2] - b[1] + 1)/(b[2] - a[1] + 1)*(a[3]/(a[2] - a[1] + 1) - b[3]/(b[2] - b[1] + 1))^2
            sub.sums[l + 1, 1] <- a[1]; sub.sums[l + 1, 3:4] <- sub.sums[l + 1, 3:4]+a[3:4]
            sub.sums <- sub.sums[-l,, drop = FALSE]
            cand <- setdiff(cand, k)
            k.ind <- which(all.cpts[, 1] == k)
            cand.ind <- setdiff(cand.ind, k.ind); pool <- setdiff(pool, k.ind); rule.seq <- setdiff(rule.seq, k.ind)
            cand_used[k.ind] <- TRUE
          }
        }
      }
    }
    
    if (doExhaustiveSearch) {
      # step 4
      # performs exhaustive search (Algorithm 2)
      out <- exhaust_sc(cand = cand, sub_sums = sub.sums, 
                        strength = pen.exp, log_penalty = log.penalty, 
                        n = n, auc = length(current), min_cost = min.cost)
      est.cpts <- c(est.cpts, out$est_cpts)
      current.est.cpts.ind <- all.cpts[all.cpts[, 1] %in% out$est_cpts, 7]
      est.cpts.ind <- c(est.cpts.ind, current.est.cpts.ind)
      
      # steps 5, 6
      # removal of candidates
      rm.set <- c(j, current.est.cpts.ind)
      if(length(current.est.cpts.ind) > 0){
        rm.set <- c(rm.set, cand.ind[cand.ind %in% min(current.est.cpts.ind):max(current.est.cpts.ind)])
        if(li_final) rm.set <- c(rm.set, cand.ind[cand.ind <= max(current.est.cpts.ind)])
        if(ri_final) rm.set <- c(rm.set, cand.ind[cand.ind >= min(current.est.cpts.ind)])
      }        
      
      # tmp <- (all.cpts[cand.ind, 1] - sub.sums[1, 1] >= all.cpts[cand.ind, 2]) & 
      #   (sub.sums[nrow(sub.sums), 2] - all.cpts[cand.ind, 1] >= all.cpts[cand.ind, 3])
      # if(li > 0) tmp <- tmp & all.cpts[cand.ind, 1] - all.cpts[li, 1] >= all.cpts[li, 3]
      # if(ri < ac + 1) tmp <- tmp & all.cpts[ri, 1] - all.cpts[cand.ind, 1] >= all.cpts[ri, 2]
      # rm.set <- c(rm.set, cand.ind[tmp])
      # 
      # if(length(current.est.cpts.ind) > 0){
      #   rm.set <- c(rm.set, cand.ind[cand.ind %in% min(current.est.cpts.ind):max(current.est.cpts.ind)])
      #   if(li_final) rm.set <- c(rm.set, cand.ind[cand.ind <= max(current.est.cpts.ind)])
      #   if(ri_final) rm.set <- c(rm.set, cand.ind[cand.ind >= min(current.est.cpts.ind)])
      # }
      # rm.set <- min(rm.set):max(rm.set)
      
      pool <- setdiff(pool, rm.set)
      cand_used[rm.set] <- TRUE
      rule.seq <- setdiff(rule.seq, rm.set)
      current <- c(pool, est.cpts.ind)
      current_cands <- is.element(cand, all.cpts[current, 1])
      ind_star <- get_comb_ind(current_cands)
      min.cost <- min.cost + adj - out$sc[nrow(out$sc), 1] + out$sc[ind_star + 1, 1]
    }
  }
  est.cpts <- sort(as.vector(est.cpts)); est.cpts.ind <- sort(as.vector(est.cpts.ind))
  
  return(list(est.cpts = est.cpts, est.cpts.ind = est.cpts.ind, min.cost = min.cost))
}

#' Identify the local environment for exhaustive search
#' @keywords internal
#' @noRd
local.env <- function(j, est.cpts.ind, all.cpts, current, ac){
  li_final <- ri_final <- TRUE
  if(sum(est.cpts.ind < j)) li <- max(est.cpts.ind[est.cpts.ind < j]) else li <- 0
  if(j > 1){
    ind <- ((li + 1):(j - 1))[(li + 1):(j - 1) %in% current]
    if(length(ind) > 0 && sum(tmp <- all.cpts[j, 1] - all.cpts[ind, 1] >= pmax(all.cpts[j, 2], all.cpts[ind, 3]))) li_tmp <- max(ind[tmp]) else li_tmp <- li
    if(li_tmp > li){ li <- li_tmp; li_final <- FALSE }
  }
  
  if(sum(est.cpts.ind > j)) ri <- min(est.cpts.ind[est.cpts.ind > j]) else ri <- ac + 1
  if(j < ac){
    ind <- ((j + 1):(ri - 1))[(j + 1):(ri - 1) %in% current]
    if(length(ind) > 0 && sum(tmp <- all.cpts[ind, 1] - all.cpts[j, 1] >= pmax(all.cpts[j, 3], all.cpts[ind, 2]))) ri_tmp <- min(ind[tmp]) else ri_tmp <- ri
    if(ri_tmp < ri){ ri <- ri_tmp; ri_final <- FALSE }
  }
  
  list(li = li, li_final = li_final, ri = ri, ri_final = ri_final)
}