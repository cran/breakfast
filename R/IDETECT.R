#' @title Solution path generation via the Isolate-Detect method
#' @description This function arranges all possible change-points in the mean of the input vector, or in its linear trend, in the order of importance, via the Isolate-Detect (ID) method.
#' It is developed to be used with the sdll and information criterion (ic) model selection rules.
#' @details
#' The Isolate-Detect method and its algorithm is described in 
#' "Detecting multiple generalized change-points by isolating single ones", A. Anastasiou & P. Fryzlewicz (2022), Metrika, https://doi.org/10.1007/s00184-021-00821-6.
#' @param x A numeric vector containing the data to be processed.
#' @param type The model type considered. \code{type = "const"}, \code{type = "lin.cont"}, \code{type = "lin.discont"} mean, respectively, that the signal (mean of \code{x}) is piecewise constant, 
#'   piecewise linear and continuous, and piecewise linear but not necessarily continuous. If not given, the default is \code{type = "const"}
#' @param thr_ic_cons A positive real number with default value equal to 0.9. It is used to create the solution path for the piecewise-constant model. The lower the value, the longer the solution path.
#' @param thr_ic_lin A positive real number with default value 1.25. Used to create the solution path if \code{type == "lin.cont"} or \code{type == "lin.discont"}
#' @param points A positive integer with default value equal to 3. It defines the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, as described in the Isolate-Detect methodology.
#' @return An S3 object of class \code{cptpath}, which contains the following fields: 
#' \item{solutions.nested}{\code{TRUE}, i.e., the change-point outputs are nested}
#' \item{solution.path}{Locations of possible change-points in the mean of \code{x}, arranged in decreasing order of change-point importance}
#' \item{solution.set}{Empty list}
#' \item{x}{Input vector \code{x}}
#' \item{type}{The input parameter \code{type}}
#' \item{cands}{Matrix of dimensions length(\code{x}) - 1 by 4. The first two columns are (start, end)-points of the detection intervals of the corresponding possible change-point location in the third column. The fourth column is a measure of strength of the corresponding possible change-point. The order of the rows is the same as the order returned in \code{solution.path}}
#' \item{method}{The method used, which has value "idetect" here}
#' @seealso \code{\link{sol.idetect_seq}}, \code{\link{sol.not}}, \code{\link{sol.wbs}}, \code{\link{sol.wbs2}}, \code{\link{sol.tguh}}
#' @references A. Anastasiou & P. Fryzlewicz (2022). Detecting multiple generalized change-points by isolating single ones. \emph{Metrika}, https://doi.org/10.1007/s00184-021-00821-6.
#' @examples
#' r3 <- rnorm(1000) + c(rep(0,300), rep(2,200), rep(-4,300), rep(0,200))
#' sol.idetect(r3)
#' @export
sol.idetect <- function(x, type = "const", thr_ic_cons = 0.9, thr_ic_lin = 1.25, points = 3) {
  #veryfing the input parameters - x
  x <- as.numeric(x)
  storage.mode(x) <- "double"
  lx <- length(x)
  check.input(x)
  stopifnot(type %in% c('const', 'lin.cont', 'lin.discont'))
  solutions.nested <- TRUE
  solution.set <- list()
  cands <- matrix(NA, 0, 4)
  if (lx < points) {solution.path <- integer()}
  else{
    points <- as.integer(points)
    if(type=="lin.cont")  {step1 <- window.idetect.th_lin(x, thr_con = thr_ic_lin, w_points = points)}
    if(type=="lin.discont"){step1 <- window.idetect.th_discont_lin(x, thr_con = thr_ic_lin, w_points = points)}
    if(type == "const"){step1 <- window.idetect.th(x, thr_con = thr_ic_cons, w_points = points)}
    s1 <- as.matrix(step1$full_information)
    #    print(s1)
    if (dim(s1)[2] == 1) {s1 <- t(s1)}
    cpt_lower <- step1[[1]]
    ord <- order(s1[,4], decreasing = T)
    cands <- s1[ord, ,drop=FALSE]
    if (length(cpt_lower) == 0){solution.path <- integer(0)
    ord <- order(s1[,4], decreasing = T)
    cands <- s1[ord, ,drop=FALSE]
    }
    else{
    lcpt_ic <- length(cpt_lower)
    seb_set <- c(0, cpt_lower, lx)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    CS <- matrix(cpt_lower,1,lcpt_ic)
    CS_start <- matrix(seb_set[1:(lseb_set - 2)] + 1,1,lcpt_ic)
    CS_end <- matrix(seb_set[3:(lseb_set)],1,lcpt_ic)
    while (lseb_set >= 3) {
      if (type == "lin.cont"){
        Rs <- IDetect_linear_contr_one(x, seb_set[1:(lseb_set - 2)] + 1, seb_set[3:(lseb_set)],
                                       seb_set[2:(lseb_set - 1)])}
      if (type == "lin.discont"){
        Rs <- IDetect_linear_discontr_one(x, seb_set[1:(lseb_set - 2)] + 1, seb_set[3:(lseb_set)],
                                          seb_set[2:(lseb_set - 1)])}
      if(type == "const"){Rs <- IDetect_cusum_one(x, seb_set[1:(lseb_set - 2)] + 1, seb_set[3:(lseb_set)],
                                                  seb_set[2:(lseb_set - 1)])}
      #      print(Rs)
      indic <- c(which.min(Rs))
      s1 <- setdiff(cpt_lower, seb_set[2:(lseb_set - 1)])
      d <- numeric(lcpt_ic)
      d_s <- numeric(lcpt_ic)
      d_e <- numeric(lcpt_ic)
      if(length(s1) == 0){d <- Rs
      d_s <- seb_set[1:(lseb_set - 2)] + 1
      d_e <- seb_set[3:(lseb_set)]}
      else{
        indic2 <- match(s1, cpt_lower)
        d[-indic2] <- Rs
        d_s[-indic2] <- seb_set[1:(lseb_set - 2)] + 1
        d_e[-indic2] <- seb_set[3:(lseb_set)]}
      CS <- rbind(CS,d)
      CS_start <- rbind(CS_start,d_s)
      CS_end <- rbind(CS_end,d_e)
      m_rs <- min(Rs)
      seb_set_temp <- seb_set[2:(lseb_set-1)]
      min_Rs <- seb_set_temp[indic]
      cands <- rbind(cands, c(seb_set[indic]+1, seb_set[indic + 2], min_Rs, m_rs))
      min_C <- c(min_C, min_Rs)
      if (seb_set_temp[indic] == 1){
        seb_set <- seb_set[-2]
      } else{
        seb_set <- seb_set[-which(seb_set == min_Rs)]}
      lseb_set <- lseb_set - 1
    }
    #   print(cands)
    solution.path <- min_C[length(min_C):1]
    #cusum_m <- apply(CS[-1,,drop = FALSE],2,max)
    #indic3 <- match(cpt_lower, cands[,3])
    #cands[indic3,4] <- cusum_m
   # ord <- order(cands[,4], decreasing = T)
    #cands <- cands[ord, ,drop=FALSE]#[-(length(solution.path)+1), ,drop = FALSE]
    #    print(cands)
    #cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
    #cands <- cands[!cands[,4] == 0,,drop = FALSE]

    #print(CS)
    #print(CS_start)
    #print(CS_end)
  cands2 <- matrix(NA, lcpt_ic,4)
  cands2[,4] <- apply(CS[2:nrow(CS),,drop = FALSE], 2, max)
  cands2[,3] <- cpt_lower
  for (i in 1:lcpt_ic){
    pos <- which(CS[,i] == cands2[,4][i])[1]
    cands2[i,1] <- CS_start[pos,i]
    cands2[i,2] <- CS_end[pos,i]
  }
  ord2 <- order(cands2[,4], decreasing = T)
  cands2 <- cands2[ord2, ,drop=FALSE]
  ord <- order(cands[,4], decreasing = T)
  cands <- cands[ord, ,drop=FALSE]
  cands <- rbind(cands2, cands)
    }
    cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
    cands <- cands[!cands[,4] == 0,,drop = FALSE]
    }
#    if(is.na(solution.path[1])){solution.path <- integer(0)}}
  ret <- list()
  ret$solutions.nested = TRUE
  ret$solution.path <- solution.path
  ret$solution.set <- solution.set
  ret$x = x
  ret$type = type
  ret$cands = as.matrix(cands)
  ret$method = "idetect"
  class(ret) <- "cptpath"
  
  return(ret)
}

idetect.th <- function(x, sigma = stats::mad(diff(x) / sqrt(2)), thr_const = 1,
                               thr_fin = sigma * thr_const * sqrt(2 * log(length(x))),
                               s = 1, e = length(x), points = 3, k_l = 1, k_r = 1) {
  y <- c(0, cumsum(x))
  l <- length(x)
  if (sigma == 0) {
    s0 <- all_shifts_are_cpts(x)
    cpt <- s0$cpts
    l.cpt <- length(cpt)
    if(l.cpt == 0){
      cpt <- integer(0)
      Res_fin <- matrix(0, 1, 4)
    }
    else{
      if(l.cpt == 1){
        CS <- IDetect_cusum_one(x, s = 1, e = l, b = cpt)
        Res_fin <- cbind(1,l,cpt,CS)
      }
      else{
    s1 <- c(1,cpt[-l.cpt]+1)
    e1 <- c(cpt[2:l.cpt], l)
    CS <- IDetect_cusum_one(x, s = s1, e = e1, b = cpt)
    Res_fin <- cbind(s1,e1,cpt,CS)
    }
    }
  }
  else{
  Res <- matrix(0, 1, 4)
  points <- as.integer(points)
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s < 2) {
    Res_fin <- matrix(0, 1, 4)
    cpt <- integer(0)
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- start_end_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ( (chp == 0) & (k_r < min(k_l, rur))) {
        ind <- c(s, right_points[k_r])
        tmp <- max_cusum(ind, y)
        pos_r[k_r] <- tmp[1]
        CUSUM_r[k_r] <- tmp[2]
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
          indic <- 0
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ( (chp == 0) & (k_l < min(k_r, lur))) {
        ind <- c(left_points[k_l], e)
        tmp <- max_cusum(ind, y)
        pos_l[k_l] <- tmp[1]
        CUSUM_l[k_l] <- tmp[2]
        Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
        if (CUSUM_l[k_l] > thr_fin) {
          chp <- pos_l[k_l]
          indic <- 1
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        ind <- c(s, right_points[k_r])
        tmp <- max_cusum(ind, y)
        pos_r[k_r] <- tmp[1]
        CUSUM_r[k_r] <- tmp[2]
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
          indic <- 0
        } else {
          ind <- c(left_points[k_l], e)
          tmp <- max_cusum(ind, y)
          pos_l[k_l] <- tmp[1]
          CUSUM_l[k_l] <- tmp[2]
          Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
            indic <- 1
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (indic == 1) {
        r <- idetect.th(x, s = s, e = chp, points = points,
                        thr_fin = thr_fin, k_r = k_r, k_l = 1)
      } else {
        r <- idetect.th(x, s = chp + 1, e = e, points = points,
                        thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
      }
      cpt <- c(chp, r[[1]])
      Res_fin <- rbind(Res, r[[2]])
    } else {
      cpt <- chp
      Res_fin <- Res
    }
  }
  cpt <- cpt[cpt != 0]
  Res_fin <- Res_fin[which(Res_fin[,3] != 0), , drop = FALSE]
  }
  return(list(changepoints = sort(cpt), full_information = Res_fin, y = y))
}

window.idetect.th <- function(xd, sigma = stats::mad(diff(xd) / sqrt(2)), thr_con = 1,
                       c_win = 5000, w_points = 3) {
  
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= c_win) {
    u <- idetect.th(x = xd, thr_const = thr_con, points = w_points)
    return(u)
  } else {
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- list()
    uaddition <- list()
    tsm[[1]] <- xd[1:c_win]
    ufin <- idetect.th(tsm[[1]], thr_fin = t, points = w_points)
    uaddition[[1]] <- list()
    uaddition[[1]] <- idetect.th(x = xd[(max(1, c_win - (10 * w_points) + 1)):min( (c_win + (10 * w_points)), lg)], thr_fin = t, points = 2)
    uaddition[[1]][[1]] <- uaddition[[1]][[1]] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,1] <- uaddition[[1]][[2]][,1] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,2] <- min(uaddition[[1]][[2]][,2] + c_win - (10 * w_points),min( (c_win + (10 * w_points)), lg))
    uaddition[[1]][[2]][,3] <- uaddition[[1]][[2]][,3] + c_win - (10 * w_points)
    ufin[[1]] <- c(ufin[[1]], uaddition[[1]][[1]])
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- list()
      u[[i]] <- idetect.th(x = tsm[[i]], thr_fin = t, points = w_points)
      u[[i]][[1]] <- u[[i]][[1]] + (i - 1) * c_win
      u[[i]][[2]][,1] <- u[[i]][[2]][,1] + (i - 1) * c_win
      u[[i]][[2]][,2] <- u[[i]][[2]][,2] + (i - 1) * c_win
      u[[i]][[2]][,3] <- u[[i]][[2]][,3] + (i - 1) * c_win
      uaddition[[i]] <- list()
      uaddition[[i]] <- idetect.th(x = xd[(max(1, i * c_win - (10 * w_points) + 1)):(min(i * c_win + (10 * w_points), lg))], thr_fin = t, points = 2)
      uaddition[[i]][[1]] <- uaddition[[i]][[1]] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,1] <- uaddition[[i]][[2]][,1] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,2] <- min(uaddition[[i]][[2]][,2] + i * c_win - (10 * w_points), min(i * c_win + (10 * w_points), lg))
      uaddition[[i]][[2]][,3] <- uaddition[[i]][[2]][,3] + i * c_win - (10 * w_points)
      ufin[[1]] <- c(ufin[[1]],u[[i]][[1]], uaddition[[i]][[1]])
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- list()
    u[[K]] <- idetect.th(tsm[[K]], thr_fin = t, points = w_points)
    u[[K]][[1]] <- u[[K]][[1]]  + (K - 1) * c_win
    u[[K]][[2]][,1] <- u[[K]][[2]][,1]  + (K - 1) * c_win
    u[[K]][[2]][,2] <- u[[K]][[2]][,2]  + (K - 1) * c_win
    u[[K]][[2]][,3] <- u[[K]][[2]][,3]  + (K - 1) * c_win
    ufin_cpt <- c(ufin[[1]], u[[K]][[1]])
    Res_fin <- matrix(0, 1, 4)
    Res_fin <- rbind(Res_fin, ufin[[2]], uaddition[[1]][[2]])
    if (K > 2){
    for (i in 2:(K-1)){Res_fin <- rbind(Res_fin,u[[i]][[2]], uaddition[[i]][[2]])}}
    Res_fin <- rbind(Res_fin, u[[K]][[2]])
    Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
    return(list(changepoints = unique(sort(ufin_cpt)), full_information = Res_fin,  y = c(0, cumsum(xd))))
    #return(list(changepoints = sort(ufin_cpt), full_information = Res_fin,  y = c(0, cumsum(xd))))
    
  }
}

idetect.th_lin <- function(x, sigma = stats::mad(diff(diff(x)))/sqrt(6), thr_const = 1.4, 
                           thr_fin = sigma * thr_const * sqrt(2 * log(length(x))), s = 1, 
                           e = length(x), points = 3, k_l = 1, k_r = 1) {
  l <- length(x)
  if (sigma == 0) {
    s0 <- all_slopechanges_are_cpts(x)
    cpt <- s0$cpts
    l.cpt <- length(cpt)
    if(l.cpt == 0){
      cpt <- integer(0)
      Res_fin <- matrix(0, 1, 4)
    }
    else{
      if(l.cpt == 1){
        CS <- IDetect_linear_contr_one(x, s = 1, e = l, b = cpt)
        Res_fin <- cbind(1,l,cpt,CS)
      }
      else{
        s1 <- c(1,cpt[-l.cpt]+1)
        e1 <- c(cpt[2:l.cpt], l)
        CS <- IDetect_linear_contr_one(x, s = s1, e = e1, b = cpt)
        Res_fin <- cbind(s1,e1,cpt,CS)
      }
    }
  }
  else{
    Res <- matrix(0, 1, 4)
    points <- as.integer(points)
    r_e_points <- seq(points, l, points)
    l_e_points <- seq(l - points + 1, 1, -points)
    chp <- 0
    if (e - s < 2) {
      Res_fin <- matrix(0, 1, 4)
      cpt <- integer(0)
    } else {
      pos_r <- numeric()
      CUSUM_r <- numeric()
      pos_l <- numeric()
      CUSUM_l <- numeric()
      moving_points <- start_end_points(r_e_points, l_e_points, s, e)
      right_points <- moving_points[[1]]
      left_points <- moving_points[[2]]
      lur <- length(left_points)
      rur <- length(right_points)
      if (k_r < k_l) {
        while ( (chp == 0) & (k_r < min(k_l, rur))) {
          x_temp_r <- x[s:right_points[k_r]]
          ipcr <- IDetect_cumsum_lin(x_temp_r)
          pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
          CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
          Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
          if (CUSUM_r[k_r] > thr_fin) {
            chp <- pos_r[k_r]
            indic <- 0
          } else {
            k_r <- k_r + 1
          }
        }
      }
      if (k_l < k_r) {
        while ( (chp == 0) & (k_l < min(k_r, lur))) {
          x_temp_l <- x[left_points[k_l]:e]
          ipcl <- IDetect_cumsum_lin(x_temp_l)
          pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
          CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
          Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
            indic <- 1
          } else {
            k_l <- k_l + 1
          }
        }
      }
      if (chp == 0) {
        while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
          x_temp_r <- x[s:right_points[k_r]]
          ipcr <- IDetect_cumsum_lin(x_temp_r)
          pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
          CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
          Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
          if (CUSUM_r[k_r] > thr_fin) {
            chp <- pos_r[k_r]
            indic <- 0
          } else {
            x_temp_l <- x[left_points[k_l]:e]
            ipcl <- IDetect_cumsum_lin(x_temp_l)
            pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
            CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
            Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
            if (CUSUM_l[k_l] > thr_fin) {
              chp <- pos_l[k_l]
              indic <- 1
            } else {
              k_r <- k_r + 1
              k_l <- k_l + 1
            }
          }
        }
      }
      if (chp != 0) {
        if (indic == 1) {
          r <- idetect.th_lin(x, s = s, e = chp, points = points,
                              thr_fin = thr_fin, k_r = k_r, k_l = 1)
        } else {
          r <- idetect.th_lin(x, s = chp + 1, e = e, points = points,
                              thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
        }
        cpt <- c(chp, r[[1]])
        Res_fin <- rbind(Res, r[[2]])
      } else {
        cpt <- chp
        Res_fin <- Res
      }
    }
    cpt <- cpt[cpt != 0]
    Res_fin <- Res_fin[which(Res_fin[,3] != 0), , drop = FALSE]
  }
  return(list(changepoints = sort(cpt), full_information = Res_fin))
}

window.idetect.th_lin <- function(xd, sigma = stats::mad(diff(diff(xd)) / sqrt(6)), thr_con = 1.4,
                                  c_win = 5000, w_points = 3) {
  
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= c_win) {
    u <- idetect.th_lin(x = xd, thr_const = thr_con, points = w_points)
    return(u)
  } else {
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- list()
    uaddition <- list()
    tsm[[1]] <- xd[1:c_win]
    ufin <- idetect.th_lin(tsm[[1]], thr_fin = t, points = w_points)
    uaddition[[1]] <- list()
    uaddition[[1]] <- idetect.th_lin(x = xd[(max(1, c_win - (10 * w_points) + 1)):min( (c_win + (10 * w_points)), lg)], thr_fin = t, points = 2)
    uaddition[[1]][[1]] <- uaddition[[1]][[1]] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,1] <- uaddition[[1]][[2]][,1] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,2] <- min(uaddition[[1]][[2]][,2] + c_win - (10 * w_points),min( (c_win + (10 * w_points)), lg))
    uaddition[[1]][[2]][,3] <- uaddition[[1]][[2]][,3] + c_win - (10 * w_points)
    ufin[[1]] <- c(ufin[[1]], uaddition[[1]][[1]])
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- list()
      u[[i]] <- idetect.th_lin(x = tsm[[i]], thr_fin = t, points = w_points)
      u[[i]][[1]] <- u[[i]][[1]] + (i - 1) * c_win
      u[[i]][[2]][,1] <- u[[i]][[2]][,1] + (i - 1) * c_win
      u[[i]][[2]][,2] <- u[[i]][[2]][,2] + (i - 1) * c_win
      u[[i]][[2]][,3] <- u[[i]][[2]][,3] + (i - 1) * c_win
      uaddition[[i]] <- list()
      uaddition[[i]] <- idetect.th_lin(x = xd[(max(1, i * c_win - (10 * w_points) + 1)):(min(i * c_win + (10 * w_points), lg))], thr_fin = t, points = 2)
      uaddition[[i]][[1]] <- uaddition[[i]][[1]] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,1] <- uaddition[[i]][[2]][,1] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,2] <- min(uaddition[[i]][[2]][,2] + i * c_win - (10 * w_points), min(i * c_win + (10 * w_points), lg))
      uaddition[[i]][[2]][,3] <- uaddition[[i]][[2]][,3] + i * c_win - (10 * w_points)
      ufin[[1]] <- c(ufin[[1]],u[[i]][[1]], uaddition[[i]][[1]])
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- list()
    u[[K]] <- idetect.th_lin(tsm[[K]], thr_fin = t, points = w_points)
    u[[K]][[1]] <- u[[K]][[1]]  + (K - 1) * c_win
    u[[K]][[2]][,1] <- u[[K]][[2]][,1]  + (K - 1) * c_win
    u[[K]][[2]][,2] <- u[[K]][[2]][,2]  + (K - 1) * c_win
    u[[K]][[2]][,3] <- u[[K]][[2]][,3]  + (K - 1) * c_win
    ufin_cpt <- c(ufin[[1]], u[[K]][[1]])
    Res_fin <- matrix(0, 1, 4)
    Res_fin <- rbind(Res_fin, ufin[[2]], uaddition[[1]][[2]])
    if (K > 2){
      for (i in 2:(K-1)){Res_fin <- rbind(Res_fin,u[[i]][[2]], uaddition[[i]][[2]])}}
    Res_fin <- rbind(Res_fin, u[[K]][[2]])
    Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
    return(list(changepoints = unique(sort(ufin_cpt)), full_information = Res_fin))
    #return(list(changepoints = sort(ufin_cpt), full_information = Res_fin,  y = c(0, cumsum(xd))))
    
  }
}

idetect_th_discont_lin <- function(x, sigma = stats::mad(diff(diff(x)))/sqrt(6), thr_const = 1.4, 
                           thr_fin = sigma * thr_const * sqrt(2 * log(length(x))), s = 1, 
                           e = length(x), points = 4, k_l = 1, k_r = 1) {
  if (points < 3) {
    stop("The input in `points' should be an integer greater than 2.")}
  l <- length(x)
  if (sigma == 0) {
    s0 <- all_slopechanges_are_cpts(x)
    cpt <- s0$cpts
    l.cpt <- length(cpt)
    if(l.cpt == 0){
      cpt <- integer(0)
      Res_fin <- matrix(0, 1, 4)
    }
    else{
      if(l.cpt == 1){
        CS <- IDetect_linear_contr_one(x, s = 1, e = l, b = cpt)
        Res_fin <- cbind(1,l,cpt,CS)
      }
      else{
        s1 <- c(1,cpt[-l.cpt]+1)
        e1 <- c(cpt[2:l.cpt], l)
        CS <- IDetect_linear_contr_one(x, s = s1, e = e1, b = cpt)
        Res_fin <- cbind(s1,e1,cpt,CS)
      }
    }
  }
  else{
    Res <- matrix(0, 1, 4)
    if (e - s < 4) {
      Res_fin <- matrix(0, 1, 4)
      cpt <- integer(0)
    }
    else {
    #print(e)
    #print(s)
    points <- as.integer(points)
    r_e_points <- seq(points, l, points)
    l_e_points <- seq(l - points + 1, 1, -points)
    chp <- 0
      pos_r <- numeric()
      CUSUM_r <- numeric()
      pos_l <- numeric()
      CUSUM_l <- numeric()
      moving_points <- start_end_points(r_e_points, l_e_points, s, e)
      right_points <- moving_points[[1]]
      left_points <- moving_points[[2]]
      lur <- length(left_points)
      rur <- length(right_points)
      #print(chp)
      #print(k_r)
      #print(k_l)
      if (k_r < k_l) {
        while ( (chp == 0) & (k_r < min(k_l, rur))) {
          ma <- max(right_points[k_r], s+points)
          x_temp_r <- x[s:ma]
          l_temp_r <- length(x_temp_r)
          ipcr <- IDetect_linear_discontr_one(x_temp_r, s = rep(1, l_temp_r - 2), e = rep(l_temp_r, l_temp_r - 2), b = 1:(l_temp_r - 2))
          l_ipcr <- length(ipcr)
          if(isTRUE(all.equal(ipcr[l_ipcr], ipcr[l_ipcr - 1]))){
          pos_r[k_r] <- which.max(ipcr) + s - 1
          }
          else{
            pos_r[k_r] <- which.max(ipcr) + s - 1
          }
          CUSUM_r[k_r] <- max(ipcr)
          Res <- rbind(Res, c(s,ma,pos_r[k_r],CUSUM_r[k_r]))
          if (CUSUM_r[k_r] > thr_fin) {
            chp <- pos_r[k_r]
            indic <- 0
          } else {
            k_r <- k_r + 1
          }
        }
      }
      if (k_l < k_r) {
        while ( (chp == 0) & (k_l < min(k_r, lur))) {
          mi <-  min(left_points[k_l], e-4)
          x_temp_l <- x[mi:e]
          l_temp_l <- length(x_temp_l)
          ipcl <- IDetect_linear_discontr_one(x_temp_l, s = rep(1, l_temp_l - 2), b = 1:(l_temp_l - 2), e = rep(l_temp_l, l_temp_l - 2))
          l_ipcl <- length(ipcl)
          if(isTRUE(all.equal(ipcl[l_ipcl], ipcl[l_ipcl - 1]))){
            pos_l[k_l] <- which.max(ipcl) + left_points[k_l] -1
          }
          else{
            pos_l[k_l] <- which.max(ipcl) + left_points[k_l] - 1
          }
          CUSUM_l[k_l] <- max(ipcl)
          Res <- rbind(Res, c(mi, e, pos_l[k_l], CUSUM_l[k_l]))
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
            indic <- 1
          } else {
            k_l <- k_l + 1
          }
        }
      }
      if (chp == 0) {
        while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
          ma <- max(right_points[k_r], s+points)
          x_temp_r <- x[s:ma]
          l_temp_r <- ma - s + 1
          ipcr <- IDetect_linear_discontr_one(x_temp_r, s = rep(1, l_temp_r - 2), b = 1:(l_temp_r - 2), e = rep(l_temp_r, l_temp_r - 2))
          l_ipcr <- length(ipcr)
          if(isTRUE(all.equal(ipcr[l_ipcr], ipcr[l_ipcr - 1]))){
            pos_r[k_r] <- which.max(ipcr) + s - 1
          }
          else{
            pos_r[k_r] <- which.max(ipcr) + s - 1
          }
          CUSUM_r[k_r] <- max(ipcr)
          Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
          if (CUSUM_r[k_r] > thr_fin) {
            chp <- pos_r[k_r]
            indic <- 0
          }
          else {
            mi <- min(left_points[k_l], e-4)
            x_temp_l <- x[mi:e]
            l_temp_l <- e - mi + 1
            ipcl <- IDetect_linear_discontr_one(x_temp_l, s = rep(1, l_temp_l - 2), b = 1:(l_temp_l - 2), e = rep(l_temp_l, l_temp_l - 2))
            l_ipcl <- length(ipcl)
            if(isTRUE(all.equal(ipcl[l_ipcl], ipcl[l_ipcl - 1]))){
              pos_l[k_l] <- which.max(ipcl) + left_points[k_l] - 1
            }
            else{
              pos_l[k_l] <- which.max(ipcl) + left_points[k_l] - 1
            }
            CUSUM_l[k_l] <- max(ipcl)
            Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
            if (CUSUM_l[k_l] > thr_fin) {
              chp <- pos_l[k_l]
              indic <- 1
            } else {
              k_r <- k_r + 1
              k_l <- k_l + 1
            }
          }
        }
      }
      #print(chp)
      #print(k_r)
      #print(k_l)
      #print(indic)
      if (chp != 0) {
        if (indic == 1) {
          r <- idetect_th_discont_lin(x, s = s, e = chp-1, points = points,
                              thr_fin = thr_fin, k_r = k_r, k_l = 1)
        } else {
          r <- idetect_th_discont_lin(x, s = chp + 2, e = e, points = points,
                              thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
        }
        cpt <- c(chp, r[[1]])
        Res_fin <- rbind(Res, r[[2]])
      } else {
        cpt <- chp
        Res_fin <- Res
      }
    }
    cpt <- cpt[cpt != 0]
    Res_fin <- Res_fin[which(Res_fin[,3] != 0), , drop = FALSE]
  }
  return(list(changepoints = sort(cpt), full_information = Res_fin))
}

window.idetect.th_discont_lin <- function(xd, sigma = stats::mad(diff(diff(xd)) / sqrt(6)), thr_con = 1.4,
                                  c_win = 5000, w_points = 4) {
  
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= c_win) {
    u <- idetect_th_discont_lin(x = xd, thr_const = thr_con, points = w_points)
    return(u)
  } else {
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- list()
    uaddition <- list()
    tsm[[1]] <- xd[1:c_win]
    ufin <- idetect_th_discont_lin(tsm[[1]], thr_fin = t, points = w_points)
    uaddition[[1]] <- list()
    uaddition[[1]] <- idetect_th_discont_lin(x = xd[(max(1, c_win - (10 * w_points) + 1)):min( (c_win + (10 * w_points)), lg)], thr_fin = t, points = 2)
    uaddition[[1]][[1]] <- uaddition[[1]][[1]] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,1] <- uaddition[[1]][[2]][,1] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,2] <- min(uaddition[[1]][[2]][,2] + c_win - (10 * w_points),min( (c_win + (10 * w_points)), lg))
    uaddition[[1]][[2]][,3] <- uaddition[[1]][[2]][,3] + c_win - (10 * w_points)
    ufin[[1]] <- c(ufin[[1]], uaddition[[1]][[1]])
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- list()
      u[[i]] <- idetect_th_discont_lin(x = tsm[[i]], thr_fin = t, points = w_points)
      u[[i]][[1]] <- u[[i]][[1]] + (i - 1) * c_win
      u[[i]][[2]][,1] <- u[[i]][[2]][,1] + (i - 1) * c_win
      u[[i]][[2]][,2] <- u[[i]][[2]][,2] + (i - 1) * c_win
      u[[i]][[2]][,3] <- u[[i]][[2]][,3] + (i - 1) * c_win
      uaddition[[i]] <- list()
      uaddition[[i]] <- idetect_th_discont_lin(x = xd[(max(1, i * c_win - (10 * w_points) + 1)):(min(i * c_win + (10 * w_points), lg))], thr_fin = t, points = 2)
      uaddition[[i]][[1]] <- uaddition[[i]][[1]] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,1] <- uaddition[[i]][[2]][,1] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,2] <- min(uaddition[[i]][[2]][,2] + i * c_win - (10 * w_points), min(i * c_win + (10 * w_points), lg))
      uaddition[[i]][[2]][,3] <- uaddition[[i]][[2]][,3] + i * c_win - (10 * w_points)
      ufin[[1]] <- c(ufin[[1]],u[[i]][[1]], uaddition[[i]][[1]])
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- list()
    u[[K]] <- idetect_th_discont_lin(tsm[[K]], thr_fin = t, points = w_points)
    u[[K]][[1]] <- u[[K]][[1]]  + (K - 1) * c_win
    u[[K]][[2]][,1] <- u[[K]][[2]][,1]  + (K - 1) * c_win
    u[[K]][[2]][,2] <- u[[K]][[2]][,2]  + (K - 1) * c_win
    u[[K]][[2]][,3] <- u[[K]][[2]][,3]  + (K - 1) * c_win
    ufin_cpt <- c(ufin[[1]], u[[K]][[1]])
    Res_fin <- matrix(0, 1, 4)
    Res_fin <- rbind(Res_fin, ufin[[2]], uaddition[[1]][[2]])
    if (K > 2){
      for (i in 2:(K-1)){Res_fin <- rbind(Res_fin,u[[i]][[2]], uaddition[[i]][[2]])}}
    Res_fin <- rbind(Res_fin, u[[K]][[2]])
    Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
    return(list(changepoints = unique(sort(ufin_cpt)), full_information = Res_fin))
    #return(list(changepoints = sort(ufin_cpt), full_information = Res_fin,  y = c(0, cumsum(xd))))
    
  }
}

