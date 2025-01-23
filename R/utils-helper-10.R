

#' An internal function to compute rcs spline basis
#' 
#' @param x A numeric vector 
#' @param hypothesis A vector 
#' 
#' @return A matrix
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
make_spline_matrix <- function(x, knots) {
  X <- x
  N <- length(X)
  nk <- length(knots)
  basis_evals <- matrix(0, N, nk - 1)
  basis_evals[, 1] <- X
  basis_evals[, 1] <- X
  Xx <- matrix(0, N, nk)
  km1 <- nk - 1
  j <- 1
  knot1 <- knots[1]
  knotnk <- knots[nk]
  knotnk1 <- knots[nk - 1]
  kd <- (knotnk - knot1) ^ (2)
  for (ia in 1:N) {
    for (ja in 1:nk) {
      Xx[ia, ja] <- ifelse(X[ia] - knots[ja] > 0, X[ia] - knots[ja], 0)
    }
  }
  while (j <= nk - 2) {
    jp1 <- j + 1
    basis_evals[, jp1] <-
      (
        Xx[, j] ^ 3 - (Xx[, km1] ^ 3) * (knots[nk] - knots[j]) /
          (knots[nk] - knots[km1]) + (Xx[, nk] ^ 3) *
          (knots[km1] - knots[j]) / (knots[nk] - knots[km1])
      ) /
      (knots[nk] - knots[1]) ^ 2
    j <- j + 1
  }
  return(basis_evals)
}



#' An internal function to b spline basis
#' 
#' @param x A numeric vector
#' @param degree An integer 
#' @param knots A vector
#' @param bknots A vector
#' @param calcderiv A real number
#' 
#' @return A matrix
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
#' 
GS_bs <- function(x, degree, knots, bknots, calcderiv) {
  Nobs       <- length(x)
  fullknots  <- c(rep(bknots[1], degree+1), knots, rep(bknots[2], degree+1))
  Nintervals <- length(fullknots) - 1
  M1         <- matrix(0, Nobs, Nintervals)
  for (i in 1:Nintervals) {
    M1[, i] <- as.numeric(fullknots[i] <= x  & x < fullknots[i + 1])
  }
  lastknot_index <- which(x == bknots[2])
  if (length(lastknot_index) > 0) {
    M1[lastknot_index, (Nintervals - degree):Nintervals] <- 1
  }
  
  for (p in 1:degree) {
    M2 <- matrix(0, Nobs, Nintervals - p)
    for (i in 1:(Nintervals - p)) {
      if (fullknots[i + p] == fullknots[i]) {
        C1 <- 0
      } else {
        C1 <- (x - fullknots[i]) / (fullknots[i + p] - fullknots[i])
      }
      
      if (fullknots[i + p + 1] == fullknots[i + 1]) {
        C2 <- 0
      } else {
        C2 <- (fullknots[i + p + 1] - x) / (fullknots[i + p + 1] - fullknots[i + 1])
      }
      
      M2[, i] <- C1 * M1[, i] + C2 * M1[, i + 1]
    }
    if(p != degree) M1 <- M2
  }
  
  splinevars <- M2 
  if(!calcderiv) {
    out <- splinevars
  }
  
  if(calcderiv) {
    deriv <- matrix(NA, Nobs, Nintervals - degree)
    for (i in 1:(Nintervals - degree)) {
      if (fullknots[i + degree] == fullknots[i]) {
        C1 <- 0
      } else {
        C1 <- degree / (fullknots[i + degree] - fullknots[i])
      }
      
      if (fullknots[i + degree + 1] == fullknots[i + 1]) {
        C2 <- 0
      } else {
        C2 <- degree / (fullknots[i + degree + 1] - fullknots[i + 1])
      }
      
      deriv[, i] <- C1 * M1[, i] - C2 * M1[, i + 1]
    }
    
    dsplinevars <- deriv # [, 2:ncol(deriv)]
    out <- dsplinevars
  }
  return(out)
}



#' An internal function to b spline basis
#' 
#' @param x A numeric vector
#' @param knots A vector
#' @param bknots A vector
#' @param intercept An integer
#' @param calcderiv A real number
#' 
#' @return A matrix
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
#' 
GS_ns <- function(x, knots, bknots, intercept, calcderiv) {
  
  Nintk     <- length(knots) 
  Nk        <- Nintk + 2
  allknots  <- c(rep(bknots[1], 4), knots,rep(bknots[2], 4))
  Nintervals <- length(allknots) - 1
  bs        <- matrix(0, length(x), Nk)
  bsderiv   <- matrix(0, length(x), Nk)
  # Compute B-spline values and derivatives
  if (calcderiv) {
    bs      <- GS_bs(x, 3, knots, bknots, 0)
    bsderiv <- GS_bs(x, 3, knots, bknots, 1)
  } else {
    bs      <- GS_bs(x, 3, knots, bknots, 0)
  }
  # Handle boundary cases (for values below and above the boundary knots)
  x_below_boundary <- sum(x < bknots[1])
  x_above_boundary <- sum(x > bknots[2])
  
  if (x_below_boundary || x_above_boundary) {
    bs_bknots <- matrix(0, 2, Nk)
    bsderiv_bknots <- matrix(0, 2, Nk)
    bs_bknots      <- GS_bs(bknots, 3, knots, bknots, 0)
    bsderiv_bknots <- GS_bs(bknots, 3, knots, bknots, 1)
   
    if (x_below_boundary) {
      xselect <- which(x < bknots[1])
      Nxselect <- length(xselect)
      below_azx1 <- matrix(rep(bs_bknots[1,], Nxselect), byrow = T, nrow = Nxselect) 
      below_azx2 <- matrix(rep(bsderiv_bknots[1,], Nxselect), byrow = T, nrow = Nxselect) 
      below_azx3 <- (x[xselect] - bknots[1]) 
      bs[xselect, ] <- below_azx1 + below_azx2 * below_azx3
      if (calcderiv) bsderiv[xselect, ] <- bsderiv_bknots[1, ]
    }
    if (x_above_boundary) {
      xselect <- which(x > bknots[2])
      Nxselect <- length(xselect)
      above_azx1 <- matrix(rep(bs_bknots[2,], Nxselect), byrow = T, nrow = Nxselect)
      above_azx2 <- matrix(rep(bsderiv_bknots[2,], Nxselect), byrow = T, nrow = Nxselect) 
      above_azx3 <- (x[xselect] - bknots[2]) 
      bs[xselect, ] <- above_azx1 + above_azx2 * above_azx3
      if (calcderiv) bsderiv[xselect, ] <- bsderiv_bknots[2, ]
    }
  }
  
  H <- GS_ns_getH(allknots, 1) # , 1 normalize T/F 1/0
  
  out <- (bs %*% H)
  if (calcderiv) out <- (bsderiv %*% H)
  
  if(intercept) out <- out else out <- out[, 2:ncol(out)]
  
  return(out)
}




#' An internal function to computer H matrix
#' 
#' @param knots A vector
#' @param normalize An integer
#' 
#' @return A matrix
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
#' 
GS_ns_getH <- function(knots, normalize) {
  # Step 1: Define constants based on the knots
  Nintk <- length(knots) - 8
  
  C11 <- 6 / ((knots[5] - knots[2]) * (knots[5] - knots[3]))
  C31 <- 6 / ((knots[6] - knots[3]) * (knots[5] - knots[3]))
  C21 <- -C11 - C31
  Cp22 <- 6 / ((knots[Nintk + 6] - knots[Nintk + 3]) * (knots[Nintk + 6] - knots[Nintk + 4]))
  Cp2  <- 6 / ((knots[Nintk + 7] - knots[Nintk + 4]) * (knots[Nintk + 6] - knots[Nintk + 4]))
  Cp12 <- -Cp22 - Cp2
  # Step 2: Build the matrix H depending on Nintk
  if (Nintk == 0) {
    H <- t(matrix(c(3, 0, 
                    2, 1, 
                    1, 2, 
                    0, 3), nrow = 4, byrow = TRUE))
    # print(H)
  } else if (Nintk == 1) {
    H <- matrix(c(-C21 / C11 ,        1, 0, 0, 0,
                  0, -C31 / C21, 1, -Cp22 / Cp12, 0,
                  # 0, C31 / C11, 1, C31 / Cp22, 0,
                  0, 0, 0, 1, -Cp12 / Cp2),
                nrow = 3, byrow = TRUE)
  } else {
    H1 <- rbind(matrix(1, nrow = 1, ncol = 3), 
                c(0, 1, -C21 / C31), 
                matrix(0, nrow = Nintk - 2, ncol = 3), 
                matrix(0, nrow = 2, ncol = 3))
    
    H2 <- rbind(matrix(0, nrow = 2, ncol = Nintk - 2),
                diag(Nintk - 2), 
                matrix(0, nrow = 2, ncol = Nintk - 2))
    
    H3 <- rbind(matrix(0, nrow = 2, ncol = 3), 
                matrix(0, nrow = Nintk - 2, ncol = 3), 
                c(-Cp12 / Cp22, 1, 0), 
                matrix(1, nrow = 1, ncol = 3))
    H <- cbind(H1, H2, H3)
  }
  
  if(normalize) {
    sumH <- rowSums(H)
    notzero <- which(sumH != 0)
    H[notzero, ] <- H[notzero, ]  / sumH[notzero]
  }
 
  return(t(H))
  
}


#' An internal function to b spline basis
#' 
#' @param x A numeric vector
#' @param knots A vector
#' @param bknots A vector
#' @param intercept An integer
#' @param derivs An integer
#' @param calcderiv A real number
#' 
#' @return A matrix
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
#' 
GS_ns_call <- function(x, knots, bknots, intercept, derivs, centerval) {
  
  if(derivs > 1) {
    stop("Second and higher order derivatives are not supported yet")
  } else {
    calcderiv <- derivs
  }
  
  if(intercept) {
    if(centerval != 0) {
      stop("centerval should be '0' when intercept = TRUE")
      centerval <- 0
    }
  }
  
  if(length(knots) > 0) {
    if(bknots[1] > knots[1]) {
      stop("Lower boundary is greater than lower internal knot")
    } else if(bknots[2] < knots[length(knots)]) {
      stop("Upper Boundary knot is less than upper internal knot")
    }
  }
  
  df     <- length(knots) + 1 + intercept
  out <- GS_ns(x, knots, bknots, intercept = intercept, calcderiv = calcderiv)
  # Centering
  if (centerval != 0) {
    cenout <- GS_ns(centerval, knots, bknots, intercept = intercept, calcderiv = calcderiv)
    if(!is.matrix(cenout)) cenout <- matrix(cenout, nrow = 1) 
    if (!calcderiv) {
      if(intercept) {
        for (i in 2:ncol(cenout)) {
          out[,i] <- out[,i] - cenout[,i]
        }
      } else if(!intercept) {
        for (i in 1:ncol(cenout)) {
          out[,i] <- out[,i] - cenout[,i]
        }
      }
    } else if (calcderiv) {
      out <- out
    }
  } # if (centerval != 0) {
  
  return(out)
}


