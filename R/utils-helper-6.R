



#' An internal function to prepare export rcs string for Stan, and auxiliary 
#' R function
#' 
#' Note that now all three functions nsp, nsk and rcs are set using 'utils-helper-7'
#' This makes 'utils-helper-6' reduntant. BUT DON'T DELTE it yet
#'
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_rcs_call_stan_get <- function() {
  GS_rcs_call_stan_str <- 
    "
  //  functions { // comment out functions{
  
  /**
   * Calculates the Restricted Cubic Spline (RCS) basis matrix.
   *
   * This Stan function provides the core logic for generating a Restricted Cubic Spline
   * basis matrix, given an input vector `x` and pre-defined fullknots.
   *
   * IMPORTANT: Knot selection logic (e.g., quantiles, `Hmisc::rcspline.eval`)
   * must be performed OUTSIDE of Stan, in your data preparation script (R, Python, etc.).
   * The inclx is now defined within the functiin and not as an argument. 
   * The arguments are exact same as GS_nsp/nsk stan
   *
   * Args:
   * x: A vector of observations for which to calculate the basis.
   * fullknots: A vector of ordered unique fullknots. Must have at least 3 fullknots.
   * derivs: Integer representing the derivative order (0 for value, 1 for 1st derivs, 2 for 2nd derivs).
   * intercept: Integer (0 or 1). If 1, adds an intercept column.
   * inclx: Integer (0 or 1). If 1, includes the 'x' column (for derivs = 0).
   *
   * Returns:
   * A matrix where rows correspond to observations in `x` and columns
   * correspond to the basis functions.
   */
   
  matrix GS_rcs_call_stan(vector x, 
                        vector knotsx, 
                        vector bknotsx,
                        int intercept, 
                        int derivs, 
                        int centerval,
                        int normalize,
                        int preH) {
                        
    int inclx = 1;                  
    int N = num_elements(x);
    
    vector[num_elements(knotsx)+2] fullknots;
    vector[num_elements(fullknots)-2] knots;
    vector[2] bknots;
    fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
    knots = segment(fullknots, 2, num_elements(fullknots)-2);
    bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
    
    int nk = num_elements(fullknots);
    
    if (nk < 3) {
      reject(\"rcs_matrix: Number of fullknots must be at least 3.\");
    }
    
    matrix[N, nk - 1] basis_evals;
    real knot1 = fullknots[1];
    real knot_nk = fullknots[nk];
    real knot_nk_minus_1 = fullknots[nk - 1];
    real denom_common = (knot_nk - knot1)^2;
    real denom_spline_diff = (knot_nk - knot_nk_minus_1);

    matrix[N, nk] Xx_pos; 
    for (n_idx in 1:N) {
      for (k_idx in 1:nk) {
        Xx_pos[n_idx, k_idx] = fmax(0.0, x[n_idx] - fullknots[k_idx]);
      }
    }

    if (derivs == 0) {
      basis_evals[:, 1] = x;
    } else if (derivs == 1) {
      basis_evals[:, 1] = rep_vector(1.0, N);
    } else if (derivs == 2) {
      basis_evals[:, 1] = rep_vector(0.0, N);
    } else {
      reject(\"rcs_matrix: derivs must be 0, 1, or 2.\");
    }

    for (j in 1:(nk - 2)) {
      int jp1 = j + 1; 
      if (derivs == 0) {
        basis_evals[:, jp1] =
          (pow(Xx_pos[:, j], 3) -
           pow(Xx_pos[:, nk - 1], 3) * (fullknots[nk] - fullknots[j]) / denom_spline_diff +
           pow(Xx_pos[:, nk], 3) * (fullknots[nk - 1] - fullknots[j]) / denom_spline_diff) / denom_common;
      } else if (derivs == 1) {
        basis_evals[:, jp1] =
          (3.0 * pow(Xx_pos[:, j], 2)) / denom_common -
          (3.0 * pow(Xx_pos[:, nk - 1], 2)) * (fullknots[nk] - fullknots[j]) / (denom_spline_diff * denom_common) +
          (3.0 * pow(Xx_pos[:, nk], 2)) * (fullknots[nk - 1] - fullknots[j]) / (denom_spline_diff * denom_common);
      } else if (derivs == 2) {
        basis_evals[:, jp1] =
          (6.0 * Xx_pos[:, j]) / denom_common -
          (6.0 * Xx_pos[:, nk - 1]) * (fullknots[nk] - fullknots[j]) / (denom_spline_diff * denom_common) +
          (6.0 * Xx_pos[:, nk]) * (fullknots[nk - 1] - fullknots[j]) / (denom_spline_diff * denom_common);
      }
    }

    matrix[N, 0] final_basis_evals_empty = rep_matrix(0.0, N, 0); 

    if (intercept == 1) {
      if (derivs == 0) {
        if (inclx == 1) {
          matrix[N, nk] temp_basis;
          temp_basis[:, 1] = rep_vector(1.0, N);
          temp_basis[:, 2:(nk)] = basis_evals; 
          return temp_basis;
        } else {
          matrix[N, nk-1] temp_basis; 
          temp_basis[:, 1] = rep_vector(1.0, N);
          temp_basis[:, 2:(nk-1)] = basis_evals[:, 2:(nk-1)];
          return temp_basis;
        }
      } else if (derivs == 1) {
        if (inclx == 1) {
            matrix[N, nk] temp_basis;
            temp_basis[:, 1] = rep_vector(0.0, N); 
            temp_basis[:, 2:(nk)] = basis_evals; 
            return temp_basis;
        } else {
            matrix[N, nk-1] temp_basis; 
            temp_basis[:, 1] = rep_vector(0.0, N); 
            temp_basis[:, 2:(nk-1)] = basis_evals[:, 2:(nk-1)]; 
            return temp_basis;
        }
      } else if (derivs == 2) {
        if (inclx == 1) {
             matrix[N, nk] temp_basis; 
             temp_basis[:, 1] = rep_vector(0.0, N); 
             temp_basis[:, 2] = rep_vector(0.0, N); 
             print(temp_basis);
             print(basis_evals);
             temp_basis[:, 3:(nk - 0)] = basis_evals[:, 2:(nk - 1)]; 
             return temp_basis;
        } else {
             matrix[N, nk -1] temp_basis; 
             temp_basis[:, 1] = rep_vector(0.0, N); 
             temp_basis[:, 2] = rep_vector(0.0, N); 
             temp_basis[:, 3:(nk-1)] = basis_evals[:, 3:(nk - 1)]; 
             return temp_basis;
        }
      }
    } else { // intercept == 0
      if (inclx == 0) {
        return basis_evals[:, 2:(nk-1)]; 
      } else {
        return basis_evals;
      }
    }
    // Should not reach here
    return final_basis_evals_empty;
  }
  
// } // comment out functions{

  "
  return(GS_rcs_call_stan_str)
} # GS_rcs_call_stan_get










#' Create rcs spline design matrix
#'
#' @param x A numeric vector representing a predictor variable (e.g., age)
#' @param df An integer. It is defined as \code{nk - 1}
#' @param deriv An integer
#' @param add_intercept A logical (default \code{FALSE}) to indicate whether to
#'   add intercept column to the design matrix. This is useful when using
#'   \code{rcs_matrix} for creating design matrix for derivatives such as
#'   \code{deriv = 1} and \code{deriv = 2} where first (\code{deriv = 1}) or,
#'   the first and second (\code{deriv = 2}) columns are automatically set as
#'   \code{'0'}.
#' @param verbose A logical (default \code{FALSE}) to indicate if infornation
#'   need to be displayed.
#'
#' @inherit Hmisc::rcspline.eval params
#' 
#' @return An object of class \code{bgmfit} 
#' @keywords internal
#' @noRd
#'
rcs_matrix <- function(x, 
                       df, 
                       knots = NULL, 
                       deriv = 0,
                       add_intercept = FALSE,
                       inclx = TRUE, 
                       knots.only = FALSE,
                       type = "ordinary", 
                       norm = 2, 
                       rpm = NULL, 
                       pc = FALSE,
                       fractied = 0.05,
                       verbose = FALSE,
                       ...) {
  
  nk <- df + 1
  ##
  # borrow from Hmisc::rcspline.eval
  
  if (!length(knots)) {
    xx <- x[!is.na(x)]
    n <- length(xx)
    if (n < 6) 
      stop("knots not specified, and < 6 non-missing observations")
    if (nk < 3) 
      stop("nk must be >= 3")
    xu <- sort(unique(xx))
    nxu <- length(xu)
    if ((nxu - 2) <= nk) {
      warning(sprintf("%s knots requested with %s unique values of x.  knots set to %s interior values.", 
                      nk, nxu, nxu - 2))
      knots <- xu[-c(1, length(xu))]
    }
    else {
      outer <- if (nk > 3) 
        0.05
      else 0.1
      if (nk > 6) 
        outer <- 0.025
      knots <- numeric(nk)
      overrideFirst <- overrideLast <- FALSE
      nke <- nk
      firstknot <- lastknot <- numeric(0)
      if (fractied > 0 && fractied < 1) {
        f <- table(xx)/n
        if (max(f[-c(1, length(f))]) < fractied) {
          if (f[1] >= fractied) {
            firstknot <- min(xx[xx > min(xx)])
            xx <- xx[xx > firstknot]
            nke <- nke - 1
            overrideFirst <- TRUE
          }
          if (f[length(f)] >= fractied) {
            lastknot <- max(xx[xx < max(xx)])
            xx <- xx[xx < lastknot]
            nke <- nke - 1
            overrideLast <- TRUE
          }
        }
      }
      if (nke == 1) 
        knots <- median(xx)
      else {
        if (nxu <= nke) 
          knots <- xu
        else {
          p <- if (nke == 2) 
            seq(0.5, 1 - outer, length = nke)
          else seq(outer, 1 - outer, length = nke)
          knots <- quantile(xx, p)
          if (length(unique(knots)) < min(nke, 3)) {
            knots <- quantile(xx, seq(outer, 1 - outer, 
                                      length = 2 * nke))
            if (length(firstknot) && length(unique(knots)) < 
                3) {
              midval <- if (length(firstknot) && length(lastknot)) 
                (firstknot + lastknot)/2
              else median(xx)
              knots <- sort(c(firstknot, 
                              midval, if (length(lastknot)) lastknot else quantile(xx, 
                                                                                   1 - outer)))
            }
            if ((nu <- length(unique(knots))) < 3) {
              cat("Fewer than 3 unique knots.  Frequency table of variable:\n")
              print(table(x))
              stop()
            }
            warning(paste("could not obtain", nke, "interior knots with default algorithm.\n", 
                          "Used alternate algorithm to obtain", nu, 
                          "knots"))
          }
        }
        if (length(xx) < 100) {
          xx <- sort(xx)
          if (!overrideFirst) 
            knots[1] <- xx[5]
          if (!overrideLast) 
            knots[nke] <- xx[length(xx) - 4]
        }
      }
      knots <- c(firstknot, knots, lastknot)
    }
  }
  knots <- sort(unique(knots))
  nk <- length(knots)
  if (nk < 3) {
    cat("fewer than 3 unique knots.  Frequency table of variable:\n")
    print(table(x))
    stop()
  }
  if (knots.only) 
    return(knots)
  if (length(rpm)) 
    x[is.na(x)] <- rpm
  # end of borrow from Hmisc::rcspline.eval
  ##
  
  X <- x
  N <- length(X)
  nk <- length(knots)
  
  basis_evals <- matrix(0, N, nk-1)
  
  if(inclx) basis_evals[,1] = X
  
  if(deriv == 0) basis_evals[,1] = X;
  if(deriv == 1) basis_evals[,1] = 1;
  if(deriv == 2) basis_evals[,1] = 0;
  
  Xx <- matrix(0, N, nk)
  km1 = nk - 1;
  j = 1;
  knot1   <- knots[1     ]
  knotnk  <- knots[nk    ]
  knotnk1 <- knots[nk - 1]
  kd <-     (knotnk - knot1) ^ (2)
  
  for(ia in 1:N) {
    for(ja in 1:nk) {
      Xx[ia,ja] = ifelse(X[ia] - knots[ja] > 0, X[ia] - knots[ja], 0)
    }
  }
  
  if(deriv == 0) {
    while (j <= nk - 2) {
      jp1 = j + 1;
      basis_evals[,jp1] = 
        (Xx[,j]^3-(Xx[,km1]^3)*(knots[nk]-knots[j])/
           (knots[nk]-knots[km1]) + (Xx[,nk]^3)*(knots[km1]-knots[j])/
           (knots[nk]-knots[km1])) / (knots[nk]-knots[1])^2;
      j = j + 1;
    }
  }
  
  if(deriv == 1) {
    while (j <= nk - 2) {
      jp1 = j + 1;
      basis_evals[,jp1] =
        (3*Xx[,j]^2) * (1/((knots[nk]-knots[1])^2))  - 
        (3*Xx[,km1]^2)*(knots[nk]-knots[j]) /
        ((knots[nk]-knots[km1]) * (knots[nk]-knots[1])^2) + 
        (3*Xx[,nk]^2)*(knots[km1]-knots[j])/
        ((knots[nk]-knots[km1]) * (knots[nk]-knots[1])^2) ;
      j = j + 1;
    }
  }
  
  if(deriv == 2) {
    while (j <= nk - 2) {
      jp1 = j + 1;
      basis_evals[,jp1] =
        (6*Xx[,j]^1) * (1/((knots[nk]-knots[1])^2))  - 
        (6*Xx[,km1]^1)*(knots[nk]-knots[j]) /
        ((knots[nk]-knots[km1]) * (knots[nk]-knots[1])^2) + 
        (6*Xx[,nk]^1)*(knots[km1]-knots[j])/
        ((knots[nk]-knots[km1]) * (knots[nk]-knots[1])^2) ;
      j = j + 1;
    }
  }
  
  if(!inclx) basis_evals <- basis_evals[,-1,drop=FALSE]
  
  
  if(add_intercept) {
    if(deriv == 0) {
      mat_intercept <- matrix(1, nrow(basis_evals), 1)
      basis_evals <- cbind(mat_intercept, basis_evals)
      if(verbose) message("Intercept column added. Please use ~0 + formula")
    }
    if(deriv == 1) {
      mat_intercept <- matrix(0, nrow(basis_evals), 1)
      basis_evals <- cbind(mat_intercept, basis_evals)
      if(verbose) message("Intercept set to '0' for deriv = 1")
    }
    if(deriv == 2) {
      mat_intercept <- matrix(0, nrow(basis_evals), 2)
      basis_evals   <- basis_evals[, -1, drop = FALSE]
      basis_evals   <- cbind(mat_intercept, basis_evals)
      if(verbose) message("Intercept and first term (x) set to '0' for deriv = 2")
    }
  } # if(add_intercept) {
  
  
  return(basis_evals)
} # end rcs_matrix





