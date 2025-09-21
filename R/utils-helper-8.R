

###############################################################################
###############################################################################
######################## bsp, msp & isp FUNCTIONS #############################
###############################################################################
###############################################################################

# GS_bsp_tuple_intermediate_stan
# GS_bsp_intermediate_stan
# GS_bsp_call_stan
# 
# GS_isp_tuple_stan
# GS_isp_stan
# GS_isp_call_stan
# 
# GS_msp_tuple_stan
# GS_msp_stan
# GS_msp_call_stan





################################################################################
# GS_bsp_tuple_intermediate_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_tuple_intermediate_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  matrix GS_bsp_intermediate_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
 if (calcderiv) {
    if (intercept) {
    return bsderiv;
  } else {
    return bsderiv[, 2:cols(bsderiv)];
  }
 } else {
    if (intercept) {
    return bs;
  } else {
    return bs[, 2:cols(bs)];
  }
 }  
}
"
  return(set_stan_str)
} # GS_bsp_tuple_intermediate_stan_R_str



################################################################################
# GS_bsp_intermediate_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_intermediate_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  matrix GS_bsp_intermediate_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  } 
 if (calcderiv) {
    if (intercept) {
    return bsderiv;
  } else {
    return bsderiv[, 2:cols(bsderiv)];
  }
 } else {
    if (intercept) {
    return bs;
  } else {
    return bs[, 2:cols(bs)];
  }
 }  
}
"
  return(set_stan_str)
} # GS_bsp_intermediate_stan_R_str



################################################################################
# GS_bsp_call_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_bsp_call_stan_R_str <- function() {
  set_stan_str <- "
   /////////////////////////////////////////////////////////////////////////
  // call GS_bsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_bsp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  int degree_in = degree + 0;
  int ord = degree_in + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;  
  // int ncolselect = Nk + intercept - 1;  
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;  
  int ncolselect = Nintk + degree_in + 1 + intercept - 1;  
  // generate basis of degree 0
  if(degree == 0) {
      matrix[N, Nintervals] out_degree0 = rep_matrix(0, N, Nintervals);
      if (calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            out_degree0[i, j] = 0.0;
          }
        }
      }
      if (!calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            real denom =  allknots[j+1] - allknots[j];
          //  real inv_denom = 1 / denom;
            real inv_denom = 1;
            // for last cloumn to match splines::msp, need this condition which 
            // differe by x[i] <= allknots[j + 1]
            if(j < Nintervals) {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] < allknots[j + 1]) ? inv_denom : 0;
            } else {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] <= allknots[j + 1]) ? inv_denom : 0;
            }
          }
        }
      }
     
     if(intercept) {
        return out_degree0[, 1:cols(out_degree0)];
     } else {
        return out_degree0[, 2:cols(out_degree0)];
     } 
  } // end if(degree == 0) {
  
  matrix[N, ncolselect] out = GS_bsp_intermediate_stan(x, knots, bknots, fullknots, allknots, N, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_bsp_intermediate_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots, allknots, 1, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }

  return out;
}     
"
  return(set_stan_str)
} # GS_bsp_call_stan_R_str




################################################################################
# GS_isp_tuple_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_isp_tuple_stan_R_str <- function() {
  set_stan_str <- "
/////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
 // Main I-spline function
matrix GS_isp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
  int Psize = degree + 1 +  size(knots);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  // Calculate rowsum of consecutive bs columns  
  if (calcderiv) {
   for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bsderiv[r, (i+0):(Psize+0)]);
     }
    }
  } else {
    for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bs[r, (i+0):(Psize+0)]);
     }
    }
  }
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_isp_tuple_stan_R_str





################################################################################
# GS_isp_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_isp_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
 // Main I-spline function
matrix GS_isp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  }
  int Psize = degree + 1 +  size(knots);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  // Calculate rowsum of consecutive bs columns  
  if (calcderiv) {
   for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bsderiv[r, (i+0):(Psize+0)]);
     }
    }
  } else {
    for (i in 1:Psize) {
     for (r in 1:N) {
      M[r, i] = sum(bs[r, (i+0):(Psize+0)]);
     }
    }
  }
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_isp_stan_R_str




################################################################################
# GS_isp_call_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_isp_call_stan_R_str <- function() {
  set_stan_str <- "
   /////////////////////////////////////////////////////////////////////////
  // call GS_bsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_isp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  // here in _call_stan, degree is renamed as degree_in because int degree part of fun
  int degree_in = degree + 1;
  int ord = degree_in + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;
  
  // int ncolselect = Nk + intercept - 1;
  int Psize = degree_in + 1 + size(knots) + intercept - 1;
  int ncolselect = Psize;
  
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;
  matrix[N, ncolselect] out = GS_isp_stan(x, knots, bknots, fullknots, allknots, N, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_isp_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots, allknots, 1, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }  
  return out[, 2:cols(out)];
}
"
  return(set_stan_str)
} # GS_isp_call_stan_R_str



################################################################################
# GS_msp_tuple_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_msp_tuple_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  // Main M-spline function
matrix GS_msp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  tuple(matrix[N, n_basis], matrix[N, n_basis]) my_tuple_main = 
  GS_bsp_tuple_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bs = my_tuple_main.1;
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = my_tuple_main.2;
  }
  int Psize = cols(bs);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  if (calcderiv) {
   for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bsderiv[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  } else {
    for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bs[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  }        
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_msp_tuple_stan_R_str



################################################################################
# GS_msp_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_msp_stan_R_str <- function() {
  set_stan_str <- "
 /////////////////////////////////////////////////////////////////////////  
  // Function to compute B-splines and their derivatives
  ///////////////////////////////////////////////////////////////////////// 
  // Main M-spline function
matrix GS_msp_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bsp_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  } 
  int Psize = cols(bs);
  matrix[N, Psize] M = rep_matrix(0.0, N, Psize);
  if (calcderiv) {
   for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bsderiv[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  } else {
    for (i in 1:Psize) {
      M[, i] = (degree+1 ) * bs[,i] / (fullknots[i+degree+1] - fullknots[i]);
    } 
  }        
  // Avoid unnecessary allocation
  if (intercept) {
    return M;
  } else {
    return M[, 2:cols(M)];
  }
}
"
  return(set_stan_str)
} # GS_msp_stan_R_str



################################################################################
# GS_msp_call_stan_R_str
################################################################################

#' An internal function 
#' 
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
GS_msp_call_stan_R_str <- function() {
  set_stan_str <- "
   /////////////////////////////////////////////////////////////////////////
  // call GS_bsp_call_stan
  /////////////////////////////////////////////////////////////////////////
  matrix GS_msp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int degree, int intercept, int derivs, real centerval, int normalize,
                        int preH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  // here in _call_stan, degree is renamed as degree_in because int degree part of fun
  int degree_in = degree + 0;
  int ord = degree_in + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;
  
  vector[2*(degree+1) + Nintk] fullknots_mspl;
  fullknots_mspl = append_row(append_row(rep_vector(bknotsx[1], (degree+1)), knotsx), rep_vector(bknotsx[2], (degree+1)));
  
  // int ncolselect = Nk + intercept - 1;
  int Psize = degree_in + 1 + size(knots) + intercept - 1;
  int ncolselect = Psize;
  
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;
  
  // generate basis of degree 0
  if(degree == 0) {
      matrix[N, Nintervals] out_degree0 = rep_matrix(0, N, Nintervals);
      if (calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            out_degree0[i, j] = 0.0;
          }
        }
      }
      if (!calcderiv) {
        for (i in 1:N) {
          for (j in 1:Nintervals) {
            real denom =  allknots[j+1] - allknots[j];
            // for last cloumn to match splines::msp, need this condition which 
            // differe by x[i] <= allknots[j + 1]
            if(j < Nintervals) {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] < allknots[j + 1]) ? 1/denom : 0;
            } else {
                out_degree0[i, j] = (x[i] >= allknots[j] && x[i] <= allknots[j + 1]) ? 1/denom : 0;
            }
          }
        }
      }
     
     if(intercept) {
        return out_degree0[, 1:cols(out_degree0)];
     } else {
        return out_degree0[, 2:cols(out_degree0)];
     } 
  } // end if(degree == 0) {
  
  
  matrix[N, ncolselect] out = GS_msp_stan(x, knots, bknots, fullknots_mspl, allknots, N, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_msp_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots_mspl, allknots, 1, degree_in, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }
  
  return out;
} 
"
  return(set_stan_str)
} # GS_msp_call_stan_R_str




###########################################################################
###########################################################################
###########################################################################



# Don't know why it is hanging here 

# get_GS_nsk_call_stan <- function() {
#   str <- 
#     "
#     
#     int VecSame(vector A, vector B) {
#   int N = num_elements(A);
#   for (i in 1:N) {
#     if (A[i] != B[i]) return 0;
#   }
#   return 1;
# }
# 
# int MatSame(matrix A, matrix B) {
#   int R = rows(A);
#   int C = cols(A);
#   for (i in 1:R) {
#     for (j in 1:C) {
#       if (A[i, j] != B[i, j]) return 0;
#     }
#   }
#   return 1;
# }
# 
# /*
# vector repeat_vector(vector input, int K) {
#   int N = num_elements(input);
#   vector[N*K] repvec;
#   for (k in 1:K)
#     for (i in 1:N)
#       repvec[i + (k - 1) * N] = input[i];
#   return repvec;
# }
# */
# 
# /*
# matrix repeat_matrix(matrix input, int K) {
#   int N = rows(input);
#   int M = cols(input);
#   matrix[N*K, M] repmat;
#   for (k in 1:K)
#     for (i in 1:N)
#       repmat[i + (k - 1) * N] = input[i];
#   return repmat;
# }
# */
# 
#     
#     /////////////////////////////////////////////////////////////////////////
#   // Function to calculate H matrix (as defined in GS_ns_getH_stan)
#   /////////////////////////////////////////////////////////////////////////
#   matrix GS_ns_getH_stan(vector knots, int normalize) {
#     int Nintk = num_elements(knots) - 8;
#     real C11 = 6 / ((knots[5] - knots[2]) * (knots[5] - knots[3]));
#     real C31 = 6 / ((knots[6] - knots[3]) * (knots[5] - knots[3]));
#     real C21 = -C11 - C31;
#     real Cp22 = 6 / ((knots[Nintk + 6] - knots[Nintk + 3]) * (knots[Nintk + 6] - knots[Nintk + 4]));
#     real Cp2 = 6 / ((knots[Nintk + 7] - knots[Nintk + 4]) * (knots[Nintk + 6] - knots[Nintk + 4]));
#     real Cp12 = -Cp22 - Cp2;
#     matrix[Nintk+1, Nintk+1] dH; // dummy
#     if (Nintk == 0) {
#       matrix[2, 4] H;
#       H = to_matrix(transpose([3, 0, 2, 1, 1, 2, 0, 3]), 2, 4);
#       if (normalize) {
#         vector[2] sumH = H * rep_vector(1.0, 4);
#         for (i in 1:2) {
#           if (sumH[i] != 0) {
#             H[i, ] = H[i, ] / sumH[i];
#           }
#         }
#       } // if (normalize) 
#       return transpose(H);
#     } else if (Nintk == 1) {
#       matrix[3, 5] H;
#       H = transpose(to_matrix([
#         -C21 / C11, 1, 0, 0, 0, 
#         0, -C31 / C21, 1, -Cp22 / Cp12, 0, 
#         0, 0, 0, 1, -Cp12 / Cp2
#       ], 5, 3));
#       if (normalize) {
#         vector[Nintk + 2] sumH = H * rep_vector(1.0, Nintk + 2 + 2);
#         for (i in 1:(Nintk + 2)) {
#           if (sumH[i] != 0) {
#             H[i, ] = H[i, ] / sumH[i];
#           }
#         }
#       } // if (normalize) 
#       return transpose(H);
#     } else if (Nintk == 2) {
#       matrix[Nintk + 2, 3] H1;
#       matrix[Nintk + 2, 3] H3;
#       matrix[Nintk + 2, Nintk + 2 + 2] H;
#       H1 = append_row(
#             append_row(
#               append_row(
#                 rep_matrix(1, 1, 3), [0, 1, -C21 / C31]
#                 ),
#                 rep_matrix(0,  Nintk - 2, 3)
#               ),
#             rep_matrix(0,  2, 3)
#           );
#       H3 = append_row(
#             append_row(
#               append_row(
#                 rep_matrix(0, 2, 3), rep_matrix(0, Nintk - 2, 3)
#                 ),
#                 [-Cp12 / Cp22, 1, 0]
#               ),
#             rep_matrix(1,  1, 3)
#           );
#       H = append_col(H1, H3);
#       if (normalize) {
#         vector[Nintk + 2] sumH = H * rep_vector(1.0, Nintk + 2 + 2);
#         for (i in 1:(Nintk + 2)) {
#           if (sumH[i] != 0) {
#             H[i, ] = H[i, ] / sumH[i];
#           }
#         }
#       } // if (normalize) 
#       return transpose(H);
#     } else if (Nintk > 2) {
#       matrix[Nintk + 2, 3] H1;
#       matrix[Nintk + 2, 3] H3;
#       matrix[Nintk + 2, Nintk - 2] H2;
#       matrix[Nintk + 2, Nintk + 2 + 2] H;
#       H1 = append_row(
#             append_row(
#               append_row(
#                 rep_matrix(1, 1, 3), [0, 1, -C21 / C31]
#                 ),
#                 rep_matrix(0,  Nintk - 2, 3)
#               ),
#             rep_matrix(0,  2, 3)
#           );
#       H2 = append_row(
#             append_row(
#               rep_matrix(0, 2, Nintk - 2), diag_matrix(rep_vector(1, Nintk - 2))
#               ),
#               rep_matrix(0, 2, Nintk - 2)
#             );
#       H3 = append_row(
#             append_row(
#               append_row(
#                 rep_matrix(0, 2, 3), rep_matrix(0, Nintk - 2, 3)
#                 ),
#                 [-Cp12 / Cp22, 1, 0]
#               ),
#             rep_matrix(1,  1, 3)
#           );
#       H = append_col(append_col(H1, H2), H3);
#       if (normalize) {
#         vector[Nintk + 2] sumH = H * rep_vector(1.0, Nintk + 2 + 2);
#         for (i in 1:(Nintk + 2)) {
#           if (sumH[i] != 0) {
#             H[i, ] = H[i, ] / sumH[i];
#           }
#         }
#       } // if (normalize) 
#       return transpose(H);
#     } 
#     return dH;
#   } // end matrix GS_ns_getH_stan
#  /////////////////////////////////////////////////////////////////////////
#   // Find elements in vector x that equal/below/above real number y
#  /////////////////////////////////////////////////////////////////////////
#   int num_matches(vector x, real y, real z) {
#     int n = 0;
#     for (i in 1:rows(x))
#       if(z == 0) {
#         if (x[i] == y) n += 1;
#       } else if(z == 1.0) {
#         if (x[i] > y) n += 1;
#       } else if(z == -1.0) {
#         if (x[i] < y) n += 1;
#       }
#     return n;
#   }
#   /////////////////////////////////////////////////////////////////////////    
#   // Find the indexes of the elements in the vector x that equal real number y
#   /////////////////////////////////////////////////////////////////////////
#   array[] int which_equal(vector x, real y, real z) {
#     array[num_matches(x, y, z)] int match_positions;
#     int pos = 1;
#     for (i in 1:size(x)) {
#       if(z == 0) {
#         if (x[i] == y) {
#           match_positions[pos] = i;
#           pos += 1;
#         }
#       } else if(z == 1.0) {
#         if (x[i] > y) {
#           match_positions[pos] = i;
#           pos += 1;
#         }
#       } else if(z == -1.0) {
#         if (x[i] < y) {
#           match_positions[pos] = i;
#           pos += 1;
#         }
#       }
#     }
#     return match_positions;
#   }
#   /////////////////////////////////////////////////////////////////////////
#   // Repeating input vector K times
#   /////////////////////////////////////////////////////////////////////////
#   vector repeat_vector(vector input, int K) {
#     int N = rows(input);
#     vector[N*K] repvec; // stack N-vector K times
#     // assign i-th value of input to i+(k-1)*N -th value of repvec
#     for (k in 1:K) {
#       for (i in 1:N) {
#         repvec[i+(k-1)*N] = input[i]; 
#       }
#     }
#     return repvec;
#   }
#   /////////////////////////////////////////////////////////////////////////
#   // Repeating input matrix K times
#   /////////////////////////////////////////////////////////////////////////
#   matrix repeat_matrix(matrix input, int K) {
#     int N = rows(input);
#     int M = cols(input);
#     matrix[N*K,M] repmat; // stack N*M matrix K times
#     // assign i-th row of input to i+(k-1)*N -th row of repmat
#     for (k in 1:K) {
#       for (i in 1:N) {
#         repmat[i+(k-1)*N] = input[i]; 
#       }
#     }
#     return repmat;
#   }  
#   /////////////////////////////////////////////////////////////////////////  
#   // Function to compute B-splines and their derivatives
#   ///////////////////////////////////////////////////////////////////////// 
#   matrix GS_bs_stan(vector x, vector knots, vector bknots, 
#                     vector fullknots, vector allknots, 
#                     int N, int degree, int ord, 
#                     int Nintk, int Nk, int Nintervals,
#                     int intercept, int calcderiv) {    
#     matrix[N, Nintervals] M1 = rep_matrix(0, N, Nintervals);    
#     for(i in 1:N) {
#       for (j in 1:Nintervals) {
#         M1[i,j] = (allknots[j] <= x[i]  && x[i] < allknots[j + 1] ? 1 : 0);
#       }
#     }
#     array[num_matches(x, bknots[2], 0)] int lastknot_index = which_equal(x, bknots[2], 0); 
#     if (size(lastknot_index) > 0) {
#       for(i in lastknot_index) {
#         M1[i, (Nintervals - degree):Nintervals] = rep_row_vector(1, (Nintervals-(Nintervals - degree))+1);
#       }
#     }
#     matrix[N, Nintervals - 0] M2x = rep_matrix(0, N, Nintervals);
#     matrix[N, Nintervals - 0] M1x = M1;
#     vector[N] C1;
#     vector[N] C2;
#     for (p in 1:degree) {
#       for (i in 1:(Nintervals - p)) {
#         if (allknots[i + p] == allknots[i]) {
#           C1 = rep_vector(0.0, N);
#         } else {
#           C1 = (x - allknots[i]) / (allknots[i + p] - allknots[i]);
#         }
#         if (allknots[i + p + 1] == allknots[i + 1]) {
#           C2 = rep_vector(0.0, N);
#         } else {
#           C2 = (allknots[i + p + 1] - x) / (allknots[i + p + 1] - allknots[i + 1]);
#         }
#         M2x[, i] = C1 .* M1x[, i] + C2 .* M1x[, i + 1];
#       }
#       if(p != degree) M1x = M2x;
#     }
#     matrix[N, Nintervals - degree]    M2 = M2x[, 1:(Nintervals - degree)];
#     matrix[N, (Nintervals - degree)+1] M1xx = M1x[, 1:(Nintervals - degree)+1];
#     if (calcderiv) {
#       matrix[N, Nintervals - degree] deriv;
#       for (i in 1:(Nintervals - degree)) {
#         real dC1 = (allknots[i + degree] == allknots[i]) ? 0 : degree / (allknots[i + degree] - allknots[i]);
#         real dC2 = (allknots[i + degree + 1] == allknots[i + 1]) ? 0 : degree / (allknots[i + degree + 1] - allknots[i + 1]);
#         deriv[, i] = dC1 * M1xx[, i] - dC2 * M1xx[, i + 1];
#       }
#       return deriv;
#     } else {
#       return M2;
#     }
#   } // end matrix GS_bs_stan  
#   /////////////////////////////////////////////////////////////////////////
#   // Main function for B-splines and their derivatives
#   /////////////////////////////////////////////////////////////////////////
#   matrix GS_ns_stan(vector x, vector knots, vector bknots, 
#                     vector fullknots, vector allknots, 
#                     int N, int degree, int ord, 
#                     int Nintk, int Nk, int Nintervals,
#                     int intercept, int calcderiv, int normalize,
#                     int preH) {                    
#     matrix[num_elements(x), Nintervals - degree] bs;
#     matrix[num_elements(x), Nintervals - degree] bsderiv;
#     if (calcderiv) {
#       bs      = GS_bs_stan(x, knots, bknots, 
#                           fullknots, allknots,  
#                           N, degree, ord,
#                           Nintk, Nk, Nintervals,
#                           intercept, 0);
#       bsderiv = GS_bs_stan(x, knots, bknots, 
#                           fullknots, allknots,  
#                           N, degree, ord,
#                           Nintk, Nk, Nintervals,
#                           intercept, 1);
#     } else {
#       bs      = GS_bs_stan(x, knots, bknots, 
#                           fullknots, allknots,
#                           N, degree, ord,
#                           Nintk, Nk, Nintervals,
#                           intercept, 0);
#     }
#     array[num_matches(x, bknots[1], -1)] int x_below_boundary  = which_equal(x, bknots[1], -1); 
#     array[num_matches(x, bknots[2],  1)] int x_above_boundary  = which_equal(x, bknots[2],  1); 
#     int size_x_below_boundary = size(x_below_boundary);
#     int size_x_above_boundary = size(x_above_boundary);
#     if (size_x_below_boundary > 0 || size_x_above_boundary > 0) {
#       matrix[2, Nintervals - degree] bs_bknots;
#       matrix[2, Nintervals - degree] bsderiv_bknots;
#       bs_bknots      = GS_bs_stan(bknots, knots, bknots, 
#                           fullknots, allknots, 
#                           2, degree, ord,
#                           Nintk, Nk, Nintervals,
#                           intercept, 0);
#       bsderiv_bknots = GS_bs_stan(bknots, knots, bknots, 
#                           fullknots, allknots,
#                           2, degree, ord,
#                           Nintk, Nk, Nintervals,
#                           intercept, 1);      
#     if (size_x_below_boundary > 0) {
#       array[size_x_below_boundary] int xselect = x_below_boundary;
#       int Nxselect = size_x_below_boundary;
#       matrix[size_x_below_boundary, (Nintervals - degree)] below_azx1;
#       matrix[size_x_below_boundary, (Nintervals - degree)] below_azx2;
#       matrix[size_x_below_boundary, (Nintervals - degree)] below_azx3;      
#       below_azx1 = transpose(to_matrix((repeat_vector(transpose(bs_bknots[1,]), 
#                     size_x_below_boundary)), Nintervals - degree, 
#                     size_x_below_boundary));
#       below_azx2 = transpose(to_matrix((repeat_vector(transpose(bsderiv_bknots[1,]), 
#                     size_x_below_boundary)), Nintervals - degree, 
#                     size_x_below_boundary));
#       below_azx3 = (to_matrix((repeat_vector(((x[xselect] - bknots[1])), 
#                     Nintervals - degree)), size_x_below_boundary, 
#                     Nintervals - degree));
#         bs[xselect, ] = below_azx1 + (below_azx2 .* below_azx3);
#         if (calcderiv) {
#           for(i in xselect) {
#             bsderiv[i, 1: Nintervals - degree] = (bsderiv_bknots[1,]);
#           }
#         }
#       } // if (size_x_below_boundary > 0)      
#       if (size_x_above_boundary > 0) {
#         array[size_x_above_boundary] int xselect = x_above_boundary;
#         int Nxselect = size_x_above_boundary; 
#         matrix[size_x_above_boundary, (Nintervals - degree)] below_azx1;
#         matrix[size_x_above_boundary, (Nintervals - degree)] below_azx2;
#         matrix[size_x_above_boundary, (Nintervals - degree)] below_azx3;
#         below_azx1 = transpose(to_matrix((repeat_vector(transpose(bs_bknots[2,]), 
#                                                         size_x_above_boundary)), Nintervals - degree, 
#                                          size_x_above_boundary));
#         below_azx2 = transpose(to_matrix((repeat_vector(transpose(bsderiv_bknots[2,]), 
#                                                         size_x_above_boundary)), Nintervals - degree, 
#                                          size_x_above_boundary));
#         below_azx3 = (to_matrix((repeat_vector(((x[xselect] - bknots[2])), 
#                                                Nintervals - degree)), size_x_above_boundary, 
#                                 Nintervals - degree));        
#         bs[xselect, ] = below_azx1 + (below_azx2 .* below_azx3);        
#         if (calcderiv) {
#             for(i in xselect) {
#               bsderiv[i,1: Nintervals - degree ] = (bsderiv_bknots[2,]);
#            }         
#         }
#       } // if (size_x_above_boundary > 0)    
#     } // if (size_x_below_boundary > 0 || size_x_above_boundary > 0) {
#     matrix[Nk+2, Nk] H;
#     if(preH) {
#    //  H = GS_ns_getH_pre(Nk+2, Nk);
#     } else {
#      H = GS_ns_getH_stan(allknots, normalize);
#     }
#     matrix[num_elements(x), Nk] out;
#     int ncolselect = Nk+intercept-1;
#     matrix[num_elements(x), ncolselect] result;
#     if (calcderiv) {
#      out = (bsderiv * H);
#       if(intercept) {
#        result = out;
#       } else {
#        result = out[, 2:cols(out)];
#       }
#     } else {
#      out = bs * H;
#      if(intercept) {
#        result = out;
#       } else {
#        result = out[, 2:cols(out)];
#       }
#     }
#     return result;
#   } // end matrix GS_ns_stan  
#   /////////////////////////////////////////////////////////////////////////
#   // call GS_nsp_call_stan
#   /////////////////////////////////////////////////////////////////////////
#   matrix GS_nsp_call_stan(vector x, vector knotsx, vector bknotsx, 
#                         int degree,int intercept, int derivs, real centerval, int normalize,
#                         int preH) {
#   int N        = num_elements(x);
#   vector[num_elements(knotsx)+2] fullknots;
#   vector[num_elements(fullknots)-2] knots;
#   vector[2] bknots;
#   fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
#   knots = segment(fullknots, 2, num_elements(fullknots)-2);
#   bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
#    int Nintk   = num_elements(knots);
#   // int degree  = 3;
#    int ord     = degree + 1;
#    int Nk      = Nintk + 2;
#    int df      = Nintk + 1 + intercept;
#    int calcderiv;
#    if(derivs > 1) {
#     } else {
#       calcderiv = derivs;
#    }
#    if(intercept) {
#     if(centerval != 0) {
#       // centerval = 0;
#     }
#   }
#   if(Nintk > 0) {
#     if(bknots[1] > knots[1]) {
#       } else if(bknots[2] < knots[Nintk]) {
#     }
#   }
#   int ncolselect = Nk+intercept-1; 
#   matrix[N, ncolselect] out;
#   vector[Nintk + 2*ord] allknots;
#    allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), 
#                           rep_vector(bknots[2], ord));    
#    int Nintervals  = num_elements(allknots) - 1;
#     out = GS_ns_stan(x, knots, bknots, fullknots, allknots, 
#                         N, degree, ord, 
#                         Nintk, Nk, Nintervals, 
#                         intercept, calcderiv, normalize, preH);                        
#     if (centerval != 0) {
#     matrix[1, df] cenout;
#     cenout = GS_ns_stan(centerval+rep_vector(0.0, 1), knots, bknots, fullknots, allknots,
#                         1, degree, ord, 
#                         Nintk, Nk, Nintervals, 
#                         intercept, calcderiv, normalize, preH);      
#       if (!calcderiv) {
#         if(intercept) {
#           for (i in 2:cols(cenout)) {
#             out[,i] = out[,i] - cenout[1,i];
#           }
#         } else if(!intercept) {
#           for (i in 1:cols(cenout)) {
#             out[,i] = out[,i] - cenout[1,i];
#           }
#         }
#       } 
#     } // if (centerval != 0) {
#   return(out);
#   } // end matrix GS_nsp_call_stan  
#   /////////////////////////////////////////////////////////////////////////
#   // call GS_nsk_call_stan
#   /////////////////////////////////////////////////////////////////////////
#   matrix GS_nsk_call_stan(vector x, vector knotsx, vector bknotsx, 
#                         int degree, int intercept, int derivs, real centerval, int normalize,
#                         int preH) {       
#   int N        = num_elements(x);
#   vector[num_elements(knotsx)+2] fullknots;
#   vector[num_elements(fullknots)-2] knots;
#   vector[2] bknots;
#   fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
#   knots = segment(fullknots, 2, num_elements(fullknots)-2);
#   bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
#    int Nintk   = num_elements(knots);
#  //  int degree  = 3;
#    int ord     = degree + 1;
#    int Nk      = Nintk + 2;
#    int df      = Nintk + 1 + intercept;
#    int calcderiv;
#    if(derivs > 1) {
#     } else {
#       calcderiv = derivs;
#    }
#    if(intercept) {
#     if(centerval != 0) {
#       // centerval = 0;
#     }
#   }
#   if(Nintk > 0) {
#     if(bknots[1] > knots[1]) {
#       } else if(bknots[2] < knots[Nintk]) {
#     }
#   }  
#   int ncolselect = Nk+intercept-1;
#   matrix[N, ncolselect] basis;
#   matrix[Nk, ncolselect] kbasis;
#   matrix[N, ncolselect+1] kout;
#   matrix[N, ncolselect] out;
#    if(!calcderiv) {
#      basis  = GS_nsp_call_stan(x, knots, bknots, degree, intercept, derivs, centerval, normalize, preH);
#      kbasis = GS_nsp_call_stan(fullknots, knots, bknots, degree,intercept, derivs, centerval, normalize, preH);
#    } else {
#      basis  = GS_nsp_call_stan(x, knots, bknots, degree,intercept, derivs, centerval, normalize, preH);
#      kbasis = GS_nsp_call_stan(fullknots, knots, bknots, degree, intercept, 0, centerval, normalize, preH);
#    }
#     if(intercept) {
#        out = basis * inverse(kbasis);
#     } else {
#        kout = append_col(rep_vector(1, N), basis) * 
#               inverse(append_col(rep_vector(1, Nk), kbasis));
#        out = kout[, 2:cols(kout)];
#     } 
#     return(out);
#   } // end matrix GS_nsk_call_stan
# 
#   
#   "
#   return(str)
# }






