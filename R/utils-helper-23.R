
#############################################################
############### GS_gps_parms_R ##################
#############################################################

GS_gps_parms_R <- function(nlp_a,
                           nlp_b,
                           nlp_c,
                           nlp_d,
                           SParMat,
                           xknots,
                           spline_eval_array,
                           xg_array,
                           xg_curve_array,
                           degree,
                           shift_indicator,
                           spline_subset_indicator,
                           spline_precomputed_indicator,
                           set_spread,
                           return_indicator,
                           drawni) {
  
  pieces_dim      <- (length(xknots)-1)
  degree_dim      <- (degree + 1L)
  deriv_root      <- 2
  solvex          <- 0
  PiecePolyCoef   <- matrix(NA, degree_dim, pieces_dim)
  rrootsmat       <- matrix(NA, 1, pieces_dim)
  
  piece_id        <- NULL
  mark_peak       <- TRUE
  
  parm_deriv      <- c(0,1)
  curve_deriv     <- c(0,1,2)
  
  n_rowni         <- length(nlp_a)
  
  degree_dim_set_spread <- degree_dim * set_spread
  
  if(0 %in%  return_indicator) {
    coef_array  <- array(NA, dim = c(degree_dim, pieces_dim, n_rowni))
  }
  
  if(1 %in%  return_indicator) {
    parm_mat_dim <- 6
    if(!is.null(drawni)) {
      parm_mat_dim <- parm_mat_dim + 1
    }
    parm_array  <- array(NA, dim = c(pieces_dim, parm_mat_dim, n_rowni))
  }
  
  
  if(2 %in%  return_indicator) {
    deriv_array_dim_1 <- pieces_dim * degree_dim * set_spread
    deriv_array_dim_2 <- (length(curve_deriv) + 1+1+1)
    if(!is.null(drawni)) {
      deriv_array_dim_2 <- deriv_array_dim_2 + 1
    }
    deriv_array <- array(NA, dim = c(deriv_array_dim_1, deriv_array_dim_2, n_rowni))
  }
  
  
  ############################################################
  # Start - over 1:n_rowni
  ############################################################
  for (rowni in 1:n_rowni) {
    if(spline_subset_indicator) {
      smat <- SParMat[1, ]
    } else {
      smat <- SParMat[rowni, ]
    }
    par_a <- nlp_a[rowni]
    par_b <- nlp_b[rowni]
    par_c <- nlp_c[rowni]
    par_d <- nlp_d[rowni]
    
    ############################################################
    # Start - coefficients and roots
    ############################################################
    i <- 1L
    while (i <= pieces_dim) {
      if(spline_precomputed_indicator == 0) {
        # xg <- seq.int(xknots[i], xknots[i + 1L], length.out = degree_dim)
        xg = seq_fun_R(xknots[i], xknots[i+1], degree_dim);
      }
      if(spline_precomputed_indicator == 1) {
        xg = xg_array[,i];
      }
      
      if(spline_precomputed_indicator == 0) {
        # eval(SplineCall) %*% smat
        SplineCall[[2]] <- quote(xg)
        eval_Spline = eval(SplineCall)
      }
      if(spline_precomputed_indicator == 1) {
        eval_Spline = spline_eval_array[,,i];
      }
      
      yg <- eval_Spline %*% smat
      # Imp, adjust for b and  c after eval(SplineCall) %*% smat
      xg <- ((xg/exp(par_c)) + par_b) - shift_indicator * xknots[i]
      
      
      Xg <- outer(xg, 0:degree, "^")
      A  <- base::crossprod(Xg)
      b  <- base::crossprod(Xg, yg)
      U  <- chol.default(A)
      PiecePolyCoef[,i] <- pc <- base::backsolve(U, base::forwardsolve(t.default(U), b))
      if (deriv_root > 0) {
        pc.i <- pc[-seq_len(deriv_root)] * choose(deriv_root:degree, deriv_root) * 
          factorial(deriv_root)
      }
      pc.i[1] <- pc.i[1] - solvex
      
      # croots <- base::polyroot(pc.i)
      # rroots <- Re(croots)[round(Im(croots), 10) == 0]
      
      rroots <- - pc.i[1] / pc.i[2]
      
      if (shift_indicator) {
        rroots <- rroots + xknots[i]
      }
      xi_  <- (xknots[i]   / exp(par_c)) + par_b
      xi_1 <- (xknots[i+1] / exp(par_c)) + par_b
      
      
      
      if((rroots >= xi_ ) & (rroots <= xi_1 )) {
        get_rroots <- rroots # [(rroots >= xi_) & (rroots <= xi_1 )]
      } else {
        get_rroots <- NA
      }
      
      # Also replace first and last piece to NA, boundry
      if(i == 1 | i == pieces_dim) {
        get_rroots <- NA
      }
      
      # if(rowni == 12 | rowni == 12) {
      #   print(xi_)
      #   print(xi_1)
      #   print(rroots)
      #   print(get_rroots)
      # }
      
      rrootsmat[,i] <- get_rroots
      i <- i + 1L
    } # while (i <= pieces_dim) {
    
    if(0 %in%  return_indicator) {
      coef_array[,,rowni]   <- PiecePolyCoef
    }
    
    
    
    #  [[-3.688,inf,inf,inf,-0.728897,inf,1.12573,2.53888,4.602]] 
    #  NA -2.585848 -1.888408 -1.616958 -0.3704589   NA   NA   NA   NA
    #  [[inf,-2.58603,-1.88834,-1.61689,-0.370454,inf,inf,inf,4.67975]] 
    ############################################################
    # End - coefficients and roots
    ############################################################
    
    
    
    ############################################################
    # Start - parameters
    ############################################################
    newx       <- rrootsmat %>% as.vector()
    if(is.null(piece_id)) {
      set_piece_id <- seq(1, pieces_dim)
      set_piece_id <- set_piece_id[!is.na(newx)]
      newx     <- newx[!is.na(newx)]
    } else {
      set_piece_id <- piece_id
      if(length(set_piece_id) != length(newx))
        stop("lengths of 'piece_id' and 'newx' must match")
    }
    ind <- split.default(seq_len(length(newx)), set_piece_id)
    unique_piece_id <- as.integer(names(ind))
    n_pieces <- length(unique_piece_id)
    
    parm_mat_dim  <- parm_deriv + 1+1+1+1
    if(mark_peak) {
      parm_mat_dim <- parm_mat_dim + 1
      if(!1 %in% parm_deriv) stop("please set parm_deriv = 1")
    }
    if(!is.null(rowni))  parm_mat_dim <- parm_mat_dim + 1
    if(!is.null(drawni)) parm_mat_dim <- parm_mat_dim + 1
    parm_mat_dim_nrows <- pieces_dim
    parm_mat           <- matrix(NA, parm_mat_dim_nrows, parm_mat_dim)
    
    i <- 1L
    while (i <= n_pieces) {
      ii      <- unique_piece_id[i]
      xg_newx <- newx[ind[[i]]]
      xg      <- xg_newx - shift_indicator * xknots[ii]
      pc      <- PiecePolyCoef[, ii]
      collectd_mat_i <- 0
      for (dxi in parm_deriv) {
        if(!is.null(parm_deriv)) {
          if (dxi >= degree) {
            stop("'parm_deriv' should be less than 'degree' i.e., ", degree)
          }
        }
        collectd_mat_i <- collectd_mat_i + 1
        if(dxi == 0) {
          pc.i <- pc
        } else if (dxi > 0) {
          pc.i <- pc[-seq_len(dxi)] * choose(dxi:degree, dxi) * factorial(dxi)
        }
        # already above shifted -> newx[ind[[i]]] - shift_indicator * x[ii]
        collectd_matx <- c(outer(xg - 0, 0:(degree - dxi), "^") %*% pc.i)
        if(dxi == 0) {
          collectd_matx <- collectd_matx + par_a
        }
        parm_mat[ii, 1]                    <- xg_newx
        parm_mat[ii, collectd_mat_i+1]     <- collectd_matx
      }
      i <- i + 1L
    } # while (i <= n_pieces) {
    add_incre_1 <- 0
    if(mark_peak) {
      ya_sitar <- parm_mat[, length(parm_deriv)+1] %>% as.vector()
      ya_sitar[is.na(ya_sitar)] <- 0
      d2_test <- diff(ya_sitar)
      
      peak_id   <- which(d2_test < 0)
      trough_id <- which(d2_test > 0)
      parm_mat[peak_id,   collectd_mat_i+1+1] <- 1
      parm_mat[trough_id, collectd_mat_i+1+1] <- 0
      # Also replace first and last piece to NA, boundary
      parm_mat[1,              collectd_mat_i+1+1] <- NA
      parm_mat[nrow(parm_mat), collectd_mat_i+1+1] <- NA
      add_incre_1 <- 1
    }
    parm_mat[, collectd_mat_i+1+1+add_incre_1]   <- seq(1, parm_mat_dim_nrows)
    parm_mat[, collectd_mat_i+1+1+1+add_incre_1] <- rowni
    if(!is.null(drawni)) parm_mat[, collectd_mat_i+1+1+1+1+add_incre_1] <- drawni
    
    if(1 %in%  return_indicator) {
      parm_array[,,rowni] <- parm_mat
    }
    
    
    
    ############################################################
    # End - parameters
    ############################################################
    
    
    
    
    ############################################################
    # Start - derivatives
    ############################################################
    
    collectd_mat_dim                      <- length(curve_deriv) + 1+1+1
    if(!is.null(drawni)) collectd_mat_dim <- collectd_mat_dim + 1
    collectd_mat_nrows <- 0
    collectd_mat  <- matrix(NA, collectd_mat_nrows, collectd_mat_dim) # + 1 piece + 1 for x + 1 rowni + 1 for draw
    
    i <- 1L
    while (i <= pieces_dim) {
      if(spline_precomputed_indicator == 0) {
        # xg_curve <- seq.int(xknots[i], xknots[i + 1L], length.out = set_spread * degree_dim)
        xg_curve = seq_fun_R(xknots[i], xknots[i+1], degree_dim_set_spread);
      }
      if(spline_precomputed_indicator == 1) {
        xg_curve = xg_curve_array[,i];
      }
      xg_curve <- seq.int(xknots[i], xknots[i + 1L], length.out = set_spread * degree_dim)
      xg_curve <- ((xg_curve/exp(par_c)) + par_b)
      pc <- PiecePolyCoef[, i]
      collectd_mat_x      <- matrix(NA, length(xg_curve), collectd_mat_dim)
      collectd_mat_i <- 0
      for (dxi in curve_deriv) {
        if(!is.null(curve_deriv)) {
          if (dxi >= degree) {
            stop("'curve_deriv' should be less than 'degree' i.e., ", degree)
          }
        }
        collectd_mat_i <- collectd_mat_i + 1
        if(dxi == 0) {
          pc.i <- pc
        } else if (dxi > 0) {
          pc.i <- pc[-seq_len(dxi)] * choose(dxi:degree, dxi) * factorial(dxi)
        }
        collectd_matx <- c(outer(xg_curve - shift_indicator * xknots[i], 0:(degree - dxi), "^") %*% pc.i)
        if(dxi == 0) {
          collectd_matx <- collectd_matx + par_a
        }
        collectd_mat_x[, 1]                    <- xg_curve
        collectd_mat_x[, collectd_mat_i+1]     <- collectd_matx
        collectd_mat_x[, collectd_mat_i+1+1]   <- i
        collectd_mat_x[, collectd_mat_i+1+1+1] <- rowni
        if(!is.null(drawni)) collectd_mat_x[, collectd_mat_i+1+1+1+1] <- drawni
      }
      collectd_mat <- rbind(collectd_mat, collectd_mat_x)
      i <- i + 1L
    } # while (i <= pieces_dim) {
    
    if(2 %in%  return_indicator) {
      deriv_array[,,rowni] <- collectd_mat
    }
    
    ############################################################
    # End - derivatives
    ############################################################
    
  } # for (rowni in 1:n_rowni) {
  
  ############################################################
  # End - over 1:n_rowni
  ############################################################
  
  out <- list()
  
  if(0 %in%  return_indicator) {
    out[['coef']]    <- coef_array
  }
  
  if(1 %in% return_indicator) {
    out[['parm']]   <- parm_array
  }
  
  if(2 %in% return_indicator) {
    out[['deriv']]    <- deriv_array
  }
  
  if(length(out) > 1) {
    return(out)
  } else if(length(out) == 1) {
    return(out[[1]])
  }
  
} # end GS_gps_parms_R


#############################################################
############ -------- my_counter -------- ##########
#############################################################


seq_fun_R <- function(start, end, N_by) { 
  h = (end - start) / (N_by - 1)
  out_c <- c()
  for (i in 1:N_by) { 
    out=start + (i - 1) * h
    out_c <- c(out_c, out)
  }
  return(out_c)
}


#############################################################
############ -------- my_counter -------- ##########
#############################################################

counter_function <- function() {
  env <- environment()
  env$counter <- 0
  increment <- function() {
    env$counter <- env$counter + 1
    return(env$counter)
  }
  get_environment <- function() {
    return(env)
  }
  reset_counter <- function(new_value = 0) {
    env$counter <- new_value
  }
  list(next_value = increment,
       get_environment = get_environment,
       reset = reset_counter)
}
my_counter <- counter_function()


#############################################################
# -------- wraper_for_drawni function for w/t future ------ #
#############################################################

wraper_for_drawni <- function(setdat_mat, 
                              drawni, 
                              callvia,
                              return_indicator, 
                              subset_data_by,
                              subset_data_by_names,
                              create_abcd_names_vector,
                              create_s_names_vector,
                              xknots,
                              degree,
                              spline_subset_indicator = spline_subset_indicator,
                              spline_precomputed_indicator = spline_precomputed_indicator,
                              shift_indicator,
                              set_spread,
                              spline_eval_array,
                              xg_array,
                              xg_curve_array,
                              call_R_stan,
                              GS_gps_parms_assign) {
  if(callvia == 'base') {
    drawniid <- drawni
  } else if(callvia == 'future') {
    # When plan multisession, my_counter will generate sequ per session
    # to get unique sequnce later assign_new_sequence
    drawniid <- my_counter$next_value() 
    pid <- Sys.getpid()
    time_us <- as.numeric(Sys.time()) * 1e6
    # session_prefix <- paste0(pid, "_", floor(time_us))
    # drawniid <- paste0(session_prefix, "_", drawniid) #%>% as.numeric()
    drawniid <- drawniid+time_us+drawniid
  }
  
  if(!is.null(subset_data_by)) {
    setdat_mat <- unique(data.table::as.data.table(setdat_mat),
                                                 by= subset_data_by_names) %>% 
      as.matrix() 
  }
  set_frame_abcd <- setdat_mat[, create_abcd_names_vector, drop = FALSE]
  set_frame_smat <- setdat_mat[, create_s_names_vector,    drop = FALSE]
  # if(call_R_stan == "Stan") {
  #   xg_array          <- xg_array_c
  #   xg_curve_array    <- xg_curve_array_c
  #   spline_eval_array <- spline_eval_array_c
  # }
  
  mat_parm <-  GS_gps_parms_assign(
    nlp_a = set_frame_abcd[, 'a'],
    nlp_b = set_frame_abcd[, 'b'],
    nlp_c = set_frame_abcd[, 'c'],
    nlp_d = set_frame_abcd[, 'd'],
    SParMat = set_frame_smat,
    xknots,
    spline_eval_array,
    xg_array,
    xg_curve_array,
    degree,
    shift_indicator,
    spline_subset_indicator,
    spline_precomputed_indicator = spline_precomputed_indicator,
    set_spread,
    return_indicator = return_indicator,
    drawni = drawniid)
}


#############################################################
############ -------- assign_new_sequence -------- ##########
#############################################################

assign_new_sequence <- function(mat, col) {
  matrix_column <- as.matrix(mat[,col])
  new_sequence <- integer(nrow(matrix_column))
  current_value <- NULL
  sequence_number <- 0
  for (i in 1:nrow(matrix_column)) {
    if (is.null(current_value) || matrix_column[i, 1] != current_value) {
      sequence_number <- sequence_number + 1
      current_value <- matrix_column[i, 1]
    }
    new_sequence[i] <- sequence_number
  } 
  mat[,col] <- new_sequence
  return(mat)
}


#############################################################
############### brms_posterior_summary #####################
#############################################################

brms_posterior_summary <- function(.x, 
                                   probs = c(0.025, 0.975), 
                                   robust = FALSE, 
                                   ...) {
  brms::posterior_summary(.x, probs = probs, robust = robust)
}


#############################################################
############### collapse_posterior_summary ##################
#############################################################

collapse_posterior_summary <- function(.x, 
                                       probs = c(0.025, 0.975), 
                                       robust = FALSE, 
                                       setcolnames = FALSE, 
                                       ...) {
  collapse_mad <- function(x, constant = 1.4826) {
    collapse::fmedian(abs(x - collapse::fmedian(x))) * constant
  }
  if(!robust) {
    out <- c(collapse::fmean(.x), 
             collapse::fsd(.x), 
             collapse::fquantile(.x, probs = probs) 
    ) 
  } else if(robust) {
    out <-  c(collapse::fmedian(.x), 
              collapse_mad(.x), 
              collapse::fquantile(.x, probs = probs) 
    ) 
  }
  # out <- matrix(out, ncol =  1 )
  if(setcolnames) {
    colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 
                                                         100))
  }
  out
} # end collapse_posterior_summary

