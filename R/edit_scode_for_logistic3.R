
#' An internal function to edit stancode for tripple logistic model
#' 
#' @param stancode A string character of stan code
#' 
#' @param set_positive_ordered A logical (default \code{TRUE}) to indicate
#'   whether to set \code{transformed parameters} to \code{positive_ordered} or
#'   \code{ordered}. If \code{TRUE}, \code{transformed parameters} are set to
#'   \code{positive_ordered} and to \code{ordered} if \code{FALSE}. if 
#'   \code{NULL}, set to \code{vector}.  
#' 
#' @param constraint A logical (default \code{TRUE}) to indicate whether to 
#' put constraints on the positive ordered vector space.
#' 
#' @param normalize A logical (default \code{TRUE}) to indicate whether to 
#' include the normalizing constant in the prior target density.
#' 
#' @keywords internal
#' @return A character string.
#' @noRd
#'
edit_scode_for_logistic3 <- function(stancode, 
                                  set_positive_ordered = TRUE,
                                  constraint = TRUE,
                                  normalize = TRUE) {
  
  setorder_d <- c(2, 3, 1)
  setorder_v <- c(3, 1, 2)
  setorder_t <- c(1, 2, 3)
  
  # Seems both set_positive_ordered and  constraint should be TRUE

  true_name_p      <- 'parameters'
  true_name_tp     <- 'transformed parameters'
  true_name_td     <- 'transformed data'
  true_name_model  <- 'model'
  
  tempt_name_p  <- 'ppppppppp'
  tempt_name_tp <- 'tptptptpt'
  tempt_name_td <- 'tdtdtdtdt'
  
  clines_tp <- get_par_names_from_stancode(stancode,
                                           section = true_name_tp,
                                           semicolan = TRUE,
                                           full = TRUE)
  
  clines_p <- get_par_names_from_stancode(stancode,
                                          section = true_name_p,
                                          semicolan = TRUE,
                                          full = TRUE)
  
  clines_m <- get_par_names_from_stancode(stancode,
                                          section = true_name_model,
                                          semicolan = TRUE,
                                          full = TRUE)
  
  
  editedcode    <- stancode 
  editedcode    <- gsub(true_name_tp, tempt_name_tp, editedcode, fixed = T)
  editedcode    <- gsub(true_name_p,  tempt_name_p,  editedcode, fixed = T)
  editedcode    <- gsub(true_name_td,  tempt_name_td,  editedcode, fixed = T)
  
  editedcode2 <- editedcode
  
  
  clines_tp2 <- c()
  for (il in clines_tp) {
    il <- gsub(pattern = "//", replacement = "//", x = il, fixed = T)
    il <- gsub(pattern = "//[^\\\n]*", replacement = "", x = il)
    if(!grepl('^lprior', gsub_space(il)) & 
       !grepl('^reallprior', gsub_space(il)) & 
       !grepl('^-', gsub_space(il))) {
      if(!is_emptyx(il)) {
        clines_tp2 <- c(clines_tp2, il) 
      }
    }
  }
  
  clines_tp <- clines_tp2
  
  
  clines_p2 <- c()
  for (il in clines_p) {
    il <- gsub(pattern = "//", replacement = "//", x = il, fixed = T)
    il <- gsub(pattern = "//[^\\\n]*", replacement = "", x = il)
    if(!grepl('^lprior', gsub_space(il)) & 
       !grepl('^reallprior', gsub_space(il)) & 
       !grepl('^-', gsub_space(il))) {
      if(!is_emptyx(il)) {
        clines_p2 <- c(clines_p2, il) 
      }
    }
  }
  
  clines_p <- clines_p2
  
  
  
  
  b_what_by_pair <- matrix(NA, 9, 3)
  K_name <- "K"
  move_to_tp <- add_move_to_tp <- c()
  for (clines_tpi in clines_p) {
    parameter_name_temp <- sub('.+](.+)', '\\1', clines_tpi)
    parameter_name_temp <- gsub(";", "", parameter_name_temp)
    parameter_name_temp <- gsub_space(parameter_name_temp)
    for (igr in 1:9) {
      if(grepl(paste0("^b", "_"), parameter_name_temp) & 
         grepl(paste0("_", letters[igr]), parameter_name_temp)) {
        move_to_tp <- c(move_to_tp, clines_tpi)
        if(grepl("<lower=", clines_tpi, fixed = T) |
           grepl("<upper=", clines_tpi, fixed = T)) {
          check_bounds_y_n <- 'y'
        } else {
          check_bounds_y_n <- 'n'
        }
        if(letters[igr] == 'a')  poparm <- paste0("raw_d[", setorder_d[1], "]")
        if(letters[igr] == 'd')  poparm <- paste0("raw_d[", setorder_d[2], "]")
        if(letters[igr] == 'g')  poparm <- paste0("raw_d[", setorder_d[3], "]")
        
        if(letters[igr] == 'b')  poparm <- paste0("raw_v[", setorder_v[1], "]")
        if(letters[igr] == 'e')  poparm <- paste0("raw_v[", setorder_v[2], "]")
        if(letters[igr] == 'h')  poparm <- paste0("raw_v[", setorder_v[3], "]")
        
        if(letters[igr] == 'c')  poparm <- paste0("raw_t[", setorder_t[1], "]")
        if(letters[igr] == 'f')  poparm <- paste0("raw_t[", setorder_t[2], "]")
        if(letters[igr] == 'i')  poparm <- paste0("raw_t[", setorder_t[3], "]")
        
        add_move_to_tp_k_dim <- paste0(K_name, "_", letters[igr])
        add_move_to_tp_ <- paste0(" = ", poparm, " + ", "rep_vector(0.0, ",
                                  add_move_to_tp_k_dim, ")")
        add_move_to_tp_2 <- gsub(";", paste0(add_move_to_tp_, ";"),  
                                 clines_tpi, fixed = T)
        
        add_move_to_tp <- c(add_move_to_tp, add_move_to_tp_2)
        parameter_name <- sub('.+](.+)', '\\1', clines_tpi)
        parameter_name <- gsub(";", "", parameter_name)
        parameter_name <- gsub_space(parameter_name)
        b_what_by_pair[igr , 1] <- parameter_name
        b_what_by_pair[igr , 2] <- poparm
        b_what_by_pair[igr , 3] <- check_bounds_y_n
      } 
    }
  } 
  
   # b_what_by_pairx <<- b_what_by_pair
  
  
  b_what_it_c <- b_what_by_c <- c()
  for (igr in 1:nrow(b_what_by_pair)) {
    set_check_bounds_y_n_cnt <- 0
    for (igri in 1:2) {
      set_check_bounds_y_n_cnt <- set_check_bounds_y_n_cnt + 1
      b_what_it <- b_what_by_pair[igr , 1]
      b_what_by <- b_what_by_pair[igr , 2]
      check_bounds_y_n <- b_what_by_pair[igr , 3]
      if(check_bounds_y_n == "n") {
        b_what_by <- gsub("[", paste0("[", igri, ","), b_what_by, fixed = T)
        b_what_it <- paste0(b_what_it, "[", igri,"]")
        b_what_it_c <- c(b_what_it_c, b_what_it)
        b_what_by_c <- c(b_what_by_c, b_what_by)
      } else if(check_bounds_y_n == "y") {
        if(set_check_bounds_y_n_cnt == 1) {
          b_what_by <- gsub("[", paste0("[", ","), b_what_by, fixed = T)
          b_what_it <- b_what_it
          b_what_it_c <- c(b_what_it_c, b_what_it)
          b_what_by_c <- c(b_what_by_c, b_what_by)
        }
      } # else if(check_bounds_y_n
    }
  }
  
  # b_what_it_cx <<- b_what_it_c
  # b_what_by_cx <<- b_what_by_c
  
  set_b_a <- paste0("b_a = to_vector(raw_dx[", ",", setorder_d[1], "]);")
  set_b_d <- paste0("b_d = to_vector(raw_dx[", ",", setorder_d[2], "]);")
  set_b_g <- paste0("b_g = to_vector(raw_dx[", ",", setorder_d[3], "]);")
  
  set_b_b <- paste0("b_b = to_vector(raw_vx[", ",", setorder_v[1], "]);")
  set_b_e <- paste0("b_e = to_vector(raw_vx[", ",", setorder_v[2], "]);")
  set_b_h <- paste0("b_h = to_vector(raw_vx[", ",", setorder_v[3], "]);")
  
  set_b_c <- paste0("b_c = to_vector(raw_tx[", ",", setorder_t[1], "]);")
  set_b_f <- paste0("b_f = to_vector(raw_tx[", ",", setorder_t[2], "]);")
  set_b_i <- paste0("b_i = to_vector(raw_tx[", ",", setorder_t[3], "]);")
  
  set_b_elements <- paste(set_b_a, set_b_b, set_b_c,
                          set_b_d, set_b_e, set_b_f,
                          set_b_g, set_b_h, set_b_i,
                          sep = "\n   ")
  
  
  if(!is.null(set_positive_ordered)) {
    if(set_positive_ordered) {
      move_to_tp_add_ordered_positive_ordered <- 
        "array[Kedit] positive_ordered[Cedit] raw_dx;
   array[Kedit] positive_ordered[Cedit] raw_vx;
   array[Kedit] positive_ordered[Cedit] raw_tx;"
    }
    
    if(!set_positive_ordered) {
      move_to_tp_add_ordered_positive_ordered <- 
        "array[Kedit] ordered[Cedit] raw_dx;
   array[Kedit] ordered[Cedit] raw_vx;
   array[Kedit] ordered[Cedit] raw_tx;"
    }
  } else if(is.null(set_positive_ordered)) {
    move_to_tp_add_ordered_positive_ordered <- 
      "array[Kedit] vector[Cedit] raw_dx;
   array[Kedit] vector[Cedit] raw_vx;
   array[Kedit] vector[Cedit] raw_tx;"
  }
  
  # print(set_positive_ordered)
  # print(move_to_tp_add_ordered_positive_ordered)
  
  if(constraint) {
    move_to_tp_add_constraint <- 
      " for(k in 1:Kedit) {
     raw_dx[k, ] = ordered_lb_ub_lp(raw_d[k,], min_d[k], max_d[k]);
     raw_vx[k, ] = ordered_lb_ub_lp(raw_v[k,], min_v[k], max_v[k]);
     raw_tx[k, ] = ordered_lb_ub_lp(raw_t[k,], min_t[k], max_t[k]);
    } "
  } else if(!constraint) {
    move_to_tp_add_constraint <- 
    " for(k in 1:Kedit) {
     raw_dx[k, ] = raw_d[k,];
     raw_vx[k, ] = raw_v[k,];
     raw_tx[k, ] = raw_t[k,];
    } "
  } 
  
  
  move_to_tp_add <- paste(move_to_tp_add_ordered_positive_ordered,
                          move_to_tp_add_constraint,
                          sep = "\n  ")
                          
  
  move_to_tp_add <- paste(move_to_tp_add, 
                          set_b_elements, 
                          sep = "\n   ")
  

  move_to_tp2 <- c()
  for (move_to_tpi in move_to_tp) {
    move_to_tp2 <- c(move_to_tp2, paste0("   ", move_to_tpi))
  }
  
  tpcode <- paste0(paste(move_to_tp2, collapse = "\n"), 
                   "\n    ", 
                   move_to_tp_add)
  
  
  pcode <- 
    "array[Kedit] vector[Cedit] raw_d;
  array[Kedit] vector[Cedit] raw_v;
  array[Kedit] vector[Cedit] raw_t;
  "
  
  
  # parameter names - remove from the parameters block
  for (il in move_to_tp) {
    editedcode2 <- gsub(pattern = "//", replacement = "//", 
                        x = editedcode2, fixed = T)
    editedcode2 <- gsub(pattern = "//[^\\\n]*", replacement = "", 
                        x = editedcode2)
    editedcode2 <- gsub(paste0(il, ""), "", editedcode2, fixed = T)
    
  }
  
  
  
  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      zz_c <- c(zz_c, zz_in)
    }
  }
  editedcode2 <- paste(zz_c, collapse = '\n')
  
  
  p_block_syb_by <- paste0("", tempt_name_tp, " {")
  p_block_syb_it <- paste0(p_block_syb_by, "\n", tpcode)
  editedcode2 <- gsub(paste0("", p_block_syb_by), p_block_syb_it, 
                      editedcode2, fixed=T, perl=F)
  

  editedcode2 <- gsub(tempt_name_tp, true_name_tp, editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_p,  true_name_p,  editedcode2, fixed = T)
  editedcode2 <- gsub(tempt_name_td,  true_name_td,  editedcode2, fixed = T)
  
  
  # Replace parameters in priors
  for (igr in 1:length(b_what_it_c)) {
    parameter_name <- b_what_it_c[igr]
    poparm         <- b_what_by_c[igr]
    editedcode2 <- gsub(paste0("(", parameter_name),
                        paste0("(", poparm),
                        editedcode2, fixed = T)
  }
  
  # Remove empty lines
  zz <- strsplit(editedcode2, "\n")[[1]]
  zz_c <- c()
  for (iz in 1:length(zz)) {
    if(!is_emptyx(gsub_space(zz[iz]))) {
      zz_in <- zz[iz]
      zz_c <- c(zz_c, zz_in)
    }
  }
  
  editedcode2 <- paste(zz_c, collapse = '\n')
  
  # https://github.com/stan-dev/math/issues/2959
  
  fcode <- 
  "vector ordered_lb_ub_lp (vector y, real lb, real ub) {
    int N = rows(y);
    vector[N] x;
    x[1] = lb + (ub - lb) * inv_logit(y[1]);
    target += log(ub - lb) + log_inv_logit(y[1]) + log1m_inv_logit(y[1]);
    for (i in 2:N) {
      x[i] = x[i - 1] + (ub - x[i - 1]) * inv_logit(y[i] - log(N+1-i));
      target += log(ub - x[i - 1]) + log_inv_logit(y[i] - log(N+1-i)) + 
        log1m_inv_logit(y[i] - log(N+1-i));
    }
    return x;
  }"
  
  # cat(editedcode2)
  # stop()
  
  out <- list(editedcode = editedcode2, 
              tpcode = tpcode, 
              pcode = pcode, 
              fcode = fcode)
  return(out)
}



