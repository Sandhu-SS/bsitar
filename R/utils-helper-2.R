


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
  
  # Seems both set_positive_ordered and constraint should be TRUE
  
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
  
  
  
  out <- list(editedcode = editedcode2, 
              tpcode = tpcode, 
              pcode = pcode, 
              fcode = fcode)
  return(out)
}






#' Identify (and remove) outliers with abnormal growth velocity
#'
#' @description The function identifies and remove the putative outliers for y
#'   in data, Codes range from 0 (normal) to 8, where 4 and 6 are conventional
#'   outliers
#'
#' @inherit sitar::velout params description
#' 
#' @inherit sitar::zapvelout params
#' 
#' @param remove A logical (default \code{FALSE}) to indicate whether identified
#'   remove to be removed. When \code{FALSE}, outliers are set as \code{NA}. If
#'   \code{TRUE}, outliers set as \code{NA} are removed.
#'   
#' @param verbose A logical (default \code{FALSE}) to show frequency of
#'   observations with different codes (see [sitar::velout()]), and the number
#'   of observations (outliers) removed when \code{remove=TRUE}.
#'
#' @return A data frame with outliers removed when \code{remove=TRUE}
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
#' @examples
#' 
#' \dontrun{
#' data <- berkeley_mdata
#' 
#' outliers(x = "age", y = "height", id = "id", data = data, remove=TRUE)
#' }
#' 
outliers <-
  function (x,
            y,
            id,
            data,
            icode = c(4:6),
            lag = 1,
            velpower = 0.5,
            limit = 5,
            linearise = FALSE,
            remove = FALSE,
            verbose = TRUE) {
    
    mcall <- match.call()
    if(is.symbol(x)) xx_ <- deparse(substitute(x)) else xx_ <- x
    if(is.symbol(y)) yy_ <- deparse(substitute(y)) else yy_ <- y
    if(is.symbol(id)) idid_ <- deparse(substitute(id)) else idid_ <- id
    
    data <- data %>% dplyr::mutate(order = dplyr::row_number())
    data <- data %>% dplyr::arrange(idid_, xx_)
    
    dc <- data %>% dplyr::select(!!as.symbol(xx_), 
                                 !!as.symbol(yy_),
                                 !!as.symbol(idid_))
    
    colnames_data_ex <- colnames(data)
    colnames_dc_ex <- colnames(dc)
    data_ex <- data %>% dplyr::select(-colnames(dc))
    
    nrow <- nrow(data)
    dc <- na.omit(cbind(dc, count = 1:nrow))
    dc <- dc[order(dc[, 3], dc[, 1]),]
    if (linearise) {
      spline.lm <- loess(dc[, 2] ~ dc[, 1])
      dc[, 2] <- residuals(spline.lm)
    }
    dt1 <- diff(dc[, 1], lag = lag)
    vel1 <- diff(dc[, 2], lag = lag) / dt1 ^ velpower
    dt2 <- diff(dc[, 1], lag = lag * 2)
    vel2 <- diff(dc[, 2], lag = lag * 2) / dt2 ^ velpower
    dt1 <- dt1 == 0
    idlev <- as.numeric(dc[, 3])
    dt1[diff(idlev, lag = lag) != 0] <- FALSE
    vel1[diff(idlev, lag = lag) != 0] <- NA
    vel2[diff(idlev, lag = lag * 2) != 0] <- NA
    vel1 <- trunc(vel1 / mad(vel1, na.rm = TRUE))
    vel2 <- trunc(vel2 / mad(vel2, na.rm = TRUE))
    vel3 <- c(rep(NA, lag), vel2, rep(NA, lag))
    vel2 <- c(vel1, rep(NA, lag))
    vel1 <- c(rep(NA, lag), vel1)
    code <- (as.numeric(abs(vel1) >= limit) + as.numeric(abs(vel2) >=
                                                           limit)) * 2 + 
      as.numeric(abs(vel3) >= limit)
    
    dt2 <- c(dt1, rep(FALSE, lag))
    dt1 <- c(rep(FALSE, lag), dt1)
    code[dt2 | dt1] <- 8
    code[dt2 & !dt1] <- 7
    t <- is.na(vel3) & !(dt1 | dt2)
    code[t] <- (as.numeric(!is.na(vel1[t]) & abs(vel1[t]) >=
                             limit) + as.numeric(!is.na(vel2[t]) &
                                                   abs(vel2[t]) >=
                                                   limit)) * 6
    dc <- cbind(dc[, c(3, 1, 2, 4)], code, vel1, vel2, vel3)
    mat <- as.data.frame(matrix(
      nrow = nrow,
      ncol = dim(dc)[2],
      dimnames = list(row.names(data), dimnames(dc)[[2]])
    ))
    attr(mat, "data") <- deparse(mcall$data)
    mat[dc$count,] <- dc
    if (is.factor(dc[, 1])) {
      mat[, 1] <- as.factor(mat[, 1])
      levels(mat[, 1]) <- levels(dc[, 1])
    }
    mat$count <- NULL
    mat$code <- factor(mat$code)
    if (verbose) {
      cat("code frequencies\n")
      print(summary(mat$code))
    }
    
    mat <- cbind(mat, data_ex)
    mat <- mat %>% dplyr::relocate(dplyr::all_of(colnames_data_ex))
    
    if (remove) {
      zap <- mat$code %in% icode
      if (verbose) {
        cat(sum(zap), yy_, "values set missing\n")
      }
      mat[mat$code %in% icode, yy_] <- NA
      mat <- mat %>% dplyr::select(dplyr::all_of(colnames_data_ex))
      mat <- mat %>% tidyr::drop_na() %>% droplevels()
    }
    mat <- mat %>% dplyr::arrange(order)
    mat <- mat %>% dplyr::select(-order)
    mat
  }



#' An internal function to set default initials for model specific parameters
#'
#' @param cargs A character string specifying the model fitted. 
#' 
#' @param fargs A character string specifying the initials. 
#' 
#' @param dargs A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param verbose A logical (default \code{FALSE}) to indicate whetehr to print
#'  \code{warnings} and \code{messages} during the function evaluation.
#'
#' @return A list.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
evaluate_call_args <- function(cargs = NULL,
                            fargs = NULL,
                            dargs = NULL,
                            verbose = FALSE) {
  for (fargsi in names(dargs)) {
    if(is.null(cargs[[fargsi]])) cargs[[fargsi]] <- fargs[[fargsi]]
  }
  for (fargsi in names(fargs)) {
    if(is.null(cargs[[fargsi]])) cargs[[fargsi]] <- fargs[[fargsi]]
  }
 
  cargs
}



#' An internal function to perform checks when calling post-processing functions
#'
#' @description The \code{post_processing_args_sanitize} perform essential 
#'   checks for the arguments passed on to the \code{brms} post-processing
#'   functions.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param xcall The \code{match.call()} from the post-processing function.
#'
#' @param resp Response variable (default \code{NULL}) specified as a character
#'   string. This is needed when processing \code{univariate_by} and
#'   \code{multivariate} models (see \code{bgmfit} function for details).
#'
#' @param deriv An integer value to specify whether to estimate distance curve
#'   (i.e., model estimated curve(s)) or the velocity curve (first derivative of
#'   the model estimated curve(s)). A value \code{0} (default) is for distance
#'   curve and  \code{1} for the velocity curve.
#'   
#' @param dots A list passing the \code{...} arguments. Default \code{NULL}.
#' 
#' @param misc A vector of character strings specifying the miscellanous
#'   arguments to be checked. Default \code{NULL}.
#'
#' @param envir Indicator to set the environment of function evaluation. The
#'   default is \code{parent.frame}.
#'   
#' @param verbose A logical (default \code{FALSE}) to indicate whetehr to print
#'  \code{warnings} and \code{messages} during the function evaluation.
#'
#' @return A list with the filtered necessary arguments required for the 
#'   successgull execition of the post-processing function
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
post_processing_args_sanitize <- function(model, 
                                          xcall, 
                                          resp = NULL, 
                                          deriv = NULL,
                                          dots = NULL,
                                          misc = NULL,
                                          envir = NULL, 
                                          verbose = FALSE) {
  
  if(is.null(envir)) envir <- parent.frame()
  if(is.null(deriv)) deriv <- 0
  
  if(!'bgmfit' %in% class(model)) {
    stop("The class of model object should be 'bgmfit' ")
  }
  
  allargs <- c(as.list(xcall), dots)
  
  excall_ <- c("plot_ppc", "loo_validation")
  excall_ <- c(excall_, paste0(excall_, ".", "bgmfit"))
  
  checkwhat_all <- c()
  if (strsplit(deparse((xcall[1])), "\\.")[[1]][1] %in% excall_) {
    # Check deriv
    checkwhat <- ""
    checkwhat <- 'deriv'
    if(!is.null(allargs[[checkwhat]])) {
      checkwhat_all <- c(checkwhat_all, checkwhat)
      if(verbose) {
        message(
          "\nargument 'deriv' is not allowed for the ",
          " post-processing function",  " '",
          strsplit(deparse((xcall[1])), "\\.")[[1]][1], "'",
          "\n ",
          "Therefore, it is set to i.e., deriv = NULL"
        )
      }
    } # if(!is.null(allargs$idata_method)) {
    
    # Check idata_method
    checkwhat <- ""
    checkwhat <- 'idata_method'
    if(!is.null(allargs[[checkwhat]])) {
      checkwhat_all <- c(checkwhat_all, checkwhat)
      if(verbose) {
        message(
          "\nargument 'idata_method' is not allowed for the ",
          " post-processing function",  " '",
          strsplit(deparse((xcall[1])), "\\.")[[1]][1], "'",
          "\n ",
          "Therefore, it is set to i.e., idata_method = NULL"
        )
      }
    } # if(!is.null(allargs$idata_method)) {
    
    # Check re_formula
    checkwhat <- ""
    checkwhat <- 're_formula'
    if(!is.null(allargs[[checkwhat]])) {
      checkwhat_all <- c(checkwhat_all, checkwhat)
      if(verbose) {
        message(
          "\nargument 'idata_method' is not allowed for the ",
          " post-processing function",  " '",
          strsplit(deparse((xcall[1])), "\\.")[[1]][1], "'",
          "\n ",
          "Therefore, it is set to i.e., idata_method = NULL"
        )
      }
    } # if(!is.null(allargs$idata_method)) {
    
    
  } # if (strsplit(deparse((xcall[1])), "\\.")[[1]][1] %in% excall_) {
  
  
  if(length(checkwhat_all) == 0) checkwhat_all <- NULL
  checkwhat_all <- c(checkwhat_all, misc)
  
  if(!is.null(checkwhat_all)) {
    allargs[which(names(allargs)%in%checkwhat_all)] <- NULL
  }
  
  allargs <- allargs[-1]
  names(allargs) <- gsub("model", "object", names(allargs))
 
  return(allargs)
}



#' An internal function to perform checks when calling post-processing functions
#'
#' @description The \code{post_processing_checks} perform essential checks (such
#'   as the validity of model class, response etc.) during post-processing of
#'   posterior draws.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param xcall The \code{match.call()} from the post-processing function.
#'
#' @param resp Response variable (default \code{NULL}) specified as a character
#'   string. This is needed when processing \code{univariate_by} and
#'   \code{multivariate} models (see \code{bgmfit} function for details).
#'
#' @param deriv An integer value to specify whether to estimate distance curve
#'   (i.e., model estimated curve(s)) or the velocity curve (first derivative of
#'   the model estimated curve(s)). A value \code{0} (default) is for distance
#'   curve and  \code{1} for the velocity curve.
#'   
#' @param all A logical (default \code{NULL}) to specify whether to return all
#'   the exposed functions.
#'
#' @param envir Indicator to set the environment of function evaluation. The
#'   default is \code{parent.frame}.
#'   
#' @param verbose A logical (default \code{FALSE}) to indicate whetehr to print
#'  \code{warnings} and \code{messages} during the function evaluation.
#'
#' @return A string with the error captured or else a list with necessary
#'   information needed when executing the post-processing function
#'   
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
post_processing_checks <- function(model, 
                                   xcall, 
                                   resp = NULL, 
                                   deriv = NULL,
                                   all = FALSE,
                                   envir = NULL, 
                                   verbose = FALSE, 
                                   check_d0 = FALSE,
                                   check_d1 = FALSE,
                                   check_d2 = FALSE) {
  
  if(is.null(envir)) envir <- parent.frame()
  if(is.null(deriv)) deriv <- 0
  
  if(!'bgmfit' %in% class(model)) {
    stop("The class of model object should be 'bgmfit' ")
  }
  
  excall_ <- c("plot_ppc", "loo_validation")

  xcall_check_it <- paste(deparse(substitute(xcall)), collapse = "")
  xcall_check_it <- gsub_space(xcall_check_it)
  check_it       <- sub(" *\\(.*", "", xcall_check_it)
 
  check_it_sss <- strsplit(check_it, "\\.")[[1]][1]
  
  # Now using checkresp_info which is also used in getmodel_info
  checkresp_info(model, resp)
 
  # if (model$model_info$nys == 1 & !is.null(resp)) {
  #   stop(
  #     "You have fit a univariate model",
  #     " but set resp option as: ",
  #     resp,
  #     ".",
  #     "\n ",
  #     " For univariate model, the resp option should be NULL",
  #     "\n ",
  #     " (i.e., resp = NULL)"
  #   )
  # }
  # if (model$model_info$nys > 1 & is.null(resp)) {
  #   if (!is.na(model$model_info$univariate_by$by)) {
  #     stop(
  #       "You have fit a univariate_by model for ",
  #       model$model_info$univariate_by$by,
  #       "\n ",
  #       " but did not correctly specified the 'resp' option",
  #       " (which is NULL at present).",
  #       "\n ",
  #       " The response options are: ",
  #       paste(model$model_info$ys, collapse = ", ")
  #     )
  #   }
  #   if (model$model_info$multivariate$mvar) {
  #     stop(
  #       "You have fit a multivariate model ",
  #       "\n ",
  #       " but dit not set the the resp options correctly",
  #       " (which is NULL at present).",
  #       "\n ",
  #       " The response options are: ",
  #       paste(model$model_info$ys, collapse = ", ")
  #     )
  #   }
  # }
  
  if (is.null(resp)) {
    resp_    <- resp
    revresp_ <- ""
  } else if (!is.null(resp)) {
    resp_    <- paste0(resp, "_")
    revresp_ <- paste0("_", resp)
  }
  
  # assign expose default funs 
  if(model$model_info[['expose_method']] == 'R') {
    assign(paste0(resp_, 
                  model$model_info[['namesexefuns']], 
                  '0'), 
           model$model_info$exefuns[[paste0(resp_, 
                                            model$model_info[['namesexefuns']], 
                                            '0')]], envir = envir)

    if(model$model_info[['select_model']] == 'sitar' |
       model$model_info[['select_model']] == 'rcs') {
      assign(paste0(resp_, 'getKnots'), 
             model$model_info$exefuns[[paste0(resp_, 'getKnots')]], 
             envir = envir)
    }
  }
 
  if(!all) {
    out <-
      list(
        paste0(resp_, model$model_info[['namesexefuns']], ''),
        paste0(resp_, model$model_info[['namesexefuns']], deriv)
      ) 
  }
  
  if(all) {
    out <- model$model_info[['exefuns']]
  } 
  
  # 6.03.2025
  # check_d1
  if(check_d1) {
    available_d0 <- available_d1 <- available_d2 <- FALSE
    for (i in names(model$model_info$exefuns)) {
      check_funds <- ifelse(grepl("\\d$", i), sub(".*?(\\d+)$", "\\1", i), "")
      if(grepl("0", check_funds)) {
        available_d0 <- TRUE
      }
      if(grepl("1", check_funds)) {
        available_d1 <- TRUE
      }
      if(grepl("2", check_funds)) {
        available_d2 <- TRUE
      }
    }
    
    if(verbose) { 
      if(!available_d0) {
        stop("No 'd0' found")
      }
      if(!available_d1) {
        message("No 'd1' found, setting 'model_deriv = FALSE', 'deriv = 0'")
      }
      if(!available_d2) {
        # message("No 'd2' found, setting 'model_deriv = FALSE', 'deriv = 0'")
      }
    }
    
    # sigma_model <- model$model_info[['sigma_model']]
    
    # When sigma_model == 'ls', then only sigma function is defined in stan
    # if(!is.null(sigma_model)) {
    #   if(sigma_model == "mu") {
    #     available_d1 <- FALSE
    #     available_d2 <- FALSE
    #   }
    # }
    
    sigma_model_      <- paste0('sigma_model', revresp_)
    sigma_model_name_ <- paste0('sigmabasicfunname', revresp_)
    sigma_model_attr_ <- paste0('sigmabasicfunattr', revresp_)
    
    sigma_model       <- model$model_info[[sigma_model_]]
    sigma_model_name  <- model$model_info[[sigma_model_name_]]
    sigma_model_attr  <- model$model_info[[sigma_model_attr_]]
    
    check_namespace_for_sigma_d1 <- c('splines2', 'bsitar')
    
    check_namespace_for_sigma_d1_test <- 
      setdiff(sigma_model_attr, check_namespace_for_sigma_d1)
    
    if(!is.null(sigma_model)) {
      if(sigma_model == "ls") {
        available_d1 <- available_d1
        available_d2 <- available_d2
      } else if(sigma_model == "basic") {
        if(length(check_namespace_for_sigma_d1_test) == 0) {
          available_d1 <- available_d1
          available_d2 <- available_d2
          if(deriv > 0) {
            for (i in sigma_model_name) {
              model$model_info$exefuns[[i]] <- 
                formals(model$model_info$exefuns[[i]])[['derivs']] <- deriv
            }
          }
        } else {
          available_d1 <- FALSE
          available_d2 <- FALSE
        }
      } else {
        available_d1 <- available_d1
        available_d2 <- available_d2
      }
    } else if(is.null(sigma_model)) {
      available_d1 <- available_d1
      available_d2 <- available_d2
    }
    
    
    
    
    
    # Force available_d1 = FALSE when model_deriv = FALSE
    if(!is.null(model$model_info[['model_deriv']])) {
      if(!model$model_info[['model_deriv']]) {
        available_d1 <- FALSE
        available_d2 <- FALSE
      }
    }
    
    
    # if(model$model_info[['dpar']] == 'sigma') {
    #   for (i in 1:length(out)) {
    #     out[[i]] <- ""
    #   }
    # }
    
    # print(out)
    # stop()
    
    out[['available_d0']] <- available_d0
    out[['available_d1']] <- available_d1
    out[['available_d2']] <- available_d2
  } # if(check_d1) {
  
  
  return(out)
}




#' An internal function to set up exposed functions and their environment
#'
#' @param o A logical (default \code{FALSE}) to indicate whether to
#' return the object as a character string.
#' @inherit growthparameters.bgmfit params
#' @param ... Additional arguments. Currently ignored.
#' @keywords internal
#' @return A list comprised of exposed functions.
#' @noRd
#'

setupfuns <- function(model,
                      resp = NULL,
                      o = NULL,
                      oall = NULL,
                      usesavedfuns = NULL,
                      deriv = NULL,
                      envir = NULL,
                      model_deriv = NULL,
                      verbose = FALSE,
                      ...) {
  
  if(is.null(envir)) {
    envir <- parent.frame()
  }
  
  if (is.null(resp)) {
    resp_ <- resp
  } else if (!is.null(resp)) {
    resp_ <- paste0(resp, "_")
  }
  
  if(is.null(model$xcall)) {
    xcall <- strsplit( deparse(sys.calls()[[sys.nframe()-1]]) , "\\(")[[1]][1]
  } else {
    xcall <- model$xcall
  }
  
  
  excall_ <- c("plot_ppc", "loo_validation")
  
  check_it <- strsplit(deparse((xcall[1])), "\\.")[[1]][1] 
  check_it <- gsub("\"",  "", check_it)
  
  if (check_it %in% excall_) {
    if(is.null(as.list(xcall)[['deriv']])) deriv <- ''
    if (!is.null(as.list(xcall)[['deriv']])) {
      deriv <- ''
      if(verbose) {
        message(
          "\nargument 'deriv' is not allowed for the ",
          " post-processing function",  " '",
          check_it, "'",
          "\n ",
          "Therefore, it is set to missing i.e., deriv = ''"
        )
      }
    } # if(!is.null(chcallls$idata_method)) {
  }
  
  if(!usesavedfuns) {
    if(is.null(check_if_functions_exists(model, o, xcall))) {
      return(invisible(NULL))
    }
  }
  
  if(usesavedfuns) {
    if(is.null(check_if_functions_exists(model, o, xcall,
                                         verbose = F))) {
      envir <- envir
    } else {
      #  envir <- getEnv(o[[1]], geteval = TRUE)
    }
    # envir <- getEnv(o[[1]], geteval = TRUE)
    
    oall <- model$model_info[['exefuns']]
    oalli_c <- names(oall)
    for (oalli in oalli_c) {
      assign(oalli, oall[[oalli]], envir = envir)
    }
  }
  
  if(!is.null(deriv)) {
    if(deriv == 0) {
      assignfun <- paste0(model$model_info[['namesexefuns']], deriv)
      assignfun <- paste0(resp_, assignfun)
      assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
    } else if(deriv > 0) {
      if(model_deriv) {
        assignfun <- paste0(model$model_info[['namesexefuns']], deriv)
        assignfun <- paste0(resp_, assignfun)
      } else if(!model_deriv) {
        assignfun <- paste0(model$model_info[['namesexefuns']], '0')
        assignfun <- paste0(resp_, assignfun)
      }
      assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
    }
  }
  
  if(is.null(deriv)) {
    assignfun <- paste0(model$model_info[['namesexefuns']], "")
    assignfun <- paste0(resp_, assignfun)
    assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
  }
  return(envir)
}







#' An internal function to set default prior for model specific parameters
#'
#' @param select_model A character string specifying the model fitted. 
#' 
#' @param prior A character string specifying the prior. 
#' 
#' @param class A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param parameter A character string specifying the parameter name. Options
#'  are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'}, \code{'e'}, \code{'f'},
#'  \code{'g'}, \code{'h'}, and \code{'i'}. Default \code{NULL} indicates that
#'  parameter name in infered automatically. 
#'
#' @return A character string.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
set_default_priors <- function(select_model,
                               prior,
                               class = NULL,
                               parameter = NULL) {
  if (is.null(parameter)) {
    parameter <- strsplit(deparse(substitute(prior)), "_")[[1]][1]
  }
  if (is.null(class)) {
    get_suffix <- strsplit(deparse(substitute(prior)), "_")[[1]]
    get_suffix <- get_suffix[length(get_suffix)]
    if (grepl("^beta", get_suffix))
      class <- 'b'
    if (grepl("^sd", get_suffix))
      class <- 'sd'
  }
  
  
  ##############################################################
  # class b
  ##############################################################
  
  # parameter a class b
  if (parameter == 'a' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(ymin, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter b class b
  if (parameter == 'b' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymaxs, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymaxs, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.5, 0.5, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter c class b
  if (parameter == 'c' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0.1, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(5, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter d class b
  if (parameter == 'd' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(ymeanxmid, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter e class b
  if (parameter == 'e' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(7, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0.2, 0.1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter f class b
  if (parameter == 'f' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(13, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(5, 3, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class b
  if (parameter == 'g' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- 
            "normal(ymeanxmidxmaxdiff, ysdxmidxmaxdiff, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter h class b
  if (parameter == 'h' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.5, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter i class b
  if (parameter == 'i' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(14, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  
  
  
  ##############################################################
  # class sd
  ##############################################################
  
  # parameter a class sd
  if (parameter == 'a' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter b class sd
  if (parameter == 'b' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter c class sd
  if (parameter == 'c' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter d class sd
  if (parameter == 'd' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter e class sd
  if (parameter == 'e' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 0.15, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter f class sd
  if (parameter == 'f' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'g' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, ysdxmidxmaxdiff, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter h class sd
  if (parameter == 'h' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter i class sd
  if (parameter == 'i' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  return(prior_out)
}







#' An internal function to set default initials for model specific parameters
#'
#' @param select_model A character string specifying the model fitted. 
#' 
#' @param init A character string specifying the initials. 
#' 
#' @param class A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param parameter A character string specifying the parameter name. Options
#'  are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'}, \code{'e'}, \code{'f'},
#'  \code{'g'}, \code{'h'}, and \code{'i'}. Default \code{NULL} indicates that
#'  parameter name in infered automatically. 
#'
#' @return A character string.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
set_default_inits <- function(select_model,
                              init,
                              class = NULL,
                              parameter = NULL) {
  if (is.null(parameter)) {
    parameter <- strsplit(deparse(substitute(init)), "_")[[1]][1]
  }
  if (is.null(class)) {
    get_suffix <- strsplit(deparse(substitute(init)), "_")[[1]]
    get_suffix <- get_suffix[length(get_suffix)]
    if (grepl("^beta", get_suffix))
      class <- 'b'
    if (grepl("^sd", get_suffix))
      class <- 'sd'
  }
  
  
  ##############################################################
  # class b
  ##############################################################
  
  # parameter a class b
  if (parameter == 'a' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "ymean"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "ymean"
      } else if (grepl('^pb', select_model)) {
        init_out <- "ymax"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "ymax"
        }
        if (select_model == 'logistic2') {
          init_out <- "ymax"
        }
        if (select_model == 'logistic3') {
          init_out <- "ymin"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter b class b
  if (parameter == 'b' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- "ymaxs"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 0.1
        }
        if (select_model == 'logistic2') {
          init_out <- "ymaxs"
        }
        if (select_model == 'logistic3') {
          init_out <- 1.5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter c class b
  if (parameter == 'c' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- 0.1
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 5.0
        }
        if (select_model == 'logistic2') {
          init_out <- 0.1
        }
        if (select_model == 'logistic3') {
          init_out <- 0.1
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter d class b
  if (parameter == 'd' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 1.0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- 1.2
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 1.0
        }
        if (select_model == 'logistic2') {
          init_out <- 1.2
        }
        if (select_model == 'logistic3') {
          init_out <- "ymeanxmid"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter e class b
  if (parameter == 'e' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        if(select_model == 'pb1') init_out <- 13
        if(select_model == 'pb2') init_out <- 13
        if(select_model == 'pb3') init_out <- 13
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- 7
        }
        if (select_model == 'logistic3') {
          init_out <- 0.15
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter f class b
  if (parameter == 'f' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        if(select_model == 'pb1') init_out <- NULL
        if(select_model == 'pb2') init_out <- 2
        if(select_model == 'pb3') init_out <- 1
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- 13
        }
        if (select_model == 'logistic3') {
          init_out <- 5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter g class b
  if (parameter == 'g' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 
            "ymeanxmidxmaxdiff"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter h class b
  if (parameter == 'h' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 1.5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter i class b
  if (parameter == 'i' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 13
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  
  
  
  ##############################################################
  # class sd
  ##############################################################
  
  # parameter a class sd
  if (parameter == 'a' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter b class sd
  if (parameter == 'b' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter c class sd
  if (parameter == 'c' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter d class sd
  if (parameter == 'd' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter e class sd
  if (parameter == 'e' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 0.15, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter f class sd
  if (parameter == 'f' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'g' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmidxmaxdiff, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter h class sd
  if (parameter == 'h' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter i class sd
  if (parameter == 'i' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  return(init_out)
}






#' An internal function to set variance covariance initials values to zero
#'
#' @param xscode A character string to specify the stancode.
#'
#' @param xsdata A character string to specify the standata.
#'
#' @param full A logical (default \code{TRUE}) to indicate whether full names
#'   should be extracted.
#'
#' @param what A character string to specify the variance covariance parameter.
#'   Default \code{L} indicating the correlation parameters. Other options are
#'   \code{sd} and \code{z}.
#'   
#' @param parameterization A character string to specify the 
#' the parameterization, CP or NCP.
#'
#' @param sd_value A numeric value to set intials for standard deviation
#'   parameters. Default \code{1} which is translated to zero initial i.e.,
#'   \code{log(1) = 0}.
#'
#' @param z_value A numeric value (default \code{0}) to set initials for
#'   \code{z} parameter which is part of the non centered parameterisation
#'   implemented in the [brms::brm()].
#'
#' @param L_value A numeric value (default \code{0}) to set initials for
#'   correlation parameter, \code{L}.
#'   
#' @return A prior object.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
set_init_gr_effects <- function(xscode,
                                xsdata,
                                full = TRUE,
                                what = 'L',
                                parameterization = 'ncp',
                                sd_value = 1,
                                z_value = 0,
                                r_value = 0,
                                L_value = 0) {
  xscode <-
    get_par_names_from_stancode(xscode, full = full, what = what)
  
  sdi_c <- c()
  parm_c <- c()
  sd_init_list_c_ <- list()
  getindexl <- 0
  for (sdi in xscode) {
    fullxscode_i <- sdi
    getparm <- tail(strsplit(sdi, split = " ")[[1]], 1)
    parm_c <- getparm # c(parm_c, getparm)
    sdi <- gsub("[", "(", sdi, fixed = T)
    sdi <- gsub("]", ")", sdi, fixed = T)
    str_d <-
      regmatches(sdi, gregexpr("(?=\\().*?(?<=\\))", sdi, perl = T))[[1]]
    # this to remove dim after par name
    parm_cm <- parm_c
    if (grepl("[", parm_cm, fixed = T)) {
      parm_cm <- gsub("[", "(", parm_cm, fixed = T)
      parm_cm <- gsub("]", ")", parm_cm, fixed = T)
      parm_cm_ <-
        regmatches(parm_cm,
                   gregexpr("(?=\\().*?(?<=\\))", parm_cm, perl = T))[[1]]
      parm_cm2 <- gsub(parm_cm_, "", parm_cm, fixed = T)
      parm_c <- parm_cm2
    } else if (!grepl("[", parm_cm, fixed = T)) {
      parm_c <- parm_c
    }
    parm_c <- parm_c
    # sdi_c <- c(sdi_c, str_d)
    sdi <- str_d
    sdi <- gsub("(", "", sdi, fixed = T)
    sdi <- gsub(")", "", sdi, fixed = T)
    str_d <- sdi
    #
    sdi <- gsub("(", "", sdi, fixed = T)
    sdi <- gsub(")", "", sdi, fixed = T)
    str_d <- gsub("[[:space:]]", "", sdi)
    str_d_ <- strsplit(str_d, ",") %>% unlist()
    
    # for student_nu distribution parameter
    if (grepl("sd_nu", parm_c, fixed = T)) {
      set_value <- 3
    }
    
    if (grepl("sd", parm_c, fixed = T)) {
      set_value <- sd_value
    }
    if (grepl("z", parm_c, fixed = T)) {
      set_value <- z_value
    }
    if (grepl("L", parm_c, fixed = T)) {
      set_value <- L_value
    }
    
    # to exclude student_nu distribution parameter
    if (grepl("sd", parm_c, fixed = T)) {
      # to exclude student_nu distribution parameter
      if (!grepl("sd_nu", parm_c, fixed = T)) {
        if (length(str_d_) == 1) {
          dim1 <- xsdata[[str_d_[1]]]
          out <- rep(set_value, dim1)
        } else if (length(str_d_) == 2) {
          dim1 <- xsdata[[str_d_[1]]]
          dim2 <- xsdata[[str_d_[2]]]
          out <- matrix(set_value, dim1, dim2)
        }
        if (is.vector(out)) {
          out <- array(out, dim = length(out))
        }
      } # if(!grepl("sd_nu", parm_c, fixed = T)) {
    } # if(grepl("sd", parm_c, fixed = T)) {
    
    if (grepl("z", parm_c, fixed = T)) {
      if (length(str_d_) == 1) {
        dim1 <- xsdata[[str_d_[1]]]
        out <- rep(set_value, dim1)
      } else if (length(str_d_) == 2) {
        dim1 <- xsdata[[str_d_[1]]]
        dim2 <- xsdata[[str_d_[2]]]
        out <- matrix(set_value, dim1, dim2)
      }
      
      if (is.vector(out)) {
        out <- array(out, dim = length(out))
      }
      
      if (ncol(out) == 1) {
        out <- t(out)
      }
    } # if(grepl("z", parm_c, fixed = T)) {
    
    
    if(parameterization == 'cp') {
      # if (grepl("r", parm_c, fixed = T)) {
      if (grepl("r", parm_c, fixed = T) & 
          !grepl("Lrescor", parm_c, fixed = T)) {
        
        if (length(str_d_) == 1) {
          dim1 <- xsdata[[str_d_[1]]]
          out <- rep(set_value, dim1)
        } else if (length(str_d_) == 2) {
          dim1 <- xsdata[[str_d_[1]]]
          dim2 <- xsdata[[str_d_[2]]]
          out <- matrix(set_value, dim1, dim2)
        }
       
        if (is.vector(out)) {
          out <- array(out, dim = length(out))
        }
        
        if (ncol(out) == 1) {
          out <- t(out)
        }
       
      } # if(grepl("r", parm_c, fixed = T)) {
    } # if(parameterization == 'cp') {
    
   
    
    if (grepl("L", parm_c, fixed = T)) {
      if (length(str_d_) == 1) {
        dim1 <- xsdata[[str_d_[1]]]
        out <- matrix(set_value, dim1, dim1)
        diag(out) <- 1
      }
      if (length(str_d_) == 2) {
        dim1 <- xsdata[[str_d_[1]]] 
        dim2 <- xsdata[[str_d_[2]]] 
        # outxx <- array(0, dim = c(dim2, dim1, dim1))
        outxx <- array(0, dim = c(dim1, dim2, dim2)) # 13 12 23
        for (i in 1:dim(outxx)[2]) {
          outxx[, i, i] <- 1
        }
        out <- outxx
      }
    } # if(grepl("L", parm_c, fixed = T)) {
    sd_init_list_c_[[parm_c]] <- out
  }
  sd_init_list_c_
}




#' An internal function to get initials from pathfinder algorithm
#'
#' @param pthf A pathfinder draw object
#' @param ndraws An integer specifying the number of paths
#' @param init_structure A list of initial to set the appropriate dimensions of 
#'  inits returned from the pathfinder
#' @param variables A character vector specifying the variable names
#' @param model An object of class \code{bgmfit}
#' @param compile_init_model_methods A logical to indicate whether to compile
#' \code{init_model_methods()}
#' @param verbose A logical
#' @return A named list.
#' @keywords internal
#' @noRd
#'
get_pathfinder_init <- function(pthf = NULL, 
                                ndraws = 1, 
                                init_structure = NULL,
                                variables = NULL,
                                model = NULL,
                                compile_init_model_methods = NULL,
                                verbose = FALSE) {
  
  if(is.null(pthf) & is.null(model)) 
    stop('Specify at least pthf or model')
  if(!is.null(pthf) & !is.null(model)) 
    stop('Specify either least pthf or model')
  
  if(!is.null(model)) {
    mod <- attr(model$fit, "CmdStanModel")
    dat <- brms::standata(model)
    ini <- model$stan_args$init
    threads <- model$threads$threads
    if(!is.null(threads)) {
      pth1 <- mod$pathfinder(data = dat, init = ini,  num_threads = threads)
    } else if(!is.null(threads)) {
      pth1 <- mod$pathfinder(data = dat, init = ini)
    }
    if(is.null(compile_init_model_methods)) compile_init_model_methods <- TRUE
  } else {
    if(is.null(compile_init_model_methods)) compile_init_model_methods <- FALSE
    pthf <- pthf
  }
  
  if(compile_init_model_methods) {
    pth1$init_model_methods()
  }
  
  lp__ <- NULL;
  lp_approx__ <- NULL;
  lw <- NULL;
  .draw <- NULL;
  lw.x <- NULL;
  lp__ <- NULL;
  
  
  as_inits <- function(draws, variable=NULL, ndraws=ndraws) {
    ndraws <- min(posterior::ndraws(draws),ndraws)
    if (is.null(draws)) {variable = variables(draws)}
    draws <- draws  %>%  posterior::as_draws_matrix()
    inits <- lapply(1:ndraws,
                    function(drawid) {
                      sapply(variable,
                             function(var) {
                               as.numeric(posterior::subset_draws(draws, 
                                                                  variable=var, 
                                                                  draw=drawid))
                             })
                    })
    if (ndraws==1) { inits[[1]] } else { inits }
  }
  
  if (is.null(variables)) {
    # set variable names to be list of parameter names
    variables <- names(pthf$variable_skeleton(transformed_parameters = FALSE,
                                              generated_quantities = FALSE))
  }
  draws <- pthf$draws(format="df")
  draws <- draws %>% 
    posterior::mutate_variables(lw = lp__ - lp_approx__)
  ndist <- dplyr::n_distinct(posterior::extract_variable(draws,"lw"))
  if (ndist < ndraws) {
    stop(paste0("Not enough distinct draws (", ndist, ") to create inits."))
  }
  if (ndist < 0.95*posterior::ndraws(draws)) {
    # Resampling has been done in Stan, compute weights for distinct draws
    #these are now non Pareto smoothed as we have lost the original information
    draws <- draws %>% 
      dplyr::group_by(lw) %>% 
      dplyr::summarise(.draw=min(.draw)) %>% 
      dplyr::left_join(draws, by = ".draw") %>% 
      posterior::as_draws_df() %>% 
      posterior::mutate_variables(lw = lw.x,
                                  w = exp(lw-max(lw)))
  } else {
    # Resampling was not done in Stan, compute Pareto smoothed weights
    draws <- draws %>% 
      posterior::mutate_variables(w=posterior::pareto_smooth(exp(lw-max(lw)), 
                                                             tail="right"))
  }
  
  out <- 
    draws %>% 
    posterior::weight_draws(weights=posterior::extract_variable(draws,"w"), 
                            log=FALSE) %>% 
    posterior::resample_draws(ndraws=ndraws, method = "simple_no_replace") %>% 
    as_inits(variable=variables, ndraws=ndraws)
  
  
  if(!is.null(init_structure)) {
    path_inits <- out # inits_pathfinder
    init_str_x <- init_structure # fit_m$stan_args$init[[1]]
    for (stri in names(init_str_x)) {
      if(is.array( init_str_x[[stri]] )) {
        if(!is.null(path_inits[[stri]])) {
          path_inits[[stri]] <- array(path_inits[[stri]], dim = dim(init_str_x[[stri]]) )
        }
      } else if(is.vector( init_str_x[[stri]] )) {
      } else if(is.numeric( init_str_x[[stri]] )) {
      }
      out <- path_inits
    }
  } # if(!is.null(init_structure)) {
  
  return(out)
}




# check if a str obj is actually numeric
# @description check if a str obj is actually numeric
# https://stackoverflow.com/questions/13638377/test-for-numeric-elements-in-a-character-string
# is.numeric.like -> changed to check_is_numeric_like -> is called in bsitar.R only
# Thid because of notes in rmd check which says 
# Mismatches for apparent methods not registered
# This perhaps because is.numeric.like sound like is.numeric
#' @param x a str vector, or a factor of str vector, or numeric vector. x will
#'   be coerced and trimws.
#' @param na.strings case sensitive strings that will be treated to NA.
#' @param naAsTrue whether NA (including actual NA and na.strings) will be
#'   treated as numeric like
#' @return a logical vector (vectorized).
#' @note Using regular expression
#' \cr TRUE for any actual numeric c(3,4,5,9.9) or c("-3","+4.4",
#' "-42","4L","9L",   "1.36e4","1.36E4",    NA, "NA", "","NaN", NaN):
#' \cr positive or negative numbers with no more than one decimal c("-3","+4.4")
#' OR
#' \cr positive or negative integers (e.g., c("-42","4L","39L")) OR
#' \cr positive or negative numbers in scientific notation c("1.36e4","1.36E4")
#' \cr NA, or na.strings
#' @keywords internal
#' @noRd
#'
check_is_numeric_like <- function(x, 
                            naAsTrue = TRUE, 
                            na.strings = 
                              c('','.','NA','na','N/A','n/a','NaN','nan')
                            ){
  x = trimws(x,'both')
  x[x %in% na.strings] = NA
  # https://stackoverflow.com/a/21154566/2292993
  result = grepl("^[\\-\\+]?[0-9]+[\\.]?[0-9]*$|^[\\-\\+]?[0-9]+[L]?$|^[\\-\\+]?[0-9]+[\\.]?[0-9]*[eE][0-9]+$",x,perl=TRUE)
  if (naAsTrue) result = result | is.na(x)
  return((result))
}




#' Title
#'
#' @param data data frame
#' @param id id
#' @param outcome outcome
#' @param time time variable such as age
#' @param timeval the value of time beyound which changes made
#' @param nset number of last time point to be flattened
#' @param inc increment for last n time point.
#'
#' @keywords internal
#' @noRd
#' 
flattten_last_time <-
  function (data,
            id,
            outcome,
            time,
            timeval = NULL,
            nset = 1,
            inc = NULL) {
    occtemp <- NULL;
    temdata <-
      data %>% dplyr::group_by_at(id) %>% 
      dplyr::mutate(`:=`("occtemp",
                         dplyr::row_number())) %>% 
      dplyr::mutate(`:=`("nocctemp",
                         max(.data[["occtemp"]])))
    setseq <- seq(1, nset, 1) - 1
    inc <- rev(inc)
    if (is.null(inc)) {
      inc <- 0
      inc <- rep(inc, length(nset))
    }
    else if (length(inc) == 1) {
      inc <- rep(inc, nset)
    }
    else if (length(inc) != nset) {
      stop("lenhth of 'inc' must be either 1 or same as the 'nset'")
    }
    if (is.null(timeval))
      settime <- min(.data[[time]])
    else
      settime <- timeval
    if (nset == 1) {
      j = 0
      for (i in setseq) {
        j <- j + 1
        addinc <- inc[j]
        temdata <- temdata %>% dplyr::group_by_at(id) %>%
          dplyr::mutate(`:=`(
            !!base::as.symbol(outcome),
            dplyr::if_else(
              .data[[time]] > settime & occtemp ==
                max(.data[["occtemp"]]) - i,
              (.data[[outcome]] +
                 addinc),
              .data[[outcome]]
            )
          ))
      }
    }
    else {
      j = 0
      for (i in setseq) {
        j <- j + 1
        addinc <- inc[j]
        temdata <- temdata %>% dplyr::group_by_at(id) %>%
          dplyr::mutate(`:=`(
            !!base::as.symbol(outcome),
            dplyr::if_else(
              .data[[time]] > settime & occtemp ==
                max(.data[["occtemp"]]) - i,
              cummax(.data[[outcome]] +
                       addinc),
              .data[[outcome]]
            )
          ))
      }
    }
    temdata2 <- temdata %>% dplyr::select(-c("occtemp",
                                             "nocctemp"))
    return(temdata2)
  }




#' Save list of ggplot2 objects to single pdf
#'
#' @param list A list of ggplot2 objects.
#' @param filename A character string to name the pdf filr.
#'
#' @return Invisible NULL.
#' @keywords internal
#' @noRd
#'
#' @examples
#' #plot histogram of each numeric variable in iris
#' list_iris = map(names(iris[-5]), ~ggplot(iris, aes_string(.)) + geom_histogram())
#' #save to a single pdf
#' GG_save_pdf(list_iris, "test.pdf")
GG_save_pdf = function(list, filename, 
                       width = 10, height = 7,
                       onefile = TRUE, compress = TRUE) {
  #start pdf
  filename <- paste0(filename, ".", "pdf")
  grDevices::pdf(filename, width = width, height = height,
      onefile = onefile, compress = compress)
  #loop
  for (p in list) {
    print(p)
  }
  #end pdf
  grDevices::dev.off()
  invisible(NULL)
}




#' An internal to check if 'bsitar' argument such as xfun is set or NULL
#'
#' @param x A symbol or a character string 
#'
#' @return A logical TRUE/FALSE
#' @keywords internal
#' @noRd
#'
check_if_arg_set <- function(x) {
  if(is.null(x)) {
    set_x <- FALSE
  } else if(is.character(x)) {
    if(x == "NULL") {
      set_x <- FALSE
    } else {
      set_x <- TRUE
    }
  } else if(is.list(x)) {
    if(length(x) == 1) {
      if(is.null(x[[1]])) {
        set_x <- FALSE
      } else {
        set_x <- TRUE
      }
    } else if(length(x) > 1) {
      if(is.null(x[[1]][1])) {
        set_x <- FALSE
      } else {
        set_x <- TRUE
      }
    }
  }
  return(set_x)
}



#' Check if 'bsitar' argument such as xfun is set or NULL
#'
#' @param x A character string 
#' @param splitat A character at which string to be split default (\code{NULL}) 
#'
#' @return A logical TRUE/FALSE
#' @keywords internal
#' @noRd
#'
remove_between_first_last_parnth <- function(x, splitat = NULL) {
  a <- sub("^.*?\\(", "", x)
  b <- sub(")\\s*$", "", a)
  if(!is.null(splitat)) {
    b <- strsplit(b, splitat)[[1]]
  }
  b
}


# borrowed from sitar::ifun -> inverse_transform

#' An internal function to set up transformation for \code{bsitar} model
#'
#' @param expr A function \code{expression}.
#'
#' @param verbose Logical to print transformation steps.
#'   
#' @param envir A logical (default \code{TRUE})
#'
#' @return A function \code{expression}.
#'
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
inverse_transform <- function (expr, verbose = FALSE, envir = NULL) {
  vars <- function(expr) {
    (av <- all.vars(expr, unique = FALSE))[grep("^pi$", av, 
                                                invert = TRUE)]
  }
  recur <- function(fun, funinv = quote(x), verbose = verbose) {
    fun <- as.expression(fun)[[1]]
    if (verbose) {
      print(fun)
      print(funinv)
      cat("---\n")
    }
    while (length(fun) > 1 && fun[[1]] == as.name("(")) fun <- fun[[2]]
    if (!is.name(fun)) {
      x1 <- which(vapply(fun, function(f) length(vars(f)) == 
                           1, TRUE))[-1]
      if (grepl("pi", fname <- as.name(fun[[1]]))) {
        fun[[1]] <- as.name(sub("pi", "", fname))
        f <- quote(x * pi)
        f[[2]] <- fun[[2]]
        fun[[2]] <- f
      }
      nf <- which(vapply(fns, function(f) f[[1]] == fun[[1]] && 
                           length(f) == length(fun) && f[[x1]] == "x", TRUE))
      if (length(nf) == 0) 
        stop(paste("unrecognised name:", deparse(fun[[1]])))
      if (length(nf) > 1 && length(fun) == 3) {
        nft <- which(vapply(fns[nf], function(f) f[[5 - 
                                                      x1]] == fun[[5 - x1]], 
                            TRUE))
        if (length(nft)) 
          nf <- nf[nft]
      }
      nf <- nf[[1]]
      fn2 <- fns[[nf - 1 + 2 * (nf%%2)]]
      x2 <- which(as.list(fn2) == "x")
      if (length(fn2) == 3) {
        f <- function(n) {
        }
        body(f) <- fn2[[5 - x2]]
        fn2[[5 - x2]] <- f(eval(fun[[5 - x1]]))
      }
      fun <- fun[[x1]]
      fn2[[x2]] <- funinv
      funinv <- fn2
      if (!is.name(fun)) {
        results <- recur(fun, funinv, verbose = verbose)
        fun <- results$fun
        funinv <- results$funinv
      }
    }
    return(list(funinv = funinv, fun = fun))
  }
  fns <- quote(c(x + n, x - n, x * n, x/n, x^n, x^(1/n), sqrt(x), 
                 x^2, exp(x), log(x), expm1(x), log1p(x), n^x, log(x, 
                                                                   n), 
                 log10(x), 10^x, log2(x), 2^x, n + x, x - n, n - 
                   x, n - x, n * x, x/n, n/x, n/x, +x, +x, -x, -x, identity(x), 
                 identity(x), I(x), I(x), cos(x), acos(x), sin(x), asin(x), 
                 tan(x), atan(x), cosh(x), acosh(x), sinh(x), asinh(x), 
                 tanh(x), atanh(x)))
  fns[[1]] <- NULL
  varname <- vars(expr)
  if (length(varname) != 1) 
    stop("expression should contain just one instance of one name")
  fn <- function(x) {
  }
  body(fn) <- with(fns, recur(expr, verbose = verbose))$funinv
  attr(fn, "varname") <- varname
  fn
}




#' Check if any variable has all NAs
#'
#' @param data A \code{data.frame}.
#' @param factor_var A character string specifying the a factor variable name
#'   (default \code{NULL}).
#' @param envir An environment of evaluation (default \code{NULL}).
#' @param return A logical to indicate if data to be returned.
#'
#' @return A \code{data.frame} when \code{return = TRUE}, or \code{NULL} if
#'   \code{return = FALSE}.
#'   
#' @keywords internal
#' @noRd
#'
check_if_any_varibale_all_NA <- function(data, 
                                         factor_var = NULL, 
                                         envir = NULL,
                                         return = FALSE) {
  if(is.null(envir)) {
    envir <- parent.frame()
  }
  if(!is.null(factor_var)) {
    if(!is.character(factor_var)) {
      if(is.symbol(factor_var)) {
        factor_var <- deparse(substitute(factor_var))
      } else {
        stop("The factor_var must be a variable names")
      }
    }
    for (l in levels(data[[factor_var]])) {
      tempdata <- data %>% dplyr::filter(!! as.name(factor_var) == l)
      for (i in names(tempdata)) {
        if(all(is.na(tempdata[[i]]))) {
          stop("The variable '", i, "' contains all NA ", 
               "for subset '", l,  "'", 
               ". ", "Please check data")
        }
      }
    } # for (l in levels(data[[uvarby]])) {
  } else if(is.null(factor_var)) {
    tempdata <- data
    for (i in names(tempdata)) {
      if(all(is.na(tempdata[[i]]))) {
        stop("The variable '", i, "' contains all NA ", 
             ". ", "Please check data")
      }
    }
  }
  
  if(return) {
    return(data)
  } else {
    return(invisible(NULL))
  }
  
}


#' Check if any variable has all NAs
#'
#' @param data A \code{data.frame}.
#' @param variables A vector character string specifying the the variable names.
#' @param envir An environment of evaluation (default \code{NULL}).
#' @param return A logical to indicate if data to be returned.
#'
#' @return A \code{data.frame} when \code{return = TRUE}, or \code{NULL} if
#'   \code{return = FALSE}.
#'   
#' @keywords internal
#' @noRd
#'
check_variable_exists <- function(data, 
                                  variables, 
                                  envir = NULL, 
                                  return = FALSE) {
  for (j in variables) {
    if(is_emptyx(data[[j]]) | is.null(data[[j]]) ) {
      stop(paste0("variable '", j, "' not found in the data"))
    } 
  }
  if(return) {
    return(data)
  } else {
    return(invisible(NULL))
  }
}


#' Check if any variable has all NAs
#'
#' @param data A \code{data.frame}.
#' @param variables A vector character string specifying the the variable names.
#' @param envir An environment of evaluation (default \code{NULL}).
#' @param return A logical to indicate if data to be returned.
#'
#' @return A \code{data.frame} when \code{return = TRUE}, or \code{NULL} if
#'   \code{return = FALSE}.
#'   
#' @keywords internal
#' @noRd
#'
check_variable_numeric_exists <- function(data, 
                                          variables, 
                                          envir = NULL, 
                                          return = FALSE) {
  for (j in variables) {
    if(is.numeric(data[[j]])) {
      if(length(data[[j]]) == 0) {
        stop("The lenght of ", j, " is zero. Check your data and formual")
      }
    } 
  }
  if(return) {
    return(data)
  } else {
    return(invisible(NULL))
  }
}


#' Check and replace arg in simple function
#'
#' @param fun A \code{function}.
#' @param checkname A character string specifying the replacement name.
#'
#' @return A \code{function}.
#'   
#' @keywords internal
#' @noRd
#'
check_and_rename_funs_args_to_x <- function(fun, checkname = "x") {
  if(!is.character(checkname)) {
    stop("'checkname' must be a character")
  }
  # https://stackoverflow.com/questions/33850219/change-
  # argument-names-inside-a-function
  rep_vars <- function(expr, keyvals) {
    if (!length(expr)) return()
    if(is.symbol(expr)) {
      expr <- deparse(expr)
      expr <- gsub("\"", "", expr)
      expr <- str2expression(expr) # Imp, must be an expression
    }
    for (i in seq_along(expr)) {
      if (is.call(expr[[i]])) expr[[i]][-1L] <- Recall(expr[[i]][-1L], keyvals)
      if (is.name(expr[[i]]) && deparse(expr[[i]]) %in% names(keyvals))
        expr[[i]] <- as.name(keyvals[[deparse(expr[[i]])]])
    }
    return( expr )
  }
  
  formalArgs_names_in <- formalArgs(args(fun))
  
  deparse_fun_str <- deparse(fun)
  deparse_fun_str <- gsub_space(paste(deparse_fun_str, collapse = ""))
  
  if(length(formalArgs_names_in) > 1) {
    stop("Function '", deparse(fun) , "' must have only one argument")
  }
  
  if(grepl("\\{", deparse_fun_str) | grepl("}", deparse_fun_str)
  ) {
    stop("'", deparse_fun_str, "' must be a simple function without",
         "curly braces '{}'.",
         "\n ",
         " Examples: 'function(x)x' 'function(x)log(x)' function(x)log(x+1)")
  }
  
  if(grepl("return\\(", deparse_fun_str)) {
    stop("Function '", deparse_fun_str, "' must be a simple function without",
         "'return()'.",
         "\n ",
         " Examples: 'function(x)x' 'function(x)log(x)' function(x)log(x+1)")
  }
  
  ##############################################################
  if(formalArgs_names_in == checkname) {
    return(fun)
  } else {
    # https://stackoverflow.com/questions/44097516/r-paste-two-
    # string-with-an-equal-sign-between-it-stringa-stringb
    # newvals <- c("z" = "x")
    newvals <- checkname
    names(newvals) <- formalArgs_names_in
    newbod <- rep_vars(body(fun), newvals)
    # formals(fun) <- pairlist(x =bquote())
    formals(fun) <- ept(paste0("pairlist(", checkname, "=bquote())" ))
    body(fun) <- newbod
    # formals(fun) <- pairlist(x=bquote())
    # body(fun)    <- rep_vars(body(fun), newvals)
    return(fun)
  }
} 



#' Assign function 
#'
#' @param fun A \code{function}.
#' @param fun A character string
#' @param envir An environment
#'
#' @return A \code{function}.
#'   
#' @keywords internal
#' @noRd
#'
assign_function_to_environment <- function(fun, funname, envir = NULL) {
  if(is.null(envir)) {
    envir <- parent.frame()
  } else {
    envir <- envir
  }
  
  if(!is.null(fun)) {
    if(is.function(fun) & !is.primitive(fun)) {
      if(!is.primitive(fun)) {
        fun <- gsub_space(paste(deparse(fun), collapse = ""))
      }
    } else if(!is.character(fun)) {
      stop(paste0("The fun argument must be either a string ('log' or 'sqrt'),", 
                  "\n  ",
                  "or a function such as function(x)log(x)"))
    }
  }
  
  set_transform_draws      <- check_if_arg_set(fun)
  
  if (!set_transform_draws) {
    fun_eval <- function(x)x
    assign(funname, fun_eval, envir = envir)
  } else if (set_transform_draws) {
    if(fun == "log") {
      fun_eval <- function(x)log(x)
    } else if(fun == "sqrt") {
      fun_eval <- function(x)sqrt(x)
    } else  if(is.function(ept(fun))) {
      fun_eval <- ept(fun)
    } else {
      stop(paste0("The fun argument must be either a string ('log' or 'sqrt'),", 
                  "\n  ",
                  "or a function such as function(x)log(x)"))
    }
    assign(funname, fun_eval, envir = envir)
  }
  
  assign(funname, check_and_rename_funs_args_to_x(fun_eval, checkname = 'x'),
         envir = envir)
  
}






#' Extracted variable names from call
#'
#' @param model An object of class \code{bgmfit}.
#' @param arg A character string.
#' @param xcall A \code{mcall} object. Evaluated only when \code{model = NULL}.
#' @param envir An environment of evaluation (default \code{NULL}).
#'
#' @return A character string.
#'   
#' @keywords internal
#' @noRd
#'
extract_names_from_call <- function(model = NULL, 
                                    arg = NULL, 
                                    xcall = NULL, 
                                    envir = NULL) {
  if(!is.null(model)) {
    if(is.null(xcall)) xcall <- model$model_info$call.full.bgmfit
  } else if(is.null(model)) {
    if(is.null(xcall)) stop("specify either 'model' or 'xcall'")
  }
  
  if(is.null(arg)) {
    stop("specify 'arg'")
  } else if(!is.null(arg)) {
    if(!is.character(arg)) stop("'arg' must be a single character")
  }
  extracted <- xcall[[arg]]
  extracted <- toString(extracted)
  extracted <- strsplit(extracted, ",")[[1]]
  if(length(extracted) > 1) {
    extracted <- extracted[-1]
  }
  extracted <- gsub_space(extracted)
  extracted
} # end extract_names_from_mcall


#' An internal function to get the inverse transformation call 
#'
#' @param itransform A character string or \code{NULL}.
#' @param ... Additional argumenst. Currently ignored.
#'
#' @return A character string.
#'   
#' @keywords internal
#' @noRd
#'
get_itransform_call <- function(itransform, ...) {
  if(is.null(itransform)) {
    itransform_set <- "x"
  } else if(!is.null(itransform)) {
    if(is.logical(itransform)) {
      if(itransform) {
        itransform_set <- c('x', 'y', 'sigma')
      }
      if(!itransform) {
        itransform_set <- ""
      }
    } else if(itransform == "") {
      itransform_set <- ""
    } else {
      itransform_set <- itransform
    }
  } # if(is.null(itransform)) {
  itransform_set
}






#' An internal function to call function via eval()
#' 
#' @details
#' https://stackoverflow.com/questions/11054208/lapply-and-do-call-running-very-slow
#' 
#' @param what A language object (\code{function()}) to be called
#' @param args A list (\code{arguments})
#' @param quote A logical
#' @param envir A calling environment
#' @keywords internal
#' @return A object default of class inheretited from \code{what}
#' @noRd
#'
CustomDoCall <- function(what, args, quote = FALSE, envir = NULL){
  
  if(is.null(envir))
    envir <-   parent.frame()
  
  if (quote)
    args <- lapply(args, enquote)
  
  if (is.null(names(args))){
    argn <- args
    args <- list()
  }else{
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }
  
  
   get_class_what <- class(what)
  
  # if(is.character(class(what))) {
   if (get_class_what == "character"){
    if(is.character(what)){
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if(length(fn)==1) {
        get(fn[[1]], envir=envir, mode="function")
      } else {
        get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
      }
    }
    call <- as.call(c(list(what), argn))
  # } else if(is.function(class(what))) { 
  } else if (get_class_what == "function"){ 
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  # } else if(is.name(class(what))) { 
  } else if (get_class_what == "name"){
    call <- as.call(c(list(what, argn)))
  }
   
  args$verbose <- eval(args$verbose)
  
  eval(call, envir = args, enclos = envir)
}

# using with CustomDoCall()

#' An internal function to work with CustomDoCall()
#' 
#' @details sanitize the argument list
#' 
#' @param what A character string
#' @param arguments A list
#' @param check_formalArgs A logical
#' @param check_trace_back A logical
#' @param envir A logical or a list
#' @keywords internal
#' @return A object default of class inherited from \code{what}
#' @noRd
#'
sanitize_CustomDoCall_args <- function(what,
                                       arguments,
                                       check_formalArgs = NULL,
                                       check_formalArgs_exceptions = NULL,
                                       check_trace_back = NULL,
                                       envir = NULL) {
  
  if(is.null(envir)) {
    envir <- parent.frame()
  } else {
    envir <- envir
  }
  
  if(is.null(check_trace_back)) {
    set_trace_back <- rlang::trace_back() 
  } else {
    set_trace_back <- check_trace_back
  }
  
  if(!is.null(check_formalArgs)) {
    if(is.function(check_formalArgs)) {
      ownargs <- formalArgs(check_formalArgs)
    } else if(is.list(check_formalArgs)) {
      ownargs <- check_formalArgs
    } else {
      stop("check_formalArgs must be a function or a list")
    }
    
    for (i in setdiff(names(arguments), ownargs)) {
      if(!is.null(check_formalArgs_exceptions)) {
        if(!i %in% check_formalArgs_exceptions) {
          arguments[[i]] <- NULL
        }
      } else if(is.null(check_formalArgs_exceptions)) {
        arguments[[i]] <- NULL
      }
    } # for (i in setdiff(names(arguments), ownargs)) {
    
  } # if(!is.null(check_formalArgs)) {
  
  eval_CustomDoCall <- FALSE
  if(any(grepl(what, set_trace_back))) {
    eval_CustomDoCall <- TRUE
  }
  
  
  # remove empty argument 
  if(eval_CustomDoCall) {
    for (i in names(arguments)) {
      if(is.symbol(arguments[[i]])) {
        if(deparse(arguments[[i]]) == "") {
          arguments[[i]] <- NULL
        }
      }
    }
    
    for (i in names(arguments)) {
      if(!is.null(arguments[[i]])) {
        arguments[[i]] <- eval(arguments[[i]], envir = envir)
      }
    }
    
  } # if(eval_CustomDoCall) {
  
  return(arguments)
} # sanitize_CustomDoCall_args



# not using

#' An internal function to work with CustomDoCall()
#' 
#' @details
#' https://stackoverflow.com/questions/11054208/lapply-and-do-call-running-very-slow
#' 
#' @param scall A language object (\code{sys.call()}) to be called
#' @param return_tf A logical
#' @param return_str A logical
#' @param return_str A logical or a list
#' @keywords internal
#' @return A object default of class inheretited from \code{what}
#' @noRd
#'
check_CustomDoCall_fun <- function(scall, 
                                   return_tf = TRUE, 
                                   return_str = FALSE ) {
  check_CustomDoCall <- gsub_space(paste(deparse(scall), collapse = ""))
  eval_CustomDoCall <- FALSE
  if(grepl("CustomDoCall\\(", check_CustomDoCall)) {
    eval_CustomDoCall <- TRUE
    check_CustomDoCall <- regmatches(check_CustomDoCall, 
                                     gregexpr("(?<=\\().*?(?=\\))", 
                                              check_CustomDoCall, perl=T))[[1]]
    check_CustomDoCall <- strsplit(check_CustomDoCall[1], "\\(")[[1]][1]
  }
  
  if(return_tf & return_str) {
    out <- list()
    out[['eval_CustomDoCall']]  <- eval_CustomDoCall
    out[['check_CustomDoCall']] <- check_CustomDoCall
  } else if(return_tf) {
    out <- eval_CustomDoCall
  } else if(return_str) {
    out <- check_CustomDoCall
  } else {
    out <- NULL
  }
  return(out)
}




#' An internal function to replace part of string
#' 
#' @details An internal function to extract / replace an exact part of string
#' 
#' @param x A string
#' @param start A string
#' @param end A string
#' @param replace A string
#' @param extract A logica to indicate whether the matched pattern
#' @param cat_str A logica
#' @param exclude_start A logica
#' @param exclude_end A logica
#' 
#' @keywords internal
#' @return A string
#' @noRd
#'
replace_string_part <- function(x, 
                                start, 
                                end, 
                                replace = "",
                                extract = FALSE,
                                cat_str = FALSE,
                                exclude_start = FALSE,
                                exclude_end = FALSE) {
  
  if(!is.character(x))       
    stop("Argument 'x' must be a character string")
  if(!is.character(start))   
    stop("Argument 'start' must be a character string")
  if(!is.character(end))     
    stop("Argument 'end' must be a character string")
  if(!is.character(replace)) 
    stop("Argument 'replace' must be a character string")
  if(!is.logical(extract))   
    stop("Argument 'extract' must be a logical (TRUE/FALSE)")
  if(!is.logical(cat_str))   
    stop("Argument 'cat_str' must be a logical (TRUE/FALSE)")
  
  original_string  <- x
  start_pattern    <- start
  end_pattern      <- end
  replacement_text <- replace
  extract_pattern  <- extract
  catit            <- cat_str
  
  
  
  start_pattern_raw <- start_pattern
  end_pattern_raw   <- end_pattern
  # Helper function to escape special regex characters
  # These are the common special regex characters that need escaping.
  # Order matters for some, e.g., escape '\' before '['
  
  special_chars <- c("\\", ".", "+", "*", "?", "^", "$", "(", ")", "[", "]", "{", "}", "|")
  escape_regex <- function(string, special_chars) {
    # Use 'fixed = TRUE' to treat the search pattern as a literal string
    # when replacing, so we don't accidentally escape the escapes themselves.
    for (char in special_chars) {
      string <- gsub(char, paste0("\\", char), string, fixed = TRUE)
    }
    return(string)
  }
  
  start_pattern_escaped <- escape_regex(start_pattern_raw, special_chars)
  end_pattern_escaped   <- escape_regex(end_pattern_raw, special_chars)
  # Construct the full regex pattern with (?s) flag for DOTALL mode
  regex_pattern <- paste0(
    start_pattern_escaped,
    # (?s) makes the dot match newlines; .*? is non-greedy
    "(?s).*?", 
    end_pattern_escaped
  )
  
  
 
  # if(extract_pattern) {
  #   export_pattern <- paste0(start_pattern_escaped, end_pattern_escaped)
  #   clean_the_export_pattern <- function(string, special_chars) {
  #     for (char in special_chars) {
  #       string <- gsub(paste0("\\", char), char, string, fixed = TRUE)
  #     }
  #     return(string)
  #   } # end clean_the_export_pattern
  #   export_pattern <- clean_the_export_pattern(export_pattern, special_chars)
  #   return(export_pattern)
  # } # if(extract_pattern) {
  
  # Perform the replacement or extract
  if(extract_pattern) {
    match_info <- regexpr(regex_pattern, original_string, perl = TRUE)
    out_str    <- regmatches(original_string, match_info)
  } else if(!extract_pattern) {
    out_str <- gsub(regex_pattern, replacement_text, original_string, perl = TRUE)
  }
  
  
  if(exclude_start) {
    if(start != "") {
      out_str <- gsub(start, "", out_str, fixed = TRUE)
    }
  }
  if(exclude_end) {
    if(end != "") {
      out_str <- gsub(end, "", out_str, fixed = TRUE)
    }
  }
  
  if(catit) {
    out_str <- cat(out_str)
  }
  return(out_str)
}



#' An internal function to extrcat R function from stan function block
#' 
#' @details An internal function to inverse matrix 
#' This is exact same as extract_r_fun_from_scode that is used in utility 7
#' Here we use it for sigma_formula_manual but can be adapted for utility 7
#' 
#' @param xstaring A string
#' @param what A string
#' @param decomp A string
#' @param SplinefunxStan A string
#' @param SplinefunxR A string
#'  
#' @return A string
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
extract_r_fun_from_scode_sigma <- function(xstaring, 
                                           what = NULL, 
                                           decomp = NULL,
                                           SplinefunxStan,
                                           SplinefunxR) {
  if(is.null(xstaring)) {
    return(xstaring)
  }
  
  dparm_part_of_SplQc <- NULL;
  smat_intercept <- NULL;
  QR_Xmat <- NULL;
  QR_center <- NULL;
  QR_complete <- NULL;
  QR_flip <- NULL;
  QR_scale <- NULL;
  XR_inv_name <- NULL;
 
  
  add_separte_getknots_fun <- FALSE
  dparm_set_fixed_or_random <- FALSE
  
  getfunnamestr <- deparse(substitute(xstaring))
  xstaring <- gsub("[[:space:]]" , "", xstaring)
  xstaring <- gsub(";" , ";\n", xstaring)
  xstaring <- gsub("\\{" , "{\n", xstaring)
  xstaring <- gsub("}" , "}\n", xstaring)
  xstaring <- gsub("vector[N]" , "", xstaring, fixed = T)
  xstaring <- gsub("vector" , "", xstaring, fixed = T)
  xstaring <- gsub("int" , "", xstaring, fixed = T)
  if(!add_separte_getknots_fun) {
    xstaring <- gsub("knots=[" , "knots=c(", xstaring, fixed = T)
    xstaring <- gsub("]';" , ");", xstaring, fixed = T)
  }
  xstaring <- gsub("smat_normalize" , "normalize", xstaring, fixed = T)
  xstaring <- gsub("ercept" , "intercept", xstaring, fixed = T)
  xstaring <- gsub("[nknots-2]iknotsx;", "", xstaring, fixed = T)
  xstaring <- gsub("[2]bknotsx;", "", xstaring, fixed = T)
  xstaring <- gsub("Spl=matrix(0,N,SbasisN);" , "", xstaring, fixed = T)
  xstaring <- gsub("Spl[,1]=rep(1.0,N);" , "", xstaring, fixed = T)
  xstaring <- gsub("Spl[,1]=rep(0.0,N);" , "", xstaring, fixed = T)
  xstaring <- gsub("segment(knots,2,nknots-2);", "knots[2:(length(knots)-1)];", xstaring, fixed = T)
  xstaring <- gsub("append_row(head(knots,1),tail(knots,1));", "c(knots[1], knots[length(knots)]);", xstaring, fixed = T)
  xstaring <- gsub(SplinefunxStan , SplinefunxR, xstaring, fixed = T)
  xstaring <- gsub("matrix[N,nknots]Spl;" , "", xstaring, fixed = T)
  xstaring <- gsub("Spl[,1]=rep(1.0,N);" , "", xstaring, fixed = T)
  xstaring <- gsub("Spl[,1]=rep(0.0,N);" , "", xstaring, fixed = T)
  xstaring <- remove_spaces_and_tabs(xstaring)
  
  xstaring <- gsub("real" , "", xstaring, fixed = T)
  xstaring <- gsub(paste0("jp1;", "\n"), "", xstaring, fixed = T)
  xstaring <- gsub("rep_vector" , "rep", xstaring, fixed = T)
  xstaring <- gsub("rep_" , "rep", xstaring, fixed = T)
  xstaring <-
    gsub(
      "Xx[ia,ja]=(X[ia]-knots[ja]>0?X[ia]-knots[ja]:0);" ,
      "Xx[ia,ja]=ifelse(X[ia]-knots[ja]>0,X[ia]-knots[ja],0);",
      xstaring,
      fixed = T
    )
  xstaring <- gsub("num_elements" , "length", xstaring, fixed = T)
  xstaring <-
    gsub("matrix[N,SbasisN]Spl" ,
         "Spl=matrix(0,N,SbasisN)",
         xstaring,
         fixed = T)
  xstaring <-
    gsub("matrix[SbasisN,N]rcs" ,
         "rcs=matrix(0,SbasisN,N)",
         xstaring,
         fixed = T)
  xstaring <-
    gsub("matrix[N,SbasisN+1]Xx" ,
         "Xx=matrix(0,N, SbasisN+1)",
         xstaring,
         fixed = T)
  
  
  xstaring <-
    gsub("for(iain1:N)" , "for(ia in 1:N)", xstaring, fixed = T)
  xstaring <- gsub("for(jain1:(SbasisN+1))" ,
                   "for(ja in 1:(SbasisN+1))",
                   xstaring,
                   fixed = T)
  xstaring <- gsub(".*" , "*", xstaring, fixed = T)
  xstaring <- gsub("./" , "/", xstaring, fixed = T)
  funame__ <- strsplit(xstaring, "\\(")[[1]][1]
  xstaring <- gsub(funame__ , paste0(funame__, "<-function"),
                   xstaring, fixed = T)
  xstaring <- sub("//[^//]+$", "", xstaring)
  # To remove stanadlon ";
  xstaring <-
    gsub(paste0(";\n;\n", ""), ";\n", xstaring, fixed = T)
  xstaring <- gsub("[nknots]knots" , "knots", xstaring, fixed = T)
  
  
  if(dparm_set_fixed_or_random) {
    if(dparm_part_of_SplQc) {
      xstaring <- gsub("vector[N]dpredictor=Xm;" ,
                       "dpredictor=Xm;",
                       xstaring, fixed = T)
      xstaring <- gsub("append_col" ,
                       "cbind",
                       xstaring, fixed = T)
      xstaring <- gsub("rep_vector" ,
                       "rep",
                       xstaring, fixed = T)
      xstaring <- gsub("transpose" ,
                       "t",
                       xstaring, fixed = T)
      
      xstaring <- gsub("[:" ,
                       "[",
                       xstaring, fixed = T)
      
      xstaring <- gsub("matrix[N,ncol(XQ)]sfull_betas_temp" ,
                       "sfull_betas_temp=matrix(NA,N,ncol(XQ))",
                       xstaring, fixed = T)
      xstaring <- gsub("matrix[N,ncol(XQ)]sfull_matrix_temp" ,
                       "sfull_matrix_temp=matrix(NA,N,ncol(XQ))",
                       xstaring, fixed = T)
      
      xstaring <- gsub("sfull_betas_temp=sfull_matrix_temp*" ,
                       "sfull_betas_temp=sfull_matrix_temp %*% ",
                       xstaring, fixed = T)
    } # if(dparm_part_of_SplQc) {
  } # if(dparm_set_fixed_or_random) {
  
  
  if(dparm_set_fixed_or_random) {
    if(dparm_part_of_SplQc) {
      if(smat_intercept == 0) {
        xstaring <- gsub("matrix[N,SbasisN+1]Spl;" ,
                         "",
                         xstaring, fixed = T)
        xstaring <- gsub("matrix[N,SbasisN+1]SplQRd0;" ,
                         "",
                         xstaring, fixed = T)
        xstaring <- gsub("matrix[N,SbasisN+1]SplQRd1;" ,
                         "",
                         xstaring, fixed = T)
      } else if(smat_intercept == 1) {
        xstaring <- gsub("matrix[N,SbasisN+1+1]Spl;" ,
                         "",
                         xstaring, fixed = T)
        xstaring <- gsub("matrix[N,SbasisN+1+1]SplQRd0;" ,
                         "",
                         xstaring, fixed = T)
        xstaring <- gsub("matrix[N,SbasisN+1+1]SplQRd1;" ,
                         "",
                         xstaring, fixed = T)
      }
    } # if(dparm_part_of_SplQc) {
  } # if(dparm_set_fixed_or_random) {
  
  
  xstaring <- gsub("vector[N]QRdbeta=d" , "QRdbeta=d", xstaring, fixed = T)
  
  # add QR
  # make QR chnages
  if (!is.null(decomp)) {
    if (decomp == 'QR') {
      if(getfunnamestr == "rcsfun_raw" |
         getfunnamestr == "spl_d0" |
         getfunnamestr == "spl_d1" |
         getfunnamestr == "spl_d2") {
        
        set_QR_decomp_R <- paste0("QR_decomp_R(", 
                                  "X=",        QR_Xmat,     "," ,
                                  "center=",   QR_center, "," ,
                                  "complete=", QR_complete, "," ,
                                  "flip=",     QR_flip,     "," ,
                                  "scale=",    QR_scale, ")")
        set_QR_Xmat <- "QRRinv"
        set_QR_decomp_R <- paste0(set_QR_Xmat, "=", set_QR_decomp_R)
        
        getQmat    <- paste0(set_QR_Xmat, "[[", "'Q'", "]]")
        getRmat    <- paste0(set_QR_Xmat, "[[", "'R'", "]]")
        getRinvmat <- paste0(set_QR_Xmat, "[[", "'Rinv'", "]]")
        
        getQmat    <- paste0("XQ", "=", getQmat)
        getRmat    <- paste0("XR", "=", getRmat)
        getRinvmat <- paste0(XR_inv_name, "=", getRinvmat)
        
        if(getfunnamestr == "rcsfun_raw") {
          Qc_str_name     <- paste0("Qc=Spl;", "\n")
          Qc_str_mat_name <- paste0("matrix[N,QK]Qc=Spl;", "\n")
        } else if(getfunnamestr == "spl_d0") {
          Qc_str_name     <- paste0("Qc=Spl;", "\n")
          Qc_str_mat_name <- paste0("matrix[N,QK]Qc=Spl;", "\n")
        } else if(getfunnamestr == "spl_d1") {
          Qc_str_name     <- paste0("Qc=SplQRd0;", "\n")
          Qc_str_mat_name <- paste0("matrix[N,QK]Qc=SplQRd0;", "\n")
        } else if(getfunnamestr == "spl_d2") {
          
        } 
        
        set_QR_decomp_str <- paste0(Qc_str_name, 
                                    set_QR_decomp_R, "\n", 
                                    getQmat, "\n", 
                                    getRmat, "\n", 
                                    getRinvmat)
        
        xstaring <- replace_string_part(x = xstaring, 
                                        start = Qc_str_mat_name,
                                        end =  "inverse(XR);", 
                                        replace = set_QR_decomp_str,
                                        extract = FALSE,
                                        cat_str = FALSE)
        
        
        xstaring <- gsub("sfull_betas=sfull_matrix*transpose(XR_inv)" ,
                         "sfull_betas=sfull_matrix %*% t(XR_inv)",
                         xstaring, fixed = T)
        
        xstaring <- gsub("Spl=matrix(0,N,SbasisN)QRd0" ,
                         "QRd0=matrix(0,N,SbasisN)",
                         xstaring, fixed = T)
        xstaring <- gsub("Spl=matrix(0,N,SbasisN)QRd1" ,
                         "QRd1=matrix(0,N,SbasisN)",
                         xstaring, fixed = T)
        xstaring <- gsub("matrix[N,SbasisN]sfull_matrix" ,
                         "sfull_matrix=matrix(0,N,SbasisN)",
                         xstaring, fixed = T)
        xstaring <- gsub("matrix[N,SbasisN]R_tall" ,
                         "R_tall=matrix(0,N,SbasisN)",
                         xstaring, fixed = T)
        
        xstaring <- gsub("matrix[N,SbasisN]sfull_betas" ,
                         "sfull_betas=matrix(0,N,SbasisN)",
                         xstaring, fixed = T)
        
        xstaring <- gsub("(SplQRd1*sfull_betas)*rep(1.0,QK)" ,
                         "Matrix::rowSums(SplQRd1 %*% sfull_betas)",
                         xstaring, fixed = T)
        
      } # if(getfunnamestr == "rcsfun_raw" |....
    } # if (decomp == 'QR') {
  } # if (!is.null(decomp)) {
  
  xstaring <- gsub("cols" , "ncol", xstaring, fixed = T)
  
  xstaring <- gsub("matrixXp" , "Xp", xstaring, fixed = T) # spfnameX
  
  # This needed because \code{bsp}, \code{msp}, and \code{isp} may have NULL
  # / numeric(0). This same approach based on the 
  # function checkgetiknotsbknots() is used in bsitar and other funs
  # Note that in Stan, when knots = NULL, the vector is empty []
  gsub_it <- "iknotsx=knots[2:(length(knots)-1)]"
  gsub_by <- "iknotsx=checkgetiknotsbknots(knots,'iknots')"
  xstaring <- gsub(gsub_it, gsub_by, xstaring, fixed = T)
  gsub_it <- "bknotsx=c(knots[1], knots[length(knots)])"
  gsub_by <- "bknotsx=checkgetiknotsbknots(knots,'bknots')"
  xstaring <- gsub(gsub_it, gsub_by, xstaring, fixed = T)
  
  # for sigma var func
  gsub_it <- "){"
  gsub_by <- ") {"
  xstaring <- gsub(gsub_it, gsub_by, xstaring, fixed = T)
  
  gsub_it <- "N=size(x);"
  gsub_by <- "N=length(x);"
  xstaring <- gsub(gsub_it, gsub_by, xstaring, fixed = T)
  
  gsub_it <- "out;"
  gsub_by <- "out=rep(NA,N);"
  xstaring <- gsub(gsub_it, gsub_by, xstaring, fixed = T)
  
  gsub_it <- "for(iin1:N) {"
  gsub_by <- "for(i in 1:N) {"
  xstaring <- gsub(gsub_it, gsub_by, xstaring, fixed = T)
  
  xstaring
} # extract_r_fun_from_scode_sigma




#' An internal function to inverse matrix 
#' 
#' @details An internal function to inverse matrix 
#' 
#' @param x A numeric matrix of betas with dim(N, N)
#'  
#' @return A matrix
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
trysolveit <- function(x) {
  enverr. <- environment()
  assign('err.', FALSE, envir = enverr.)
  tryCatch(
    expr = {
      solvedit <- solve(x)
    },
    error = function(e) {
      assign('err.', TRUE, envir = enverr.)
    }
  )
  err. <- get('err.', envir = enverr.)
  if (err.) {
    solvedit <- MASS::ginv(x)
  } else {
    solvedit <- solvedit
  }
  return(solvedit)
}






#' An internal function to perform QR  decomposition
#' 
#' @details An internal function to perform QR decomposition 
#' 
#' @param X A numeric matrix
#' @param center A logical to indicate whether to center each column of matrix
#' before applying the QR decomposition 
#' @param complete A logical to indicate whether to get the full matrix or not
#' @param flip A logical to indicate whether to flip negative sign or not. 
#' @param scale A string to pass on the scaling factor. Typically it is 
#' \code{'sqrt(N-1)'} that is the default in Stan. 
#' 
#' @keywords internal
#' @return A list comprise of Q, R AND Rinv matrices
#' @noRd
#'
QR_decomp_R <- function(X, center = FALSE, complete = FALSE, 
                        flip = TRUE, scale = NULL) {
  
  if(is.null(scale)) {
    QR_scale_str <- sqrt(nrow(X) - 1)
  } else {
    if(is.numeric(scale)) {
      QR_scale_str <- scale
    } else if(is.character(scale)) {
      QR_scale_str <- scale
      QR_scale_str <- gsub("N", "nrow(X)", QR_scale_str, fixed = T)
      QR_scale_str <- ept(QR_scale_str)
    } else {
      stop("'scale' must be a numeric or string such as sqrt(N-1)")
    }
  }
  scale   <- round(QR_scale_str, 2)
  if(center) {
    for(i in 1:ncol(X)) {
      X[,i] <- X[,i] - mean(X[,i])
    }
  }
  qr_result_R <- qr(X)
  Q_X_R <- qr.Q(qr_result_R, complete = complete) * scale
  R_X_R <- qr.R(qr_result_R, complete = complete) / scale
  # Stan ensures diagonal elements of R are positive. R's qr() does not.
  # We need to check if any diagonal elements of R_X_R are negative and flip signs.
  if(flip) {
    diag_R_R <- diag(R_X_R)
    sign_flips <- sign(diag_R_R)
    sign_flips[sign_flips == 0] <- 1 # Treat zero as positive (no flip)
    # Apply sign flips to both Q and R from R
    # Multiply columns of Q by sign_flips and rows of R by sign_flips
    Q <- Q_X_R %*% diag(sign_flips)
    R <- diag(sign_flips) %*% R_X_R
  } else {
    Q <- Q_X_R
    R <- R_X_R
  }
  
  # if(complete) {
  #  # stop("use complete = 'FALSE' for Rinv")
  #   Rinv <- NULL
  # } else {
  #   Rinv <- solve(R)
  # }
  # gemini https://gemini.google.com/app/1af6d967c1f4e5b9
  Rinv <- trysolveit(R)
  list(Q = Q, R = R, Rinv = Rinv)
}



#' An internal function to add default options to nlf() and lf() functions
#' 
#' @details The \code{'add_default_args_to_nlf_lf'} function is used when 
#' using sigma_formula_manual for modelling sigma by mu and location scale
#' model. This \code{'add_default_args_to_nlf_lf'} is can also be used to 
#' extract co variate from the population part of the lf() form. The 
#' \code{'add_default_args_to_nlf_lf'} can be used for dpar_formula.
#' 
#' @param str A string
#' @param nys An integer to know whether to set the resp variable (if nys > 1) 
#' @param ysi A string to indicate the response variable name. Should be of 
#' length one because \code{'add_default_args_to_nlf_lf'} opetates within the
#' loop i.e., at the level of sigma_formula_manualsi
#' @param check A logical (\code{FALSE}). It should be set as \code{TRUE} when
#' using \code{'add_default_args_to_nlf_lf'} for dpar_formual
#' @param extract_nlpar A logical (\code{FALSE}). Set it as \code{TRUE} if 
#' want to extract \code{nlpar} from the nlf() form.
#' @param extract_covar A logical (\code{FALSE}). Set it as \code{TRUE} if 
#' want to extract co variate from the population part of the lf() form.
#' @param verbose A logical (\code{FALSE}).
#' @keywords internal
#' @return A string \code{extract_covar = FALSE} or a vector of string when 
#' \code{extract_covar = TRUE}
#' @noRd
#'
add_default_args_to_nlf_lf <- function(str, 
                                       nys, 
                                       ysi, 
                                       check = FALSE,
                                       extract_covar = FALSE, 
                                       extract_nlpar = FALSE, 
                                       data_varnames = NULL,
                                       verbose = FALSE) {
  
  if(extract_covar & extract_nlpar) {
    stop("specify either 'extract_covar' or 'extract_nlpar', not both")
  }
  
  if(length(ysi) > 1) {
    if(verbose)
      message("The length of ysi is > 1, only the first element is used")
  }
  str <- paste(gsub_space(str), collapse = "")
  str <- gsub("\"" , "'", str, fixed = T)
  temp_str <- gsub("\\+(nlf|lf)\\(", "###SPLIT###\\1(", str, perl = TRUE)
  split_result <- strsplit(temp_str, split = "###SPLIT###", fixed = TRUE)[[1]]
  split_result_c <- c()
  get_lf_covars_c <- c()
  get_nlf_nlpars_c <- c()
  for (i in 1:length(split_result)) {
    split_result_ith <- split_result[i]
    if (!is.null(split_result_ith)) {
      if(check) {
        if (grepl("^1$", split_result_ith)) {
          split_result_ith <- paste0("lf(", "sigma", "~", split_result_ith, ")")
        } else if (grepl("^~1", split_result_ith)) {
          split_result_ith <- paste0("lf(", "sigma", split_result_ith, ")")
        } else if (grepl("^sigma~1", split_result_ith)) {
          split_result_ith <- paste0("lf(", "", split_result_ith, ")")
        } else {
          split_result_ith <- split_result_ith
        }
      } # if(check) {
      if (grepl("lf\\(", split_result_ith) |
          grepl("nlf\\(", split_result_ith)) {
        if (!grepl("^lf\\(", split_result_ith) &
            grepl("nlf\\(", split_result_ith)) {
          lf_list <- c('flist',
                       # 'dpar',
                       'resp',
                       'loop')
          getnlpar <- ept(split_result_ith) [[1]][-2] %>% all.vars()
          getnlpar <- setdiff(getnlpar, data_varnames)
          get_nlf_nlpars_c <- c(get_nlf_nlpars_c, getnlpar)
        } else if (grepl("^lf\\(", split_result_ith) &
                   !grepl("^nlf\\(", split_result_ith)) {
          lf_list <- c('flist',
                       # 'dpar',
                       'resp',
                       'center',
                       'cmc',
                       'sparse',
                       'decomp')
          # Start extraction of covariates from lf() 
          get_lf_part_sigma_formula_manualsi <- ept(split_result_ith)[[1]]
          nthtau <- length(get_lf_part_sigma_formula_manualsi)
          get_lf_part_sigma_formula_manualsi_form <- 
            get_lf_part_sigma_formula_manualsi[[nthtau]]
          sigma_formulasi <- get_lf_part_sigma_formula_manualsi_form %>% deparse()
          sigma_formulasi <- paste0(gsub_space(sigma_formulasi), collapse = "")
          sigma_formulasi_check <- strsplit(sigma_formulasi, "+(", fixed = T)[[1]]
          sigma_formulasi <- sigma_formulasi_check[1]
          if(length(strsplit(sigma_formulasi, "~", fixed = T)[[1]]) > 1) {
            sigma_formulasi <- strsplit(sigma_formulasi, "~", fixed = T)[[1]][-1]
          } else {
            sigma_formulasi <- sigma_formulasi
          }
          sigma_formulasi <- paste0("~", sigma_formulasi)
          sigma_formulasi_covar <- all.vars(ept(sigma_formulasi))
          if(length(sigma_formulasi_covar) == 0) {
            sigma_formulasi_covar <- NULL
          }
          get_lf_covars_c <- c(get_lf_covars_c, sigma_formulasi_covar)
          # End extraction of covariates from lf() 
        } # end else if (grepl("^lf\\(", split_result_ith) &
        lf_list_c <- c()
        for (lf_listi in lf_list) {
          if (!grepl(lf_listi, split_result_ith)) {
            if (lf_listi == 'flist') {
              o. <- paste0(lf_listi, "=", 'NULL')
            } else if (lf_listi == 'dpar') {
              o. <- paste0(lf_listi, "=", paste0("'", 'sigma', "'"))
            } else if (lf_listi == 'center') {
              o. <- paste0(lf_listi, "=", 'NULL')
            } else if (lf_listi == 'cmc') {
              o. <- paste0(lf_listi, "=", 'NULL')
            } else if (lf_listi == 'sparse') {
              o. <- paste0(lf_listi, "=", 'NULL')
            } else if (lf_listi == 'decomp') {
              o. <- paste0(lf_listi, "=", 'NULL')
            } else if (lf_listi == 'resp') {
              # Imp: even when nys > 1, ysi should be a single resp varible 
              # That's why this is implemented at the level of _si
              if (nys > 1) {
                o. <- paste0(lf_listi, "=", paste0("'", ysi[1], "'"))
                # o. <- paste0(lf_listi, "=", 'NULL')
              } else {
                o. <- paste0(lf_listi, "=", 'NULL')
              }
            } else if (lf_listi == 'loop') {
              o. <- paste0(lf_listi, "=", 'FALSE')
            } else {
              o. <- o.
            }
            lf_list_c <- c(lf_list_c, o.)
          }
        }
        lf_list_c <- paste(lf_list_c, collapse = ",")
        if (lf_list_c != "")
          lf_list_c <- paste0(",", lf_list_c)
        split_result_ith <- gsub(")$", lf_list_c, split_result_ith)
        split_result_ith <- paste0(split_result_ith, ")")
      }
    } # if (!is.null(split_result_ith)) {
    split_result_ith <- paste(gsub_space(split_result_ith), collapse = "")
    split_result_c <- c(split_result_c, split_result_ith)
  } # for (i in 1:length(split_result)) {
  out <- paste(split_result_c, collapse = "+")
  out <- paste(gsub_space(out), collapse = "")
  if(extract_covar) {
    return(get_lf_covars_c)
  } else if(extract_nlpar) {
    return(get_nlf_nlpars_c)
  }
  return(out)
} # end add_default_args_to_nlf_lf <- function(str) {



#' An internal function to extract portion of string between pair of characters
#' 
#' @details An internal function to compute average prediction
#' 
#' @param str A string
#' @param start A string
#' @param end A string
#' @param verbose A logical (\code{FALSE}).
#' @keywords internal
#' @return A string
#' @noRd
#'
extract_between_specl_chars <- function(str, start, end, verbose = FALSE) {
  if(is.null(str) | length(str) == 0) {
    if(verbose)  print("String is NULL")
    return(invisible(NULL))
  } else if(!is.null(str)) {
    if(is.character(str)) {
      if(str == "NULL") return(invisible(NULL))
    }
    if(verbose)  print("String is NULL")
  } else if(is.null(str)) {
    if(verbose)  print("String is NULL")
    return(invisible(NULL))
  } 
  pattern <- paste0(start, "(.*?)\\", end)
  match_info <- regexpr(pattern, str, perl = TRUE)
  if (match_info[1] != -1) { # Check if a match was found
    full_match <- regmatches(str, match_info)
    extracted_string <- gsub("^~|\\*$", "", full_match, perl = TRUE)
    return(extracted_string)
  } else {
    if(verbose)  print("No match found.")
    return(invisible(NULL))
  }
} 


#' An internal function to check if all inner lists have the same length
#' 
#' @details An internal function
#' 
#' @param x A list
#' @keywords internal
#' @return A string
#' @noRd
#'
all_inner_lengths_equal_in_list <- function(x) {
  # Handle empty or single-element lists gracefully
  if (length(x) < 2) {
    return(TRUE)
  }
  # Get lengths of all inner lists and check for uniqueness
  length(unique(sapply(x, length))) == 1
}



#' An internal function to check if all inner lists elements are identical
#' 
#' @details An internal function
#' 
#' @param x A list
#' @keywords internal
#' @return A string
#' @noRd
#'
all_elements_identical_in_list <- function(x) {
  if (length(x) < 2) {
    return(TRUE)
  }
  # Compare every inner list to the first inner list
  all(sapply(x, identical, x[[1]]))
}



#' An internal function to Remove empty lines from code strings
#' 
#' @details Remove empty lines from code strings
#' 
#' @param x A string
#' @keywords internal
#' @return A string
#' @noRd
#'
remove_spaces_and_tabs <- function(x) {
  if(!is.null(x)) {
    x <- gsub("^ *|(?<= ) | *$", "", x, perl = TRUE)
    # '\\L\\1' converts first letter beyoind .* to lower
    # x <- gsub("(\\..*?[A-Z]|^[A-Z])", '\\L\\1', x, perl=T)
    x <- gsub("(\\..*?[A-Z]|^[A-Z])", '\\1', x, perl=T)
    x <- x[x != ""]
    x <- gsub("\\s*\n\\s*","\n",x) 
    xx <- x
  } else {
    xx <- x
  }
  return(xx)
}



#' An internal function to compute average prediction from individual curves
#' 
#' @details An internal function to compute average prediction
#' 
#' @param data A data frame
#' @param idvar A string
#' @param xvar A string
#' @param yvar NULL
#' @param name A string
#' @param xvar_orig TRUE
#' @param min_xvar NULL
#' @param max_xvar NULL
#' @param length.out NULL
#' @keywords internal
#' @return A list comprise of Q, R AND Rinv matrices
#' @noRd
#'
mean_curve_over_ids <- function(data, 
                                idvar, 
                                xvar, 
                                yvar, 
                                name = NULL,
                                xvar_orig = TRUE,
                                min_xvar = NULL, 
                                max_xvar = NULL, 
                                length.out = NULL) {
  
  
  if(is.null(length.out)) {
    set.length.out <- nrow(data)
  } else {
    set.length.out <- 100
  }
  
  if(is.null(min_xvar)) {
    set.min_xvar <- min(data[[xvar]])
  } else {
    set.min_xvar <- min_xvar
  }
  
  if(is.null(max_xvar)) {
    set.max_xvar <- max(data[[xvar]])
  } else {
    set.max_xvar <- max_xvar
  }
  
  min_age <- set.min_xvar
  max_age <- set.max_xvar
  
  if(xvar_orig) {
    common_age_x <- data[[xvar]] 
  } else {
    common_age_x <- seq(min_age, max_age, length.out = set.length.out) 
  }
  
  if(is.null(name)) {
    yvar_name <- "avg_pred"
  } else {
    yvar_name <- name
  }
  
  num_individuals <- length(unique(data[[idvar]]))
  
  interpolated_measurements <- matrix(NA, nrow = length(common_age_x), ncol = num_individuals)
  
  
  individual_data_list <- split(data, data[[idvar]])
  
  for (i in 1:num_individuals) {
    current_individual_data <- individual_data_list[[i]]
    # Use approx() for linear interpolation
    # This function naturally handles uneven 'x' (age) points and different 'length' of data
    interpolated_measurements[, i] <- stats::approx(
      x = current_individual_data[[xvar]],
      y = current_individual_data[[yvar]],
      xout = common_age_x,
      rule = 2 # Extrapolate with the closest known value if common_age_x goes beyond individual's observed age range
    )$y
  }
  
  mean_measurement_y <- rowMeans(interpolated_measurements, na.rm = TRUE)
  
  if(nrow(data) == length(mean_measurement_y)) {
    mean_curve_df <- data
    mean_curve_df[[yvar_name]] <- mean_measurement_y
  } else {
    # Create a data frame for the mean curve
    mean_curve_df <- data.frame(
      xvar = common_age_x,
      yvar_name = mean_measurement_y,
      idvar = "Mean Curve" # Assign a distinct ID for the mean curve
    )
  }
  return(mean_curve_df) 
}



#' An internal function to clean individual curves
#' 
#' @details An internal function to clean individual curves
#' 
#' @param data A data frame
#' @param idvar A character string specifying the individual identifier
#' @param xvar A character string specifying the age variable
#' @param yvar A character string specifying the outcome
#' @param carryforward A logical to indicate whether or not to carryforward
#' @param rate A numeric value to set the growth rate (mm/year)
#' 
#' @keywords internal
#' @return A data frame
#' @noRd
#'
get_clean_data <- function(data, 
                           idvar, 
                           xvar, 
                           yvar, 
                           carryforward = FALSE,
                           rate = 10.0) {
  names_in        <- colnames(data)
  max_growth_rate <- rate
  idvar           <- deparse(dplyr::ensym(idvar))
  xvar            <-  deparse(dplyr::ensym(xvar))
  yvar            <- deparse(dplyr::ensym(yvar))
  ept <- function (x) {
    eval(parse(text = x), envir = parent.frame())
  }
  z________idvar     <- z________xvar      <- z________yvar     <- NULL;
  y1_diff_prev       <- age_diff_prev      <- growth_rate_prev  <- NULL;
  y1_diff_next       <- age_diff_next      <- growth_rate_next  <- NULL;
  rate_extreme_prev  <- is_decreasing_prev <- rate_extreme_next <- NULL;
  is_decreasing_next <- is_outlier         <- copod_cleaned     <- NULL;
  data_with_checks <- data %>%
    dplyr::mutate(z________idvar := ept(idvar)) %>%
    dplyr::mutate(z________xvar := ept(xvar)) %>% 
    dplyr::mutate(z________yvar := ept(yvar)) %>% 
    dplyr::arrange(z________idvar, z________xvar) %>%
    dplyr::group_by(z________idvar) %>%
    dplyr::mutate(
      # Interval to previous
      y1_diff_prev = z________yvar - dplyr::lag(z________yvar),
      age_diff_prev = z________xvar - dplyr::lag(z________xvar),
      growth_rate_prev = y1_diff_prev / age_diff_prev,
      rate_extreme_prev = !is.na(growth_rate_prev) & 
        growth_rate_prev > max_growth_rate,
      is_decreasing_prev = !is.na(y1_diff_prev) & y1_diff_prev < 0,
      # Interval to next
      y1_diff_next = dplyr::lead(z________yvar) - z________yvar,
      age_diff_next = dplyr::lead(z________xvar) - z________xvar,
      growth_rate_next = y1_diff_next / age_diff_next,
      rate_extreme_next = !is.na(growth_rate_next) & 
        growth_rate_next > max_growth_rate,
      is_decreasing_next = !is.na(y1_diff_next) & y1_diff_next < 0,
      # Combine: if growth rate or decrease is outlier before or after
      is_outlier = rate_extreme_prev | is_decreasing_prev |
        rate_extreme_next | is_decreasing_next
    ) %>% dplyr::ungroup()
  
  ###################################################################
  zoo_na.locf.default_func <- function (object, 
                                        na.rm = TRUE, 
                                        fromLast, 
                                        rev, 
                                        maxgap = Inf, 
                                        rule = 2, ...) {
    
    zoo_na.locf.default_func <- na.approx <- na.trim <- NULL;
    ###
    zoo_na.locf0_fun <- function (object, 
                                  fromLast = FALSE, 
                                  maxgap = Inf, 
                                  coredata = NULL) {
      zoo_na.locf.default_func <- na.approx <- na.trim <- NULL;
      .fill_short_gaps <- function (x, fill, maxgap) {
        if (maxgap <= 0) 
          return(x)
        if (maxgap >= length(x)) 
          return(fill)
        naruns <- rle(is.na(x))
        naruns$values[naruns$lengths > maxgap] <- FALSE
        naok <- inverse.rle(naruns)
        x[naok] <- fill[naok]
        return(x)
      }
      
      if (is.null(coredata)) 
        coredata <- inherits(object, "ts") || inherits(object, 
                                                       "zoo") || inherits(object, "its") || inherits(object, 
                                                                                                     "irts")
      if (coredata) {
        x <- object
        object <- if (fromLast) 
          rev(coredata(object))
        else coredata(object)
      }
      else {
        if (fromLast) 
          object <- rev(object)
      }
      ok <- which(!is.na(object))
      if (is.na(object[1L])) 
        ok <- c(1L, ok)
      gaps <- diff(c(ok, length(object) + 1L))
      object <- if (any(gaps > maxgap)) {
        .fill_short_gaps(object, rep(object[ok], gaps), maxgap = maxgap)
      }
      else {
        rep(object[ok], gaps)
      }
      if (fromLast) 
        object <- rev(object)
      if (coredata) {
        x[] <- object
        return(x)
      }
      else {
        return(object)
      }
    } # end of zoo_na.locf0_fun
    
    ###
    
    L <- list(...)
    if ("x" %in% names(L) || "xout" %in% names(L)) {
      if (!missing(fromLast)) {
        stop("fromLast not supported if x or xout is specified")
      }
      return(na.approx(object, na.rm = na.rm, maxgap = maxgap, 
                       method = "constant", rule = rule, ...))
    }
    if (!missing(rev)) {
      warning("na.locf.default: rev= deprecated. Use fromLast= instead.")
      if (missing(fromLast)) 
        fromLast <- rev
    }
    else if (missing(fromLast)) 
      fromLast <- FALSE
    rev <- base::rev
    # na.locf0
    object[] <- if (length(dim(object)) == 0) 
      zoo_na.locf0_fun(object, fromLast = fromLast, maxgap = maxgap)
    else apply(object, length(dim(object)), zoo_na.locf0_fun, fromLast = fromLast, 
               maxgap = maxgap)
    if (na.rm) 
      na.trim(object, is.na = "all")
    else object
  } # end of zoo_na.locf.default_func
  
  ###################################################################
  
  if (carryforward) stop("carryforward = TRUE not working yet, set it as FALSE")
 
  if (carryforward) {
    data_cleaned <- data_with_checks %>%
      dplyr::group_by(z________idvar) %>%
      # Create a new column for the cleaned data
      dplyr::mutate(
        # Temporarily set outliers to NA to prepare for imputation
        copod_cleaned = dplyr::if_else(is_outlier, NA_real_, z________yvar),
        # Pass 1: Last Observation Carried Forward - zoo::na.locf
        copod_cleaned = zoo_na.locf.default_func(copod_cleaned, na.rm = FALSE),
        # Pass 2: Next Observation Carried Backward (to fix any leading NAs)
        copod_cleaned = zoo_na.locf.default_func(copod_cleaned, fromLast = TRUE, 
                                                 na.rm = FALSE)
      ) %>% 
      dplyr::mutate(z________yvar = copod_cleaned) %>% 
      dplyr::select(-copod_cleaned) %>% 
      dplyr::ungroup()
  } else if (!carryforward) {
    # Inspect flagged outliers
    flagged_obs <- data_with_checks %>% dplyr::filter(is_outlier == TRUE)
    data_cleaned <- data_with_checks %>% dplyr::filter(is_outlier == FALSE)
  }
  names_unwanted <- colnames(data_cleaned)
  names_unwanted <- setdiff(names_unwanted, names_in)
  data_cleaned <- data_cleaned %>% dplyr::select(-dplyr::all_of(names_unwanted))
  return(data_cleaned)
}





#' An internal function to clean individual curves using growthcleanr
#' 
#' @details An internal function to clean individual curves
#' 
#' @param data A data frame
#' @param idvar A character string specifying the individual identifier
#' @param xvar A character string specifying the age variable
#' @param yvar A character string specifying the outcome
#' @param xvar_unit A character string specifying the unit of xvar that need
#' to be converted in age days 
#' @param adjustcarryforward Adjust carry forward values in cleaned data
#' @param exclude_opt passed on to \code{adjustcarryforward} function. The 
#' \code{adjustcarryforward} uses absolute height velocity to identify values
#' excluded as carried forward values for re inclusion. \code{exclude_opt}
#' is a number from 0 to 3 indicating which option to use to handle strings
#' of carried-forwards: 
#'    0. no change.
#'    1. when deciding to exclude values, if we have a string of carried forwards,
#'    drop the most deviant value, and all CFs in the same string, and move on as
#'    normal.
#'    2. when deciding to exclude values, if the most deviant in a
#'    string of carried forwards is flagged, check all the CFs in that
#'    string from 1:N. Exclude all after the first that is flagged for
#'    exclusion when comparing to the Include before and after. Do not
#'    remove things designated as include.
#'    3. when deciding to exclude values, if the most deviant in a
#'    string of carried forwards is flagged, check all the CFs in that
#'    string from 1:N. Exclude all after the first that is flagged for
#'    exclusion when comparing to the Include before and after. Make sure
#'    remove things designated as include.
#' 
#' @param gender sex = 0 = male, sex = 1 = female, if NULL, then set as 0
#' 
#' @inherit growthcleanr::cleangrowth params
#' 
#' @keywords internal
#' @return A data frame
#' @noRd
#'
get_clean_data_growthcleanr <- function(data, 
                                        idvar, 
                                        xvar, 
                                        yvar,
                                        xvar_unit = 'year',
                                        param = "LENGTHCM",
                                        gender = NULL,
                                        recover.unit.error = FALSE,
                                        sd.extreme = 25,
                                        z.extreme = 25,
                                        lt3.exclude.mode = "default",
                                        height.tolerance.cm = 2.0,
                                        error.load.mincount = 2,
                                        error.load.threshold = 0.5,
                                        sd.recenter = NA,
                                        sdmedian.filename = "",
                                        sdrecentered.filename = "",
                                        include.carryforward = T,
                                        adjustcarryforward = T,
                                        exclude_opt = 0,
                                        ewma.exp = -1.5,
                                        ref.data.path = "",
                                        log.path = NA,
                                        parallel = FALSE,
                                        num.batches = NA,
                                        quietly = TRUE,
                                        adult_cutpoint = 20,
                                        weight_cap = Inf,
                                        adult_columns_filename = "",
                                        prelim_infants = FALSE) {
  # ‘HEIGHTCM’, ‘LENGTHCM’, or ‘WEIGHTKG’
  
  try(insight::check_if_installed(c("growthcleanr"), stop = FALSE, 
                                  prompt = FALSE))
  
  names_in        <- colnames(data)
  idvar           <- deparse(dplyr::ensym(idvar))
  xvar            <-  deparse(dplyr::ensym(xvar))
  yvar            <- deparse(dplyr::ensym(yvar))
  ept <- function (x) {
    eval(parse(text = x), envir = parent.frame())
  }
  z________idvar     <- z________xvar      <- z________yvar     <- NULL;
  z________param       <- z________sexvar      <- gcr_result  <- NULL;
  locf0_fun       <- fill_short_gaps      <- na.approx  <- NULL;
  na.trim  <- is_decreasing_prev <- rate_extreme_next <- NULL;
  # is_decreasing_next <- is_outlier         <- NULL;
  
  if(xvar_unit == 'year') {
    xvar_days_multi = 12 * 30.4375
  } else if(xvar_unit == 'month') {
    xvar_days_multi = 1 * 30.4375
  } else if(xvar_unit == 'day') {
    xvar_days_multi = 1 * 1
  }
  
  if(is.null(gender)) {
    gender <- 0
  }
  
  datafor_growthcleanr <- data %>%
    dplyr::mutate(z________idvar := ept(idvar)) %>%
    dplyr::mutate(z________xvar := ept(xvar)) %>% 
    dplyr::mutate(z________yvar := ept(yvar)) %>% 
    dplyr::mutate(z________param = param) %>% 
    dplyr::mutate(z________xvar = z________xvar * xvar_days_multi) %>%
    dplyr::mutate(z________sexvar = gender) 
  
  # prepare data as a data.table
  datafor_growthcleanr <- data.table:: as.data.table(datafor_growthcleanr)
  
  # set the data.table key for better indexing
  data.table::setkey(datafor_growthcleanr, 
                     z________idvar, 
                     z________xvar,
                     z________param)
  
  # generate new exclusion flag field using function
  cleaned_data <- datafor_growthcleanr[, 
                                       gcr_result := growthcleanr::cleangrowth(
                                         subjid = z________idvar, 
                                         param = z________param, 
                                         agedays = z________xvar, 
                                         sex = z________sexvar, 
                                         measurement = z________yvar,
                                         recover.unit.error = recover.unit.error,
                                         sd.extreme = sd.extreme,
                                         z.extreme = z.extreme,
                                         lt3.exclude.mode = lt3.exclude.mode,
                                         height.tolerance.cm = height.tolerance.cm,
                                         error.load.mincount = error.load.mincount,
                                         error.load.threshold = error.load.threshold,
                                         sd.recenter = sd.recenter,
                                         sdmedian.filename = sdmedian.filename,
                                         sdrecentered.filename = sdrecentered.filename,
                                         include.carryforward = include.carryforward,
                                         ewma.exp = ewma.exp,
                                         ref.data.path = ref.data.path,
                                         log.path = log.path,
                                         parallel = parallel,
                                         num.batches = num.batches,
                                         quietly = quietly,
                                         adult_cutpoint = adult_cutpoint,
                                         weight_cap = weight_cap,
                                         adult_columns_filename = adult_columns_filename,
                                         prelim_infants = prelim_infants
                                       )]
  
  
  if(adjustcarryforward) {
    cleaned_data <- growthcleanr::adjustcarryforward(subjid = cleaned_data$z________idvar,
                                                     param = cleaned_data$z________param,
                                                     agedays = cleaned_data$z________xvar,
                                                     sex = cleaned_data$z________sexvar,
                                                     measurement = cleaned_data$z________yvar,
                                                     orig.exclude = cleaned_data$gcr_result,
                                                     exclude_opt = exclude_opt,
                                                     sd.recenter = NA,
                                                     ewma.exp = -1.5,
                                                     ref.data.path = "",
                                                     quietly = TRUE,
                                                     minfactor = 0.5,
                                                     maxfactor = 2,
                                                     banddiff = 3,
                                                     banddiff_plus = 5.5,
                                                     min_ht.exp_under = 2,
                                                     min_ht.exp_over = 0,
                                                     max_ht.exp_under = 0.33,
                                                     max_ht.exp_over = 1.5)
  }
  
  # extract data limited only to values flagged for inclusion:
  only_included_data <- cleaned_data[gcr_result == "Include"]
  
  only_included_data <- only_included_data %>% data.frame()
  
  names_unwanted <- colnames(only_included_data)
  names_unwanted <- setdiff(names_unwanted, names_in)
  only_included_data <- only_included_data %>% 
    dplyr::select(-dplyr::all_of(names_unwanted))
  
  return(only_included_data)
}





#' Create rcs spline design matrix. 
#' @details used in bsitar
#' 
#' @keywords internal
#' @noRd
#' 
check_for_nan_inf <- function(x) {
  suppressWarnings({
    nan_inf <- FALSE
    if(is.infinite(x)) {
      nan_inf <- TRUE
    } else if(is.na(x)) {
      nan_inf <- TRUE
    } else if(is.nan(x)) {
      nan_inf <- TRUE
    }
  })
  nan_inf
}



#' Create rcs spline design matrix. 
#' @details used in bsitar
#' 
#' @keywords internal
#' @noRd
#' 
gkn <- function(x, df, bounds) {
  c(min(x) - bounds * (max(x) - min(x)),
    quantile(x, (1:(df - 1)) / df, na.rm = TRUE), # 28 01 2024
    max(x) +
      bounds * (max(x) - min(x)))
}



#' Create rcs spline design matrix. 
#' @details used in bsitar
#' 
#' @keywords internal
#' @noRd
#' 
eval_xoffset_bstart_args <- function(x, 
                                     y, 
                                     knots, 
                                     data, 
                                     eval_arg, 
                                     xfunsi, 
                                     arg = 'xoffset',
                                     smat,
                                     degree, 
                                     intercept, 
                                     derivs, 
                                     centerval, 
                                     normalize,
                                     preH, 
                                     sfirst, 
                                     sparse) {
  iknots <- checkgetiknotsbknots(knots, 'iknots')
  bknots <- checkgetiknotsbknots(knots, 'bknots')
    if (eval_arg == "mean") {
      eval_arg.o <- mean(data[[x]])
    } else if (eval_arg == "min") {
      eval_arg.o <- min(data[[x]])
    } else if (eval_arg == "max") {
      eval_arg.o <- max(data[[x]])
    } else if (eval_arg == "apv") {
      if(smat == 'rcs') {
        # mat_s <- make_spline_matrix(data[[x]], knots)
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_rcs_call(x = data[[x]], knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'nsp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_nsp_call(x = data[[x]], knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'nsk') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_nsk_call(x = data[[x]], knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'bsp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_bsp_call(x = data[[x]], knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'msp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_msp_call(x = data[[x]], knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'isp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_isp_call(x = data[[x]], knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      }
      lmform <- as.formula(paste0(y, "~1+", "mat_s"))
      lmfit <- lm(lmform, data = data)
      eval_arg.o <- sitar::getPeak(data[[x]],
                                   predict(smooth.spline(data[[x]],
                                                         fitted(lmfit)),
                                           data[[x]], deriv = 1)$y)[1]
      if(is.na(eval_arg.o)) {
        stop(arg, " specified as '", eval_arg, "' returned NA.",
             "\n ",
             " Please change ", arg,
             " argument to 'mean' or a numeric value.")
      }
    } else {
      eval_arg.o <- ept(eval_arg)
    }
    out <- as.numeric(eval_arg.o)
    out <- round(out, 3)
    return(out)
  }



#' Create rcs spline design matrix. 
#' @details used in bsitar
#' 
#' @keywords internal
#' @noRd
#' 
eval_xoffset_cstart_args <- function(x, 
                                     y, 
                                     knots, 
                                     data, 
                                     eval_arg, 
                                     xfunsi,
                                     smat,
                                     degree, 
                                     intercept, 
                                     derivs, 
                                     centerval, 
                                     normalize,
                                     preH, 
                                     sfirst, 
                                     sparse) {
  iknots <- checkgetiknotsbknots(knots, 'iknots')
  bknots <- checkgetiknotsbknots(knots, 'bknots')
    if (eval_arg == "pv") {
      # mat_s <- make_spline_matrix(data[[x]], knots)
      if(smat == 'rcs') {
        # mat_s <- make_spline_matrix(data[[x]], knots)
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_nsp_call(x = data[[x]], 
                             knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'nsp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_nsp_call(x = data[[x]], 
                             knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'nsk') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_nsk_call(x = data[[x]], 
                             knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'bsp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_bsp_call(x = data[[x]], 
                             knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'msp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_msp_call(x = data[[x]], 
                             knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      } else if(smat == 'isp') {
        # iknots <- knots[2:(length(knots)-1)]
        # bknots <- c(knots[1], knots[length(knots)])
        mat_s <- GS_nsk_call(x = data[[x]], 
                             knots = iknots, bknots = bknots, 
                             degree = degree,
                             intercept = intercept, 
                             derivs = derivs, 
                             centerval = centerval, 
                             normalize = normalize,
                             preH = preH,
                             sfirst = sfirst, 
                             sparse = sparse)
      }
      lmform <- as.formula(paste0(y, "~1+", "mat_s"))
      lmfit <- lm(lmform, data = data)
      eval_arg.o <- sitar::getPeak(data[[x]],
                                   predict(smooth.spline(data[[x]],
                                                         fitted(lmfit)),
                                           data[[x]], deriv = 1)$y)[2]
      if(is.na(eval_arg.o)) {
        stop("cstart specified as '", eval_arg, "' returned NA.",
             "\n ",
             " Please change cstart argument to 'mean' or a numeric value.")
      }
    } else {
      eval_arg.o <- ept(eval_arg)
    }
    out <- as.numeric(eval_arg.o)
    out <- round(out, 3)
    return(out)
  }




