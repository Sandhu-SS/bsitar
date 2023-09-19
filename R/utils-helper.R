

# load("z.RData")

expose_optimize_fit <- function(optimize_fit, 
                                subset_list = NULL, 
                                expose_function = T) {
  optimize_fit_models <-  optimize_fit
  
  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop("models must a numeric vector")
    optimize_fit_models <- optimize_fit_models[subset_list]
  } else {
    optimize_fit_models <- optimize_fit_models
  }
  
  m_list <- list()
  for (il in 1:length(optimize_fit_models)) {
    if(is.null(expose_function)) {
      m_list[[il]] <- optimize_fit_models[[il]]
    } else if(!is.null(expose_function)) {
      if(expose_function) {
        message("Exposing Stan function...", " (for model no. ", il, ")")
        m_list[[il]] <- 
          expose_bsitar_functions(optimize_fit_models[[il]], 
                                  optimize_fit_models[[il]]$bmodel)
      } else if(!expose_function) {
        m_list[[il]] <- optimize_fit_models[[il]]
      }
    }
  }
  m_list <- m_list[!sapply(m_list, is.null)]
  m_list
}

# elist <- expose_optimize_fit(optimize_fit2$models, subset_list = c(1, 6))



plot_optimize_fit <- function(optimize_fit, subset_list = NULL, what = "plot", 
                              expose_function = F, print= T, ...) {
  
  optimize_fit_models <-  optimize_fit # optimize_fit$models
  
  dots <- list(...)
  
  for (i in names(dots)) {
    if(!i %in% formalArgs(plot_bsitar.bsitar)) 
      stop("arguments must be be one of the following",
           "\n ",
           formalArgs(plot_bsitar.bsitar))
  }
  
  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop("models must a numeric vector")
    optimize_fit_models <- optimize_fit_models[subset_list]
  } else {
    optimize_fit_models <- optimize_fit_models
  }
  
  m_list <- list()
  for (il in 1:length(optimize_fit_models)) {
    if(is.null(expose_function)) {
      m_list[[il]] <- optimize_fit_models[[il]]
    } else if(!is.null(expose_function)) {
      if(expose_function) {
        message("Exposing Stan function...", " (for model no. ", il, ")")
        m_list[[il]] <- 
          expose_bsitar_functions(optimize_fit_models[[il]], 
                                  optimize_fit_models[[il]]$bmodel)
      } else if(!expose_function) {
        m_list[[il]] <- optimize_fit_models[[il]]
      }
    }
  }
  
  m_list <- m_list[!sapply(m_list, is.null)]
  
  nx <- function(.x, bx, args_) {
    message("Working on model no. ", .x)
    
    dots$model <- bx[[.x]]
    dots$... <- NULL
    
    if(is.null(what)) what <- 'plot'
    
    if(what == "plot") {
      out_ <- do.call(plot_bsitar, dots)
      title_ <- bx[[.x]]$model_info$optimization_info
      out_ <- out_ + ggplot2::labs(title = title_)
    }
    
    if(what == "gparameters") {
      out_ <- do.call(gparameters, dots)
    }
    
    if(!is.null(print)) {
      if(print) {
        print(out_)
      }
    }
    return(out_)
  }
  
  out <- purrr::map(1:length(m_list), ~nx(.x, m_list, args_))
  return(out)
}




transform.sec.axis <- function(primary, secondary, na.rm = TRUE) {
  from <- range(secondary, na.rm = na.rm)
  to   <- range(primary, na.rm = na.rm)
  zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
    if (length(x) == 1) {
      return(TRUE)
    }
    if (length(x) != 2)
      stop("x must be length 1 or 2")
    if (any(is.na(x))) {
      return(NA)
    }
    if (x[1] == x[2]) {
      return(TRUE)
    }
    if (all(is.infinite(x))) {
      return(FALSE)
    }
    m <- min(abs(x))
    if (m == 0) {
      return(FALSE)
    }
    abs((x[1] - x[2]) / m) < tol
  }
  rescale.numeric_ <-
    function(x,
             to = c(0, 1),
             from = range(x, na.rm = TRUE, finite = TRUE),
             ...) {
      if (zero_range(from) || zero_range(to)) {
        return(ifelse(is.na(x), NA, mean(to)))
      }
      (x - from[1]) / diff(from) * diff(to) + to[1]
    }
  forward <- function(x) {
    rescale.numeric_(x, from = from, to = to)
  }
  reverse <- function(x) {
    rescale.numeric_(x, from = to, to = from)
  }
  list(fwd = forward, rev = reverse)
}











get_gr_str_coef_id <- function(tsx, data) {
  tsx <- strsplit(tsx, "+(", fixed = T)[[1]] 
  tsx_id_w_or_wo_gr <- c()
  for (tsx_id_w_or_wo_gri in 1:length(tsx)) {
    tsx_id_w_or_wo_gr_get <- get_x_random2_asitis(tsx[tsx_id_w_or_wo_gri])
    tsx_id_w_or_wo_gr <- c(tsx_id_w_or_wo_gr, tsx_id_w_or_wo_gr_get)
  }
  tsx <- gsub("(", "", tsx, fixed = T)
  tsx <- gsub(")", "", tsx, fixed = T)
  # tsx <- gsub("0+", "", tsx, fixed = T)
  # tsx <- gsub("1+", "", tsx, fixed = T)
  tsx_c_coef  <- tsx_c_id    <- set_form_gr_it      <- list()
  set_ncov_it <- set_corr_it <- set_corr_true_false <- list()
  for (i in 1:length(tsx)) {
    tsx_c <- strsplit(tsx[i], "|", fixed = T)[[1]]
    set_corr_it_get <- tsx_c[2]
    tsx_c1 <- tsx_c[1]
    tsx_c3 <- tsx_c[3]
    if(!grepl("^~", tsx_c1)) tsx_c1 <- paste0("~", tsx_c1)
    if(grepl("^~0", tsx_c1)) set_form_0_gr <- TRUE
    if(grepl("^~1", tsx_c1)) set_form_0_gr <- FALSE
    set_form_gr <- tsx_c1
    # tsx_c1 <- strsplit(tsx_c1, "+", fixed = T)[[1]]
    tsx_c1_mat <- eval(parse(text = paste0(
      "model.matrix(",
      tsx_c1, ",data = data)"
    )))
    if (ncol(tsx_c1_mat) == 1)
      nlcov <- NULL
    else
      nlcov <- ncol(tsx_c1_mat) - 1
    nlcovoefnames <- colnames(tsx_c1_mat)
    nlcovoefnames <- gsub("\\(|)", "", nlcovoefnames)
    if(length(nlcovoefnames) > 1) tsx_c3 <- rep(tsx_c3, length(nlcovoefnames))
    tsx_c_coef[[i]] <- nlcovoefnames
    tsx_c_id[[i]]   <- tsx_id_w_or_wo_gr[i] # tsx_c3[1]
    set_form_gr_it[[i]]   <- set_form_gr
    if(set_form_0_gr) {
      set_ncov_it_get <- length(nlcovoefnames)
    }
    if(!set_form_0_gr) {
      if(length(nlcovoefnames) == 1) set_ncov_it_get <- NULL
      if(length(nlcovoefnames) > 1) set_ncov_it_get <- length(nlcovoefnames) - 1
    }
    if(is.null(set_ncov_it_get)) {
      set_corr_true_false[[i]] <- FALSE
    } else if(!is.null(set_ncov_it_get)) {
      if(set_corr_it_get == "") set_corr_true_false[[i]] <- FALSE
      if(set_corr_it_get != "") set_corr_true_false[[i]] <- TRUE
    }
    set_corr_it[[i]] <- set_corr_it_get
    set_ncov_it[[i]] <- set_ncov_it_get
  } # for (i in 1:length(tsx)) {
  
  # print(set_corr_it) %>% unlist()
  
  if(length(tsx_c_coef) != length(tsx_c_id)) 
    stop("coef and id length should be same")
  list(tsx_c_coef = tsx_c_coef, tsx_c_id = tsx_c_id, 
       set_form_gr_it = set_form_gr_it, set_ncov_it = set_ncov_it,
       set_corr_it = set_corr_it, set_corr_true_false = set_corr_true_false)
}




# This function get corr true false from | | syntax in _str 
# Used in prepare_formual 

get_str_corr_tf_function_new_better <- function(str_id_all_list, 
                                                str_corr_all_list, 
                                                str_corr_tf_all_list) {
  
  if(length(str_id_all_list) > 0 ) {
    id_corr_tf_bind <- cbind(unlist(str_id_all_list), 
                             unlist(str_corr_all_list), 
                             unlist(str_corr_tf_all_list))
    
    
    checkdi_c <- group_id_unique <- str_corr_tf <- c()
    
    for (checkdi in 1:length(str_id_all_list)) {
      checkdi_c <- c(checkdi_c,  length(str_corr_all_list[[checkdi]]) )
    }
    
    for (id_corr_tf_bind_1i in unique(id_corr_tf_bind[ , 1])) {
      temp_mat <- id_corr_tf_bind[which(id_corr_tf_bind == id_corr_tf_bind_1i),]
      if(!is.matrix(temp_mat)) temp_mat <- matrix(temp_mat) %>% t()
      get_check_id_mat      <- temp_mat[ , 1]
      get_check_corr_mat    <- temp_mat[ , 1]
      get_check_corr_tf_mat <- temp_mat[ , 1]
      if(any(grepl("TRUE", get_check_corr_tf_mat))) 
        str_corr_tf_single <- TRUE else str_corr_tf_single <- FALSE
      if(max(table(get_check_corr_mat)) > 1) {
        if(!str_corr_tf_single) str_corr_tf_single <- TRUE
      }
      str_corr_tf <- c(str_corr_tf, str_corr_tf_single)
      group_id_unique <- c(group_id_unique, id_corr_tf_bind_1i)
    }
    list(str_corr_tf = str_corr_tf, group_id_unique = group_id_unique)
  } else {
    list(str_corr_tf = NULL, group_id_unique = NULL)
  }
  
}




# This function will append priors to the above bpriors
# And, will out stanvar and inits to be added to stanvar_priors and initials

extract_prior_str_lv <- function(tempx) {
  if(!is.list(tempx) & !is.vector(tempx)) {
    out_prior_str <- tempx
  } else if(is.vector(tempx) & length(tempx) == 1 & 
            !grepl("c\\(", tempx) & 
            !grepl("list\\(", tempx) ) {
    out_prior_str <- tempx
  } else if(is.vector(tempx) & length(tempx) == 1 & 
            grepl("list\\(", tempx) ) {
    tempx2 <- str2lang(tempx)
    tempx2[[1]] <- NULL
    out_prior_str <- c()
    for (tempx2i in 1:length(tempx2)) {
      get_it_ <- tempx2[[tempx2i]] %>% deparse()
      get_it_ <- gsub("\"", "", get_it_)
      out_prior_str <- c(out_prior_str, get_it_)
    }
  } else if(is.list(tempx) | is.vector(tempx)  ) {
    tempx2 <- str2lang(tempx)
    tempx2[[1]] <- NULL
    out_prior_str <- c()
    for (tempx2i in 1:length(tempx2)) {
      get_it_ <- tempx2[[tempx2i]] %>% deparse()
      get_it_ <- gsub("\"", "", get_it_)
      out_prior_str <- c(out_prior_str, get_it_)
    }
  }
  out_prior_str
}







restore_paranthese_grgr_str_form <- function(strx) {
  restore_paranthese_grgr_str <- function(strx2) {
    if(!grepl("gr", strx2, fixed = T)) {
      strx_ <- strx2
      if(grepl("|", strx2, fixed = T) & !grepl("^\\(", strx2, fixed = F)) {
        strx_ <- paste0("(", strx2, "")
      }
    } else if(grepl("|gr", strx2)) {
      if(!grepl("|gr(", strx2, fixed = T)) {
        strx_ <- gsub("|gr" , "|gr(", strx2, fixed = T)
        strx_ <- paste0("(", strx_, ")")
      } else if(grepl("^|gr\\(", strx2)) {
        strx_ <- paste0("(", strx2, "")
      } 
    } else if(grepl("|gr", strx2)) {
      strx_ <-  paste0("(", strx2, "")
    }
    strx_
  }
  abx_c <- c()
  strx <- gsub("+(", "_xxxx_", strx, fixed = T)
  abxs <- strsplit(strx, "_xxxx_", fixed = T)[[1]]
  for (abxi in 1:length(abxs)) {
    abx_c[abxi] <- restore_paranthese_grgr_str(abxs[abxi])
    abx_c[abxi]
  }
  abx_c <- paste0(abx_c, collapse = "+")
  abx_c
}



get_x_random2 <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- strsplit(x, ")+" )[[1]]
  x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
  if(any(grepl("^|gr", x))) {
    x <- sub(".*gr", "", x) 
    x <- strsplit(x, ",")[[1]][1]
  } 
  x <- sub(".*\\|", "", x) 
  x <- unique(unlist(strsplit(x, ":")) )
  x
}

get_x_random2_asitis <- function(x) {
  x <- gsub("[[:space:]]", "", x)
  x <- strsplit(x, ")+" )[[1]]
  x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
  if(any(grepl("^|gr", x))) {
    x <- sub(".*gr", "", x) 
    x <- strsplit(x, ",")[[1]][1]
  } 
  x <- sub(".*\\|", "", x) 
  x
}




get_o_paranthesis <- function(x) {
  if(!grepl("lf\\(", x)) {
    x <- gsub("^lf\\(", "", x)
    x <- gsub(")$", "", x)
  }
  if(!grepl("nlf\\(", x)) {
    x <- gsub("^nlf\\(", "", x)
    x <- gsub(")$", "", x)
  }
  x <- strsplit(x, "~")[[1]][2]
  x
}


get_o_paranthesis2 <- function(x) {
  x <- gsub("^\\(", "", x)
  x <- gsub(")$", "", x)
  x <- strsplit(x, "~")[[1]][2]
  x
}



getcovlist <- function(x) {
  if (is.character(x))
    x <- x
  else
    x <- deparse(x)
  x <- gsub("~", "", gsub("\\s", "", x))
  x <- strsplit(x, "+", fixed = T)[[1]]
  if (length(x) == 1)
    x <- NULL
  else
    x <- x[-1]
  return(x)
}


ept <- function(x)
  eval(parse(text = x), envir = parent.frame())






########################
########################


get_args_ <- function(arguments, xcall) {
  `%!in%` <- Negate(`%in%`)
  f_bsitar_arg <- formals(paste0(xcall, '.bsitar'))
  nf_bsitar_arg_names <-
    intersect(names(arguments), names(f_bsitar_arg))
  arguments <-
    c(arguments, f_bsitar_arg[names(f_bsitar_arg) %!in% nf_bsitar_arg_names])
  arguments
}


########################
########################

get.cores <- function(cores.arg) {
  cores_ <- eval(cores.arg)
  if (!is.null(cores_)) {
    if (cores_ == "maximise") {
      max.cores <-
        as.numeric(future::availableCores(methods = "system", omit = 0))
      if (max.cores < 1)
        max.cores <- 1
    } else if (cores_ == "optimize") {
      max.cores <-
        as.numeric(future::availableCores(methods = "system", omit = 1))
      if (max.cores < 1)
        max.cores <- 1
    } else {
      max.cores <- eval(cores_)
    }
  } else if (is.null(cores_)) {
    max.cores <- NULL
  }
  
  if (!is.null(cores_)) {
    if (Sys.info()["sysname"] == "Windows") {
      .cores_ps <- 1
    } else {
      .cores_ps <- max.cores
    }
  } else if (is.null(cores_)) {
    .cores_ps <- 1
  }
  
  list(max.cores = max.cores, .cores_ps = .cores_ps)
}




########################
########################

validate_response <- function(model, resp = NULL) {
  if (model$model_info$nys == 1 & !is.null(resp)) {
    stop(
      "You have fit a univariate model",
      " but set resp option as ",
      resp,
      ".",
      "\n ",
      " The resp option should be appropriately set to NULL",
      "\n ",
      " (i.e., resp = NULL)"
    )
  }
  if (model$model_info$nys > 1 & is.null(resp)) {
    if (!is.na(model$model_info$univariate_by)) {
      stop(
        "You have fit a univariate-by-subset model for ",
        model$model_info$univariate_by,
        "\n ",
        " but dit not set the the resp options appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
    if (model$model_info$multivariate) {
      stop(
        "You have fit a multivariate model ",
        "\n ",
        " but dit not set the the resp options appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
  }
  
  if (!is.null(resp)) {
    if (!resp %in% model$model_info[['ys']]) {
      stop(
        "Response should be one of the following: ",
        paste(model$model_info[['ys']], collapse = " "),
        "\n ",
        " but you have specified: ",
        resp
      )
    }
  }
  
}





########################
########################

setup_higher_priors <- function(new_prior_list) {
  o_l <- list()
  ixi = 0
  group_ <- class_ <- nlpar_ <- resp_ <- cor_check <- sd_check <- NA
  group_ <-
    class_ <- nlpar_ <- resp_ <- cor_check <- sd_check <- c()
  for (new_prior_listi in 1:length(new_prior_list)) {
    ixi <- ixi + 1
    pstr <- new_prior_list[[new_prior_listi]]
    if (is.null(pstr[['prior']]))
      prior_i <- ''
    else
      prior_i <- pstr[['prior']]
    if (is.null(pstr[['group']]))
      group_i <- ''
    else
      group_i <- pstr[['group']]
    if (is.null(pstr[['class']]))
      class_i <- ''
    else
      class_i <- pstr[['class']]
    if (is.null(pstr[['nlpar']]))
      nlpar_i <- ''
    else
      nlpar_i <- pstr[['nlpar']]
    if (is.null(pstr[['coef']]))
      coef_i  <- ''
    else
      coef_i  <- pstr[['coef']]
    if (is.null(pstr[['resp']]))
      resp_i  <- ''
    else
      resp_i  <- pstr[['resp']]
    if (is.null(pstr[['lb']]))
      lb_i    <- ''
    else
      lb_i    <- pstr[['lb']]
    if (is.null(pstr[['ub']]))
      ub_i    <- ''
    else
      ub_i <- pstr[['ub']]
    if (is.null(pstr[['dpar']]))
      dpar_i  <- ''
    else
      dpar_i  <- pstr[['dpar']]

    prior_i <- gsub("[[:space:]]", "", prior_i)
    group_i <- gsub("[[:space:]]", "", group_i)
    class_i <- gsub("[[:space:]]", "", class_i)
    nlpar_i <- gsub("[[:space:]]", "", nlpar_i)
    coef_i  <- gsub("[[:space:]]", "", coef_i)
    resp_i  <- gsub("[[:space:]]", "", resp_i)
    lb_i    <- gsub("[[:space:]]", "", lb_i)
    ub_i    <- gsub("[[:space:]]", "", ub_i)
    dpar_i  <- gsub("[[:space:]]", "", dpar_i)
    
    group_ <- c(group_, group_i)
    class_ <- c(class_, class_i)
    nlpar_ <- c(nlpar_, nlpar_i)
    resp_ <- c(nlpar_, resp_i)
    if (class_i == 'sd')
      sd_check <- c(sd_check, group_i)
    if (class_i == 'cor')
      cor_check <- c(cor_check, group_i)
    
    if (lb_i == '' & ub_i == '') {
      o_l[[ixi]] <- prior_string(
        prior_i,
        group = group_i,
        class = class_i,
        nlpar = nlpar_i,
        coef = coef_i,
        resp = resp_i,
        # lb = lb_i,
        # ub = ub_i,
        dpar = dpar_i
      )
    } else if (lb_i != '' | ub_i != '') {
      o_l[[ixi]] <- prior_string(
        prior_i,
        group = group_i,
        class = class_i,
        nlpar = nlpar_i,
        coef = coef_i,
        resp = resp_i,
        # lb = lb_i,
        # ub = ub_i,
        dpar = dpar_i
      )
    } # if(lb_i == '' & ub_i == '' ) {
  } # for (new_prior_listi in 1:length(new_prior_list)) {
  o_l %>%  do.call(rbind, .)
}





# From brms

# rename specified patterns in a character vector
# @param x a character vector to be renamed
# @param pattern the regular expressions in x to be replaced
# @param replacement the replacements
# @param fixed same as for 'gsub'
# @param check_dup: logical; check for duplications in x after renaming
# @param ... passed to 'gsub'
# @return renamed character vector of the same length as x
rename <- function(x, pattern = NULL, replacement = NULL,
                   fixed = TRUE, check_dup = FALSE, ...) {
  pattern <- as.character(pattern)
  replacement <- as.character(replacement)
  if (!length(pattern) && !length(replacement)) {
    # default renaming to avoid special characters in coeffcient names
    pattern <- c(
      " ", "(", ")", "[", "]", ",", "\"", "'",
      "?", "+", "-", "*", "/", "^", "="
    )
    replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ")
  }
  if (length(replacement) == 1L) {
    replacement <- rep(replacement, length(pattern))
  }
  stopifnot(length(pattern) == length(replacement))
  # avoid zero-length pattern error
  has_chars <- nzchar(pattern)
  pattern <- pattern[has_chars]
  replacement <- replacement[has_chars]
  out <- x
  for (i in seq_along(pattern)) {
    out <- gsub(pattern[i], replacement[i], out, fixed = fixed, ...)
  }
  dup <- duplicated(out)
  if (check_dup && any(dup)) {
    dup <- x[out %in% out[dup]]
    stop2("Internal renaming led to duplicated names. \n",
          "Occured for: ", collapse_comma(dup))
  }
  out
}



# split vector at factor indices

splitAt2 <- function(x, pos) {
  x <- droplevels(x)
  out <- c()
  pos2 <- c(1, pos, length(x)+1)
  for (i in seq_along(pos2[-1])) {
    out[[i]] <- x[pos2[i]:(pos2[i+1]-1)]
  }
  return(out)
}



'%||%' <- function(x, y) {
  if (is.null(x)) x <- y
  x
}


####

# Find 
# (.)$
#   replace
# \1 <- NULL;
# 
# Both regrex and wrap ticked ues

