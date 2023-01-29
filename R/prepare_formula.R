

#' Prepare formula for fitting Bayesian SITAR growth curve model
#' 
#' The \code{prepare_formula}) prepares \code{brms::brmsformual} which is  
#' passed on to the [bsitar::bsitar] function. For univariate-by-
#' subgroup model (specified by using the \code{univariate_by}) and 
#' multivariate model (specified by using the \code{multivariate}),
#' the \code{x}, \code{y}, \code{id}, \code{knots}, \code{nknots}, are 
#' automatically set to match the sub-model(s). See \code{brms::brmsformual} 
#' for details. 
#'
#' @param x vector of predictor (typically age in years).
#' @param y vector of outcome (i.e., repeated growth measurements). 
#' @param id a factor variable identifying the groups (typically individuals).
#' @param knots vector of values for knots.
#' @param nknots an integer specifying the number of knots.
#' @param data data frame containing variables \code{x}, \code{y} and \code{id}.
#' @param internal_formula_args Other internal arguments passed from the 
#' [bsitar::bsitar] to the \code{prepare_formula}).
#'
#' @return An object of class \code{brmsformula}, which is a \code{list} 
#'   containing formulas
#'   
#' @seealso [brms::brm] [brms::brmsformula]
#'   
#' @importFrom nlme lme
#'   
#' @export
#' 
prepare_formula <- function(x,
                            y,
                            id,
                            knots,
                            nknots,
                            data,
                            internal_formula_args) {
  if (!is.null(internal_formula_args)) {
    eout <- list2env(internal_formula_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  
  
  
  if (!is.null(group_arg$cor)) {
    if (group_arg$cor == "un")
      abccorr <- TRUE
    if (group_arg$cor == "diagonal")
      abccorr <- FALSE
  } else {
    group_arg$cor <- "un"
    abccorr <- TRUE
  }
  
  
  if (!is.null(group_arg$by)) {
    group_arg$by <- group_arg$by
  } else {
    group_arg$by <- NULL
  }
  
  if (!is.null(group_arg$cov)) {
    group_arg$cov <- group_arg$cov
  } else {
    group_arg$cov <- NULL
  }
  
  if (!is.null(group_arg$dist)) {
    group_arg$dist <- group_arg$dist
  } else {
    group_arg$dist <- 'gaussian'
  }
  
  if (!is.null(group_arg$verbose)) {
    group_arg$verbose <- group_arg$verbose
  } else {
    group_arg$verbose <- FALSE
  }
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (!is.null(univariate_by$cor)) {
      if (univariate_by$cor == "un")
        uvarabccorr <- TRUE
      if (univariate_by$cor == "diagonal")
        uvarabccorr <- FALSE
    } else {
      univariate_by$cor <- "un"
      uvarabccorr <- TRUE
    }
    if (is.null(univariate_by$verbose))
      univariate_by$verbose <- FALSE
  }
  
  if ((is.na(univariate_by$by) | univariate_by$by == "NA")) {
    univariate_by$cor <- "un"
    uvarabccorr <- FALSE
    univariate_by$verbose <- FALSE
  }
  
  if (multivariate$mvar) {
    if (!is.null(multivariate$cor)) {
      if (multivariate$cor == "un")
        mvarccorr <- "WB"
      if (multivariate$cor == "un_s")
        mvarccorr <- "W"
      if (multivariate$cor == "diagonal")
        mvarccorr <- "none"
    } else {
      multivariate$cor <- "un"
      mvarccorr <- "WB"
    }
    if (!is.null(multivariate$rescor)) {
      multivariate$rescor <- multivariate$rescor
    } else {
      multivariate$rescor <- TRUE
    }
    if (is.null(multivariate$verbose))
      multivariate$verbose <- FALSE
  }
  
  if (!multivariate$mvar) {
    multivariate$cor <- "none"
    multivariate$rescor <- FALSE
    multivariate$verbose <- FALSE
  }
  
  abcnames <-
    paste0(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]], sep = ",")
  
  snames <- c()
  for (i in 1:(nknots - 1)) {
    if (i < (nknots - 1)) {
      name1 <- paste0("s", i, sep = ",")
    }
    else {
      name1 <- paste0("s", i, sep = "")
    }
    snames[i] <- name1
  }
  
  ###########
  # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
  # In fact for df > 1, it forces 'd' to be random parameter only
  if (match_sitar_d_form) {
    if (!grepl("d", fixedsi, fixed = F) &
        grepl("d", randomsi, fixed = F)) {
      abcnames <- c(abcnames, "d,")
    }
  }
  ###########
  
  fullabcsnames <- c(abcnames, snames)
  
  abcselements  <- paste0(fullabcsnames, collapse = "")
  abcselements  <- paste0(spfncname, "(", x, ",", abcselements, ")")
  abcsformfit   <- (paste0(y, " ~ ", abcselements)) # as.formula
  
  if (!(is.na(univariate_by$by) |
        univariate_by$by == "NA") &
      !is.null(subindicatorsi)) {
    abcsformfit <- (paste0(y,  "| subset(", subindicatorsi, ")",
                           " ~ ", abcselements))
  }
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") &
      !multivariate$mvar &  abccorr) {
    coridv <- "C"
  }
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") &
      !multivariate$mvar & !abccorr) {
    coridv <- ""
  }
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (uvarabccorr) {
      coridv <- y
    } else {
      coridv <- ""
    }
  }
  
  if (multivariate$mvar && mvarccorr == "none")
    coridv <- ""
  if (multivariate$mvar && mvarccorr == "W")
    coridv <- y
  if (multivariate$mvar && mvarccorr == "WB")
    coridv <- "MVC"
  
  dmatnames <- NULL
  
  if (grepl("a", fixedsi, fixed = T)) {
    afixed <- a_formulasi
  } else {
    afixed <- NULL
  }
  if (grepl("b", fixedsi, fixed = T)) {
    bfixed <- b_formulasi
  } else {
    bfixed <- NULL
  }
  if (grepl("c", fixedsi, fixed = T)) {
    cfixed <- c_formulasi
  } else {
    cfixed <- NULL
  }
  if (grepl("d", fixedsi, fixed = T)) {
    dfixed <- d_formulasi
  } else {
    dfixed <- NULL
  }
  
  sfixed <- s_formulasi
  

  if (grepl("a", randomsi, fixed = T)) {
    arandom <- a_formula_grsi
  } else {
    arandom <- NULL
  }
  if (grepl("b", randomsi, fixed = T)) {
    brandom <- b_formula_grsi
  } else {
    brandom <- NULL
  }
  if (grepl("c", randomsi, fixed = T)) {
    crandom <- c_formula_grsi
  } else {
    crandom <- NULL
  }
  if (grepl("d", randomsi, fixed = T)) {
    drandom <- d_formula_grsi
  } else {
    drandom <- NULL
  }
  
  arandom_wb <- NULL
  arandom_wb_ <- FALSE
  brandom_wb <- crandom_wb <- drandom_wb <- arandom_wb
  brandom_wb_ <- crandom_wb_ <- drandom_wb_ <- arandom_wb_
  
  get_x_random <- function(x) {
    x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
    x <- gsub("~1+1", "~1", x, fixed = T)
    if(!grepl("~", x)) x <- paste0("~", x)
    x <- sub("\\|.*", "", x)
    x
  }
  
  if(!is.null(arandom) & grepl("|", arandom, fixed = TRUE)) {
    arandom_wb <- gsub("~", "", arandom) # with bar
    arandom_wb_ <- TRUE
    arandom <- get_x_random(arandom)
  }
  if(!is.null(brandom) & grepl("|", brandom, fixed = TRUE)) {
    brandom_wb <- gsub("~", "", brandom)
    brandom_wb_ <- TRUE
    brandom <- get_x_random(brandom)
  }
  if(!is.null(crandom) & grepl("|", crandom, fixed = TRUE)) {
    crandom_wb <- gsub("~", "", crandom)
    crandom_wb_ <- TRUE
    crandom <- get_x_random(crandom)
  }
  
  if(!is.null(drandom)) {
    if(grepl("|", drandom, fixed = TRUE)) {
      drandom_wb <- gsub("~", "", drandom)
      drandom_wb_ <- TRUE
      drandom <- get_x_random(drandom)
    }
  }
  
  
  arandom_wb <- gsub("1+1", "1", arandom_wb, fixed = T)
  brandom_wb <- gsub("1+1", "1", brandom_wb, fixed = T)
  crandom_wb <- gsub("1+1", "1", crandom_wb, fixed = T)
  drandom_wb <- gsub("1+1", "1", drandom_wb, fixed = T)
  
  
  if (!is.null(afixed)) {
    acovmat <- eval(parse(text = paste0(
      "model.matrix(",
      afixed, ",data = data)"
    )))
    if (ncol(acovmat) == 1)
      nacov <- NULL
    else
      nacov <- ncol(acovmat) - 1
    acovcoefnames <- colnames(acovmat)
    acovcoefnames <- gsub("\\(|)", "", acovcoefnames)
  } else if (is.null(afixed)) {
    nacov <- acovcoefnames <- NULL
  }
  
  if (!is.null(bfixed)) {
    bcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      bfixed, ",data = data)"
    )))
    if (ncol(bcovmat) == 1)
      nbcov <- NULL
    else
      nbcov <- ncol(bcovmat) - 1
    bcovcoefnames <- colnames(bcovmat)
    bcovcoefnames <- gsub("\\(|)", "", bcovcoefnames)
  } else if (is.null(bfixed)) {
    nbcov <- bcovcoefnames <- NULL
  }
  
  if (!is.null(cfixed)) {
    ccovmat <- eval(parse(text = paste0(
      "model.matrix(",
      cfixed, ",data = data)"
    )))
    if (ncol(ccovmat) == 1)
      nccov <- NULL
    else
      nccov <- ncol(ccovmat) - 1
    ccovcoefnames <- colnames(ccovmat)
    ccovcoefnames <- gsub("\\(|)", "", ccovcoefnames)
  } else if (is.null(cfixed)) {
    nccov <- ccovcoefnames <- NULL
  }
  
  if (!is.null(dfixed)) {
    dcovmat <- eval(parse(text = paste0(
      "model.matrix(",
      dfixed, ",data = data)"
    )))
    if (ncol(dcovmat) == 1)
      ndcov <- NULL
    else
      ndcov <- ncol(dcovmat) - 1
    dcovcoefnames <- colnames(dcovmat)
    dcovcoefnames <- gsub("\\(|)", "", dcovcoefnames)
  } else if (is.null(dfixed)) {
    ndcov <- dcovcoefnames <- NULL
  }
  
  
  if (!is.null(sfixed)) {
    scovmat <- eval(parse(text = paste0(
      "model.matrix(",
      sfixed, ",data = data)"
    )))
    if (ncol(scovmat) == 1) {
      nscov <- NULL
    } else {
      nscov <- ncol(scovmat) - 1
    }
    scovcoefnames <- colnames(scovmat)
    scovcoefnames <- gsub("\\(|)", "", scovcoefnames)
  }
  
  if (!is.null(arandom)) {
    acovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      arandom, ",data = data)"
    )))
    if (ncol(acovmat_gr) == 1) {
      nacov_gr <- NULL
    } else {
      nacov_gr <- ncol(acovmat_gr) - 1
    }
    acovcoefnames_gr <- colnames(acovmat_gr)
    acovcoefnames_gr <- gsub("\\(|)", "", acovcoefnames_gr)
  } else if (is.null(arandom)) {
    nacov_gr <- acovcoefnames_gr <- NULL
  }
  
  if (!is.null(brandom)) {
    bcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      brandom, ",data = data)"
    )))
    if (ncol(bcovmat_gr) == 1) {
      nbcov_gr <- NULL
    } else {
      nbcov_gr <- ncol(bcovmat_gr) - 1
    }
    bcovcoefnames_gr <- colnames(bcovmat_gr)
    bcovcoefnames_gr <- gsub("\\(|)", "", bcovcoefnames_gr)
  } else if (is.null(brandom)) {
    nbcov_gr <- bcovcoefnames_gr <- NULL
  }
  
  if (!is.null(crandom)) {
    ccovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      crandom, ",data = data)"
    )))
    if (ncol(ccovmat_gr) == 1) {
      nccov_gr <- NULL
    } else {
      nccov_gr <- ncol(ccovmat_gr) - 1
    }
    ccovcoefnames_gr <- colnames(ccovmat_gr)
    ccovcoefnames_gr <- gsub("\\(|)", "", ccovcoefnames_gr)
  } else if (is.null(crandom)) {
    nccov_gr <- ccovcoefnames_gr <- NULL
  }
  
  if (!is.null(drandom)) {
    dcovmat_gr <- eval(parse(text = paste0(
      "model.matrix(",
      drandom, ",data = data)"
    )))
    if (ncol(dcovmat_gr) == 1) {
      ndcov_gr <- NULL
    } else {
      ndcov_gr <- ncol(dcovmat_gr) - 1
    }
    dcovcoefnames_gr <- colnames(dcovmat_gr)
    dcovcoefnames_gr <- gsub("\\(|)", "", dcovcoefnames_gr)
  } else if (is.null(drandom)) {
    ndcov_gr <- dcovcoefnames_gr <- NULL
  }
  
  
  if (!is.null(afixed)) {
    aform <- paste0("a", afixed)
  } else {
    aform <- NULL
  }
  if (!is.null(bfixed)) {
    bform <- paste0("b", bfixed)
  } else {
    bform <- NULL
  }
  if (!is.null(cfixed)) {
    cform <- paste0("c", cfixed)
  } else {
    cform <- NULL
  }
  if (!is.null(dfixed)) {
    dform <- paste0("d", dfixed)
  } else {
    dform <- NULL
  }
  
  
  sinterceptelements <- paste0(snames, collapse = "+")
  sinterceptelements <- gsub("," , "" , sinterceptelements)
  sform     <- paste0(sinterceptelements, sfixed)
  
  if (!is.null(arandom)) {
    aform_gr <- gsub("^~", "", arandom, fixed = F)
  } else {
    aform_gr <- NULL
  }
  if (!is.null(brandom)) {
    bform_gr <- gsub("^~", "", brandom, fixed = F)
  } else {
    bform_gr <- NULL
  }
  if (!is.null(crandom)) {
    cform_gr <- gsub("^~", "", crandom, fixed = F)
  } else {
    cform_gr <- NULL
  }
  if (!is.null(drandom)) {
    dform_gr <- gsub("^~", "", drandom, fixed = F)
  } else {
    dform_gr <- NULL
  }
  
  
  if (!is.null(group_arg)) {
    gr_prefixss <- "gr"
    gr_varss <- group_arg$groupvar
    gr_by = group_arg$by
    gr_cov = group_arg$cov
    gr_dist = group_arg$dist
  } else {
    gr_prefixss <- NULL
    gr_varss    <- NULL
  }
  
  if (!is.null(randomsi)) {
    # these two arguments are set automaticaly if missing
    if (is.null(gr_prefixss))
      gr_prefixss <- 'gr'
    if (is.null(gr_varss))
      gr_varss <- id
    if (!is.null(group_arg)) {
      if (!is.null(gr_by))  {
        gr_byss2  <- paste0(",by=", gr_by)
      } else {
        gr_byss2 <- NULL
      }
      if (!is.null(gr_cov)) {
        gr_covss2  <- paste0(",cov=", gr_cov)
      } else {
        gr_covss2 <- NULL
      }
      if (!is.null(gr_dist)) {
        gr_distss2 <- paste0(",dist=", paste0("'", gr_dist, "'"))
      } else {
        gr_distss2 <- ''
      }
      gr__args <- paste0(gr_prefixss,
                         "(",
                         gr_varss,
                         gr_byss2,
                         gr_covss2,
                         gr_distss2,
                         ")")
      
      
      
      if (!is.null(aform) & grepl("a", randomsi, fixed = T)) {
        if(!arandom_wb_) {
          aform <- paste0(aform,
                          " + (", aform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(arandom_wb_) {
          aform <- paste0(aform,
                          " + (", arandom_wb, ")")
        }
      }
      
      if (!is.null(bform) & grepl("b", randomsi, fixed = T)) {
        if(!brandom_wb_) {
          bform <- paste0(bform,
                          " + (", bform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(brandom_wb_) {
          bform <- paste0(bform,
                          " + (", brandom_wb, ")")
        }
      }
      
      if (!is.null(cform) & grepl("c", randomsi, fixed = T)) {
        if(!crandom_wb_) {
          cform <- paste0(cform,
                          " + (", cform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(crandom_wb_) {
          cform <- paste0(cform,
                          " + (", crandom_wb, ")")
        }
      }
      
      if (!is.null(dform) & grepl("d", randomsi, fixed = T)) {
        if(!drandom_wb_) {
          dform <- paste0(dform,
                          " + (", dform_gr, "|", coridv, "|" , gr__args, ")")
        } else if(drandom_wb_) {
          dform <- paste0(dform,
                          " + (", drandom_wb, ")")
        }
      }
      
      ###########
      # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
      # In fact for df > 1, it forces 'd' to be random parameter only
      # allow random effect for 'd' even corresponding fixed effect is missing
      # only allowing for 'd'
      if (match_sitar_d_form) {
        if (is.null(dform) & grepl("d", randomsi, fixed = T)) {
          if(!drandom_wb_) {
            dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , gr__args, ")")
          } else if(drandom_wb_) {
            dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , gr__args, ")")
          }
        }
      }
      ###########
    }
    
    
    if (is.null(group_arg)) {
      if (!is.null(aform) & grepl("a", randomsi, fixed = T)) {
        if(!arandom_wb_) { 
          aform <- paste0(aform, " + (", aform_gr, "|", coridv, "|" , id, ")")
        } else if(arandom_wb_) {
          aform <- paste0(aform, " + (", arandom_wb, ")")
        }
      }
      
      if (!is.null(bform) & grepl("b", randomsi, fixed = T)) {
        if(!brandom_wb_) { 
          bform <- paste0(bform, " + (", bform_gr, "|", coridv, "|" , id, ")")
        } else if(brandom_wb_) {
          bform <- paste0(bform, " + (", brandom_wb, ")")
        }
      }
      
      if (!is.null(cform) & grepl("c", randomsi, fixed = T)) {
        if(!crandom_wb_) { 
          cform <- paste0(cform, " + (", cform_gr, "|", coridv, "|" , id, ")")
        } else if(crandom_wb_) {
          cform <- paste0(cform, " + (", crandom_wb, ")")
        }
      }
      
      if (!is.null(cform) & grepl("d", randomsi, fixed = T)) {
        if(!drandom_wb_) { 
          dform <- paste0(dform, " + (", dform_gr, "|", coridv, "|" , id, ")")
        } else if(drandom_wb_) {
          dform <- paste0(dform, " + (", drandom_wb, ")")
        }
      }
      
      ###########
      # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
      # In fact for df > 1, it forces 'd' to be random parameter only
      # allow random effect for 'd' even corresponding fixed effect is missing
      # only allowing for 'd'
      if (match_sitar_d_form) {
        if (is.null(dform) & grepl("d", randomsi, fixed = T)) {
          if(!drandom_wb_) { 
            dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , id, ")")
          } else if(drandom_wb_) {
            dform <- paste0(dform, "d ~ 0 + (1|", coridv , "|" , id, ")")
          }
        }
      }
      ###########
    }
  }
  
  

   
   add_higher_level_str <- function(form, str) {
     if(!is.null(str[[1]])) {
       get_n_str <- strsplit(str, ")+(", fixed = T)[[1]][-1]
       get_n_str_length <- length(get_n_str)
     } else {
       get_n_str_length <- 0
     }
     if(get_n_str_length != 0) {
       str <- gsub(")+(", ")_(", str, fixed = T)
       str_ <- sub("^[^_]*_", "", str)
       str_ <- gsub(")_(", ")+(", str_, fixed = T)
       form <- paste0(form, "+", str_)
     } else {
       form <- form
     }
     form <- gsub("[[:space:]]", "", form)
     form
   }
   
   
  
  
  if(set_higher_levels) {
    aform <- add_higher_level_str(aform, a_formula_gr_strsi)
    bform <- add_higher_level_str(bform, b_formula_gr_strsi)
    cform <- add_higher_level_str(cform, c_formula_gr_strsi)
    if(!is.null(dform)) {
      dform <- add_higher_level_str(dform, d_formula_gr_strsi)
    }
    
    extract_xx <- NULL
    if(!is.null(a_formula_gr_strsi[[1]])) {
      extract_xx <- a_formula_gr_strsi
    } else if(!is.null(b_formula_gr_strsi[[1]])) {
      extract_xx <- b_formula_gr_strsi
    } else if(!is.null(c_formula_gr_strsi[[1]])) {
      extract_xx <- c_formula_gr_strsi
    } else if(!is.null(d_formula_gr_strsi[[1]])) {
      extract_xx <- d_formula_gr_strsi
    }
    
    extract_xx <- gsub("[[:space:]]", "", extract_xx)
    extract_xx <- strsplit(extract_xx, ")+(", fixed = T)[[1]][1]
    extract_xx <- gsub("\\)", "", extract_xx)
    gr_varss <- sub(".*\\|", "", extract_xx) 
    
    
    get_x_random2 <- function(x) {
      x <- gsub("[[:space:]]", "", x)
      x <- strsplit(x, ")+" )[[1]]
      x <- gsub("[[:space:]]", "", gsub("[()]", "", x))
      x <- sub(".*\\|", "", x) 
      x <- unique(unlist(strsplit(x, ":")) )
      x
    }
    aform_gr_names <- lapply(aform, get_x_random2)[[1]]
    bform_gr_names <- lapply(bform, get_x_random2)[[1]]
    cform_gr_names <- lapply(cform, get_x_random2)[[1]]
    if(!is.null(dform)) {
      dform_gr_names <- lapply(dform, get_x_random2)[[1]]
    } else {
      dform_gr_names <- NULL
    }
    hierarchical_gr_names <- c(aform_gr_names, bform_gr_names, 
                               cform_gr_names, dform_gr_names)
    hierarchical_gr_names <- unique(hierarchical_gr_names)
    
  } # if(set_higher_levels) {
  
   
   if(!set_higher_levels) hierarchical_gr_names <- NULL
   
   # print(gr_varss)
   # stop()
  
  # print(a_formula_gr_strsi_)
  
  
  
  
  
  
  if (!is.null(dpar_formulasi)) {
    if (!grepl("lf\\(", dpar_formulasi) &
        !grepl("nlf\\(", dpar_formulasi)) {
      dpar_covi_mat_form <- dpar_formulasi
    } else {
      dpar_covi_mat_form <- gsub("\\(|)", "",
                                 strsplit(dpar_formulasi, "~")[[1]][2])
      dpar_covi_mat_form <- paste0("~", dpar_covi_mat_form)
    }
  }
  
  
  
  if (!is.null(dpar_formulasi)) {
    dparcovmat <- eval(parse(
      text =
        paste0("model.matrix(",
               dpar_covi_mat_form, ",data = data)")
    ))
    if (ncol(dparcovmat) == 1) {
      ndparcov <- NULL
    } else {
      ndparcov <- ncol(dparcovmat) - 1
    }
    dparcovcoefnames <- colnames(dparcovmat)
    dparcovcoefnames <- gsub("\\(|)", "", dparcovcoefnames)
  }
  
  if (is.null(dpar_formulasi)) {
    ndparcov <- NULL
    dparcovcoefnames <- NULL
  }
  
  abcform <-
    paste(cbind(aform, bform, cform, dform), collapse = ",")
  
  
  
  
  
  
  
  # Imp that beyond this point the bform will be changed to combine
  
  if (!is.null(dpar_formulasi)) {
    if (!grepl("lf\\(", dpar_formulasi) &
        !grepl("nlf\\(", dpar_formulasi)) {
      sform <- paste0(sform, ",", paste0("sigma", dpar_formulasi))
    } else {
      sform <- sform
    }
  }
  
  if (!is.null(autocor_formi)) {
    sform <- paste0(sform,
                    ",",
                    paste0("autocor=", autocor_formi))
  }
  
  bform <- paste0("bf(",
                  abcsformfit,
                  ", " ,
                  abcform,
                  ", " ,
                  sform,
                  ", " ,
                  "nl = T ,loop = F)")
  
  bform <- gsub("\\s", "", bform)
  
  
  if (!is.null(dpar_formulasi)) {
    if (grepl("lf\\(", dpar_formulasi) |
        grepl("nlf\\(", dpar_formulasi)) {
      if (!grepl("\\(~", dpar_formulasi) &
          strsplit(strsplit(dpar_formulasi, "~")[[1]][1],
                   "\\(")[[1]][2] != "sigma") {
        stop(
          "The distributional parameter name on the left hand side of ",
          "\n ",
          " lf/nlf formula (i.e., before ~) should be 'sigma'"
        )
      } else if (grepl("sigma~", dpar_formulasi, fixed = T)) {
        dpar_formulasi <- dpar_formulasi
      } else {
        dpar_formulasi <- gsub("~", "sigma~", dpar_formulasi)
      }
      bform <- paste0(bform, "+", dpar_formulasi)
    } else {
      bform <- bform
    }
  }
  
  
  if (!is.null(familysi)) {
    bform <- paste0(bform, "+", familysi)
  }
  
  
  group_arg_groupvar <- gr_varss
  multivariate_rescor <-  multivariate$rescor
  univariate_by_by <- univariate_by$by
  
  # fit lm model
  
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
  
  a_covariate <- getcovlist(a_formulasi)
  b_covariate <- getcovlist(b_formulasi)
  c_covariate <- getcovlist(c_formulasi)
  d_covariate <- getcovlist(d_formulasi)
  
  s_covariate <- getcovlist(s_formulasi)
  
  covariates <- c(a_covariate, b_covariate, c_covariate, d_covariate, s_covariate)
  covariates_ <- unique(covariates)
  
  
  
  a_covariate_i <- c()
  if (grepl("\\*", a_formulasi)) {
    for (a_covariatei in a_covariate) {
      if (grepl("\\*", a_covariatei)) {
        t <- strsplit(a_covariatei, "\\*")[[1]]
        t <- c(paste0(t, collapse = "+"), paste0(t, collapse = ":"))
      } else {
        t <- a_covariatei
      }
      a_covariate_i <- c(a_covariate_i, t)
    }
  } else {
    a_covariate_i <- a_covariate
  }
  
  a_covariate <- a_covariate_i
  
  
  
  
  s_covariate_i <- c()
  if (grepl("\\*", s_formulasi)) {
    for (s_covariatei in s_covariate) {
      if (grepl("\\*", s_covariatei)) {
        t <- strsplit(s_covariatei, "\\*")[[1]]
        t <- c(paste0(t, collapse = "+"), paste0(t, collapse = ":"))
      } else {
        t <- s_covariatei
      }
      s_covariate_i <- c(s_covariate_i, t)
    }
  } else {
    s_covariate_i <- s_covariate
  }
  
  s_covariate <- s_covariate_i
  
  
  
  if (!grepl("^~1$", s_formulasi)) {
    if (!identical(strsplit(a_formulasi, "+")[[1]][1:2],
                   strsplit(s_formulasi, "+")[[1]][1:2])) {
      stop(
        "a_formula and s_formula should have the identical ",
        "\n ",
        " intercept structure i.e., ~ 1 or ~ 0"
      )
    }
    
    if (length(a_covariate) != length(s_covariate)) {
      stop(
        "s_formula formula should be intercept only or must have ",
        "\n ",
        " the same number of covariates as a_formula"
      )
    }
    
    if (length(a_covariate) == 1) {
      if (grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~0+",  a_covariate, "*", "mat_s"))
      } else if (!grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~1+",  a_covariate, "*", "mat_s"))
      }
    } else if (length(a_covariate) > 1) {
      main_cov <- a_covariate
      main_cov <-
        paste(unlist(strsplit(main_cov, "+", fixed = T)), sep = " ")
      inte_cov <- paste0(main_cov, ":", "mat_s")
      main_inte_cov <- c(main_cov, "mat_s", inte_cov)
      main_inte_cov <- paste(main_inte_cov, collapse = "+")
      if (grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~0+",  main_inte_cov))
      } else if (!grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~1+",  main_inte_cov))
      }
    }
  }
  
  if (grepl("^~1$", s_formulasi)) {
    if (length(a_covariate) == 1) {
      if (grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~0+",  a_covariate, "+", "mat_s"))
      } else if (!grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~1+",  a_covariate, "+", "mat_s"))
      }
    } else if (length(a_covariate) > 1) {
      main_cov <- a_covariate
      main_inte_cov <- c(main_cov, "mat_s")
      main_inte_cov <- paste(main_inte_cov, collapse = "+")
      if (grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~0+",  main_inte_cov))
      } else if (!grepl("~0", s_formulasi, fixed = T)) {
        lmform  <- as.formula(paste0(y, "~1+",  main_inte_cov))
      }
    }
  }
  
  
  if (grepl("^~1$", s_formulasi))
    s_covariate <- NULL
  
  
  if (grepl("^~1$", a_formulasi)) {
    if (grepl("~0", a_formulasi, fixed = T)) {
      lmform  <- as.formula(paste0(y, "~0+", "mat_s"))
    } else if (!grepl("~0", a_formulasi, fixed = T)) {
      lmform  <- as.formula(paste0(y, "~1+", "mat_s"))
    }
  }
  
  
  
  lm_fit  <- lm(lmform, data = data)
  lm_coef <- coef(lm_fit)
  
  lm_rsd  <- summary(lm_fit)$sigma
  
  # lme
  err. <- FALSE
  tryCatch(
    expr = {
      datalme <- data
      datalme[['mat_s']] <- eval(parse(text = 'mat_s'))
      randomlmer <- "Intercept"
      if (randomlmer == "slope")
        randomform <- paste0("~ 1 + ", x , " | ", id)
      if (randomlmer == "Intercept")
        randomform <- paste0("~ 1 ", " | ", id)
      randomform <- as.formula(randomform)
      lme_fit <-
        nlme::lme(fixed = lmform,
                  random = randomform,
                  data = datalme)
    },
    error = function(e) {
      err. <<- TRUE
    }
  )
  if (err.) {
    lme_coef <- lm_coef
    lme_sd_a <- sd(predict(lm_fit))
    lme_rsd <- lm_rsd
  } else if (!err.) {
    lme_coef <- unname(nlme::fixed.effects(lme_fit))
    VarCorrnumeric <- VarCorr(lme_fit)[, 2] %>% as.numeric()
    lme_sd_a <- VarCorrnumeric[1]
    if (randomlmer == "Intercept") {
      lme_rsd <- VarCorrnumeric[2]
    } else if (randomlmer == "slope") {
      lme_rsd <- VarCorrnumeric[3]
    }
  }
  
  
  
  if (grepl("\\*", a_formulasi) & grepl("^~1$", s_formulasi)) {
    intercept_ <- lm_coef[!grepl("^mat_s", names(lm_coef))]
    spls_ <- lm_coef[grepl("^mat_s", names(lm_coef))]
    lm_coef <- c(intercept_, spls_)
  }
  
  
  if (grepl("\\*", a_formulasi) & !grepl("^~1$", s_formulasi)) {
    intercept_ <- lm_coef[!grepl("mat_s", names(lm_coef))]
    spls_ <- lm_coef[grepl("mat_s", names(lm_coef))]
    lm_coef <- c(intercept_, spls_)
  }
  
  
  
  
  
  
  lm_a_all <- lm_coef[1:ncol(acovmat)]
  if (!is.null(bfixed)) {
    lm_b_all <- rep(0, ncol(bcovmat))
  } else {
    lm_b_all <- NULL
  }
  if (!is.null(cfixed)) {
    lm_c_all <- rep(0, ncol(ccovmat))
  } else {
    lm_c_all <- NULL
  }
  if (!is.null(dfixed)) {
    lm_d_all <- rep(0, ncol(dcovmat))
  } else {
    lm_d_all <- NULL
  }
  
  lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
  
  
  
  if (grepl("~1", a_formulasi, fixed = T)) {
    lm_a_all[1] <- lm_a_all[1] + lm_s_all[1] * min(knots)
  }
  
  
  names(lm_a_all) <- acovcoefnames
  names(lm_b_all) <- bcovcoefnames
  names(lm_c_all) <- ccovcoefnames
  names(lm_d_all) <- dcovcoefnames
  
  lm_a <- lm_a_all[1]
  lm_b <- lm_b_all[1]
  lm_c <- lm_c_all[1]
  lm_d <- lm_d_all[1]
  
  if (!is.null(nacov)) {
    lm_a_cov <- lm_a_all[2:length(lm_a_all)]
  } else {
    lm_a_cov <- NULL
  }
  if (!is.null(nbcov)) {
    lm_b_cov <- lm_b_all[2:length(lm_b_all)]
  } else {
    lm_b_cov <- NULL
  }
  if (!is.null(nccov)) {
    lm_c_cov <- lm_c_all[2:length(lm_c_all)]
  } else {
    lm_c_cov <- NULL
  }
  if (!is.null(ndcov)) {
    lm_d_cov <- lm_d_all[2:length(lm_d_all)]
  } else {
    lm_d_cov <- NULL
  }
  
  
  
  
  xnames <- names(lm_s_all)
  names_mat_s_scovmat <- c()
  for (x in xnames) {
    x.a <- rev(strsplit(x, "_")[[1]])
    x.a <- gsub("^mat", "Intercept", x.a)
    x.a <- gsub(":mat", "", x.a)
    x.a <- paste(x.a, collapse = "_")
    names_mat_s_scovmat <- c(names_mat_s_scovmat, x.a)
  }
  
  if (grepl("^~1$", s_formulasi)) {
    mat_s_scovmat <- mat_s
  } else if (!grepl("^~1$", s_formulasi)) {
    mat_s_scovmat <- model.matrix(lm_fit)
    mat_s_scovmat <-
      mat_s_scovmat[, (ncol(acovmat) + 1):ncol(mat_s_scovmat)]
  }
  sds_X <- rep(0, ncol(mat_s_scovmat))
  
  for (i in 1:ncol(mat_s_scovmat)) {
    sds_X[i] = sd(mat_s_scovmat[, i])
  }
  lm_sdx_all <- sd(data[[y]]) / sds_X
  
  
  names(lm_s_all) <- names_mat_s_scovmat
  names(lm_sdx_all) <- names_mat_s_scovmat
  
  
  lm_s   <- lm_s_all[1:(nknots - 1)]
  lm_sdx <- lm_sdx_all[1:(nknots - 1)]
  
  
  if (!is.null(s_covariate) & length(s_covariate) > 1) {
    lm_s_cov <- lm_s_all[nknots:length(lm_s_all)]
    lm_sdx_cov <- lm_sdx_all[nknots:length(lm_sdx_all)]
    
    inname_c_all <- c()
    for (inname in paste0("s", 1:df)) {
      t <- names(lm_a_all)[2:length(names(lm_a_all))]
      inname_c_all <- c(inname_c_all, paste0(inname, "_", t))
    }
    lm_s_cov   <- lm_s_cov[order(factor(names(lm_s_cov),
                                        levels = inname_c_all))]
    lm_sdx_cov <- lm_sdx_cov[order(factor(names(lm_sdx_cov),
                                          levels = inname_c_all))]
    
    lm_s_all <- lm_sdx_all <- c()
    for (idfi in 1:df) {
      lm_s_all <- c(lm_s_all, c(lm_s[idfi],
                                lm_s_cov[grep(paste0("s", idfi, "_"),
                                              names(lm_s_cov))]))
      lm_sdx_all <- c(lm_sdx_all, c(lm_sdx[idfi],
                                    lm_sdx_cov[grep(paste0("s", idfi, "_"),
                                                    names(lm_sdx_cov))]))
    }
  } else if ((!is.null(s_covariate) &
              length(s_covariate) == 1) |
             is.null(s_covariate)) {
    if (!grepl("~0", s_formulasi, fixed = T)) {
      lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
      names(lm_s_all) <- names_mat_s_scovmat
      names(lm_sdx_all) <- names_mat_s_scovmat
      lm_s <- lm_s_all[1:(nknots - 1)]
      lm_sdx <- lm_sdx_all[1:(nknots - 1)]
      if (length(lm_s_all) > (nknots - 1)) {
        lm_s_cov <- lm_s_all[nknots:length(lm_s_all)]
        lm_sdx_cov <- lm_sdx_all[nknots:length(lm_sdx_all)]
        tnames_s <-
          names_mat_s_scovmat[nknots:length(names_mat_s_scovmat)]
        names(lm_s_cov) <- tnames_s
        tnames_sdx <-
          names_mat_s_scovmat[nknots:length(names_mat_s_scovmat)]
        names(lm_sdx_cov) <- tnames_sdx
      } else {
        lm_s_cov <- NULL
        lm_sdx_cov <- NULL
      }
    }
    
    
    if (grepl("~0", s_formulasi, fixed = T)) {
      lm_s_all <- lm_coef[(ncol(acovmat) + 1):length(lm_coef)]
      names(lm_s_all) <- names_mat_s_scovmat
      names(lm_sdx_all) <- names_mat_s_scovmat
      inname_c_all <- c()
      for (inname in paste0("s", 1:df)) {
        t <- names(lm_a_all)[1:length(names(lm_a_all))]
        inname_c_all <- c(inname_c_all, paste0(inname, "_", t))
      }
      lm_s <- NULL
      lm_s_cov <- NULL
      lm_sdx <- NULL
      lm_sdx_cov <- NULL
    }
  }
  
  if (any(is.na(lm_coef))) {
    stop(
      "Inclusion of covariates resulted in rank-deficient design  matrix",
      "\n ",
      "(with some NA coefficients for the 'lm' model fit)",
      "\n ",
      "Please simplyfy the model"
    )
  }
  
  # brms removes while spaces from the coefficient names
  # mimicking that behaviors but keeping it seperate here as 
  # brms may later change it to underscore or something else
  
  gsubitbt <- ""
  acovcoefnames <- gsub("[[:space:]]", gsubitbt, acovcoefnames)
  bcovcoefnames <- gsub("[[:space:]]", gsubitbt, bcovcoefnames)
  ccovcoefnames <- gsub("[[:space:]]", gsubitbt, ccovcoefnames)
  dcovcoefnames <- gsub("[[:space:]]", gsubitbt, dcovcoefnames)
  acovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, acovcoefnames_gr)
  bcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, bcovcoefnames_gr)
  ccovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, ccovcoefnames_gr)
  dcovcoefnames_gr <- gsub("[[:space:]]", gsubitbt, dcovcoefnames_gr)
  dparcovcoefnames <- gsub("[[:space:]]", gsubitbt, dparcovcoefnames)
  
  
  
  list_out <- list(
    nacov = nacov,
    nbcov = nbcov,
    nccov = nccov,
    ndcov = ndcov,
    nscov = nscov,
    acovcoefnames = acovcoefnames,
    bcovcoefnames = bcovcoefnames,
    ccovcoefnames = ccovcoefnames,
    dcovcoefnames = dcovcoefnames,
    scovcoefnames = scovcoefnames,
    nacov_gr = nacov_gr,
    nbcov_gr = nbcov_gr,
    nccov_gr = nccov_gr,
    ndcov_gr = ndcov_gr,
    acovcoefnames_gr = acovcoefnames_gr,
    bcovcoefnames_gr = bcovcoefnames_gr,
    ccovcoefnames_gr = ccovcoefnames_gr,
    dcovcoefnames_gr = dcovcoefnames_gr,
    ndparcov = ndparcov,
    dparcovcoefnames = dparcovcoefnames,
    group_arg_groupvar = group_arg_groupvar,
    multivariate_rescor = multivariate_rescor,
    univariate_by_by = univariate_by_by,
    set_higher_levels = set_higher_levels,
    hierarchical_gr_names = hierarchical_gr_names,
    covariates_ = covariates_,
    lm_a_all = lm_a_all,
    lm_b_all = lm_b_all,
    lm_c_all = lm_c_all,
    lm_d_all = lm_d_all,
    lm_a = lm_a,
    lm_b = lm_b,
    lm_c = lm_c,
    lm_d = lm_d,
    lm_a_cov = lm_a_cov,
    lm_b_cov = lm_b_cov,
    lm_c_cov = lm_c_cov,
    lm_d_cov = lm_d_cov,
    lm_s = lm_s,
    lm_s_cov = lm_s_cov,
    lm_s_all = lm_s_all,
    lm_sdx = lm_sdx,
    lm_sdx_cov = lm_sdx_cov,
    lm_sdx_all = lm_sdx_all,
    lm_rsd = lm_rsd,
    lme_sd_a = lme_sd_a,
    lme_rsd = lme_rsd
  )
  
  # list_out <<- list_out
  # print(a_formula_grsi)
   # print(aform)
   # print(bform); stop()
  
  attr(bform, "list_out") <- as.list(list_out)
  
  return(bform)
}
