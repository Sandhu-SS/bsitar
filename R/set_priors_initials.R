

#' Set priors and initials for fitting Bayesian SITAR growth curve model
#'
#' @description The \code{set_priors_initials}) sets priors and initials values
#'   which are passed from the [bsitar::bsitar] function to
#'   \code{set_priors_initials}. For univariate-by- subgroup model (specified by
#'   using the \code{univariate_by}) and multivariate model (specified by using
#'   the \code{multivariate}), each argument is automatically matched with the
#'   sub-model(s).
#'
#'
#' @param a_prior_beta Set priors on the the fixed effect \code{a} parameter.
#'   See [bsitar::bsitar()] function, \code{a_prior_beta} for details.
#' @param b_prior_beta Set priors on the the fixed effect \code{b} parameter.
#'   See [bsitar::bsitar()] function, \code{b_prior_beta} for details.
#' @param c_prior_beta Set priors on the the fixed effect \code{c} parameter.
#'   See [bsitar::bsitar()] function, \code{c_prior_beta} for details.
#' @param d_prior_beta Set priors on the the fixed effect \code{d} parameter.
#'   See [bsitar::bsitar()] function, \code{d_prior_beta} for details.
#' @param s_prior_beta Set priors on the the fixed effect \code{s} parameter
#'   (i.e., spline coeficients). See [bsitar::bsitar()] function,
#'   \code{s_prior_beta} for details.
#' @param a_cov_prior_beta Set priors on the the covariate(s) for the fixed
#'   effect \code{a} parameter. See [bsitar::bsitar()] function,
#'   \code{a_cov_prior_beta} for details.
#' @param b_cov_prior_beta Set priors on the the covariate(s) for the fixed
#'   effect \code{b} parameter. See [bsitar::bsitar()] function,
#'   \code{b_cov_prior_beta} for details.
#' @param c_cov_prior_beta Set priors on the the covariate(s) for the fixed
#'   effect \code{c} parameter. See [bsitar::bsitar()] function,
#'   \code{c_cov_prior_beta} for details.
#' @param d_cov_prior_beta Set priors on the the covariate(s) for the fixed
#'   effect \code{d} parameter. See [bsitar::bsitar()] function,
#'   \code{d_cov_prior_beta} for details.
#' @param s_cov_prior_beta Set priors on the the covariate(s) for the fixed
#'   effect \code{s} parameter (i.e., spline coeficients). See
#'   [bsitar::bsitar()] function, \code{s_cov_prior_beta} for details.
#' @param a_prior_sd Set priors on the the group effect (i.e., random) effect
#'   for \code{a} parameter. See [bsitar::bsitar()] function, \code{a_prior_sd}
#'   for details.
#' @param b_prior_sd Set priors on the the group effect (i.e., random) effect
#'   for \code{b} parameter. See [bsitar::bsitar()] function, \code{b_prior_sd}
#'   for details.
#' @param c_prior_sd Set priors on the the group effect (i.e., random) effect
#'   for \code{c} parameter. See [bsitar::bsitar()] function, \code{c_prior_sd}
#'   for details.
#' @param d_prior_sd Set priors on the the group effect (i.e., random) effect
#'   for \code{d} parameter. See [bsitar::bsitar()] function, \code{d_prior_sd}
#'   for details.
#' @param a_cov_prior_sd Set priors on the the covariate(s) for the group effect
#'   (i.e., random) effect for \code{a} parameter. See [bsitar::bsitar()]
#'   function, \code{a_cov_prior_sd} for details.
#' @param b_cov_prior_sd Set priors on the the covariate(s) for the group effect
#'   (i.e., random) effect for \code{b} parameter. See [bsitar::bsitar()]
#'   function, \code{b_cov_prior_sd} for details.
#' @param c_cov_prior_sd Set priors on the the covariate(s) for the group effect
#'   (i.e., random) effect for \code{c} parameter. See [bsitar::bsitar()]
#'   function, \code{c_cov_prior_sd} for details.
#' @param d_cov_prior_sd Set priors on the the covariate(s) for the group effect
#'   (i.e., random) effect for \code{d} parameter. See [bsitar::bsitar()]
#'   function, \code{a_cov_prior_sd} for details.
#' @param gr_prior_cor Set priors on the the group level correlation parameter.
#'   See [bsitar::bsitar()] function, \code{d_prior_beta} for details.
#' @param rsd_prior_sigma Set priors on the the residual standared deviation
#'   \code{sigma} parameter. See [bsitar::bsitar()] function,
#'   \code{a_prior_beta} for details.
#' @param dpar_prior_sigma Set priors on the distributional \code{sigma}
#'   parameter (which is same as residual standared deviation for Gaussian
#'   distribution). See [bsitar::bsitar()] function, \code{a_prior_beta} for
#'   details.
#' @param dpar_cov_prior_sigma Set priors on the covariates for the
#'   distributional \code{sigma} parameter (which is same as residual standard
#'   deviation for Gaussian distribution). See [bsitar::bsitar()] function,
#'   \code{a_prior_beta} for details.
#' @param autocor_prior_acor Set priors on the the autocorrelation parameters
#'   \code{ar}, \code{ma} and \code{arma}. See [bsitar::bsitar()] function,
#'   \code{a_prior_beta} for details.
#' @param mvr_prior_rescor Set priors on the the residual correlation parameter
#'   for multivariate model. See [bsitar::bsitar()] function,
#'   \code{d_prior_beta} for details.
#' @param prior_data An optional argument as named list to pass value for prior.
#'   See [bsitar::bsitar()] function, \code{prior_data} for details.
#' @param prior_data_internal An internal data frame (named list) to pass on the
#'   relevant information on priors from the [bsitar::bsitar()] function to the
#'   \code{set_priors_initials}.
#' @param prior_args_internal An internal argument list that is passed from the
#'   [bsitar::bsitar()] function to the \code{set_priors_initials} and is used
#'   for setting the priors.
#' @param init_arguments A list containing all the init arguments specified in
#'   the [bsitar::bsitar()] function and now passed on to the
#'   \code{set_priors_initials}.
#' @param init_data An optional data argument (named list) used to pass initial
#'   values. See [bsitar::bsitar()] function, \code{prior_data} for details.
#' @param init_data_internal An internal data frame (named list) to pass on the
#'   relevant information on initials from the [bsitar::bsitar()] function to
#'   the \code{set_priors_initials}.
#' @param init_args_internal An internal argument list that is passed from the
#'   [bsitar::bsitar()] function to the \code{set_priors_initials} and is used
#'   for setting the initials.
#'
#' @return An object of class \code{brmsprior} (See \code{brmsprior}). In
#'   addition to the priors, the returned object contains a list of initial
#'   values.
#'
#' @export
#' 
set_priors_initials <- function(a_prior_beta,
                                b_prior_beta,
                                c_prior_beta,
                                d_prior_beta,
                                s_prior_beta,
                                a_cov_prior_beta,
                                b_cov_prior_beta,
                                c_cov_prior_beta,
                                d_cov_prior_beta,
                                s_cov_prior_beta,
                                a_prior_sd,
                                b_prior_sd,
                                c_prior_sd,
                                d_prior_sd,
                                a_cov_prior_sd,
                                b_cov_prior_sd,
                                c_cov_prior_sd,
                                d_cov_prior_sd,
                                gr_prior_cor,
                                rsd_prior_sigma,
                                dpar_prior_sigma,
                                dpar_cov_prior_sigma,
                                autocor_prior_acor,
                                mvr_prior_rescor,
                                prior_data,
                                prior_data_internal,
                                prior_args_internal,
                                init_arguments,
                                init_data,
                                init_data_internal,
                                init_args_internal) {
  eout <- list2env(prior_data_internal)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  
  eout <- list2env(prior_args_internal)
  for (eoutii in names(eout)) {
    assign(eoutii, eout[[eoutii]])
  }
  
  group <- group_arg_groupvar
  
  ept <- function(x)
    eval(parse(text = x), envir = parent.frame())
  
  normalize <- ept(normalize)
  
  if (resp == "")
    resp_ <- ""
  if (resp != "")
    resp_ <- paste0("_", resp)
  
  
  if (is.null(autocor_formi)) {
    autocor_prior_acor <- NULL
  }
  
  getArgNames <-
    function(value)
      formalArgs(deparse(substitute(value)[[1]]))
  
  mvar <- multivariate$mvar
  
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
  
  
  if (!multivariate$rescor) {
    mvr_prior_rescor <- NULL
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
  
  if (is.null(getcovlist(a_formulasi)))
    a_cov_prior_beta <- NULL
  if (is.null(getcovlist(b_formulasi)))
    b_cov_prior_beta <- NULL
  if (is.null(getcovlist(c_formulasi)))
    c_cov_prior_beta <- NULL
  if (is.null(getcovlist(d_formulasi)))
    d_cov_prior_beta <- NULL
  if (is.null(getcovlist(s_formulasi)))
    s_cov_prior_beta <- NULL
  
  
  if (!grepl("a", fixedsi, fixed = T))
    a_prior_beta <- NULL
  if (!grepl("b", fixedsi, fixed = T))
    b_prior_beta <- NULL
  if (!grepl("c", fixedsi, fixed = T))
    c_prior_beta <- NULL
  if (!grepl("d", fixedsi, fixed = T))
    d_prior_beta <- NULL
  
  if (!grepl("a", randomsi, fixed = T))
    a_prior_sd <- NULL
  if (!grepl("b", randomsi, fixed = T))
    b_prior_sd <- NULL
  if (!grepl("c", randomsi, fixed = T))
    c_prior_sd <- NULL
  if (!grepl("d", randomsi, fixed = T))
    d_prior_sd <- NULL
  
  if (is.null(getcovlist(a_formula_grsi)))
    a_cov_prior_sd <- NULL
  if (is.null(getcovlist(b_formula_grsi)))
    b_cov_prior_sd <- NULL
  if (is.null(getcovlist(c_formula_grsi)))
    c_cov_prior_sd <- NULL
  if (is.null(getcovlist(d_formula_grsi)))
    d_cov_prior_sd <- NULL
  
  if (!is.null(a_cov_prior_beta))
    nacov <- length(acovcoefnames)
  else
    nacov <- NULL
  if (!is.null(b_cov_prior_beta))
    nbcov <- length(bcovcoefnames)
  else
    nbcov <- NULL
  if (!is.null(c_cov_prior_beta))
    nccov <- length(ccovcoefnames)
  else
    nccov <- NULL
  if (!is.null(d_cov_prior_beta))
    ndcov <- length(dcovcoefnames)
  else
    ndcov <- NULL
  if (!is.null(s_cov_prior_beta))
    nscov <- length(scovcoefnames)
  else
    nscov <- NULL
  
  if (!is.null(a_cov_prior_sd))
    nacov_gr <- length(acovcoefnames_gr)
  else
    nacov_gr <- NULL
  if (!is.null(b_cov_prior_sd))
    nbcov_gr <- length(bcovcoefnames_gr)
  else
    nbcov_gr <- NULL
  if (!is.null(c_cov_prior_sd))
    nccov_gr <- length(ccovcoefnames_gr)
  else
    nccov_gr <- NULL
  if (!is.null(d_cov_prior_sd))
    ndcov_gr <- length(dcovcoefnames_gr)
  else
    ndcov_gr <- NULL
  
  
  if (!is.null(dpar_formulasi)) {
    if (!grepl("lf\\(", dpar_formulasi) &
        !grepl("nlf\\(", dpar_formulasi)) {
      dpar_covi_mat_form <- dpar_formulasi
    } else {
      dpar_covi_mat_form <-
        gsub("\\(|)", "", strsplit(dpar_formulasi, "~")[[1]][2])
      dpar_covi_mat_form <- paste0("~", dpar_covi_mat_form)
    }
    dpar_covi_mat_form <- strsplit(dpar_covi_mat_form, ",")[[1]][1]
  }
  
  if (is.null(dpar_formulasi)) {
    dpar_covi_mat_form <- NULL
  }
  
  # rsd_prior_sigma and dpar_prior are mutually exclusive
  
  if (!is.null(dpar_formulasi)) {
    rsd_prior_sigma <- NULL
  } else {
    dpar_prior_sigma <- dpar_cov_prior_sigma <- NULL
  }
  
  
  # no prior if no lf | nlf(sigma ~
  if (!is.null(dpar_formulasi)) {
    if (!grepl("^lf\\(", dpar_formulasi) &
        !grepl("^nlf\\(", dpar_formulasi)) {
      dpar_prior_sigma <- dpar_cov_prior_sigma <- NULL
    }
  }
  
  
  
  
  if (grepl("~0", a_formulasi, fixed = T)) {
    a_form_0 <- TRUE
    a_cov_prior_beta <- NULL
  } else {
    a_form_0 <- FALSE
  }
  
  if (grepl("~0", b_formulasi, fixed = T)) {
    b_form_0 <- TRUE
    b_cov_prior_beta <- NULL
  } else {
    b_form_0 <- FALSE
  }
  
  if (grepl("~0", c_formulasi, fixed = T)) {
    c_form_0 <- TRUE
    c_cov_prior_beta <- NULL
  } else {
    c_form_0 <- FALSE
  }
  
  if (grepl("~0", d_formulasi, fixed = T)) {
    d_form_0 <- TRUE
    d_cov_prior_beta <- NULL
  } else {
    d_form_0 <- FALSE
  }
  
  
  if (grepl("~0", s_formulasi, fixed = T)) {
    s_form_0 <- TRUE
    s_cov_prior_beta <- NULL
  } else {
    s_form_0 <- FALSE
  }
  
  
  if (grepl("~0", a_formula_grsi, fixed = T)) {
    a_form_0_gr <- TRUE
    a_cov_prior_sd <- NULL
  } else {
    a_form_0_gr <- FALSE
  }
  
  if (grepl("~0", b_formula_grsi, fixed = T)) {
    b_form_0_gr <- TRUE
    b_cov_prior_sd <- NULL
  } else {
    b_form_0_gr <- FALSE
  }
  
  if (grepl("~0", c_formula_grsi, fixed = T)) {
    c_form_0_gr <- TRUE
    c_cov_prior_sd <- NULL
  } else {
    c_form_0_gr <- FALSE
  }
  
  if (grepl("~0", d_formula_grsi, fixed = T)) {
    d_form_0_gr <- TRUE
    d_cov_prior_sd <- NULL
  } else {
    d_form_0_gr <- FALSE
  }
  
  
  if (!is.null(dpar_formulasi)) {
    if (grepl("^~1$", dpar_covi_mat_form, fixed = F)) {
      dpar_intercept_only <- TRUE
    } else {
      dpar_intercept_only <- FALSE
    }
    if (!is.null(dpar_covi_mat_form) &
        grepl("~0", dpar_covi_mat_form, fixed = T)) {
      sigma_form_0 <- TRUE
      dpar_cov_prior_sigma <- NULL
    } else {
      sigma_form_0 <- FALSE
      if (grepl("^~1$", dpar_covi_mat_form, fixed = F)) {
        dpar_cov_prior_sigma <- NULL
      }
    }
  }
  
  if (is.null(dpar_formulasi)) {
    sigma_form_0 <- FALSE
    dpar_cov_prior_sigma <- NULL
  }
  
  
  if ((is.na(univariate_by$by) |
       univariate_by$by == "NA") & !mvar & !abccorr) {
    gr_prior_cor <- NULL
  }
  
  if (is.null(randomsi[[1]]))
    gr_prior_cor <- NULL
  if (nabcrei == 1)
    gr_prior_cor <- NULL
  
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (uvarabccorr) {
      gr_prior_cor <- gr_prior_cor
    } else {
      gr_prior_cor <- NULL
    }
  }
  
  if (mvar && mvarccorr == "none") {
    gr_prior_cor <- NULL
  }
  
  
  # evaluate prior arguments
  eval_prior_args <- function(x, ...) {
    x_org <- x
    if (grepl("_beta", x))
      class <- 'b'
    if (grepl("_sd", x))
      class <- 'sd'
    if (grepl("rsd_", x) & grepl("_sigma", x))
      class <- 'sigma'
    if (grepl("_cor$", x))
      class <- 'cor'
    if (grepl("dpar", x) &
        grepl("_sigma", x))
      class <- '' # no class sigma if sigma ~ .
    
    if (grepl("_rescor$", x))
      class <- 'rescor'
    
    nlpar <- ""
    if (grepl("a_", x))
      nlpar <- 'a'
    if (grepl("b_", x))
      nlpar <- 'b'
    if (grepl("c_", x))
      nlpar <- 'c'
    if (grepl("d_", x))
      nlpar <- 'd'
    if (grepl("s_", x))
      nlpar <- 's'
    dpar <- ""
    if (grepl("dpar", x) &
        grepl("_sigma", x) & !grepl("rsd_", x))
      dpar <- "dpar"
    
    cov_dpar <- ""
    if (grepl("dpar_cov", x))
      dpar <- paste0(dpar, "_cov")
    
    dpar_cov <- dpar
    
    cov_nlpar <- ""
    if (grepl("a_cov", x))
      cov_nlpar <- 'a'
    if (grepl("b_cov", x))
      cov_nlpar <- 'b'
    if (grepl("c_cov", x))
      cov_nlpar <- 'c'
    if (grepl("d_cov", x))
      cov_nlpar <- 'd'
    if (grepl("s_cov", x))
      cov_nlpar <- 's'
    
    if (!is.null(dpar_prior_sigma) |
        !is.null(dpar_cov_prior_sigma)) {
      ndparcov <- length(dparcovcoefnames) - 1
    } else {
      ndparcov <- NULL
    }
    
    setautocorr <- FALSE
    if (grepl("autocor_", x) & grepl("_acor", x)) {
      setautocorr <- TRUE
      acorclass <- gsub("~", "", autocor_formi, fixed = T)
      acorclass <- strsplit(acorclass, "\\(")[[1]][1]
      if (acorclass == "arma")
        class <- 'ma'
      if (acorclass == "ar")
        class <- 'ar'
      if (acorclass == "ma")
        class <- 'ma'
      if (acorclass == "cosy")
        class <- 'cosy'
      if (acorclass == "car")
        class <- 'car'
      if (acorclass == "lagsar")
        class <- 'lagsar'
      if (acorclass == "errorsar")
        class <- 'errorsar'
    }
    
    
    if (class == "b" & nlpar == 'a') {
      if (a_form_0) {
        nrep_of_parms <- length(acovcoefnames)
      } else {
        if (!grepl("a_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("a_cov", x)) {
          nrep_of_parms <- length(acovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'b') {
      if (b_form_0) {
        nrep_of_parms <- length(bcovcoefnames)
      } else {
        if (!grepl("b_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("b_cov", x)) {
          nrep_of_parms <- length(bcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'c') {
      if (c_form_0) {
        nrep_of_parms <- length(ccovcoefnames)
      } else {
        if (!grepl("c_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("c_cov", x)) {
          nrep_of_parms <- length(ccovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 'd') {
      if (d_form_0) {
        nrep_of_parms <- length(dcovcoefnames)
      } else {
        if (!grepl("d_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("d_cov", x)) {
          nrep_of_parms <- length(dcovcoefnames) - 1
        }
      }
    } else if (class == "b" & nlpar == 's') {
      if (s_form_0) {
        nrep_of_parms <- df * length(scovcoefnames)
      } else {
        if (!grepl("s_cov", x)) {
          nrep_of_parms <- df
        } else if (grepl("s_cov", x)) {
          nrep_of_parms <- df * (length(scovcoefnames) - 1)
        }
      }
    } else if (class == "sd" & nlpar == 'a') {
      if (a_form_0_gr) {
        nrep_of_parms <- length(acovcoefnames_gr)
      } else {
        if (!grepl("a_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("a_cov", x)) {
          nrep_of_parms <- length(acovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'b') {
      if (b_form_0_gr) {
        nrep_of_parms <- length(bcovcoefnames_gr)
      } else {
        if (!grepl("b_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("b_cov", x)) {
          nrep_of_parms <- length(bcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'c') {
      if (c_form_0_gr) {
        nrep_of_parms <- length(ccovcoefnames_gr)
      } else {
        if (!grepl("c_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("c_cov", x)) {
          nrep_of_parms <- length(ccovcoefnames_gr) - 1
        }
      }
    } else if (class == "sd" & nlpar == 'd') {
      if (d_form_0_gr) {
        nrep_of_parms <- length(dcovcoefnames_gr)
      } else {
        if (!grepl("d_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("d_cov", x)) {
          nrep_of_parms <- length(dcovcoefnames_gr) - 1
        }
      }
    } else if (class == "sigma" &
               (class != 'b' |
                class != 'cor') & is.null(ndparcov)) {
      nrep_of_parms <- 1
    } else if (class == "" &
               (class != 'b' |
                class != 'cor') & !is.null(ndparcov)) {
      if (sigma_form_0) {
        nrep_of_parms <- length(dparcovcoefnames)
      } else {
        if (!grepl("dpar_cov", x)) {
          nrep_of_parms <- 1
        } else if (grepl("dpar_cov", x)) {
          nrep_of_parms <- length(dparcovcoefnames) - 1
        }
      }
    } else if (setautocorr &
               class == "" &
               (class != 'b' |
                class != 'cor') & is.null(ndparcov)) {
      nrep_of_parms <- 1
    } else if (class == "cor" &
               (class != 'b' |
                class != 'sigma') & is.null(ndparcov)) {
      nrep_of_parms <- 1
    } else {
      nrep_of_parms <- 1
    }
    
    
    get_priors_parms <- function(x,
                                 prior_data,
                                 prior_data_internal,
                                 resp,
                                 nlpar,
                                 dpar,
                                 class,
                                 acorclass,
                                 get_priors_parms_args) {
      if (!is.null(get_priors_parms_args)) {
        eout <- list2env(get_priors_parms_args)
        for (eoutii in names(eout)) {
          assign(eoutii, eout[[eoutii]])
        }
      }
      
      x <- gsub("\"", "", gsub("\\s", "", x))
      
      if (resp == "") {
        resp_ <- ""
      } else if (resp != "") {
        resp_ <- paste0("_", resp)
      }
      
      prior_argument <- x
      
      zz <- prior_str_arg <- eval(parse(text = x))
      zz <- strsplit(zz, "\\(")[[1]]
      dist <- zz[1]
      
      
      #################
      
      list_names <-
        c(
          'prior_str_arg',
          'dist',
          'resp',
          'resp_',
          'nlpar',
          'dpar',
          'class',
          'cov_nlpar',
          'cov_dpar',
          'nacov',
          'nbcov',
          'nccov',
          'ndcov',
          'nscov',
          'nacov_gr',
          'nbcov_gr',
          'nccov_gr',
          'ndcov_gr',
          'ndparcov',
          'nabci',
          'nabcrei',
          'fixedsi',
          'randomsi',
          'df',
          'setautocorr',
          'nrep_of_parms',
          'N_J_all',
          'ii',
          'nys',
          'a_form_0',
          'b_form_0',
          'c_form_0',
          'd_form_0',
          's_form_0',
          'a_form_0_gr',
          'b_form_0_gr',
          'c_form_0_gr',
          'd_form_0_gr',
          'sigma_form_0',
          'dpar_covi_mat_form',
          'dpar_formulasi',
          'univariate_by',
          'multivariate',
          'group_arg',
          'setautocorr',
          'acorclass',
          'initsi',
          'normalize',
          'seed',
          'verbose'
        )
      
      
      prior_internal_args <- mget(list_names)
      
      # set to NULL as appropriate to match priors
      for (ip in names(init_arguments)) {
        init_ip <- ip
        if (grepl("_init_", init_ip) &
            !grepl("r_init_z", init_ip)) {
          prior_ip <- gsub("_init_", "_prior_", init_ip)
          # if(!is.null(eval(parse(text = prior_ip)))) {
          xx   <- deparse(prior_ip)
          if (is.null(ept(prior_ip)))
            init_arguments[[ip]] <- 'NULL'
          # }
        }
      }
      
      if (ept(nabcrei) == 0)
        init_arguments[['r_init_z']] <- 'NULL'
      
      out_p_str <- prepare_priors(
        prior_argument,
        prior_data,
        prior_data_internal,
        prior_internal_args,
        init_arguments,
        init_data,
        init_data_internal,
        init_args_internal
      )
      
      
      stanvars_data_in <- out_p_str$stanvars_data
      prior_str_arg    <- out_p_str$prior_str_arg
      lowerbound       <- out_p_str$lowerbound
      upperbound       <- out_p_str$upperbound
      initial_in       <- out_p_str$initial_out
      
      return(
        list(
          dist = dist,
          lowerbound = lowerbound,
          upperbound = upperbound,
          define_ = prior_str_arg,
          stanvars_data_in = stanvars_data_in,
          initial_in = initial_in
        )
      )
    }
    
    
    if (resp != "") {
      for (i in names(prior_data_internal)) {
        names(prior_data_internal)[names(prior_data_internal) == i] <-
          paste0(i, "_", resp)
      }
    }
    
    x_name <- deparse(substitute(x_org))
    
    
    get_priors_parms_args <- list(
      df = df,
      nys = nys,
      univariate_by = univariate_by,
      multivariate = multivariate,
      group_arg = group_arg,
      fixedsi = fixedsi,
      randomsi = randomsi,
      nabci = nabci,
      nabcrei = nabcrei,
      ii = ii,
      N_J_all = N_J_all,
      cov_nlpar = cov_nlpar,
      cov_dpar = cov_dpar,
      nrep_of_parms = nrep_of_parms,
      nacov = nacov,
      nbcov = nbcov,
      nccov = nccov,
      ndcov = ndcov,
      nscov = nscov,
      nacov_gr = nacov_gr,
      nbcov_gr = nbcov_gr,
      nccov_gr = nccov_gr,
      ndcov_gr = ndcov_gr,
      ndparcov = ndparcov,
      a_form_0 = a_form_0,
      b_form_0 = b_form_0,
      c_form_0 = c_form_0,
      d_form_0 = d_form_0,
      s_form_0 = s_form_0,
      a_form_0_gr = a_form_0_gr,
      b_form_0_gr = b_form_0_gr,
      c_form_0_gr = c_form_0_gr,
      d_form_0_gr = d_form_0_gr,
      sigma_form_0 = sigma_form_0,
      dpar_covi_mat_form = dpar_covi_mat_form,
      dpar_formulasi = dpar_formulasi,
      setautocorr = setautocorr,
      initsi = initsi,
      init_arguments = init_arguments,
      init_data = init_data,
      init_data_internal = init_data_internal,
      init_args_internal = init_args_internal,
      normalize = normalize,
      seed = seed,
      verbose = verbose
    )
    
    
    if (setautocorr) {
      if (acorclass == 'arma')
        acorclassclasses <- c("ar", "ma")
      if (acorclass == 'ar')
        acorclassclasses <- c("ar")
      if (acorclass == 'ma')
        acorclassclasses <- c("ma")
      
      priors_arma_c_define <- list()
      for (acorclassi in acorclassclasses) {
        priors_parms <- get_priors_parms(
          x_name,
          prior_data = prior_data,
          prior_data_internal = prior_data_internal,
          resp = resp,
          nlpar = nlpar,
          dpar = dpar,
          class = acorclassi,
          acorclass = acorclass,
          get_priors_parms_args = get_priors_parms_args
        )
        
        priors_arma_c_define[[acorclassi]] <- priors_parms
      }
    } else {
      priors_parms <- get_priors_parms(
        x_name,
        prior_data = prior_data,
        prior_data_internal = prior_data_internal,
        resp = resp,
        nlpar = nlpar,
        dpar = dpar,
        class = class,
        get_priors_parms_args = get_priors_parms_args
      )
    }
    
    
    define_ <- priors_parms$define_
    
    dist <- priors_parms$dist
    
    lowerbound <- priors_parms$lowerbound
    upperbound <- priors_parms$upperbound
    
    stanvars_data_in <- priors_parms$stanvars_data_in
    initial_in       <- priors_parms$initial_in
    
    
    if (class == 'b') {
      # need to remove lb and ub if specifying coef, otherwise
      # Error: Argument 'coef' may not be specified when using boundaries.
      
      mnf <- paste0(nlpar, "_form_0")
      mnc <- paste0("cov_nlpar")
      if (dist != 'uniform' &
          dist != 'lognormal' &
          dist != 'gamma' &
          dist != 'inv_gamma' &
          dist != 'exponential') {
        if (all(is.na(lowerbound)) &
            all(is.na(upperbound))) {
          if (nlpar == 'a')
            coef <- acovcoefnames
          if (nlpar == 'b')
            coef <- bcovcoefnames
          if (nlpar == 'c')
            coef <- ccovcoefnames
          if (nlpar == 'd')
            coef <- dcovcoefnames
          if (nlpar == 's') {
            coef <- scovcoefnames
            if (!s_form_0 & !is.null(nscov))
              coef <- coef[1]
          }
        } else {
          if (nlpar == 'a')
            coef <- rep("", length(acovcoefnames))
          if (nlpar == 'b')
            coef <- rep("", length(bcovcoefnames))
          if (nlpar == 'c')
            coef <- rep("", length(ccovcoefnames))
          if (nlpar == 'd')
            coef <- rep("", length(dcovcoefnames))
          if (nlpar == 's') {
            coef <- rep("", length(scovcoefnames))
            if (!s_form_0 & !is.null(nscov))
              coef <- coef[1]
          }
        }
      }
      
      # nlpar a b c - betas
      if (ept(mnf) & cov_nlpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          setcoef <- coef
        }
        priors_ <-
          prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      if (!ept(mnf) & cov_nlpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          if (ept(mnc) == "") {
            setcoef <- coef[1]
          } else {
            setcoef <- coef
          }
        }
        priors_ <-
          prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      
      
      # nlpar s - betas
      
      if (nlpar == 's' & cov_nlpar == "") {
        nlpar <- paste0(nlpar, 1:df)
        if (grepl("~1", s_formulasi, fixed = T)) {
          if (all(coef == "")) {
            priors_ <- prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
          } else {
            priors_ <-   prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          }
        }
        if (grepl("~0", s_formulasi, fixed = T)) {
          nlpar <- rep(nlpar ,
                       times = 1,
                       each = length(scovcoefnames))
          if (all(coef == "")) {
            priors_ <- prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
          } else {
            priors_ <-   prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          }
        }
      }
      
      
      
      # nlpar cov a - betas
      if (!grepl("~0", a_formulasi, fixed = T)) {
        if (class == 'b' & grepl("a_cov", x) & !is.null(a_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- acovcoefnames
          } else {
            coef <- acovcoefnames[2:length(acovcoefnames)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          
        }
      }
      
      # nlpar cov b - betas
      if (!grepl("~0", b_formulasi, fixed = T)) {
        if (class == 'b' & grepl("b_cov", x) & !is.null(b_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- bcovcoefnames
          } else {
            coef <- bcovcoefnames[2:length(bcovcoefnames)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # nlpar cov c - betas
      if (!grepl("~0", c_formulasi, fixed = T)) {
        if (class == 'b' & grepl("c_cov", x) & !is.null(c_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- ccovcoefnames
          } else {
            coef <- ccovcoefnames[2:length(ccovcoefnames)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov d - betas
      if (!grepl("~0", d_formulasi, fixed = T)) {
        if (class == 'b' & grepl("d_cov", x) & !is.null(d_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- dcovcoefnames
          } else {
            coef <- dcovcoefnames[2:length(dcovcoefnames)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov s - betas
      if (!grepl("~0", s_formulasi, fixed = T)) {
        if (class == 'b' & grepl("s_cov", x) & !is.null(s_cov_prior_beta)) {
          if (ept(mnf)) {
            coef <- scovcoefnames
          } else {
            coef <- scovcoefnames[2:length(scovcoefnames)]
          }
          nlpar <- paste0(nlpar, 1:df)
          nlpar <- rep(nlpar, times = length(coef), each = 1)
          coef <- rep(coef , times = 1, each = df)
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
          
        }
      }
    } # if(class == 'b')
    
    
    
    
    
    if (class == 'sd') {
      mnf <- paste0(nlpar, "_form_0_gr")
      mnc <- paste0("cov_nlpar")
      
      if (nlpar == 'a')
        coef <- acovcoefnames_gr
      if (nlpar == 'b')
        coef <- bcovcoefnames_gr
      if (nlpar == 'c')
        coef <- ccovcoefnames_gr
      if (nlpar == 'd')
        coef <- dcovcoefnames_gr
      
      # nlpar a b c - sd
      
      if (ept(mnf) & cov_nlpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          setcoef <- coef
        }
        priors_ <-
          prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            group = group,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      if (!ept(mnf) & cov_nlpar == "") {
        if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
          define_ <- unique(define_)
          lowerbound <- unique(lowerbound)
          upperbound <- unique(upperbound)
          setcoef <- ""
        } else {
          if (ept(mnc) == "") {
            setcoef <- coef[1]
          } else {
            setcoef <- coef[-1]
          }
        }
        priors_ <-
          prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = setcoef,
            group = group,
            resp = resp,
            dpar = dpar,
            lb = lowerbound,
            ub = upperbound
          )
      }
      
      
      # nlpar cov a - sd
      if (!grepl("~0", a_formula_grsi, fixed = T)) {
        if (class == 'sd' & grepl("a_cov", x) & !is.null(a_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- acovcoefnames_gr
          } else {
            coef <- acovcoefnames_gr[2:length(acovcoefnames_gr)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # nlpar cov b - sd
      if (!grepl("~0", b_formula_grsi, fixed = T)) {
        if (class == 'sd' & grepl("b_cov", x) & !is.null(b_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- bcovcoefnames_gr
          } else {
            coef <- bcovcoefnames_gr[2:length(bcovcoefnames_gr)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      # nlpar cov c - sd
      if (!grepl("~0", c_formula_grsi, fixed = T)) {
        if (class == 'sd' & grepl("c_cov", x) & !is.null(c_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- ccovcoefnames_gr
          } else {
            coef <- ccovcoefnames_gr[2:length(ccovcoefnames_gr)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      # nlpar cov d - sd
      if (!grepl("~0", d_formula_grsi, fixed = T)) {
        if (class == 'sd' & grepl("d_cov", x) & !is.null(d_cov_prior_sd)) {
          if (ept(mnf)) {
            coef <- dcovcoefnames_gr
          } else {
            coef <- dcovcoefnames_gr[2:length(dcovcoefnames_gr)]
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              group = group,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
    }
    
    
    
    # correlation priors (lkj)
    
    # currently brms does not allow setting separate ljk prior for
    # subset and multivariate
    # removing resp leads to duplicate priors, so need to set only once
    
    if (class == 'cor') {
      if (ii == 1) {
        priors_ <-  prior_string(define_,
                                 class = class,
                                 group = group,
                                 dpar = dpar) # resp = resp,
      } else {
        priors_ <- ""
      }
    }
    
    
    if (class == 'rescor') {
      if (ii == 1) {
        priors_ <-  prior_string(define_,
                                 class = class,
                                 group = "",
                                 dpar = "")
      } else {
        priors_ <- ""
      }
    }
    
    
    
    # residual standard deviation (sigma) prior
    
    if (class == 'sigma' & dpar == "") {
      priors_ <-  prior_string(
        define_,
        class = class,
        lb = lowerbound,
        ub = upperbound,
        resp = resp,
        dpar = dpar
      )
      
    }
    
    # residual standard deviation (sigma) prior - dpar_formula formulation
    
    if (class == "" & !is.null(dpar_formulasi)) {
      # need to remove lb and ub if specifying coef, otherwise
      # Error: Argument 'coef' may not be specified when using boundaries.
      
      dpar <- 'sigma'
      
      if (!is.null(dpar_covi_mat_form) &
          !grepl("~1$", dpar_covi_mat_form, fixed = F)) {
        class <- 'b'
        mnf <- paste0(dpar, "_form_0")
        mnc <- paste0("dpar_cov")
        
        if (dist != 'uniform' &
            dist != 'lognormal' &
            dist != 'gamma' &
            dist != 'inv_gamma' &
            dist != 'exponential') {
          if (all(is.na(lowerbound)) &
              all(is.na(upperbound))) {
            if (grepl("^lf\\(", dpar_formulasi)) {
              if (grepl("cmc=F", dpar_formulasi) |
                  grepl("cmc=FALSE", dpar_formulasi)) {
                coef <- c("", dparcovcoefnames[2:length(dparcovcoefnames)])
                define_ <- c("", define_[2:length(define_)])
              } else {
                coef <- dparcovcoefnames
              }
            }
            if (!grepl("^lf\\(", dpar_formulasi) |
                !grepl("^nlf\\(", dpar_formulasi)) {
              coef <- dparcovcoefnames
            }
          } else {
            coef <- rep("", length(dparcovcoefnames))
          }
        }
        
        if (ept(mnf)) {
          if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
            define_ <- unique(define_)
            lowerbound <- unique(lowerbound)
            upperbound <- unique(upperbound)
            setcoef <- ""
          } else {
            setcoef <- coef
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = setcoef,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
        }
        
        if (!ept(mnf)) {
          if (!any(is.na(lowerbound)) | !any(is.na(upperbound))) {
            define_ <- unique(define_)
            lowerbound <- unique(lowerbound)
            upperbound <- unique(upperbound)
            setcoef <- ""
          } else {
            if (ept(mnc) != "") {
              setcoef <- coef[1]
            } else {
              setcoef <- coef
            }
          }
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = setcoef,
              resp = resp,
              dpar = dpar,
              lb = lowerbound,
              ub = upperbound
            )
        }
        
        
        
        
        # residual standard deviation (sigma) covariate prior -
        # dpar_formula formulation
        
        if (!is.null(dpar_covi_mat_form) &
            grepl("~1", dpar_covi_mat_form, fixed = T) &
            !grepl("~1$", dpar_covi_mat_form, fixed = T) &
            !is.null(dpar_cov_prior_sigma)) {
          if (grepl("dpar", x) & !grepl("dpar_cov", x)) {
            if (grepl("^lf\\(", dpar_formulasi)) {
              if (grepl("center=T", dpar_formulasi) |
                  grepl("center=TRUE", dpar_formulasi)) {
                class <- dparcovcoefnames[1]
                coef  <- ""
              } else {
                class <- 'b'
                coef  <- dparcovcoefnames[1]
              }
            } else if (!grepl("^lf\\(", dpar_formulasi) |
                       !grepl("^nlf\\(", dpar_formulasi)) {
              class <- dparcovcoefnames[1]
              coef  <- ""
            }
          }
          
          if (!grepl("dpar", x) & grepl("dpar_cov", x)) {
            class <- 'b'
            coef <- dparcovcoefnames
          }
          
          if (class == 'b') {
            if (grepl("center=T", dpar_formulasi) |
                grepl("center=TRUE", dpar_formulasi)) {
              coef <- coef[-1]
            } else {
              coef <- coef
            }
          }
          if (!grepl("center=T", dpar_formulasi) &
              !grepl("center=TRUE", dpar_formulasi)) {
            if (grepl("dpar_cov", x))
              coef <- coef[-1]
          }
          
          priors_ <-
            prior_string(
              define_,
              class = class,
              nlpar = nlpar,
              coef = coef,
              resp = resp,
              dpar = dpar
            )
        }
      }
      
      
      if (!is.null(dpar_covi_mat_form) &
          grepl("~1$", dpar_covi_mat_form, fixed = F)) {
        if (grepl("center=T", dpar_formulasi) |
            grepl("center=TRUE", dpar_formulasi)) {
          class <- dparcovcoefnames[1]
          coef  <- ""
        } else {
          class <- 'b'
          coef  <- dparcovcoefnames[1]
        }
        priors_ <-
          prior_string(
            define_,
            class = class,
            nlpar = nlpar,
            coef = coef,
            resp = resp,
            dpar = dpar
          )
      }
    }
    
    
    # autocorrelation priors
    
    if (setautocorr) {
      coef <- ""
      if (acorclass == 'arma') {
        acorclassclasses <- c("ar", "ma")
        priors_arma_c <- list()
        stanvars_data_in_c <- list()
        for (acorclassi in acorclassclasses) {
          define_ <- priors_arma_c_define[[acorclassi]]$define_
          lowerbound <-
            priors_arma_c_define[[acorclassi]]$lowerbound
          upperbound <-
            priors_arma_c_define[[acorclassi]]$upperbound
          priors_temp <-  prior_string(
            define_,
            class = acorclassi,
            lb = lowerbound,
            ub = upperbound,
            coef = coef,
            resp = resp,
            dpar = dpar
          )
          priors_arma_c[[acorclassi]] <- priors_temp
          stanvars_data_in_c[[acorclassi]] <-
            priors_arma_c_define[[acorclassi]]$stanvars_data_in
        }
        priors_ <- priors_arma_c %>% do.call(rbind, .)
        stanvars_data_in <- stanvars_data_in_c %>% do.call(rbind, .)
      } else {
        priors_ <-  prior_string(
          define_,
          class = class,
          lb = lowerbound,
          ub = upperbound,
          coef = coef,
          resp = resp,
          dpar = dpar
        )
      }
    }
    out_pr <-
      list(
        priors_ = priors_,
        stanvars_data_in = stanvars_data_in,
        initial_in = initial_in
      )
    return(out_pr)
  } 
  
  
  
  # use following custom order
  # This order ensures that corresponding initial arguments are matched
  # with the sequence of prior argument evaluation
  
  custom_order_prior <- c(
    'a_prior_beta',
    'a_cov_prior_beta',
    'b_prior_beta',
    'b_cov_prior_beta',
    'c_prior_beta',
    'c_cov_prior_beta',
    'd_prior_beta',
    'd_cov_prior_beta',
    's_prior_beta',
    's_cov_prior_beta',
    'a_prior_sd',
    'a_cov_prior_sd',
    'b_prior_sd',
    'b_cov_prior_sd',
    'c_prior_sd',
    'c_cov_prior_sd',
    'd_prior_sd',
    'd_cov_prior_sd',
    'gr_prior_cor',
    'rsd_prior_sigma',
    'dpar_prior_sigma',
    'dpar_cov_prior_sigma',
    'autocor_prior_acor',
    'mvr_prior_rescor'
  )
  
  
  
  stanvars_data_5 <- list()
  initial_in_data <- list()
  c_priors <- list()
  for (ip in custom_order_prior) {
    if (grepl("_prior_", ip)) {
      if (!is.null(eval(parse(text = ip)))) {
        xx   <- deparse(ip)
        pget <- eval_prior_args(eval(parse(text = xx)))
        zz   <- pget$priors_
        # print(zz)
        stanvars_data_5[[ip]] <- pget$stanvars_data_in
        initial_in_data[[ip]] <- pget$initial_in
        c_priors[[ip]] <- zz
      }
    }
  }
  
  evaluated_priors <- c_priors %>% do.call(rbind, .)
  
  newlist <- c()
  for (i in 1:length(stanvars_data_5)) {
    ttt <- stanvars_data_5[[i]][1:length(stanvars_data_5[[i]])]
    newlist <- c(newlist, unname(ttt))
  }
  
  svardatalistlist <- c()
  for (istanvardata in 1:length(newlist)) {
    svardatalistlist[istanvardata] <-
      paste0("newlist[[", istanvardata, "]]")
  }
  
  stanvars <-
    eval(parse(text = paste(svardatalistlist, collapse = "+")))
  
  attr(evaluated_priors, 'stanvars') <- stanvars
  
  initial_in_datazz <- initial_in_data
  
  if (is.list(initial_in_datazz) & length(initial_in_datazz) == 0) {
    initial_in_datazz <- NULL
  }
  
  if (!is.null(initial_in_datazz)) {
    if (!is.null(gr_prior_cor)) {
      list_ck <- list_ck_ <- list()
      list_ck_rescor <- list()
      ik_j <- ik_j_ <- 0
      for (ik in 1:length(initial_in_datazz)) {
        ik_j <- ik_j + 1
        ik_names <- names(initial_in_datazz[[ik]])
        if (!any(grepl("^L_|^z_|Intercept_sigma|Lrescor", ik_names))) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (any(grepl("^L_|^z_", ik_names))) {
          mn <- 0
          for (ikl in 1:length(grepl("^L_|^z_", ik_names))) {
            mn <- mn + 1
            ik_j_ <- ik_j_ + 1
            list_ck_[[ik_j_]] <- initial_in_datazz[[ik_j]][[mn]]
          }
          names(list_ck_) <- ik_names
        } else if (grepl("^Intercept_sigma", ik_names)) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (multivariate$mvar & multivariate$rescor &
                   grepl("^Lrescor", ik_names)) {
          list_ck_rescor[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck_rescor[[ik]]) <- ik_names
        }
      }
      list_ck <- list_ck[lengths(list_ck) != 0]
      keys    <- unique(unlist(lapply(list_ck, names)))
      list_ck <-
        setNames(do.call(mapply, c(FUN = c, lapply(
          list_ck, `[`, keys
        ))), keys)
      combined_inits <- c(list_ck, list_ck_)
    }
    
    
    if (is.null(gr_prior_cor)) {
      list_ck <- list_ck_z <- list_ck_sd <- list()
      list_ck_rescor <- list()
      ik_j <- ik_j_ <- 0
      for (ik in 1:length(initial_in_datazz)) {
        ik_j <- ik_j + 1
        ik_names <- names(initial_in_datazz[[ik]])
        if (!any(grepl("^L_|^z_|Intercept_sigma|Lrescor", ik_names))) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (any(grepl("^L_|^z_", ik_names))) {
          mn <- 0
          for (ikl in 1:length(grepl("^L_|^z_", ik_names))) {
            mn <- mn + 1
            ik_j_ <- ik_j_ + 1
            if (is.matrix(initial_in_datazz[[ik_j]][[mn]]))
              list_ck_z[[ik_j_]] <- initial_in_datazz[[ik_j]][[mn]]
            if (!is.matrix(initial_in_datazz[[ik_j]][[mn]]))
              list_ck_sd[[ik_j_]] <- initial_in_datazz[[ik_j]][[mn]]
          }
          
          names(list_ck_z) <- ik_names[2]
          names(list_ck_sd) <- ik_names[1]
        } else if (grepl("^Intercept_sigma", ik_names)) {
          list_ck[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck[[ik]]) <- ik_names
        } else if (multivariate$mvar & multivariate$rescor &
                   grepl("^Lrescor", ik_names)) {
          list_ck_rescor[[ik]] <- initial_in_datazz[[ik_j]]
          names(list_ck_rescor[[ik]]) <- ik_names
        }
      }
      list_ck <- list_ck[lengths(list_ck) != 0]
      keys    <- unique(unlist(lapply(list_ck, names)))
      list_ck <-
        setNames(do.call(mapply, c(FUN = c, lapply(
          list_ck, `[`, keys
        ))), keys)
      list_ck_sd <- list_ck_sd[lengths(list_ck_sd) != 0]
      for (list_ck_sd_i in 1:length(list_ck_sd)) {
        if (length(list_ck_sd[[list_ck_sd_i]]) > 1) {
          nami_ <-
            paste0(names(list_ck_sd[[list_ck_sd_i]][1]),
                   "cov",
                   2:length(list_ck_sd[[list_ck_sd_i]]) - 1)
          names(list_ck_sd[[list_ck_sd_i]]) <-
            c(names(list_ck_sd[[list_ck_sd_i]][1]), nami_)
        }
      }
      names(list_ck_sd) <-
        rep(names(list_ck_sd[1]), length(list_ck_sd))
      list_ck_sd2 <- list_ck_sd
      list_ck_z <- list_ck_z[lengths(list_ck_z) != 0]
      list_ck_z2 <- list()
      for (list_ck_i in 1:length(list_ck_z)) {
        addelemnt <-
          strsplit(gsub("\\+", " ", randomsi), " ")[[1]][list_ck_i]
        list_ck_z2[[paste0("z", "_", addelemnt, resp_, list_ck_i)]] <-
          list_ck_z[[list_ck_i]]
        attr(list_ck_z2[[paste0("z", "_", addelemnt, list_ck_i)]], "names") <-
          NULL
      }
      combined_inits <- c(list_ck, list_ck_sd2, list_ck_z2)
    }
    
    
    if (multivariate$mvar & multivariate$rescor) {
      list_ck_rescor <- list_ck_rescor[lengths(list_ck_rescor) != 0]
      list_ck_rescor <- list_ck_rescor[[1]]
      combined_inits <- c(combined_inits, list_ck_rescor)
    }
    
    # convert vector of 's' initials to named individual (s1, s2)

    nlpar_s_init <- paste0('_s', 1:df)
    if (grepl("~0", s_formulasi, fixed = T)) {
      nlpar_s_init <-
        rep(nlpar_s_init ,
            times = 1,
            each = length(scovcoefnames))
    } else if (!grepl("~0", s_formulasi, fixed = T)) {
      nlpar_s_init <- rep(nlpar_s_init , times = length(scovcoefnames))
    }
    
    subset_sparms <-
      combined_inits[grepl(".*_s$", names(combined_inits))]
    subset_sparms_name <- names(subset_sparms)
    subset_sparms_numeric <- subset_sparms[[1]]
    subset_sparms2 <- list()
    subset_sparms2names <- c()
    for (subset_sparmsi in 1:length(subset_sparms_numeric)) {
      subset_sparms_namei <-
        gsub("_s", nlpar_s_init[subset_sparmsi], subset_sparms_name)
      subset_sparms2[[subset_sparms_namei]] <-
        subset_sparms_numeric[subset_sparmsi]
      subset_sparms2names <-
        c(subset_sparms2names, subset_sparms_namei)
    }
    names(subset_sparms_numeric) <- subset_sparms2names
    subset_sparms3 <- list()
    for (isi in 1:df) {
      subset_sparms3[[paste0("b", resp_, "_s", isi)]] <-
        subset_sparms_numeric[grep(paste0("b", resp_, "_s", isi),
                                   names(subset_sparms_numeric))]
    }
    subset_sparms <- subset_sparms3
    subset_sparms <-
      subset_sparms[!names(subset_sparms) %in% subset_sparms_name]
    combined_inits <-
      append(combined_inits, subset_sparms, after = grep(
        paste0("^", subset_sparms_name, "$"),
        names(combined_inits)
      ))
    combined_inits <-
      combined_inits[!names(combined_inits) %in% paste0("",
                                                        subset_sparms_name, "")]
    initials <- combined_inits
  }
  
  
  if (is.null(initial_in_datazz)) {
    initials <- NULL
  }
  
  ###################3

  stanvar_priors_names <- names(stanvars)
  getaux <- "tau"
  stanvar_priors_names_c <- c()
  for (stanvar_priors_namesi in stanvar_priors_names) {
    t <-
      stanvar_priors_namesi[grep(paste0(getaux, '_scale', resp_), 
                                 stanvar_priors_namesi)]
    t <- gsub(paste0('_scale', resp_), "", t, fixed = T)
    stanvar_priors_names_c <- c(stanvar_priors_names_c, t)
  }
  
  add_tau <- list()
  for (stanvar_priors_names_ci in stanvar_priors_names_c) {
    fstandat <-
      unlist(stanvars)[grep(paste0(
        stanvar_priors_names_ci,
        paste0('_scale', resp_, ".sdata")
      ), names(unlist(stanvars)))] %>% as.numeric()
    add_tau[[paste0(stanvar_priors_names_ci, resp_)]] <-
      rep(1, length(fstandat))
  }
  if (length(add_tau) == 0)
    add_tau <- NULL
  
  getaux <- "nu"
  stanvar_priors_names_c <- c()
  for (stanvar_priors_namesi in stanvar_priors_names) {
    t <-
      stanvar_priors_namesi[grep(paste0(getaux, '_scale', resp_), 
                                 stanvar_priors_namesi)]
    t <- gsub(paste0('_scale', resp_), "", t, fixed = T)
    stanvar_priors_names_c <- c(stanvar_priors_names_c, t)
  }
  add_nu <- list()
  for (stanvar_priors_names_ci in stanvar_priors_names_c) {
    add_nu[[paste0(stanvar_priors_names_ci, resp_)]] <-  5
  }
  if (length(add_nu) == 0)
    add_nu <- NULL

  initials <- c(initials, add_tau, add_nu)

  ################
  revSubstr <- function(x_) {
    x__ <- substr(x_, start = 1, stop = 3)
    x___ <- paste0(rev(strsplit(x__, "_")[[1]]), collapse = "_")
    x___ <- gsub(x__, x___, x_)
    x___
  }
  
  tau_nu_init_list <- c(add_tau, add_nu)
  
  if (length(tau_nu_init_list) != 0) {
    names_tau_nu_parms <- names(tau_nu_init_list)
    names_tau_nu_parmsi_c <- c()
    for (names_tau_nu_parmsi in names_tau_nu_parms) {
      plength <- length(tau_nu_init_list[[names_tau_nu_parmsi]])
      revstr <- revSubstr(names_tau_nu_parmsi)
      if (!grepl("^b_b", names_tau_nu_parmsi, fixed = F)) {
        o <-
          paste0("vector[",
                 plength,
                 "]",
                 " ",
                 revstr,
                 " = ",
                 names_tau_nu_parmsi,
                 ";")
        names_tau_nu_parmsi_c <- c(names_tau_nu_parmsi_c, o)
      }
    }
    names_tau_nu_parmsi_c <- names_tau_nu_parmsi_c
    names_tau_nu_parmsi_cc <-
      paste(names_tau_nu_parmsi_c, collapse = "\n")
    
    scode_auxillary <-
      stanvar(scode = names_tau_nu_parmsi_cc,
              block = "genquant",
              position = 'end')
  } else if (length(tau_nu_init_list) == 0) {
    scode_auxillary <- NULL
  }
  
  ##################
  
  attr(evaluated_priors, 'initials') <- initials
  attr(evaluated_priors, 'scode_auxillary') <- scode_auxillary
  
  return(evaluated_priors)
}
