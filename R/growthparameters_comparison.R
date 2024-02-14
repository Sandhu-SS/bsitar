

#' Compare growth parameters
#' 
#'@description The \strong{growthparameters_comparison()} function estimates and
#'  compare growth parameters such as peak growth velocity and the age at peak
#'  growth velocity. This function is a wrapper around the
#'  [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()].
#'  The [marginaleffects::comparisons()] computes unit-level (conditional)
#'  estimates whereas [marginaleffects::avg_comparisons()] return average
#'  (marginal) estimates. A detailed explanation is available
#'  [here](https://marginaleffects.com). Note that for the current use case,
#'  i.e., to estimate and compare growth parameters, the arguments
#'  \code{variables} and \code{comparion} of [marginaleffects::comparisons()]
#'  and [marginaleffects::avg_comparisons()] are modified (see below).
#'
#' @details The \code{growthparameters_comparison} function estimates and
#'   returns the following growth parameters:
#' \itemize{
#'   \item pgv  - peak growth velocity
#'   \item apgv - age at peak growth velocity
#'   \item tgv  - takeoff growth velocity
#'   \item atgv - age at takeoff growth velocity
#'   \item cgv  - cessation growth velocity
#'   \item acgv - age at cessation growth velocity
#' }
#' 
#' The takeoff growth velocity is the lowest velocity just before the peak
#' starts and it indicates the beginning of the pubertal growth spurt. The
#' cessation growth velocity indicates the end of the active pubertal growth
#' spurt and is calculated as some percentage of the peak velocity (\code{pgv}).
#' Typically, a 10 percent of the \code{pgv} is considered as a good indicator
#' of the cessation of the active pubertal growth spurt
#' \insertCite{Anna2022}{bsitar}. The percentage is controlled via the
#' \code{acg_velocity} argument which takes a positive real value bounded
#' between 0 and 1 (default \code{0.1} implying 10 percent). 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param parameter A single character string, or a character vector specifying
#'   the growth parameter(s) to be estimated. Options are \code{'tgv'} (takeoff
#'   growth velocity), \code{'atgv'} (age at takeoff growth velocity),
#'   \code{'pgv'} (peak growth velocity), \code{'apgv'} (age at peak growth
#'   velocity), \code{'cgv'} (cessation growth velocity), and \code{'acgv'} (age
#'   at cessation growth velocity). For \code{parameter = NULL} (default), age
#'   at peak growth velocity (\code{'apgv'}) is estimated.
#' 
#' @param acg_velocity A real number to set the percentage of peak growth growth
#'   velocity as the cessation velocity when estimating the \code{cgv} and
#'   \code{acgv} growth parameters. The \code{acg_velocity} should be greater
#'   than \code{0} and less than \code{1}. The default \code{acg_velocity =
#'   0.10} indicates that a 10 per cent of the peak growth velocity will be used
#'   to get the cessation velocity and the corresponding age at the cessation
#'   velocity. For example if peak growth velocity estimate is \code{10
#'   mm/year}, then cessation growth velocity is \code{1 mm/year}.
#' 
#' @param digits An integer (default \code{2}) to set the decimal places for the
#'   estimated growth parameters. The \code{digits} is passed on to the
#'   [base::round()] function.
#'   
#' @param average A logical to indicate whether to internally call the
#'    [marginaleffects::comparisons()] or the
#'    [marginaleffects::avg_comparisons()] function. If \code{FALSE} (default),
#'    [marginaleffects::comparisons()] is called otherwise
#'    [marginaleffects::avg_comparisons()] when \code{average = TRUE}.
#'
#' @param variables For estimating growth parameters in the current use case,
#'   the \code{variables} is the level 1 predictor such as
#'   \code{age}/\code{time}. The \code{variables} is a named list where value is
#'   set via the \code{esp} argument (default 1e-6). If \code{NULL}, the
#'   \code{variables} is set internally by retrieving the relevant information
#'   from the \code{model}. Otherwise, user can define it as follows:
#'   \code{variables = list('age' = 1e-6)}.
#' 
#' @param comparison For estimating growth parameters in the current use case,
#'   options allowed for the \code{comparison} are \code{'difference'} and
#'   \code{'differenceavg'}. Note that \code{comparison} is a placeholder and is
#'   only used to setup the the internal function that estimates
#'   \code{'parameter'} via [sitar::getPeak()], [sitar::getTakeoff()] and
#'   [sitar::getTrough()] functions to estimate various growth parameters.
#'   Options \code{'difference'} and \code{'differenceavg'} are internally
#'   restructured according to the user specified \code{hypothesis} argument.
#'   
#' @param reformat A logical (default \code{TRUE}) to reformat the  output
#'   returned by the \code{marginaleffects} as a data.frame with column names
#'   re-defined as follows: \code{conf.low} as \code{Q2.5}, and \code{conf.high}
#'   as \code{Q97.5} (assuming that \code{conf_int = 0.95}). Also, following
#'   columns are dropped from the data frame: \code{term}, \code{contrast},
#'   \code{tmp_idx}, \code{predicted_lo}, \code{predicted_hi}, \code{predicted}.
#' 
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  marginaleffects::comparisons
#' @inheritParams  marginaleffects::avg_comparisons
#'
#' @return A data frame objects with estimates and CIs for computed parameter(s)
#' 
#' @export growthparameters_comparison.bgmfit
#' @export
#' 
#' @references
#' \insertAllCited{}
#' 
#' @inherit berkeley author
#' 
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' growthparameters_comparison(model, parameter = 'apgv', draw_ids = 1)
#' 
#' 
growthparameters_comparison.bgmfit <- function(model,
                                   resp = NULL,
                                   ndraws = NULL,
                                   draw_ids = NULL,
                                   newdata = NULL,
                                   re_formula = NA,
                                   parameter = NULL,
                                   xrange = 1,
                                   acg_velocity = 0.10,
                                   digits = 2,
                                   numeric_cov_at = NULL,
                                   aux_variables = NULL,
                                   levels_id = NULL,
                                   avg_reffects = NULL,
                                   deriv_model = NULL,
                                   idata_method = NULL,
                                   ipts = NULL,
                                   seed = 123,
                                   future = FALSE,
                                   future_session = 'multisession',
                                   cores = NULL,
                                   average = FALSE, 
                                   variables = NULL,
                                   comparison = "difference",
                                   type = NULL,
                                   by = FALSE,
                                   conf_level = 0.95,
                                   transform = NULL,
                                   cross = FALSE,
                                   wts = NULL,
                                   hypothesis = NULL,
                                   equivalence = NULL,
                                   eps = NULL,
                                   reformat = NULL,
                                   dummy_to_factor = NULL, 
                                   verbose = FALSE,
                                   expose_function = FALSE,
                                   usesavedfuns = NULL,
                                   clearenvfuns = NULL,
                                   envir = NULL,
                                   ...) {
  
  
  if(is.null(envir)) {
    envir <- model$model_info$envir
  } else {
    envir <- parent.frame()
  }
  
  if(is.null(usesavedfuns)) {
    if(!is.null(model$model_info$exefuns[[1]])) {
      usesavedfuns <- TRUE
    } else if(is.null(model$model_info$exefuns[[1]])) {
      if(expose_function) {
        model <- expose_model_functions(model, envir = envir)
        usesavedfuns <- TRUE
      } else if(!expose_function) {
        usesavedfuns <- FALSE
      }
    }
  } else { # if(!is.null(usesavedfuns)) {
    if(!usesavedfuns) {
      if(expose_function) {
        model <- expose_model_functions(model, envir = envir)
        usesavedfuns <- TRUE
      }
    } else if(usesavedfuns) {
      check_if_functions_exists(model, checks = TRUE, 
                                usesavedfuns = usesavedfuns)
    }
  }
  
  
  
  if(is.null(ndraws)) {
    ndraws <- brms::ndraws(model)
  }
   
  if(is.null(deriv_model)) {
    deriv_model <- TRUE
  }
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  
  
  o <- post_processing_checks(model = model,
                              xcall = match.call(),
                              resp = resp,
                              envir = envir,
                              deriv = 0)
  
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  scall <- sys.calls()
  
  get_xcall <- function(xcall, scall) {
    scall <- scall[[length(scall)]]
    if(any(grepl("growthparameters_comparison", scall, fixed = T)) |
       any(grepl("growthparameters_comparison.bgmfit", scall, fixed = T))) {
      xcall <- "growthparameters_comparison"
    } else {
      xcall <- xcall
    } 
  }
  
  if(!is.null(model$xcall)) {
    if(model$xcall == "growthparameters_comparison") {
      xcall <- "growthparameters_comparison"
    }
  } else {
    scall <- sys.calls()
    xcall <- get_xcall(xcall, scall)
  }
  
  
  xcall <- xcall
  
  check_if_package_installed(model, xcall = xcall)
  
  model$xcall <- xcall
  
  
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  arguments$model <- model
  arguments$usesavedfuns <- usesavedfuns
  
 
  if (is.null(eps)) eps <- 1e-6

  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  # set_names_  <- c('Estimate', 'Est.Error', probtitles)
  set_names_  <- c('Estimate', probtitles)
  
  
  
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- setincores <-  get.cores_[['max.cores']]
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if (future_session == 'multisession') {
      future::plan('multisession', workers = setincores)
    } else if (future_session == 'multicore') {
      future::plan('multicore', workers = setincores)
    }
  }
  
  
  
  full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                  fargs = formals(), 
                                  dargs = list(...), 
                                  verbose = verbose)
  
  full.args$model <- model
  full.args$deriv_model <- deriv_model
  newdata <- do.call(get.newdata, full.args)
  
  
  
  arguments$newdata  <- newdata
  arguments[["..."]] <- NULL
  
  
  ##############################################
  # Initiate non formalArgs()
  ##############################################
  term <- NULL;
  contrast <- NULL;
  tmp_idx <- NULL;
  predicted_lo <- NULL;
  predicted_hi <- NULL;
  predicted <- NULL;
  conf.high <- NULL;
  conf.low <- NULL;
  estimate <- NULL;
  `:=` <- NULL;
  `.` <- NULL;
  

  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  xvar_ <- paste0('xvar', resp_rev_)
  xvar  <- model$model_info[[xvar_]]
  
  cov_ <- paste0('cov', resp_rev_)
  cov  <- model$model_info[[cov_]]
  
  
  #####################
  
  # allowed_parms <- c(
  #   'ato',
  #   'dto',
  #   'vto',
  #   'apv',
  #   'dpv',
  #   'vpv',
  #   'afo',
  #   'dfo',
  #   'vfo',
  #   'acg',
  #   'dcg',
  #   'vcg',
  #   'dgs',
  #   'dgain',
  #   'vgain'
  # )
  
  
  
  allowed_parms <- c(
    'atgv',
    'tgv',
    'apgv',
    'pgv',
    'acgv',
    'cgv')
  
  if (is.null(parameter)) {
    parm <- 'apgv'
  } else {
    parm <- parameter
  }
  
  comparisons_arguments <- arguments
  exclude_args <- as.character(quote(
    c(
      parameter,
      xrange,
      acg_velocity,
      digits,
      numeric_cov_at,
      acg_asymptote,
      levels_id,
      avg_reffects,
      ipts,
      seed,
      future,
      future_session
    )
  ))[-1]
  
  for (exclude_argsi in exclude_args) {
    comparisons_arguments[[exclude_argsi]] <- NULL
  }
  
  
  if (!is.null(variables)) {
    if (!is.list(variables)) {
      stop("'variables' argument must be a named list and the first ", 
           "element should be ", xvar, 
           "\n ",
           " specified as follows ",
           "\n ",
          " variables = list(", xvar, "=", "1e-6",")",
          "\n ",
          " where 1e-6 is the default value for the argument 'eps'"
           )
    } else if (is.list(variables)) {
      set_variables <- variables
    }
  }
  
  
  if (is.null(variables)) {
    set_variables <- list(eps)
    names(set_variables) <- xvar
  }
  
  allowed_comparison <- c('difference', 'differenceavg')
  
  if(!comparison %in% allowed_comparison) {
    stop("Allowed comparison options are ", 
         paste(paste0("'", allowed_comparison, "'"), collapse = ", ")
         )
  }
  
  
  if(comparison == 'differenceavg') {
    if(!average) {
      stop("For comparison = 'differenceavg' ", 
           ", the argument 'average' should be TRUE")
    }
    if(is.null(hypothesis)) {
      stop("For comparison = 'differenceavg' ", 
           ", the argument 'hypothesis' is required.",
           " \n",
           "An example of Non-linear hypothesis testing via hypothesis",
           " argument is as follows:",
           " \n",
           "hypothesis = 'b2 - b1 = 0.2'",
           " where b2 and b1 are row indices",
           " \n",
           "(see https://marginaleffects.com/vignettes/comparisons.html)",
           " for more details",
           " \n",
           "Note that user need to set comparison = 'differenceavg'",
           " and average = TRUE when ",
           " \n", 
           "testing non-linear hypothesis as described above"
           )
    }
  }
  hypothesis = "b2 + b1 = 1"
  
  if(is.null(by)) {
    if(is.null(cov)) {
      set_group <- FALSE
    } else if(!is.null(cov)) {
      set_group <- cov
      if (!set_group %in% cov) {
        stop('by must be one of the ', cov)
      } 
    }
  } else if(!is.null(by)) {
    if (!isFALSE(by)) {
      set_group <- by
    } else if (isFALSE(by)) {
      set_group <- FALSE
    }
  }
  
 
  if(!isFALSE(set_group)) {
    if (length(parm) > 1) stop("For 'by' estimates/comparisons, please ",
                               "\n ",
                               " specificy only one parameter. ")
  }
  

  if (acg_velocity >= 1 | acg_velocity <= 0) {
    stop("The acg_velocity should be set between 0.01 and 0.99")
  }
  
  
  
  call_comparison_gparms_fun <- function(parm, eps, ...) {
    gparms_fun = function(hi, lo, x, ...) {
      y <- (hi - lo) / eps
      if (parm == 'apgv') {
        out <- sitar::getPeak(x = x, y = y)[1]
      } else if (parm == 'pgv') {
        out <- sitar::getPeak(x = x, y = y)[2]
      } else if (parm == 'atgv') {
        out <- sitar::getTakeoff(x = x, y = y)[1]
      } else if (parm == 'tgv') {
        out <- sitar::getTakeoff(x = x, y = y)[2]
      } else if (parm == 'acgv') {
        cgv  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
        vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
        out <-  x[vcgi]
      } else if (parm == 'cgv') {
        cgv  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
        vcgi <- which(abs(y - cgv) == min(abs(y - cgv)))[1]
        out <-  y[vcgi]
      } else if (parm == 'xxxx') {
        
      } else {
        stop('parm not valid')
      }
      out <- round(out, digits = digits)
      out
    } # gparms_fun
    
    comparisons_arguments$variables  <- set_variables
    comparisons_arguments$by         <- set_group
    comparisons_arguments$comparison <- gparms_fun
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    suppressWarnings({
      if(!average) {
        out <- do.call(marginaleffects::comparisons, comparisons_arguments)
      } else if(average) {
        out <- do.call(marginaleffects::avg_comparisons, comparisons_arguments)
      }
    })
    return(out)
  } # call_comparison_gparms_fun
  

  if (length(parm) == 1) {
    out_sf <- call_comparison_gparms_fun(parm = parm, eps = eps) %>% 
      data.frame() %>% 
      dplyr::mutate(!!as.symbol('parameter') := parm) %>% 
      dplyr::relocate(!!as.symbol('parameter'))
    
  } else if (length(parm) > 1) {
    list_cout <- list()
    list_name <- list()
    for (allowed_parmsi in parm) {
      list_cout[[allowed_parmsi]] <-
        call_comparison_gparms_fun(parm = allowed_parmsi, eps = eps)
      list_name[[allowed_parmsi]] <- allowed_parmsi
    }
    list_name2 <- do.call(rbind, list_name)
    out_sf <- do.call(rbind, list_cout) %>% data.frame() %>% 
      dplyr::mutate(!!as.symbol('parameter') := list_name2) %>% 
      dplyr::relocate(!!as.symbol('parameter'))
  }
  
  out_sf <- out_sf %>% 
    dplyr::rename(!!as.symbol('Parameter') := parameter) %>% 
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                         ~ round(., digits = digits))) %>% 
    data.frame()
  
  if(is.null(reformat)) {
    if(is.null(hypothesis))  reformat <- TRUE else reformat <- FALSE
    if(is.null(equivalence)) reformat <- TRUE else reformat <- FALSE
  }
  
  if (reformat) {
    out_sf <- out_sf %>% 
      dplyr::rename(!!as.symbol(set_names_[1]) := estimate) %>% 
      dplyr::rename(!!as.symbol(set_names_[2]) := conf.low) %>% 
      dplyr::rename(!!as.symbol(set_names_[3]) := conf.high) 
      data.frame()
    
    remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                      'predicted_hi', 'predicted')
    
    out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
  }
   
  out_sf <- out_sf %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of('Parameter'), toupper))
  
  return(out_sf)
}




#' @rdname growthparameters_comparison.bgmfit
#' @export
growthparameters_comparison <- function(model, ...) {
  UseMethod("growthparameters_comparison")
}


