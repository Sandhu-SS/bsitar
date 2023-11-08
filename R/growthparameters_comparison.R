




#' Compute growth parameters for the \code{bgmfit} model
#' 
#'@description The \code{growthparameters_comparison} function provides
#'  estimates of growth parameters (such as peak velocity and the age at peak
#'  velocity) derived from the growth velocity curve. This function is a wrapper
#'  around the [marginaleffects::comparisons()] and
#'  [marginaleffects::avg_comparisons()]. These function estimate and return
#'  predicted values at t different regressor values and and compare those
#'  predictions by computing a difference between them. The
#'  [marginaleffects::comparisons()] computes unit-level (conditional) estimates
#'  whereas [marginaleffects::avg_comparisons()] return average (marginal)
#'  estimates. A detailed explanation is available at
#'  <https://marginaleffects.com/articles/comparisons.html> and
#'  <https://marginaleffects.com/>. Note that for the current use case, i.e.,
#'  for estimating growth parameters, the arguments, especially \code{variables}
#'  and \code{comparion} are modified via custom functions (see below).
#'  
#' @details The \code{growthparameters_comparison} function estimates and
#'   returns the following growth paramaeters:
#' \itemize{
#'   \item pv  - peak velocity
#'   \item apv - age at peak velocity
#'   \item tv  - take off velocity
#'   \item atv - age at take off velocity
#'   \item cv  - cessation velocity
#'   \item acv - age at cessation velocity
#' }
#' 
#' The take off velocity is the minimum velocity before the peak velocity and it
#' indicates the beginning of the pubertal growth spurt. The cessation velocity
#' indicates the end of the active pubertal growth spurt and is calculated as
#' the percentage of the peak velocity (\code{pv}). The percentage of the peak
#' velocity used in the calculation of the the cessation velocity (\code{'cv'})
#' is controlled via the \code{acg_velocity} argument. The \code{acg_velocity}
#' takes a real value as an input (greater than 0 and less than 1). Typically, a
#' 10 percent of \code{pv} (i.e., \code{acg_velocity = 0.1}) is considered as a
#' good indicator of the cessation of the active pubertal growth spurt
#' \insertCite{Anna2022}{bsitar}.
#' 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param parameter A single character string, or a vector of character strings 
#' indicating the parameter to be estimated. Options available are \code{ato}, 
#' \code{vto}, \code{vpv}, \code{apv}, \code{vcg}, and \code{acg}. For 
#' \code{parameter = NULL} (default), age at peak growth velocity is estimated.
#' 
#' @param acg_velocity A real number to set the percentage of peak growth growth
#'   velocity as the cessation velocity when estimating the \code{vcg} and
#'   \code{acg} growth parameters. The \code{acg_velocity} should be greater
#'   than \code{0} and less than \code{1}. The default \code{acg_velocity =
#'   0.10} indicates that a 10 per cent of the peak growth velocity will be used
#'   to get the cessation velocity and the corresponding age at the cessation
#'   velocity. For example if peak growth velocity estimate is \code{10
#'   mm/year}, then cessation velocity is \code{1 mm/year}.
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
#' the \code{variables} is the level 1 predictor such as \code{age}/\code{time}
#' and is is a named list with value of list set via the \code{esp} 
#' (default 1e-6). If \code{NULL}, the \code{variables} is set internally by
#' retrieving the relevant information from the \code{model}. Otherwise, 
#' can define it as \code{variables = list('age' = 1e-6)}.
#' 
#' @param comparison For estimating growth parameters in the current use case,
#' the \code{comparison} is an internal function that uses [sitar::getPeak()], 
#' [sitar::getTakeoff()] and [sitar::getTrough()] functions to estimate
#' various growth parameters. 
#' 
#' @param reformat A logical (default \code{TRUE}) to reformat returned output
#' as a data.frame with colnames renamed as follows: \code{estimate} as 
#' \code{Estimate}, \code{conf.low} as \code{Q2.5}, and \code{conf.high} as 
#' \code{Q97.5} (assuming that \code{conf_int = 0.95}). Also, following 
#' columns are dropped from the data frame: \code{term}, \code{contrast},
#' \code{tmp_idx}, \code{predicted_lo}, \code{predicted_hi}, \code{predicted}.
#' 
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
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @examples
#' 
#' # To avoid running the model which takes some time, model fit to the
#' # \code{berkeley_mdata} has already been saved as berkeley_mfit.rda object.
#' # Please see \code{bgm} examples.
#' 
#' model <- berkeley_mfit
#' 
#' \donttest{
#' growthparameters_comparison(model, parameter = 'apv')
#' }
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
                                   reformat = TRUE,
                                   envir = globalenv(),
                                   ...) {
  
  
  # if(system.file(package='collapse') == "") {
  #   stop("Please install 'collapse' package before 
  #        calling the 'growthparameters_comparison'")
  # }
  
  required_packages <- c('tidyr', 'collapse')
  check_and_install_if_not_installed(required_packages, 
                                     'growthparameters_comparison')
  
  if (is.null(ndraws))
    ndraws  <- brms::ndraws(model)
  else
    ndraws <- ndraws
  
  o <- post_processing_checks(model = model,
                              xcall = match.call(),
                              resp = resp,
                              envir = envir,
                              deriv = 0)
  
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  
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
  
  
  model$xcall <- xcall
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  
  # This arguments$model <- model required when using pipe %>% to use gparameter
  arguments$model <- model
  
 
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
  
  
  if (is.null(newdata)) {
    newdata <- model$data
  }
  
  # if (is.null(newdata)) {
  #   newdata <- get.newdata(
  #     model,
  #     newdata = newdata,
  #     resp = resp,
  #     numeric_cov_at = NULL,
  #     aux_variables = NULL,
  #     levels_id = NULL,
  #     ipts = ipts,
  #     xrange = xrange
  #   )
  # }
  
  
  
  arguments$newdata <- newdata
  
  
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
  
  allowed_parms <- c(
    'ato',
    'dto',
    'vto',
    'apv',
    'dpv',
    'vpv',
    'afo',
    'dfo',
    'vfo',
    'acg',
    'dcg',
    'vcg',
    'dgs',
    'dgain',
    'vgain'
  )
  
  
  if (is.null(parameter)) {
    parm <- 'apv'
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
      stop('variables must be a named list')
    } else if (is.list(variables)) {
      set_variables <- variables
    }
  }
  
  
  if (is.null(variables)) {
    set_variables <- list(eps)
    names(set_variables) <- xvar
  }
  
  
  
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
  
  
  
  call_comparison_gparms_fun <- function(parm, eps...) {
    gparms_fun = function(hi, lo, x, ...) {
      y <- (hi - lo) / eps
      # xy <- xy.coords(x, y)
      # xy <- unique(as.data.frame(xy[1:2])[order(xy$x), ])
      # x  <- xy$x
      # y  <- xy$y
      # out <- x[y == max(y)][1]
      if (parm == 'apv') {
        out <- sitar::getPeak(x = x, y = y)[1]
      } else if (parm == 'vpv') {
        out <- sitar::getPeak(x = x, y = y)[2]
      } else if (parm == 'ato') {
        out <- sitar::getTakeoff(x = x, y = y)[1]
      } else if (parm == 'vto') {
        out <- sitar::getTakeoff(x = x, y = y)[2]
      } else if (parm == 'afo') {
        vcg  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
        vcgi <- which(abs(y - vcg) == min(abs(y - vcg)))[1]
        out <-  x[vcgi]
      } else if (parm == 'vfo') {
        vcg  <- acg_velocity * sitar::getPeak(x = x, y = y)[2]
        vcgi <- which(abs(y - vcg) == min(abs(y - vcg)))[1]
        out <-  y[vcgi]
      } else if (parm == 'xxxx') {
        
      } else {
        stop('parm not valid')
      }
      out <- round(out, digits = digits)
      
      out
    }
    
    comparisons_arguments$variables  <- set_variables
    comparisons_arguments$by         <- set_group
    comparisons_arguments$comparison <- gparms_fun
    

    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    
    suppressWarnings({
      if(!average) {
        out <- do.call(marginaleffects::comparisons, comparisons_arguments)
       # out <- do.call(comparisons, comparisons_arguments)
      } else if(average) {
        out <- do.call(marginaleffects::avg_comparisons, comparisons_arguments)
      }
    })
    
    out
  }
  

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
  
  
  return(out_sf)
}




#' @rdname growthparameters_comparison.bgmfit
#' @export
growthparameters_comparison <- function(model, ...) {
  UseMethod("growthparameters_comparison")
}


