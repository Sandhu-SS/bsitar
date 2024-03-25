

#' Compare growth curves
#' 
#'@description The \strong{marginal_comparison()} function estimates and
#'  compare growth curves such as distance and velocity. This function is a wrapper around the
#'  [marginaleffects::comparisons()] and [marginaleffects::avg_comparisons()].
#'  The [marginaleffects::comparisons()] computes unit-level (conditional)
#'  estimates whereas [marginaleffects::avg_comparisons()] return average
#'  (marginal) estimates. A detailed explanation is available
#'  [here](https://marginaleffects.com). 
#' 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param datagrid Generate a grid of user-specified values for use in the
#'   \code{newdata} argument in various functions of the \pkg{marginaleffects}
#'   package. This is useful to define where in the predictor space we want to
#'   evaluate the quantities of interest. See [marginaleffects::datagrid()] for
#'   details. The default value for the \code{datagrid} is \code{NULL} implying
#'   that no custom grid is constructed. To set a data grid, the argument should
#'   be a data.frame constructed by using the [marginaleffects::datagrid()]
#'   function, or else a named list which are internally used for setting up the
#'   grid. For the user convenience, we also allow setting an empty list
#'   \code{datagrid = list()} in which case essential arguments such as
#'   \code{model}, \code{newdata} are taken up from the respective arguments
#'   specified elsewhere. Further, the level 1 predictor (such as age) and any
#'   covariate included in the model fit (e.g., gender) are also automatically
#'   inferred from the \code{model} object.
#' 
#' @param digits An integer (default \code{2}) to set the decimal places for the
#'   estimates. The \code{digits} is passed on to the
#'   [base::round()] function.
#'   
#' @param average A logical to indicate whether to internally call the
#'    [marginaleffects::comparisons()] or the
#'    [marginaleffects::avg_comparisons()] function. If \code{FALSE} (default),
#'    [marginaleffects::comparisons()] is called otherwise
#'    [marginaleffects::avg_comparisons()] when \code{average = TRUE}.
#'
#' @param plot A logical to specify whether to plot comparisons by calling the
#'   [marginaleffects::plot_comparisons()] function (\code{FALSE}) or not
#'   (\code{FALSE}). If \code{FALSE} (default), then
#'   [marginaleffects::comparisons()] or [marginaleffects::avg_comparisons()]
#'   are called to compute predictions (see \code{average} for details).
#'
#' @param method A character string to specify whether to make computation at
#'   post draw stage by using the \code{'marginaleffects'} machinery i.e.,
#'   [marginaleffects::comparisons()] (\code{method = 'pkg'}, default) or via
#'   the custom functions written for efficiency and speed (\code{method =
#'   'custom'}). Note that \code{method = 'custom'} is on experimental basis
#'   and should be used cautiously. A particular use case is if user wants to
#'   compute estimates and comparisons together for a factor co variate i.e.,
#'   \code{by = 'cov'}, or at each value of predictor (typically 'age') by 
#'   setting \code{by = c('age', 'cov')}. The \code{method} is ignored when 
#'   \code{by = FALSE}.
#'   
#' @param constrats_at A named list or a logical (default \code{FALSE}) to
#'   specify the values at which estimates and contrast (post draw stage) via
#'   \code{hypothesis} argument should be computed. The \code{constrats_at} can
#'   be specified as a one of the following
#'   strings \code{'max', 'min', 'unique', 'range'} (e.g., \code{constrats_at =
#'   list(age = 'min')}) or else as a numeric value or a numeric vector (e.g.,
#'   \code{constrats_at = list(age = c(6, 7))}). When \code{constrats_at =
#'   NULL}, any level-1 predictor (such as \code{age}) is be automatically set
#'   up as \code{constrats_at = list(age = 'unique')}. Note that
#'   \code{constrats_at} only subsets the data that has been set up the
#'   [marginaleffects::datagrid()] or specified as the \code{newdata} argument.
#'   In case no match is found, an error will be triggered. The
#'   \code{constrats_at} argument is only evaluated when \code{method =
#'   'custom'} and the \code{hypothesis} is not \code{NULL}.
#'   
#' @param deriv A numeric to specify whether to estimate parameters based on the
#'   differentiation of the distance curve or the model based first derivative.
#'   Please see argument \code{variables} for more details.  
#'   
#' @param reformat A logical (default \code{TRUE}) to reformat the  output
#'   returned by the \code{marginaleffects} as a data.frame with column names
#'   re-defined as follows: \code{conf.low} as \code{Q2.5}, and \code{conf.high}
#'   as \code{Q97.5} (assuming that \code{conf_int = 0.95}). Also, following
#'   columns are dropped from the data frame: \code{term}, \code{contrast},
#'   \code{tmp_idx}, \code{predicted_lo}, \code{predicted_hi}, \code{predicted}.
#' 
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  growthparameters_comparison.bgmfit
#' @inheritParams  marginal_draws.bgmfit
#' @inheritParams  marginaleffects::comparisons
#' @inheritParams  marginaleffects::avg_comparisons
#' @inheritParams  marginaleffects::plot_comparisons
#' @inheritParams  marginaleffects::datagrid
#' @inheritParams  brms::fitted.brmsfit
#'
#' @return A data frame objects with estimates and CIs for computed parameter(s)
#' 
#' @export marginal_comparison.bgmfit
#' @export
#' 
#' @seealso [marginaleffects::comparisons()]
#'   [marginaleffects::avg_comparisons()]
#'   [marginaleffects::plot_comparisons()]
#' 
#' @references
#' \insertAllCited{}
#' 
#' @inherit berkeley author
#' 
#' @examples
#' 
#' \donttest{
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
#' marginal_comparison(model, parameter = 'apgv', draw_ids = 1)
#' }
#' 
marginal_comparison.bgmfit <- function(model,
                                   resp = NULL,
                                   ndraws = NULL,
                                   draw_ids = NULL,
                                   newdata = NULL,
                                   datagrid = NULL,
                                   re_formula = NA,
                                   allow_new_levels = FALSE,
                                   sample_new_levels = "gaussian",
                                   xrange = 1,
                                   digits = 2,
                                   numeric_cov_at = NULL,
                                   aux_variables = NULL,
                                   levels_id = NULL,
                                   avg_reffects = NULL,
                                   idata_method = NULL,
                                   ipts = NULL,
                                   seed = 123,
                                   future = FALSE,
                                   future_session = 'multisession',
                                   cores = NULL,
                                   average = FALSE, 
                                   plot = FALSE, 
                                   showlegends = NULL, 
                                   variables = NULL,
                                   deriv = NULL,
                                   deriv_model = NULL,
                                   method = 'pkg',
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
                                   constrats_at = FALSE,
                                   reformat = NULL,
                                   estimate_center = NULL,
                                   estimate_interval = NULL,
                                   dummy_to_factor = NULL, 
                                   verbose = FALSE,
                                   expose_function = FALSE,
                                   usesavedfuns = NULL,
                                   clearenvfuns = NULL,
                                   envir = NULL,
                                   ...) {
  
  
  if(!is.null(estimate_center)) {
    ec_ <- getOption("marginaleffects_posterior_center")
    options("marginaleffects_posterior_center" = estimate_center)
    on.exit(options("marginaleffects_posterior_center" = ec_), add = TRUE)
  }
  if(!is.null(estimate_interval)) {
    ei_ <- getOption("marginaleffects_posterior_interval")
    options("marginaleffects_posterior_interval" = estimate_interval)
    on.exit(options("marginaleffects_posterior_interval" = ei_), add = TRUE)
  }
  ec_agg <- getOption("marginaleffects_posterior_center")
  ei_agg <- getOption("marginaleffects_posterior_interval")
  if(is.null(ec_agg)) ec_agg <- "mean"
  if(is.null(ei_agg)) ei_agg <- "eti"
  
  try(zz <- insight::check_if_installed(c("marginaleffects"), 
                                        minversion = 
                                          get_package_minversion(
                                            'marginaleffects'
                                          ), 
                                        prompt = FALSE,
                                        stop = FALSE))
  
  if(!isTRUE(zz)) {
    message("Please install the latest version of the 'marginaleffects' package",
            "\n ",
            "remotes::install_github('vincentarelbundock/marginaleffects')")
    return(invisible(NULL))
  }
  
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
  
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  xvar_  <- paste0('xvar', resp_rev_)
  xvar   <- model$model_info[[xvar_]]
  cov_   <- paste0('cov', resp_rev_)
  cov    <- model$model_info[[cov_]]
  uvarby <- model$model_info$univariate_by
  
  # Note here, newdata is not model$data but rather model$model_info$bgmfit.data
  # This was must for univariate_by
  if(is.null(newdata)) {
    #newdata <- model$model_info$bgmfit.data
  }
  
  if(!is.na(uvarby)) {
    uvarby_ind <- paste0(uvarby, resp)
    varne <- paste0(uvarby, resp)
   # newdata <- newdata %>% dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
  }
  
  
  if(is.null(deriv) & is.null(deriv_model)) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(deriv == 0 & is.null(deriv_model)) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(deriv == 1 & is.null(deriv_model)) {
    deriv <- 1
    deriv_model <- TRUE
  } else if(is.null(deriv) & !deriv_model) {
    deriv <- 0
    deriv_model <- FALSE
  } else if(is.null(deriv) & deriv_model) {
    deriv <- 1
    deriv_model <- TRUE
  }
   
  # The deriv_model is a placeholder in marginaleffects
  # if(is.null(deriv_model)) {
  #   deriv_model <- TRUE
  # }
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  
  if(idata_method == 'm1') {
    stop("For marginaleffects based functions, the " ,
         " \n",
         " 'idata_method' argument must be either NULL or 'm2'" )
  }
  
  
  if (is.null(eps)) eps <- 1e-6
  
  
  # Initiate non formalArgs()
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
  
  
 
  
  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  # set_names_  <- c('Estimate', 'Est.Error', probtitles)
  set_names_  <- c('Estimate', probtitles)
  
  if(!is.null(model$model_info$decomp)) {
    if(model$model_info$decomp == "QR") deriv_model<- FALSE
  }
  
  expose_method_set <- model$model_info[['expose_method']]
  
  model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
  
  o <- post_processing_checks(model = model,
                              xcall = match.call(),
                              resp = resp,
                              envir = envir,
                              deriv = deriv, 
                              all = FALSE,
                              verbose = verbose)
  
  oall <- post_processing_checks(model = model,
                                 xcall = match.call(),
                                 resp = resp,
                                 envir = envir,
                                 deriv = deriv, 
                                 all = TRUE,
                                 verbose = FALSE)
  
  
  test <- setupfuns(model = model, resp = resp,
                    o = o, oall = oall, 
                    usesavedfuns = usesavedfuns, 
                    deriv = deriv, envir = envir, 
                    deriv_model = deriv_model, 
                    ...)
  
  if(is.null(test)) return(invisible(NULL))
  
  if(!isTRUE(
    check_pkg_version_exists('brms', 
                             minversion = get_package_minversion('brms'),
                             prompt = FALSE,
                             stop = FALSE,
                             verbose = FALSE))) {
    if(is.null(check_if_functions_exists(model, o, model$xcall,
                                         usesavedfuns = usesavedfuns))) {
      return(invisible(NULL))
    }
  }
  
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  scall <- sys.calls()
  
  get_xcall <- function(xcall, scall) {
    scall <- scall[[length(scall)]]
    if(any(grepl("marginal_comparison", scall, fixed = T)) |
       any(grepl("marginal_comparison.bgmfit", scall, fixed = T))) {
      xcall <- "marginal_comparison"
    } else {
      xcall <- xcall
    } 
  }
  
  if(!is.null(model$xcall)) {
    if(model$xcall == "marginal_comparison") {
      xcall <- "marginal_comparison"
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
                                  # fargs = formals(), 
                                  fargs = arguments, 
                                  dargs = list(...), 
                                  verbose = verbose)
  
  full.args$model <- model
  full.args$deriv_model <- deriv_model
  
  full.args$newdata <- newdata
  newdata           <- do.call(get.newdata, full.args)
  
  if(!is.na(uvarby)) {
    uvarby_ind <- paste0(uvarby, resp)
    varne <- paste0(uvarby, resp)
     newdata <- newdata %>% dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
  }
  full.args$newdata <- newdata
  
  # keeping ... cause marginaleffects:: argument is missing, with no default
  full.args[["..."]] <- NULL
  
  
  comparisons_arguments <- full.args
  
  # Drop that not required for marginaleffects::
  exclude_args <- as.character(quote(
    c(
      xrange,
      digits,
      numeric_cov_at,
      acg_asymptote,
      levels_id,
      avg_reffects,
      deriv,
      deriv_model,
      ipts,
      seed,
      future,
      future_session,
      dummy_to_factor,
      verbose,
      expose_function,
      usesavedfuns,
      clearenvfuns,
      envir,
      plot,
      showlegends,
      average,
      parameter,
      estimate_center,
      estimate_interval,
      reformat,
      method,
      constrats_at
    )
  ))[-1]
  
  if(plot) {
    exclude_args <- c(exclude_args, "cross")
  }
  
  for (exclude_argsi in exclude_args) {
    comparisons_arguments[[exclude_argsi]] <- NULL
  }
  
  
  if (!is.null(variables)) {
    if (!is.list(variables)) {
      set_variables <- variables
    } else if (is.list(variables)) {
      set_variables <- variables
      if(is.null(set_variables[[xvar]])) {
        if(deriv == 0) set_variables[[xvar]] <- eps
        if(deriv > 0 | deriv_model)  set_variables[[xvar]] <- 0
        if(method == 'custom') set_variables <- NULL
      } else if(!is.null(set_variables[[xvar]])) {
        if(eval(set_variables[[xvar]]) !=0) {
          if(verbose) {
            message("The value of ", xvar, " is not same as used in the ",
                    " \n", 
                    " model fit. Please check if this is intended")
          }
        }
      }
    }
  } else if (is.null(variables)) {
    if(deriv == 0) set_variables <- list(eps)
    if(deriv > 0 | deriv_model)  set_variables <- list(0)
    if(method == 'custom') set_variables <- NULL
    if(!is.null(set_variables)) names(set_variables) <- xvar
  } 
  
  
  
  if(is.null(by)) {
    if(is.null(cov)) {
      set_group <- FALSE
    } else if(!is.null(cov)) {
      set_group <- cov
      if (!set_group %in% cov) {
        stop("Argument 'by' must be one of the ", cov)
      } 
    }
  } else if(!is.null(by)) {
    if (!isFALSE(by)) {
      set_group <- by
    } else if (isFALSE(by)) {
      set_group <- FALSE
    }
  }
  
 
  
  
  
  
    
    comparisons_arguments$variables  <- set_variables
    comparisons_arguments$by         <- set_group

    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    
    
    if(is.null(showlegends)) {
      if(is.null(comparisons_arguments$re_formula)) {
        showlegends <- FALSE
      } else {
        showlegends <- TRUE
      }
    }
    
    
    # Set up datagrid
    
    if(!is.null(datagrid)) {
      if(is.data.frame(datagrid)) {
        set_datagrid <- datagrid
        comparisons_arguments$newdata <- set_datagrid
      } else if(is.list(datagrid)) {
        if(is.null(datagrid[['model']]))
          setmodel <- model
        else
          setmodel <- datagrid$model
        if (is.null(datagrid[['newdata']]))
          setnewdata <- newdata
        else
          setnewdata <- datagrid$newdata
        if (is.null(datagrid[['grid_type']]))
          setgrid_type <- "mean_or_mode"
        else
          setgrid_type <- datagrid$grid_type
        if (is.null(datagrid[['by']]))
          setgrid_by <- NULL
        else
          setgrid_by <- datagrid$by
        
        if (is.null(datagrid[['FUN_character']]))
          setgrid_FUN_character <- NULL
        else
          setgrid_FUN_character <- datagrid$FUN_character
        if (is.null(datagrid[['FUN_factor']]))
          setgrid_FUN_factor <- NULL
        else
          setgrid_FUN_factor <- datagrid$FUN_factor
        if (is.null(datagrid[['FUN_logical']]))
          setgrid_FUN_logical <- NULL
        else
          setgrid_FUN_logical <- datagrid$FUN_logical
        if (is.null(datagrid[['FUN_numeric']]))
          setgrid_FUN_numeric <- NULL
        else
          setgrid_FUN_numeric <- datagrid$FUN_numeric
        if (is.null(datagrid[['FUN_integer']]))
          setgrid_FUN_integer <- NULL
        else
          setgrid_FUN_integer <- datagrid$FUN_integer
        if (is.null(datagrid[['FUN_binary']]))
          setgrid_FUN_binary <- NULL
        else
          setgrid_FUN_binary <- datagrid$FUN_binary
        if (is.null(datagrid[['FUN_other']]))
          setgrid_FUN_other <- NULL
        else
          setgrid_FUN_other <- datagrid$FUN_other
        
        if(is.null(datagrid[[xvar]])) 
          setxvar <- newdata[[xvar]] 
        else 
          setxvar <- datagrid$newdata[[xvar]]
        
        datagrid_arguments <- list(model = setmodel,
                                   newdata = setnewdata,
                                   by = setgrid_by,
                                   grid_type = setgrid_type,
                                   FUN_character = setgrid_FUN_character,
                                   FUN_factor = setgrid_FUN_factor,
                                   FUN_logical = setgrid_FUN_logical,
                                   FUN_numeric = setgrid_FUN_numeric,
                                   FUN_integer = setgrid_FUN_integer,
                                   FUN_binary = setgrid_FUN_binary,
                                   FUN_other = setgrid_FUN_other
                                   )
        datagrid_arguments[[xvar]] <- setxvar
        if(setgrid_type == "mean_or_mode") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- set_group
        } else if(setgrid_type == "balanced") {
          if(!isFALSE(set_group)) datagrid_arguments[['by']] <- NULL
          # correctly set comparisons_arguments[['by']]  too 
          comparisons_arguments[['by']] <- NULL
        }
        set_datagrid <- do.call(marginaleffects::datagrid, datagrid_arguments)
        comparisons_arguments$newdata <- set_datagrid
      } else {
        stop("datagrid should be a data frame or named list")
      }
    } else if(is.null(datagrid)) {
      comparisons_arguments$newdata <- comparisons_arguments$newdata
    }
    
    # The datagrid argument is not allowed. It served its purpose by defining 
    # the newdata. So remove it from the arguments
    
    comparisons_arguments[['datagrid']] <- NULL
    
    # Somehow draw_ids not passed correctly if not specified explicitly as arg
    get_draw_ids <- comparisons_arguments[['draw_ids']]
    if(is.null(eval(get_draw_ids))) {
      set_draw_ids <- NULL
    } else if(is.numeric(eval(get_draw_ids))) {
      set_draw_ids <- get_draw_ids
    } else if(!eval(get_draw_ids)) {
      set_draw_ids <- NULL
    }
    comparisons_arguments[['draw_ids']] <- set_draw_ids
    
    # comparisons_arguments$average <- NULL
    # comparisons_arguments$parameter <- NULL
   # comparisons_arguments$... <- NULL
    
    out_sf_hy <- NULL
    allowed_methods <- c('pkg', 'custom')
    if(!method %in% allowed_methods) 
      stop("Argument 'method' should be one of the following:",
           "\n ", 
           collapse_comma(allowed_methods)
      )
    if(method == 'pkg') parm_via <- 'comparisons'
    if(method == 'custom') parm_via <- 'predictions'
    
    # if(is.null(method)) {
    #   if(is.null(predictions_arguments$re_formula)) {
    #     parm_via <- 'predictions'
    #   } else if(is.na(predictions_arguments$re_formula)) {
    #     parm_via <- 'comparisons'
    #   }
    # } else {
    #   allowed_methods <- c('comparisons', 'predictions')
    #   if(!method %in% allowed_methods) 
    #     stop("Argument 'method' should be one of the following:",
    #          "\n ", 
    #          collapse_comma(allowed_methods))
    #   parm_via <- method
    # }
    
    
    if(!is.null(comparisons_arguments[['by']])) {
      checbyx <- comparisons_arguments[['by']]
      if(all(checbyx == "")) parm_via <- 'comparisons'
      if(is.logical(checbyx)) {
        if(!checbyx) parm_via <- 'comparisons'
      }
    }
    
    
    # In this case parm_via == 'predictions' only relevant later for hypothesis
    if(parm_via == 'comparisons') {
      # suppressWarnings({
        if(!plot) {
          if(!average) {
            out <- do.call(marginaleffects::comparisons, comparisons_arguments)
          } else if(average) {
            out <- do.call(marginaleffects::avg_comparisons, comparisons_arguments)
          }
        }
        if(plot) {
          if(isFALSE(set_group)) comparisons_arguments$by <- NULL
          out <- do.call(marginaleffects::plot_comparisons, comparisons_arguments)
          outp <- out
          if(!showlegends) outp <- outp + ggplot2::theme(legend.position = 'none')
          return(outp)
        }
      # })
      out_sf <- out
    } # if(parm_via == 'comparisons') {
    
    # let probs be passed directly via ...
    # probs = c(0.25, 0.75), 
    get_etix <- utils::getFromNamespace("get_eti", "marginaleffects")
    get_hdix <- utils::getFromNamespace("get_hdi", "marginaleffects")
    get_pe_ci <- function(x, na.rm = TRUE, ...) {
      if(ec_agg == "mean") estimate <- mean(x, na.rm = na.rm)
      if(ec_agg == "median") estimate <- median(x, na.rm = na.rm)
      if(ei_agg == "eti") luci = get_etix(x, credMass = conf)
      if(ei_agg == "hdi") luci = get_hdix(x, credMass = conf)
      tibble::tibble(
        estimate = estimate, conf.low = luci[1],conf.high = luci[2]
      )
    }
  
    if(parm_via == 'predictions') {
      predictions_arguments <- comparisons_arguments
      predictions_arguments[['cross']] <- NULL
      predictions_arguments[['method']] <- NULL
      predictions_arguments[['hypothesis']] <- NULL # hypothesis evaluated later
      # predictions_arguments[['by']] <- xcby
      
      xcby <- predictions_arguments[['by']] 
      xvar <- intersect(xvar, xcby)         
      # cby <- setdiff(xcby, xvar)            
      # parm <- 'xz'
      
      if(!average) {
        out <- do.call(marginaleffects::predictions, predictions_arguments)
      } else if(average) {
        out <- do.call(marginaleffects::predictions, predictions_arguments)
      }
      
      onex0 <- out %>% marginaleffects::posterior_draws()
      
      
      if(isFALSE(constrats_at)) {
        constrats_at <- NULL
      } else if(!isFALSE(constrats_at)) {
        if(is.null(constrats_at)) {
          if(length(xvar) > 0) {
            constrats_at <- list()
            constrats_at[[xvar]] <- 'unique'
          }
        }
      }
      
      # if(is.null(constrats_at)) {
      #   if(length(xvar) > 0) {
      #     constrats_at <- list()
      #     constrats_at[[xvar]] <- 'unique'
      #   }
      # }
      
      # For hypothesis
      groupvarshyp1 <- c('drawid')
      groupvarshyp2 <- c('term')
      if(!is.null(constrats_at)) {
        for (caxi in names(constrats_at)) {
          if(!caxi %in% names(onex0)) {
            stop("Variable '", caxi, ". specified in 'constrats_at' is not in",
                 "\n ", 
                 " the estimates. Note that ", caxi, " should also be included",
                 "\n ", 
                 " in the 'by' argument. The current 'by' argument includes:", 
                 "\n ",
                 collapse_comma(xcby)
            )
          }
          allowed_char_constrats_at <- c('max', 'min', 'unique', 'range')
          if(is.character(constrats_at[[caxi]])) {
            if(length(constrats_at[[caxi]]) > 1) {
              stop(caxi, " specified in 'constrats_at' as character should be a single character:",
                   "\n ", 
                   collapse_comma(allowed_char_constrats_at)
              )
            }
            if(!constrats_at[[caxi]] %in% allowed_char_constrats_at) {
              stop(constrats_at[[caxi]], " specified in 'constrats_at' as character should be one of",
                   "\n ", 
                   collapse_comma(allowed_char_constrats_at)
              )
            }
            getatval<- paste0(constrats_at[[caxi]], "(", "onex0", "[['", caxi, "']]", ")" )
            getatval <- eval(parse(text = getatval))                
          } # if(is.character(constrats_at[[caxi]])) {
          
          if(!is.character(constrats_at[[caxi]])) {
            getatval <- constrats_at[[caxi]]
          }
          constrats_at[[caxi]] <- getatval
          groupvarshyp1 <- c(caxi, groupvarshyp1)
          groupvarshyp2 <- c(caxi, groupvarshyp2)
        } # for (caxi in names(constrats_at)) {
        
        
        for (caxi in names(constrats_at)) {
          onex1 <- base::subset(onex0, onex0[[caxi]] %in% constrats_at[[caxi]])
          
          # onex1 %>% 
          #   onex0 %>% dplyr::filter(eval(parse(text = caxi)) %in%  constrats_at[[caxi]])
            # onex0 %>% dplyr::filter(!! as.name(caxi) %in%  constrats_at[[caxi]])
            
          
          if(nrow(onex1) == 0) {
            stop(caxi, " specified in 'constrats_at' has resulted in zero rows:",
                 "\n ", 
                 ""
            )
          }
        }
      } # if(!is.null(constrats_at)) {
      
      
      
      if(is.null(constrats_at)) {
        onex1 <- onex0
      }
      
      # zcx <- out %>% marginaleffects::posterior_draws() %>%
      #   dplyr::mutate(!! as.name(parm) :=  eval(parse(text = 'draw'))) %>%
      #   dplyr::select(-estimate)
      # drawid_c <-  split(zcx , f = zcx[['drawid']] )
      # drawid_ci_c <- list()
      # for (drawid_ci in 1:length(drawid_c)) {
      #   outhy <- drawid_c[[drawid_ci]] %>%
      #     dplyr::rename(estimate = dplyr::all_of('draw'))
      #   if(length(xvar) > 0) {
      #     outhy <- outhy %>%
      #       dplyr::mutate(!! 'xf' :=  as.factor(eval(parse(text = xvar))))
      #   } else {
      #     outhy <- outhy %>%
      #       dplyr::mutate(!! 'xf' :=  as.factor(1))
      #   }
      #   parmi_estimate <- 'estimate'
      #   x_ci_c <- list()
      #   for (x in levels(outhy[['xf']])) {
      #     outhy2 <- outhy %>% dplyr::filter(!! as.name('xf') == x)
      #     x_ci_c[[x]] <- outhy2 %>% dplyr::reframe(
      #       dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, .unpack = TRUE),
      #       .by = dplyr::all_of(!! cby)
      #     ) %>%
      #       dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, fixed = TRUE)) %>%
      #       dplyr::mutate(!! 'at' := x) %>%
      #       dplyr::relocate(dplyr::all_of('at'))
      #   }
      #   drawid_ci_c[[drawid_ci]] <- do.call(rbind, x_ci_c)
      # }
      # out_sf <- do.call(rbind, drawid_ci_c) %>% data.frame
      # out_sf <- out_sf %>% dplyr::select(-'at')
      # row.names(out_sf) <- NULL
      
      
      # zcx <- out %>% marginaleffects::posterior_draws()
      
      parmi_estimate <- 'estimate'
      measurevar     <- 'draw'
      groupvars2 <- predictions_arguments[['by']] 
      groupvars1 <- c('drawid', groupvars2)
    
      aggform1 <- as.formula(paste(measurevar, paste(groupvars1, collapse=" + "), 
                                   sep=" ~ "))
      aggform2 <- as.formula(paste(measurevar, paste(groupvars2, collapse=" + "), 
                                   sep=" ~ "))
      
      if(ec_agg == "mean") 
        onex1 <- stats::aggregate(aggform1, data = onex1, mean) 
      
      if(ec_agg == "median") 
        onex1 <- stats::aggregate(aggform1, data = onex1, median) 
      
      
      out_sf <- onex1 %>% 
        dplyr::mutate(!! parmi_estimate := eval(parse(text = measurevar))) %>% 
        dplyr::reframe(
          dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, .unpack = TRUE),
          .by = dplyr::all_of(!! groupvars2)
        ) %>%
        dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, fixed = TRUE))
      
     
      if(!is.null(hypothesis)) {
        get_hypothesis_x_modify <- function(.x, hypothesis, by, draws, ...) {
          get_hypothesis_x(x = .x, hypothesis = hypothesis, by = by, draws = draws)
        }
        
        # Moved from here to up to even get estimates at constrats_at
        
        # groupvarshyp1 <- c('drawid')
        # groupvarshyp2 <- c('term')
        # if(!is.null(constrats_at)) {
        #   for (caxi in names(constrats_at)) {
        #     if(!caxi %in% names(onex1)) {
        #       stop(caxi, " specified in 'constrats_at' is not in the estimates",
        #            "\n ",
        #            " Note that ", caxi, " should also be included in the 'by' argument",
        #            "\n ",
        #            " Avilable names are:",
        #            "\n ",
        #            collapse_comma(xcby)
        #       )
        #     }
        #     allowed_char_constrats_at <- c('max', 'min', 'unique', 'range')
        #     if(is.character(constrats_at[[caxi]])) {
        #       if(length(constrats_at[[caxi]]) > 1) {
        #         stop(caxi, " specified in 'constrats_at' as character should be a single character:",
        #              "\n ",
        #              collapse_comma(allowed_char_constrats_at)
        #         )
        #       }
        #       if(!constrats_at[[caxi]] %in% allowed_char_constrats_at) {
        #         stop(constrats_at[[caxi]], " specified in 'constrats_at' as character should be one of",
        #              "\n ",
        #              collapse_comma(allowed_char_constrats_at)
        #         )
        #       }
        #       getatval<- paste0(constrats_at[[caxi]], "(", "onex1", "[['", caxi, "']]", ")" )
        #       getatval <- eval(parse(text = getatval))
        #     } # if(is.character(constrats_at[[caxi]])) {
        # 
        #     if(!is.character(constrats_at[[caxi]])) {
        #       getatval <- constrats_at[[caxi]]
        #     }
        #     constrats_at[[caxi]] <- getatval
        #     groupvarshyp1 <- c(caxi, groupvarshyp1)
        #     groupvarshyp2 <- c(caxi, groupvarshyp2)
        #   } # for (caxi in names(constrats_at)) {
        # 
        #   for (caxi in names(constrats_at)) {
        #     onex12 <-
        #       onex1 %>% filter(!! as.name(caxi) %in%  constrats_at[[caxi]])
        #     if(nrow(onex12) == 0) {
        #       stop(caxi, " specified in 'constrats_at' has resulted in zero rows:",
        #            "\n ",
        #            ""
        #       )
        #     }
        #   }
        #   onex1 <- onex12
        # } # if(!is.null(constrats_at)) {
        
        if(nrow(onex1) > 25) {
          if(is.null(constrats_at)) {
            message("Note that the 'marginaleffects' package does not allow" ,
                    "\n",
                    "hypothesis argument when estimates rows are more than 25",
                    "\n",
                    "To avoid this issue, you can use 'constrats_at' argument",
                    "\n"
            )
          }
        }
        
        
        out_sf_hy <-
          onex1 %>% dplyr::group_by_at(groupvarshyp1) %>% 
          dplyr::mutate(!! parmi_estimate := eval(parse(text = measurevar))) %>% 
          dplyr::group_modify(., ~get_hypothesis_x_modify(.x,
                                                          hypothesis = 'pairwise', 
                                                          by = 'class', 
                                                          draws = estimate
          ), .keep = F) %>% 
          dplyr::ungroup() %>% 
          dplyr::reframe(
            dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, .unpack = TRUE),
            .by = dplyr::all_of(!! groupvarshyp2)
          ) %>%
          dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, fixed = TRUE))
        
      }
      
      # if(!is.null(hypothesis)) {
      #   parmi_ci_c <- list()
      #   for (parmi in parm) {
      #     hypthesis_drawid_ci_c <- list()
      #     for (drawid_ci in 1:length(drawid_c)) {
      #       outhy <- drawid_c[[drawid_ci]] %>% 
      #         dplyr::rename(estimate = dplyr::all_of(parmi)) 
      #       if(length(xvar) > 0) {
      #         outhy <- outhy %>% 
      #           dplyr::mutate(!! 'xf' :=  as.factor(eval(parse(text = xvar))))
      #       } else {
      #         outhy <- outhy %>% 
      #           dplyr::mutate(!! 'xf' :=  as.factor(1))
      #       }
      #       x_ci_c <- list()
      #       for (x in levels(outhy[['xf']])) {
      #         outhy2 <- outhy %>% dplyr::filter(!! as.name('xf') == x) 
      #         x_ci_c[[x]] <- get_hypothesis_x(x = outhy2, 
      #                                         hypothesis = hypothesis, 
      #                                         by = cby, 
      #                                         draws = estimate) %>% 
      #           dplyr::mutate(!! 'at' := as.factor(x)) %>% 
      #           dplyr::mutate(!! 'drawid' := drawid_ci) %>% 
      #           dplyr::relocate(dplyr::all_of(c('drawid', 'at'))) 
      #       }
      #       hypthesis_drawid_ci_c[[drawid_ci]] <- do.call(rbind, x_ci_c)
      #     }
      #     parmi_ci_c[[parmi]] <- do.call(rbind, hypthesis_drawid_ci_c)
      #   }
      #   out4 <- hypthesis_drawid_ci_c %>% do.call(rbind, .) %>% data.frame()
      #   
      #   parm <- levels(out4[['at']])
      #   summary_c <- list()
      #   for (parmi in parm) {
      #     cby_term <- 'term'
      #     parmi_estimate <- 'estimate'
      #     summary_c[[parmi]] <- out4 %>% dplyr::filter(!! as.name('at') == parmi) %>% 
      #       dplyr::reframe(
      #         dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, .unpack = TRUE),
      #         .by = dplyr::all_of(!! cby_term) 
      #       ) %>% 
      #       dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, fixed = TRUE)) %>% 
      #       dplyr::mutate(!! 'at' := parmi) %>% 
      #       dplyr::relocate(dplyr::all_of('at'))
      #   }
      #   out5 <- summary_c %>% do.call(rbind, .) %>% data.frame()
      #   row.names(out5) <- NULL
      #   if(length(xvar) > 0) {
      #     out5 <- out5 %>% dplyr::mutate(!!xvar := as.numeric(eval(parse(text = 'at')))) %>% 
      #       dplyr::relocate(dplyr::all_of(xvar), .before = 'at')
      #   } else {
      #     out5 <- out5 %>%
      #       dplyr::select(-dplyr::all_of('at'))
      #   }
      #   out_sf_hy <- out5
      # } # if(!is.null(hypothesis)) {
    } # if(parm_via == 'predictions') {
  

  
  out_sf <- out_sf %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                         ~ round(., digits = digits))) %>%
    data.frame()
  
 
  if(is.null(reformat)) {
    if(is.null(hypothesis) && is.null(equivalence)) {
      reformat <- TRUE
    } else {
      reformat <- FALSE
    }
  }
  
  if (reformat) {
    out_sf <- out_sf %>% 
      dplyr::rename(!!as.symbol(set_names_[1]) := dplyr::all_of('estimate')) %>% 
      dplyr::rename(!!as.symbol(set_names_[2]) := dplyr::all_of('conf.low')) %>% 
      dplyr::rename(!!as.symbol(set_names_[3]) := dplyr::all_of('conf.high')) %>% 
      data.frame()
    
    remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                      'predicted_hi', 'predicted', 'rowid')
    
    out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
    row.names(out_sf) <- NULL
    
    
    
    if(!is.null(out_sf_hy)) {
      out_sf_hy <- out_sf_hy %>% 
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                    ~ round(., digits = digits))) %>% 
        # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
        dplyr::rename(!!as.symbol(set_names_[1]) := dplyr::all_of('estimate')) %>% 
        dplyr::rename(!!as.symbol(set_names_[2]) := dplyr::all_of('conf.low')) %>% 
        dplyr::rename(!!as.symbol(set_names_[3]) := dplyr::all_of('conf.high')) %>% 
        # dplyr::rename_with(., ~ tools::toTitleCase(.x)) %>% 
        dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
        data.frame()
    } # if(!is.null(out_sf_hy)) {
  } # if (reformat) {
  
  
  if(!is.null(out_sf_hy)) {
    out_sf <- out_sf %>% dplyr::ungroup()
    out_sf_hy <- out_sf_hy %>% dplyr::ungroup()
    out_sf <- list(estimate = out_sf, contrast = out_sf_hy)
  } 
  
  return(out_sf)
}




#' @rdname marginal_comparison.bgmfit
#' @export
marginal_comparison <- function(model, ...) {
  UseMethod("marginal_comparison")
}


