

#' Estimate growth curves
#' 
#' @description The \strong{marginal_draws()} function estimates and plots
#'   growth curves (distance and velocity) by using \pkg{marginaleffects}
#'   package as back-end. This function can compute growth curves (via
#'   [marginaleffects::predictions()]), average growth curves (via
#'   [marginaleffects::avg_predictions()]) or plot growth curves (via
#'   [marginaleffects::plot_predictions()]). Please see
#'   [here](https://marginaleffects.com/) for details.
#'
#' @details The \strong{marginal_draws()} estimates fitted values (via
#'   [brms::fitted.brmsfit()]) or the posterior draws from the posterior
#'   distribution (via [brms::predict.brmsfit()]) depending on the \code{type}
#'   argument.
#'   
#' @param average A logical to indicate whether to internally call the
#'    [marginaleffects::predictions()] or the
#'    [marginaleffects::avg_predictions()] function. If \code{FALSE} (default),
#'    [marginaleffects::predictions()] is called otherwise
#'    [marginaleffects::avg_predictions()] when \code{average = TRUE}.
#'
#' @param plot A logical to specify whether to plot predictions by calling the
#'   [marginaleffects::plot_predictions()] function (\code{FALSE}) or not
#'   (\code{FALSE}). If \code{FALSE} (default), then
#'   [marginaleffects::predictions()] or [marginaleffects::avg_predictions()]
#'   are called to compute predictions (see \code{average} for details).
#'   
#' @param method A character string to specify whether to make computation at
#'   post draw stage by using the \code{'marginaleffects'} machinery i.e.,
#'   [marginaleffects::predictions()] (\code{method = 'pkg'}, default) or via
#'   the custom functions written for efficiency and speed (\code{method =
#'   'custom'}). Note that \code{method = 'custom'} is on experimental basis
#'   and should be used cautiously. A particular use case is if user wants to
#'   compute estimates and comparisons together for a factor co variate i.e.,
#'   \code{by = 'cov'}, or at each value of predictor (typically 'age') by 
#'   setting \code{by = c('age', 'cov')}. The \code{method} is ignored when 
#'   \code{by = FALSE}.
#' 
#' @param deriv An integer to indicate whether to estimate distance curve or its
#'   derivative (i.e., velocity curve). The \code{deriv = 0} (default) is for
#'   the distance curve whereas \code{deriv = 1} for the velocity curve. 
#' 
#' @inheritParams growthparameters.bgmfit
#' @inheritParams growthparameters_comparison.bgmfit
#' @inheritParams marginal_comparison.bgmfit
#' @inheritParams marginaleffects::predictions
#' @inheritParams marginaleffects::plot_predictions
#' @inheritParams brms::fitted.brmsfit
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit()] 
#' function. Please see \code{brms::fitted.brmsfit()} for details on 
#' various options available.
#' 
#' @return An array of predicted mean response values. See
#'   [brms::fitted.brmsfit] for details.
#' 
#' @export marginal_draws.bgmfit
#' @export
#' 
#' @seealso [marginaleffects::predictions()]
#'   [marginaleffects::avg_predictions()]
#'   [marginaleffects::plot_predictions()]
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
#' # Population average distance curve
#' marginal_draws(model, deriv = 0, re_formula = NA)
#' 
#' 
#' # Individual-specific distance curves
#' marginal_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' marginal_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' marginal_draws(model, deriv = 1, re_formula = NULL)
#' }
#' 
marginal_draws.bgmfit <-
  function(model,
           resp = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           newdata = NULL,
           datagrid = NULL,
           re_formula = NA,
           allow_new_levels = FALSE,
           sample_new_levels = "gaussian", 
           parameter = NULL,
           xrange = 1,
           acg_velocity = 0.10,
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
           fullframe = FALSE, 
           average = FALSE, 
           plot = FALSE, 
           showlegends = NULL, 
           variables = NULL,
           condition = NULL,
           deriv = 0,
           deriv_model = TRUE,
           method = 'pkg',
           type = NULL,
           by = NULL,
           conf_level = 0.95,
           transform = NULL,
           byfun = NULL,
           wts = NULL,
           hypothesis = NULL,
           equivalence = NULL,
           constrats_by = NULL,
           constrats_at = NULL,
           constrats_subset = NULL,
           reformat = NULL,
           estimate_center = NULL,
           estimate_interval = NULL,
           usedtplyr = FALSE,
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
      message("Please install the latest version of the 'marginaleffects' 
              package",
              "\n ",
              "remotes::install_github('vincentarelbundock/marginaleffects')")
      return(invisible(NULL))
    }
    
    
    if(usedtplyr) {
      try(zz <- insight::check_if_installed(c("dtplyr"), 
                                            minversion =
                                              get_package_minversion(
                                                'dtplyr'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
      
      try(zz <- insight::check_if_installed(c("data.table"), 
                                            minversion =
                                              get_package_minversion(
                                                'data.table'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
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
   
    
    # If default marginal effects 'dydx', then 
    call_predictions <- TRUE
    call_slopes      <- FALSE
    if(!deriv_model) {
      if(deriv > 0) {
        deriv <- 0
        call_predictions <- FALSE
        call_slopes      <- TRUE
      }
    } # if(!deriv_model) {
    
    
    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
    
    if(idata_method == 'm1') {
      stop("For marginaleffects based functions, the " ,
           " \n",
           " 'idata_method' argument must be either NULL or 'm2'" )
    }
    
    
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
    
    model$model_info[['expose_method']] <- 'NA' # Over ride 'R'
    
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
      if(any(grepl("marginal_draws", scall, fixed = T)) |
         any(grepl("marginal_draws.bgmfit", scall, fixed = T))) {
        xcall <- "marginal_draws"
      } else {
        xcall <- xcall
      } 
    }
    
    if(!is.null(model$xcall)) {
      if(model$xcall == "marginal_draws") {
        xcall <- "marginal_draws"
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
      if(usedtplyr) {
        newdata <- newdata %>% dtplyr::lazy_dt() %>% 
          dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
      } else if(!usedtplyr) {
        newdata <- newdata %>% 
          dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
      }
    }
    
    full.args$newdata <- newdata
    full.args[["..."]] <- NULL
    
    predictions_arguments <- full.args
    

    # Drop that not required for marginaleffects::
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
        deriv,
        deriv_model,
        aux_variables, 
        idata_method, 
        condition, 
        reformat, 
        dummy_to_factor, 
        expose_function,
        ipts,
        seed,
        future,
        future_session,
        verbose,
        usesavedfuns,
        clearenvfuns,
        envir, 
        fullframe,
        average,
        plot,
        showlegends,
        average,
        estimate_center,
        estimate_interval, 
        reformat,
        method,
        constrats_by,
        constrats_at,
        constrats_subset,
        usedtplyr
      )
    ))[-1]
    
    if(call_predictions) exclude_args <- c(exclude_args, 'variables')
    
    if(call_slopes) exclude_args <- c(exclude_args, 'transform', 'byfun')
    
    for (exclude_argsi in exclude_args) {
      predictions_arguments[[exclude_argsi]] <- NULL
    }
    
    
    
    
    
    if(call_slopes) {
      if (!is.null(variables)) {
        if (!is.character(variables)) {
          stop("'variables' argument must be a character string such as", 
               "\n ",
               " variables = ", "'", xvar, "'"
          )
        } else {
          set_variables <- variables
          if(!grepl(xvar, variables)) {
            set_variables <- xvar
          } else if(!is.null(set_variables[[xvar]])) {
            
          }
        }
      } else if (is.null(variables)) {
        set_variables <- xvar
      } 
    } # if(call_slopes) {
    
    
    
    # Decide if set by = NULL and then here pick and replace 'by' set_group 
    
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
    
    
    if(call_slopes) predictions_arguments$variables  <- set_variables
    predictions_arguments$by         <- set_group
    
    if(is.null(predictions_arguments$by)) predictions_arguments$by < 'NULL'
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)

    if(is.null(showlegends)) {
      if(is.null(predictions_arguments$re_formula)) {
        showlegends <- FALSE
      } else {
        showlegends <- TRUE
      }
    }
    
    
    # Set up datagrid
    
    if(!is.null(datagrid)) {
      if(is.data.frame(datagrid)) {
        set_datagrid <- datagrid
        predictions_arguments$newdata <- set_datagrid
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
          # correctly set predictions_arguments[['by']]  too 
          predictions_arguments[['by']] <- 
            setdiff(predictions_arguments[['by']], cov)
        }
        set_datagrid <- do.call(marginaleffects::datagrid, datagrid_arguments)
        predictions_arguments$newdata <- set_datagrid
      } else {
        stop("datagrid should be a data frame or named list")
      }
    } else if(is.null(datagrid)) {
      predictions_arguments$newdata <- predictions_arguments$newdata
    }

    predictions_arguments[['datagrid']] <- NULL
    
    # Somehow draw_ids not passed correctly if not specified explicitly as arg
    get_draw_ids <- predictions_arguments[['draw_ids']]
    if(is.null(eval(get_draw_ids))) {
      set_draw_ids <- NULL
    } else if(is.numeric(eval(get_draw_ids))) {
      set_draw_ids <- get_draw_ids
    } else if(!eval(get_draw_ids)) {
      set_draw_ids <- NULL
    }
   predictions_arguments[['draw_ids']] <- set_draw_ids
   
   
   out_sf_hy <- NULL
   allowed_methods <- c('pkg', 'custom')
   if(!method %in% allowed_methods) 
     stop("Argument 'method' should be one of the following:",
          "\n ", 
          collapse_comma(allowed_methods)
     )
   if(!is.null(predictions_arguments[['by']])) {
     checbyx <- predictions_arguments[['by']]
     if(all(checbyx == "")) method <- 'pkg'
     if(is.logical(checbyx)) {
       if(!checbyx) method <- 'pkg'
     }
   }
   
   if(method == 'pkg') {
     if(call_predictions) {
       if(!plot) {
         if(!average) {
           out_sf <- do.call(marginaleffects::predictions, 
                             predictions_arguments)
         } else if(average) {
           out_sf <- do.call(marginaleffects::avg_predictions, 
                             predictions_arguments)
         }
       } else if(plot) {
         . <- do.call(marginaleffects::plot_predictions, predictions_arguments)
         outp <- .
         if(!showlegends) outp <- outp + 
           ggplot2::theme(legend.position = 'none')
         return(outp)
       }
     } # if(call_predictions) {
     
     if(call_slopes) {
       if(!plot) {
         if(!average) {
           out_sf <- do.call(marginaleffects::slopes, predictions_arguments)
         } else if(average) {
           out_sf <- do.call(marginaleffects::avg_slopes, predictions_arguments)
         }
       } else if(plot) {
         . <- do.call(marginaleffects::plot_slopes, predictions_arguments)
         outp <- .
         outp <- outp + ggplot2::theme(legend.position = 'none')
         return(outp)
       }
     } # if(call_slopes) {
   } # if(method == 'pkg') {

    
  
   # This from marginaleffects does not allow na.rm
   # So use stats::quantile instead
   get_etix <- utils::getFromNamespace("get_eti", "marginaleffects")
   get_etix <- stats::quantile
   get_hdix <- utils::getFromNamespace("get_hdi", "marginaleffects")
   get_pe_ci <- function(x, draw = NULL, na.rm = TRUE, ...) {
     if(data.table::is.data.table(x) | is.data.frame(x)) {
       if(is.null(draw)) {
         stop("please specify the 'draw' argument")
       }
       x <- x %>% dplyr::select(dplyr::all_of(draw)) %>% 
         unlist() %>% as.numeric()
     }
     if(ec_agg == "mean") estimate <- mean(x, na.rm = na.rm)
     if(ec_agg == "median") estimate <- median(x, na.rm = na.rm)
     if(ei_agg == "eti") luci = get_etix(x, probs = probs, na.rm = na.rm)
     if(ei_agg == "hdi") luci = get_hdix(x, credMass = conf)
     tibble::tibble(
       estimate = estimate, conf.low = luci[1],conf.high = luci[2]
     )
   }
    
   if(method == 'custom') {
     predictions_arguments[['cross']] <- NULL
     predictions_arguments[['method']] <- NULL
     predictions_arguments[['hypothesis']] <- NULL # hypothesis evaluated later

     by <- predictions_arguments[['by']] 
     # xcby <- predictions_arguments[['by']] 
     # xvar <- intersect(xvar, xcby)           
    
     if(call_predictions) {
       if(!plot) {
         if(!average) {
           out <- do.call(marginaleffects::predictions, predictions_arguments)
         } else if(average) {
           out <- do.call(marginaleffects::avg_predictions, predictions_arguments)
         }
       } else if(plot) {
         out <- do.call(marginaleffects::plot_predictions, 
                        predictions_arguments)
         outp <- out
         if(!showlegends) outp <- outp + 
           ggplot2::theme(legend.position = 'none')
         return(outp)
       }
     } # if(call_predictions) {
     
     if(call_slopes) {
       if(!plot) {
         if(!average) {
           out <- do.call(marginaleffects::slopes, predictions_arguments)
         } else if(average) {
           out <- do.call(marginaleffects::avg_slopes, predictions_arguments)
         }
       } else if(plot) {
         out <- do.call(marginaleffects::plot_slopes, predictions_arguments)
         outp <- out
         outp <- outp + ggplot2::theme(legend.position = 'none')
         return(outp)
       }
     } # if(call_slopes) {
     
     onex0 <- out %>% marginaleffects::posterior_draws()
     
     if(isFALSE(constrats_by)) {
       constrats_by <- NULL
     } else if(!isFALSE(constrats_by)) {
       if(is.null(constrats_by)) {
         if(is.null(hypothesis)) {
           constrats_by <- NULL
         } else if(!is.null(hypothesis)) {
           if(!isFALSE(by)) constrats_by <- setdiff(by, xvar) 
         }
       } else if(!is.null(constrats_by)) {
         constrats_by <- constrats_by
       }
     }
     
     
     if(isFALSE(constrats_at)) {
       constrats_at <- NULL
     } else if(!isFALSE(constrats_at)) {
       if(is.null(constrats_at)) {
         constrats_at <- list()
         for(byi in by) {
           if(is.numeric(constrats_at[[byi]])) {
             constrats_at[[byi]] <- 'unique'
           }
         }
         if(length(constrats_at) == 0) constrats_at <- NULL
       }
     }
     
    
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
                collapse_comma(by)
           )
         }
         allowed_char_constrats_at <- c('max', 'min', 'unique', 'range')
         if(is.character(constrats_at[[caxi]])) {
           if(length(constrats_at[[caxi]]) > 1) {
             stop(caxi, " specified in 'constrats_at' as character should be a 
                  single character:",
                  "\n ", 
                  collapse_comma(allowed_char_constrats_at)
             )
           }
           if(!constrats_at[[caxi]] %in% allowed_char_constrats_at) {
             stop(constrats_at[[caxi]], " specified in 'constrats_at' as 
                  character should be one of",
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
     

     if(!is.null(constrats_subset)) {
       for (caxi in names(constrats_subset)) {
         if(!caxi %in% names(onex0)) {
           stop("Variable '", caxi, ". specified in 'constrats_subset' is not ",
                " avaialble in the 'by' argument.",
                "\n ", 
                " Please include ", caxi, " in the 'by' argument.",
                "\n ", 
                " The current 'by' argument includes:", 
                "\n ",
                collapse_comma(by)
           )
         }
         getatval <- constrats_subset[[caxi]]
         constrats_subset[[caxi]] <- getatval
         groupvarshyp1 <- c(caxi, groupvarshyp1)
         groupvarshyp2 <- c(caxi, groupvarshyp2)
       } # for (caxi in names(constrats_subset)) {
       
       
       for (caxi in names(constrats_subset)) {
         onex1 <- 
           base::subset(onex1, onex1[[caxi]] %in% constrats_subset[[caxi]])
         if(nrow(onex1) == 0) {
           stop(caxi, " specified in 'constrats_subset' resulted in zero rows:",
                "\n ", 
                ""
           )
         }
       }
     } # if(!is.null(constrats_subset)) {
     
     
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
     
     
     
     if(usedtplyr) {
       get_pe_ci2 <- get_pe_ci
       hypothesisargs <- formals(get_pe_ci2)
       hypothesisargs[['x']]          <- as.name('x')
       hypothesisargs[['draw']]       <- 'estimate'
       hypothesisargs[['...']]        <- NULL
       hypothesisargs[['na.rm']]      <- TRUE
       hypothesisargs <- base::append(hypothesisargs, as.name('...') , after = 1)
       names(hypothesisargs)[2] <- '...'
       formals(get_pe_ci2) <- hypothesisargs
       out_sf <- onex1 %>% dtplyr::lazy_dt() %>%
         dplyr::mutate(!! parmi_estimate := eval(parse(text = measurevar))) %>% 
         dplyr::group_by_at(groupvars2) %>% 
         dplyr::group_modify(., get_pe_ci2, .keep = F) %>% 
         dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                      fixed = TRUE))
     } else if(!usedtplyr) {
       out_sf <- onex1 %>% 
         dplyr::mutate(!! parmi_estimate := eval(parse(text = measurevar))) %>% 
         dplyr::reframe(
           dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, 
                         .unpack = TRUE),
           .by = dplyr::all_of(!! groupvars2)
         ) %>%
         dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                      fixed = TRUE))
     }
     
    
     
     if(!is.null(hypothesis)) {
       get_hypothesis_x_modify <- function(.x, hypothesis, by, draws, ...) {
         get_hypothesis_x(x = .x, hypothesis = hypothesis, by = by, 
                          draws = draws)
       }
       
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
       
       if(usedtplyr) {
         get_hypothesis_x_modifyx2 <- get_hypothesis_x
         hypothesisargs <- formals(get_hypothesis_x_modifyx2)
         hypothesisargs[['x']]          <- as.name('x')
         hypothesisargs[['hypothesis']] <- hypothesis
         hypothesisargs[['draws']]      <- as.name(parmi_estimate)
         hypothesisargs[['by']]         <- constrats_by
         hypothesisargs[['newdata']]    <- NULL
         hypothesisargs <- base::append(hypothesisargs, as.name('...') , 
                                        after = 1)
         names(hypothesisargs)[2] <- '...'
         formals(get_hypothesis_x_modifyx2) <- hypothesisargs
         
         get_pe_ci2 <- get_pe_ci
         hypothesisargs <- formals(get_pe_ci2)
         hypothesisargs[['x']]          <- as.name('x')
         hypothesisargs[['draw']]       <- 'estimate'
         hypothesisargs[['...']]        <- NULL
         hypothesisargs[['na.rm']]      <- TRUE
         hypothesisargs <- base::append(hypothesisargs, as.name('...') , 
                                        after = 1)
         names(hypothesisargs)[2] <- '...'
         formals(get_pe_ci2) <- hypothesisargs
         
         out_sf_hy <-
           onex1 %>% dtplyr::lazy_dt() %>%
           dplyr::group_by_at(groupvarshyp1) %>% 
           dplyr::mutate(!! parmi_estimate := eval(parse(text = measurevar))) %>% 
           dplyr::group_modify(., get_hypothesis_x_modifyx2, .keep = F) %>% 
           dtplyr::lazy_dt(.) %>%
           dplyr::group_by_at(groupvarshyp2) %>% 
           dplyr::group_modify(., get_pe_ci2, .keep = F) %>% 
           dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                        fixed = TRUE))
       } else if(!usedtplyr) {
         out_sf_hy <-
           onex1 %>% 
           dplyr::group_by_at(groupvarshyp1) %>% 
           dplyr::mutate(!! parmi_estimate := eval(parse(text = measurevar))) %>% 
           dplyr::group_modify(., ~get_hypothesis_x_modify(.x,
                                                           hypothesis = 
                                                             hypothesis, 
                                                           by = constrats_by, 
                                                           draws = estimate
           ), .keep = F) %>% 
           dplyr::ungroup() %>% 
           dplyr::reframe(
             dplyr::across(c(dplyr::all_of(parmi_estimate)), get_pe_ci, 
                           .unpack = TRUE),
             .by = dplyr::all_of(!! groupvarshyp2)
           ) %>%
           dplyr::rename_with(., ~ gsub(paste0(parmi_estimate, "_"), "", .x, 
                                        fixed = TRUE))
       }
     } # if(!is.null(hypothesis)) {
   } # if(method == 'custom') {
   
   
    # Restore function(s)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    
    if(!is.null(eval(predictions_arguments$clearenvfuns))) {
      if(!is.logical(eval(predictions_arguments$clearenvfuns))) {
        stop('clearenvfuns must be NULL or a logical')
      } else {
        setcleanup <- eval(predictions_arguments$clearenvfuns)
      }
    }
    
    if(is.null(eval(predictions_arguments$clearenvfuns))) {
      if(is.null(eval(predictions_arguments$usesavedfuns))) {
        predictions_arguments$usesavedfuns <- usesavedfuns
      }
      if(eval(predictions_arguments$usesavedfuns)) {
        setcleanup <- TRUE 
      } else {
        setcleanup <- FALSE
      }
    }
    
    # Cleanup environment if requested
    if(setcleanup) {
      # suppressWarnings({
        tempgenv <- envir
        for (oalli in names(oall)) {
          if(exists(oalli, envir = tempgenv )) {
            remove(list=oalli, envir = tempgenv)
          }
        }
        tempgenv <- test
        for (oalli in names(oall)) {
          if(exists(oalli, envir = tempgenv )) {
            remove(list=oalli, envir = tempgenv)
          }
        }
      # })
    } # if(setcleanup) {
    
    
    # fullframe - for marginal_draws, already fullframe returned 
    # if(!isTRUE(eval(fullframe))) setfullframe <- FALSE
    # if( isTRUE(eval(fullframe))) setfullframe <- TRUE
    setfullframe <- FALSE
    # fullframe
    full.args$idata_method <- idata_method
    full.args$fullframe <- eval(full.args$fullframe)
    full.args$summary   <- TRUE
    if(!is.null(eval(full.args$fullframe))) {
      if(eval(full.args$fullframe)) {
        if(!eval(full.args$summary)) {
          stop("fullframe can not be combined with summary = FALSE")
        }
        if(full.args$idata_method == 'm1') {
          stop("fullframe can not be combined with idata_method = 'm1'")
        }
      }
    }
    if(is.null(eval(full.args$fullframe))) {
      if (!is.na(model$model_info$univariate_by)) {
        if(full.args$idata_method == 'm1') setfullframe <- FALSE
        if(full.args$idata_method == 'm2') setfullframe <- TRUE
      } else {
        setfullframe <- FALSE
      }
    }
    if (!is.na(model$model_info$univariate_by)) {
      if(is.null(full.args$fullframe)) 
        full.args$fullframe <- fullframe <- FALSE
      if(full.args$fullframe & full.args$idata_method == 'm1') 
        setfullframe <- FALSE
      if(full.args$fullframe & full.args$idata_method == 'm2') 
        setfullframe <- TRUE
      if(!full.args$fullframe) 
        setfullframe <- FALSE
      if(setfullframe) {
        uvarby <- model$model_info$univariate_by
        uvarbyresp <- paste0(uvarby, resp)
        uvarbynewdata <- eval(full.args$newdata) %>% 
          dplyr::filter(!!dplyr::sym(uvarbyresp) == 1)
        if(setfullframe) out_sf <- cbind(out_sf, uvarbynewdata)
      }
    }
    
    
    
    if(is.null(reformat)) {
      if(is.null(hypothesis) && is.null(equivalence)) {
        reformat <- TRUE
      } else {
        reformat <- FALSE
      }
    }
    
    out_sf <- out_sf %>% data.frame() %>% 
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                  ~ round(., digits = digits)))
    
   
    if (reformat) {
      out_sf <- out_sf %>% 
        dplyr::rename(!!as.symbol(set_names_[1]) := 
                        dplyr::all_of('estimate')) %>% 
        dplyr::rename(!!as.symbol(set_names_[2]) := 
                        dplyr::all_of('conf.low')) %>% 
        dplyr::rename(!!as.symbol(set_names_[3]) := 
                        dplyr::all_of('conf.high')) %>% 
        data.frame()
      
      remove_cols_ <- c('term',  'tmp_idx', 'predicted_lo', 
                        'predicted_hi', 'predicted')
      
      out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]

      if(!is.null(out_sf_hy)) {
        out_sf_hy <- out_sf_hy %>% data.frame() %>% 
          dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                      ~ round(., digits = digits))) %>% 
          dplyr::rename(!!as.symbol(set_names_[1]) := 
                          dplyr::all_of('estimate')) %>% 
          dplyr::rename(!!as.symbol(set_names_[2]) := 
                          dplyr::all_of('conf.low')) %>% 
          dplyr::rename(!!as.symbol(set_names_[3]) := 
                          dplyr::all_of('conf.high')) %>% 
          dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
          data.frame()
      } # if(!is.null(out_sf_hy)) {
    }
    
    
    if(!is.null(out_sf_hy)) {
      row.names(out_sf) <- NULL
      row.names(out_sf_hy) <- NULL
      out_sf <- list(estimate = out_sf, contrast = out_sf_hy)
    } 
    
    return(out_sf)
  }


#' @rdname marginal_draws.bgmfit
#' @export
marginal_draws <- function(model, ...) {
  UseMethod("marginal_draws")
}


