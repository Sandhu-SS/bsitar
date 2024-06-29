

#' Estimate growth curves
#' 
#' @description The \strong{marginal_draws()} function estimates and plots
#'   growth curves (distance and velocity) by using \pkg{marginaleffects}
#'   package as back-end. This function can compute growth curves (via
#'   [marginaleffects::predictions()]), average growth curves (via
#'   [marginaleffects::avg_predictions()]) or plot growth curves (via
#'   [marginaleffects::plot_predictions()]). Please see
#'   [here](https://marginaleffects.com/) for details.
#'  Note that \pkg{marginaleffects} package is highly flexible and therefore it
#'  is expected that user has a strong understanding of its working.
#'  Furthermore, since \pkg{marginaleffects} package is rapidly evolving, the
#'  results obtained from the current implementation should be considered
#'  experimental.
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
#' @param deriv An integer to indicate whether to estimate distance curve or its
#'   derivative (i.e., velocity curve). The \code{deriv = 0} (default) is for
#'   the distance curve whereas \code{deriv = 1} for the velocity curve. 
#'
#' @param method A character string to specify whether to make computation at
#'   post draw stage by using the \code{'marginaleffects'} machinery i.e.,
#'   [marginaleffects::comparisons()] (\code{method = 'pkg'}) or via
#'   the custom functions written for efficiency and speed (\code{method =
#'   'custom'}, default). Note that \code{method = 'custom'} is useful and
#'   rather needed when testing hypotheses. Note that when \code{method =
#'   'custom'}, [marginaleffects::predictions()] and not the
#'   [marginaleffects::comparisons()] is used internally.
#' 
#' @param constrats_by A character vector (default \code{NULL}) specifying the
#'   variable(s) by which hypotheses (post draw stage) should be tested
#'   (see \code{hypothesis} argument). Note that variable(s)
#'   specified in the \code{constrats_by} should be sub set of the variables
#'   included in the \code{'by'} argument. 
#'   
#' @param constrats_at A character vector (default \code{NULL}) specifying the
#'   variable(s) at which hypotheses (post draw stage) should be tested (see
#'   \code{hypothesis} argument). Note that variable(s) specified in the
#'   \code{constrats_at} should be sub set of the variables included in the
#'   \code{'by'} argument. The \code{constrats_at} is particularly useful when
#'   number of rows in the estimates is large. This is because the
#'   \pkg{marginaleffects} does not allow hypotheses testing when the number of
#'   rows in the estimates is more that 25.  
#'  
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
           usedtplyr = FALSE,
           usecollapse = TRUE,
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
           pdrawsp = FALSE, 
           pdrawsh = FALSE, 
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
                                          minimum_version = 
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
                                            minimum_version =
                                              get_package_minversion(
                                                'dtplyr'
                                              ),
                                            prompt = FALSE,
                                            stop = FALSE))
      
      try(zz <- insight::check_if_installed(c("data.table"), 
                                            minimum_version =
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
    drawid <- NULL;
    
    
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
                               minimum_version = get_package_minversion('brms'), 
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
        usedtplyr,
        usecollapse,
        pdrawsp,
        pdrawsh
      )
    ))[-1]
    
    # don't exclude variables
   # if(call_predictions) exclude_args <- c(exclude_args, 'variables')
    
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
   # if(!is.null(predictions_arguments[['by']])) {
   #   checbyx <- predictions_arguments[['by']]
   #   if(all(checbyx == "")) method <- 'pkg'
   #   if(is.logical(checbyx)) {
   #     if(!checbyx) method <- 'pkg'
   #   }
   # }
   
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
   
   
   pdrawsp_est <- NULL
   pdrawsh_est <- NULL
   # pdraws_est <- NULL
    
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
     
     if(!isFALSE(pdrawsp)) {
       selectchoicesr <- c("return", 'add') 
       checkmate::assert_choice(pdrawsp, choices = selectchoicesr)
       if(pdrawsp == 'return') {
         return(onex0)
       } else if(pdrawsp == 'add') {
         pdrawsp_est <- onex0
       } else {
         
       }
     }
     
     
     setdrawidparm <- by
     namesx <- c('estimate', 'conf.low', 'conf.high')
     setdrawidparm_ <- c(setdrawidparm, namesx)
     
     out_sf <- onex0 %>% collapse::fsubset(., drawid == 1) %>%
       collapse::fselect(., setdrawidparm_)
     
     
     
     if(!is.null(hypothesis)) {
       # For hypothesis
       groupvarshyp1 <- c('drawid')
       groupvarshyp2 <- c('term')
       if(!is.null(constrats_at)) {
         if(!is.character(constrats_at)) 
           stop("'constrats_at' must be a character vector")
         for (caxi in constrats_at) {
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
           groupvarshyp1 <- c(caxi, groupvarshyp1)
           groupvarshyp2 <- c(caxi, groupvarshyp2)
         } # for (caxi in names(constrats_at)) {
       } # if(!is.null(constrats_at)) {
       
       
       if(is.null(constrats_by)) {
         stop("Please specify 'constrats_by' argument when testing 'hypothesis'",
              "\n ",
              " The available options are: ",
              collapse_comma(by)
         )
       }
       
       if(!is.null(constrats_by)) {
         if(!is.character(constrats_by)) 
           stop("'constrats_by' must be a character vector")
         for (caxi in constrats_by) {
           if(!caxi %in% names(onex0)) {
             stop("Variable '", caxi, ". specified in 'constrats_by' is not in",
                  "\n ", 
                  " the estimates. Note that ", caxi, " should also be included",
                  "\n ", 
                  " in the 'by' argument. The current 'by' argument includes:", 
                  "\n ",
                  collapse_comma(by)
             )
           }
         } # for (caxi in names(constrats_by)) {
       } # if(!is.null(constrats_by)) {
       
       
       if(nrow(out_sf) > 25) {
         if(is.null(constrats_at)) {
           cat(" Note that the 'marginaleffects' package does not allow" ,
               "\n",
               "'hypothesis' argument when estimates rows are more than 25",
               "\n",
               "To avoid this issue, you can use 'constrats_at' argument",
               "\n"
           )
         }
       }
       
       
       get_pe_ci_collapse <- function(x, na.rm = TRUE,...) {
         if(ec_agg == "mean")  estimate <- 
             collapse::fmean(x, 
                             na.rm = na.rm, 
                             nthreads = arguments$cores) 
         
         if(ec_agg == "median") estimate <- 
             collapse::fmedian(x, 
                               na.rm = na.rm, 
                               nthreads = arguments$cores)
         
         if(ei_agg == "eti") luci = collapse::fquantile(x, probs = probs, 
                                                        na.rm = na.rm)
         if(ei_agg == "hdi") luci = get_hdix(x, credMass = conf)
         cbind(estimate, luci[1], luci[2]) 
       }
       
       set_constrats_by <- c(constrats_by, 'draw')
       namesx <- c('estimate', 'conf.low', 'conf.high')
       
       setdrawidparm_at <- c(constrats_at, 'term')
       setdrawidparm_at_ <- c(setdrawidparm_at, namesx)
       
       
       temhyy <-
         onex0 %>% 
         collapse::fgroup_by(groupvarshyp1) %>%
         collapse::fselect(set_constrats_by) %>% 
         collapse::frename('estimate' = 'draw') %>% 
         collapse::fsummarise(collapse::qDF(
           get_hypothesis_x(.data,
                            by = constrats_by,
                            hypothesis = hypothesis,
                            draws = 'estimate'))) 
       
       out_sf_hy <- 
       temhyy %>%
         collapse::fgroup_by(groupvarshyp2) %>%
         collapse::fsummarise(collapse::mctl(
           get_pe_ci_collapse(.data[['estimate']]))
         ) %>%
         collapse::frename(., setdrawidparm_at_)
       
       
       if(!isFALSE(pdrawsh)) {
         selectchoicesr <- c("return", 'add') 
         checkmate::assert_choice(pdrawsh, choices = selectchoicesr)
         if(pdrawsh == 'return') {
           return(temhyy)
         } else if(pdrawsh == 'add') {
           pdrawsh_est <- temhyy
         } else {
           
         }
       }
       
       # out_sf_hy <-
       #   onex0 %>% 
       #   collapse::fgroup_by(groupvarshyp1) %>%
       #   collapse::fselect(set_constrats_by) %>% 
       #   collapse::frename('estimate' = 'draw') %>% 
       #   collapse::fsummarise(collapse::qDF(
       #     get_hypothesis_x(.data,
       #                           by = constrats_by,
       #                           hypothesis = hypothesis,
       #                           draws = 'estimate'))) %>%
       #   collapse::fgroup_by(groupvarshyp2) %>%
       #   collapse::fsummarise(collapse::mctl(
       #     get_pe_ci_collapse(.data[['estimate']]))
       #   ) %>%
       #   collapse::frename(., setdrawidparm_at_)
       
       
       
       
     } # if(!is.null(hypothesis)) {
   } # if(method == 'custom') {
   
   
   
   out_sf <- out_sf %>% data.frame() %>% 
     dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                 ~ round(., digits = digits))) %>%
     data.frame()
   
   
   
   if(!is.null(pdrawsh_est)) {
     if(usecollapse) {
       pdrawsh_est <- pdrawsh_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     } else {
       pdrawsh_est <- pdrawsh_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     }
   }
   
   
   
   if(!is.null(pdrawsp_est)) {
     if(usecollapse) {
       pdrawsp_est <- pdrawsp_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     } else {
       pdrawsp_est <- pdrawsp_est %>% data.frame() %>% 
         dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                     ~ round(., digits = digits)))
     }
   }
   
   
   if(is.null(reformat)) {
     if(is.null(hypothesis) && is.null(equivalence)) {
       reformat <- TRUE
     } else {
       reformat <- FALSE
     }
   }
   
   if (reformat) {
     out_sf <- out_sf %>% 
       dplyr::rename(!!as.symbol(set_names_[1]) := 
                       dplyr::all_of('estimate')) %>% 
       dplyr::rename(!!as.symbol(set_names_[2]) := 
                       dplyr::all_of('conf.low')) %>% 
       dplyr::rename(!!as.symbol(set_names_[3]) := 
                       dplyr::all_of('conf.high')) %>% 
       data.frame()
     
     if(method == 'pkg') {
       remove_cols_ <- c('tmp_idx', 'predicted_lo', 
                         'predicted_hi', 'predicted', 'rowid')
     } else if(method == 'custom') {
       remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                         'predicted_hi', 'predicted', 'rowid')
     }
     
     
     out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
     row.names(out_sf) <- NULL
     
     
     
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
     
     
     if(!is.null(pdrawsp_est)) {
       pdrawsp_est <- pdrawsp_est %>% 
         # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         dplyr::rename(!!as.symbol(set_names_[2]) := 
                         dplyr::all_of('conf.low')) %>% 
         dplyr::rename(!!as.symbol(set_names_[3]) := 
                         dplyr::all_of('conf.high')) %>% 
         dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
         data.frame()
     } # if(!is.null(pdrawsp_est)) {
     
     if(!is.null(pdrawsh_est)) {
       pdrawsh_est <- pdrawsh_est %>% 
         # dplyr::mutate(dplyr::across(dplyr::all_of('parameter'), toupper)) %>% 
         dplyr::rename(!!as.symbol(set_names_[1]) := 
                         dplyr::all_of('estimate')) %>% 
         dplyr::rename(!!as.symbol(set_names_[2]) := 
                         dplyr::all_of('conf.low')) %>% 
         dplyr::rename(!!as.symbol(set_names_[3]) := 
                         dplyr::all_of('conf.high')) %>% 
         dplyr::rename_with(., ~ sub("(.)", "\\U\\1", .x, perl = TRUE)) %>% 
         data.frame()
     } # if(!is.null(pdrawsh_est)) {
     
     
   } # if (reformat) {
   
   
   if(!is.null(out_sf_hy)) {
     out_sf <- out_sf %>% dplyr::ungroup()
     out_sf_hy <- out_sf_hy %>% dplyr::ungroup()
     out_sf <- list(estimate = out_sf, contrast = out_sf_hy)
   } 
   
   
   if(!is.null(pdrawsp_est)) {
     pdrawsp_est <- pdrawsp_est %>% dplyr::ungroup() 
     out_sf <- base::append(out_sf, list(pdrawsp_est = pdrawsp_est), after = 0)
   }
   
   if(!is.null(pdrawsh_est)) {
     pdrawsh_est <- pdrawsh_est %>% dplyr::ungroup() 
     out_sf <- base::append(out_sf, list(pdrawsh_est = pdrawsh_est), after = 0)
   }
   
   return(out_sf)
  }


#' @rdname marginal_draws.bgmfit
#' @export
marginal_draws <- function(model, ...) {
  UseMethod("marginal_draws")
}


