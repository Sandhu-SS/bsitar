

#' Visualize conditional effects of predictor
#'
#' @details The \strong{plot_conditional_effects()} is a wrapper around the
#'   [brms::conditional_effects()]. The [brms::conditional_effects()] function
#'   from the \pkg{brms} package can used to plot the fitted (distance) curve
#'   when response (e.g., height) is not transformed. However, when the outcome
#'   is log or square root transformed, the [brms::conditional_effects()] will
#'   return the fitted curve on the log or square root scale whereas the
#'   \strong{plot_conditional_effects()} will return the fitted curve on the
#'   original scale. Furthermore, the \strong{plot_conditional_effects()} also
#'   plots the velocity curve on the original scale after making required
#'   back-transformation. Apart from these differences, both these functions
#'   ([brms::conditional_effects] and \strong{plot_conditional_effects()} work
#'   in the same manner. In other words, user can specify all the arguments
#'   which are available in the [brms::conditional_effects()].
#' 
#' @inherit brms::conditional_effects params description
#' @inherit growthparameters.bgmfit params
#' @inherit plot_curves.bgmfit params
#' @inherit fitted_draws.bgmfit params
#'
#' @param ... Additional arguments passed to the [brms::conditional_effects()]
#'   function. Please see [brms::conditional_effects()] for details.
#'
#' @return An object of class 'brms_conditional_effects' which is a named list
#'   with one data.frame per effect containing all information required to
#'   generate conditional effects plots. See brms::conditional_effects for
#'   details.
#'
#' @export plot_conditional_effects.bgmfit
#' @export
#' 
#' @seealso [brms::conditional_effects()]
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
#' # Population average distance curve
#' plot_conditional_effects(model, deriv = 0, re_formula = NA)
#' 
#' \donttest{
#' # Individual-specific distance curves
#' plot_conditional_effects(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' plot_conditional_effects(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' plot_conditional_effects(model, deriv = 1, re_formula = NULL)
#' }
#' 
plot_conditional_effects.bgmfit <-
  function(model,
           effects = NULL,
           conditions = NULL,
           int_conditions = NULL,
           re_formula = NA,
           spaghetti = FALSE,
           surface = FALSE,
           categorical = FALSE,
           ordinal = FALSE,
           method = 'posterior_epred',
           transform = NULL,
           resolution = 100,
           select_points = 0,
           too_far = 0,
           prob = 0.95,
           robust = TRUE,
           newdata = NULL,
           ndraws = NULL,
           dpar = NULL,
           draw_ids = NULL,
           levels_id = NULL,
           resp = NULL,
           ipts = 10,
           deriv = 0,
           deriv_model = NULL,
           idata_method = NULL,
           verbose = FALSE,
           dummy_to_factor = NULL, 
           expose_function = FALSE,
           usesavedfuns = NULL,
           clearenvfuns = NULL,
           funlist = NULL,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- parent.frame()
    }
    
    
    # Depending on dpar 'mu' or 'sigma', subset model_info
    model <- getmodel_info(model = model, dpar = dpar)
    

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
    } else { 
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
    
    check_if_package_installed(model, xcall = NULL)
    
    if(is.null(ndraws)) {
      ndraws <- brms::ndraws(model)
    }
    
    if(is.null(deriv_model)) {
      deriv_model <- TRUE
    }
    
    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
    
    
    
    full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                    fargs = formals(), 
                                    dargs = list(...), 
                                    verbose = verbose)
    
    full.args$model <- model
    full.args$deriv_model <- deriv_model

    
    if(!is.null(model$xcall)) {
      arguments <- get_args_(as.list(match.call())[-1], model$xcall)
      newdata <- newdata
    } else {
      newdata <- do.call(get.newdata, full.args)
    }
    full.args$newdata <- newdata
    
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
    
    
    if(!is.null(funlist)) {
      if(!is.list(funlist)) {
        stop("funlist must be a list")
      } else {
        o <- funlist
      }
    }
    
    
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
    
    
    misc <- c("verbose", "usesavedfuns", "clearenvfuns", 
              "envir", "fullframe")
    
    calling.args <- post_processing_args_sanitize(model = model,
                                                  xcall = match.call(),
                                                  resp = resp,
                                                  envir = envir,
                                                  deriv = deriv, 
                                                  dots = list(...),
                                                  misc = misc,
                                                  verbose = verbose)
    
    
    calling.args$object <- full.args$model
    if(is.null(calling.args$newdata)) {
      if(!is.null(newdata)) calling.args$newdata <- newdata
    }
    
    
    if(!eval(full.args$deriv_model)) {
      if (is.null(resp)) {
        resp_rev_ <- resp
      } else if (!is.null(resp)) {
        resp_rev_ <- paste0("_", resp)
      }
      xvar_ <- paste0('xvar', resp_rev_)
      yvar_ <- paste0('yvar', resp_rev_)
      groupvar_ <- paste0('groupvar', resp_rev_)
      xvar <- model$model_info[[xvar_]]
      yvar <- model$model_info[[yvar_]]
      hierarchical_ <- paste0('hierarchical', resp_rev_)
      if (is.null(levels_id)) {
        IDvar <- model$model_info[[groupvar_]]
        if (!is.null(model$model_info[[hierarchical_]])) {
          IDvar <- model$model_info[[hierarchical_]]
        }
      } else if (!is.null(levels_id)) {
        IDvar <- levels_id
      }
      xvar  <- xvar
      idvar <- IDvar
      if(length(idvar) > 1) idvar <- idvar[1]
      yvar  <- 'yvar'
      calling.args_ce <- calling.args_cefd <- calling.args
      calling.args_ce$newdata <- NULL
      calling.args_ce$x <- calling.args_ce$object
      calling.args_ce$object <- NULL
      out_    <- do.call(brms::conditional_effects, calling.args_ce)
      datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
      datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
      calling.args_cefd$newdata <- datace
      calling.args_cefd$model <- model
      calling.args_cefd$object <- NULL
      outx <-  do.call(fitted_draws, calling.args_cefd)
      out_[[1]][['estimate__']] <- outx[, 1]
      out_[[1]][['se__']] <- outx[, 2]
      out_[[1]][['lower__']] <- outx[, 3]
      out_[[1]][['upper__']] <- outx[, 4]
      . <- out_
    }
    
    if(eval(full.args$deriv_model)) {
      calling.args_ce <- calling.args
      calling.args_ce$newdata <- NULL
      calling.args_ce$x <- calling.args_ce$object
      calling.args_ce$object <- NULL
      .   <- do.call(brms::conditional_effects, calling.args_ce)
    }
    
    
    # Restore function(s)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    
    if(!is.null(eval(full.args$clearenvfuns))) {
      if(!is.logical(eval(full.args$clearenvfuns))) {
        stop('clearenvfuns must be NULL or a logical')
      } else {
        setcleanup <- eval(full.args$clearenvfuns)
      }
    }
    
    if(is.null(eval(full.args$clearenvfuns))) {
      if(is.null(eval(full.args$usesavedfuns))) {
        full.args$usesavedfuns <- usesavedfuns
      }
      if(eval(full.args$usesavedfuns)) {
        setcleanup <- TRUE 
      } else {
        setcleanup <- FALSE
      }
    }
    
    # Cleanup environment if requested
    if(setcleanup) {
      suppressWarnings({
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
      })
    } # if(setcleanup) {
    
    .
  }


#' @rdname plot_conditional_effects.bgmfit
#' @export
plot_conditional_effects <- function(model, ...) {
  UseMethod("plot_conditional_effects")
}





