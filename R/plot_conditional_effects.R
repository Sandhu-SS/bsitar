

#' @title Visualize conditional effects of predictor
#'
#' @details The \strong{plot_conditional_effects()} is a wrapper around the
#'   [brms::conditional_effects()]. The [brms::conditional_effects()] function
#'   from the \pkg{brms} package can be used to plot the fitted (distance) curve
#'   when response (e.g., height) is not transformed. However, when the outcome
#'   is log or square root transformed, the [brms::conditional_effects()] will
#'   return the fitted curve on the log or square root scale, whereas the
#'   \strong{plot_conditional_effects()} will return the fitted curve on the
#'   original scale. Furthermore, the \strong{plot_conditional_effects()} also
#'   plots the velocity curve on the original scale after making the required
#'   back-transformation. Apart from these differences, both these functions
#'   ([brms::conditional_effects] and \strong{plot_conditional_effects()}) work
#'   in the same manner. In other words, the user can specify all the arguments
#'   which are available in the [brms::conditional_effects()]. An alternative
#'   approach is to [marginal_draws()] function (with \code{plot = TRUE}) which
#'   is based on the \pkg{marginaleffects}.
#' 
#' @inherit brms::conditional_effects params description
#' @inherit growthparameters.bgmfit params
#' @inherit plot_curves.bgmfit params
#' @inherit fitted_draws.bgmfit params
#'
#' @param ... Additional arguments passed to the [brms::conditional_effects()]
#'   function. Please see [brms::conditional_effects()] for details.
#'
#' @return An object of class \code{'brms_conditional_effects'}, which is a
#'   named list with one data.frame per effect containing all information
#'   required to generate conditional effects plots. See
#'   [brms::conditional_effects()] for details.
#'
#' @rdname plot_conditional_effects
#' @export
#' 
#' @seealso [brms::conditional_effects()]
#'
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
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
           probs = c(0.025, 0.975),
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
           itransform = NULL,
           newdata_fixed = FALSE,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- parent.frame()
    }
    
    # conf <- conf_level
    # probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    # set_names_  <- c('Estimate', probtitles)
    
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
    
    # 6.03.2025
    # moved here up in the call
    ########################################################
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
    
    # 6.03.2025
    ########################################################
    # prepare_data2
    ifunx_ <- paste0('ixfuntransform2', resp_rev_)
    ifunx_ <- model$model_info[[ifunx_]]
    ########################################################
    
    ########################################################
    # Define lables fun for x- axis
    labels_ggfunx <- function(...) {
      out <- ifunx_(list(...)[[1]]) 
      out <- scales::number(
        out,
        accuracy = 1,
        scale = 1,
        prefix = "",
        suffix = "",
        big.mark = " ",
        decimal.mark = ".",
        style_positive = c("none", "plus", "space"),
        style_negative = c("hyphen", "minus", "parens"),
        scale_cut = NULL,
        trim = TRUE
      )
      return(out)
    }
    
    labels_ggfunx_str <- 
      "ggplot2::scale_x_continuous(labels = labels_ggfunx)"
    
    ########################################################
    
    
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
    
    
    setxcall_   <- match.call()
    post_processing_checks_args <- list()
    post_processing_checks_args[['model']]    <- model
    post_processing_checks_args[['xcall']]    <- setxcall_
    post_processing_checks_args[['resp']]     <- resp
    post_processing_checks_args[['envir']]    <- envir
    post_processing_checks_args[['deriv']]    <- deriv
    post_processing_checks_args[['all']]      <- FALSE
    post_processing_checks_args[['verbose']]  <- verbose
    post_processing_checks_args[['check_d0']] <- FALSE
    post_processing_checks_args[['check_d1']] <- TRUE
    post_processing_checks_args[['check_d2']] <- FALSE
    
    o    <- do.call(post_processing_checks, post_processing_checks_args)
    
    post_processing_checks_args[['all']]      <- TRUE
    oall <- do.call(post_processing_checks, post_processing_checks_args)
    post_processing_checks_args[['all']]      <- FALSE
    
    
    # o <- post_processing_checks(model = model,
    #                             xcall = match.call(),
    #                             resp = resp,
    #                             envir = envir,
    #                             deriv = deriv, 
    #                             all = FALSE,
    #                             verbose = verbose)
    # 
    # oall <- post_processing_checks(model = model,
    #                                xcall = match.call(),
    #                                resp = resp,
    #                                envir = envir,
    #                                deriv = deriv, 
    #                                all = TRUE,
    #                                verbose = FALSE)
    
    
    if(!is.null(funlist)) {
      if(!is.list(funlist)) {
        stop("funlist must be a list")
      } else {
        o <- funlist
      }
    }
    
    
    # 6.03.2025
    # see slopes will be mandatory
    check_fun <- FALSE
    if(deriv > 0) {
      available_d1 <- o[['available_d1']]
      if(!available_d1) {
        deriv_model <- FALSE
        call_slopes <- TRUE
        post_processing_checks_args[['deriv']]    <- 0
        o    <- do.call(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    full.args$deriv_model <- deriv_model
    
    # post_processing_checks_args[['deriv']]    <- deriv
    
    
    # Unlike marginal_... functions where assign() works, for brms::fitted...
    # the envrionment is assigned to d0/d1 functions via setupfuns()
    # Therefore check_fun block is moved above 
    # and deriv = 0 shoud also reflect in  setupfuns() 
    
    # test <- setupfuns(model = model, resp = resp,
    #                   o = o, oall = oall, 
    #                   usesavedfuns = usesavedfuns, 
    #                   deriv = deriv, envir = envir, 
    #                   deriv_model = deriv_model, 
    #                   ...)
    
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall, 
                      usesavedfuns = usesavedfuns, 
                      deriv = post_processing_checks_args[['deriv']], 
                      envir = envir, 
                      deriv_model = deriv_model, 
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
    
    # 6.03.2025
    # see slopes will be mandatory
    check_fun <- FALSE
    if(deriv > 0) {
      available_d1 <- o[['available_d1']]
      if(!available_d1) {
        deriv_model <- FALSE
        call_slopes <- TRUE
        post_processing_checks_args[['deriv']]    <- 0
        o    <- do.call(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    # post_processing_checks_args[['deriv']]    <- deriv
    
    if(check_fun) {
      if(!available_d1) {
        stop("For 'deriv = 1', use 'marginal_draws(..., deriv = 1)'",
             "\n ",
             " instead of 'plot_conditional_effects()'") 
      }
    }
    
    
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
    
    # 6.03.2025
    xcallz <- match.call()
    
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
    

    calling.args_ce <- calling.args_cefd <- calling.args
    calling.args_ce$newdata <- NULL
    calling.args_ce$x       <- calling.args_ce$object
    calling.args_ce$object  <- NULL
    
    if(!eval(full.args$deriv_model)) {
      if(is.null(calling.args_ce$re_formula)) {
        out_    <- do.call(brms::conditional_effects, calling.args_ce)
        datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
        datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
        calling.args_cefd$newdata <- datace
        calling.args_cefd$model <- model
        calling.args_cefd$object <- NULL
        calling.args_cefd$xcall_str <- paste(deparse(xcallz), collapse = "")
        outx <-  do.call(fitted_draws, calling.args_cefd)
        out_ < out_ + ggplot2::theme(legend.position = 'none')
        outx <- outx %>% dplyr::select(dplyr::all_of(set_names_))
        out_[[1]][['estimate__']] <- outx[, 1]
        out_[[1]][['se__']]       <- outx[, 2]
        out_[[1]][['lower__']]    <- outx[, 3]
        out_[[1]][['upper__']]    <- outx[, 4]
        . <- out_ 
      } else if(is.na(calling.args_ce$re_formula)) {
        out_    <- do.call(brms::conditional_effects, calling.args_ce)
        datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
        datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
        calling.args_cefd$newdata <- datace
        calling.args_cefd$model <- model
        calling.args_cefd$object <- NULL
        calling.args_cefd$xcall_str <- paste(deparse(xcallz), collapse = "")
        outx <-  do.call(fitted_draws, calling.args_cefd)
        outx <- outx %>% dplyr::select(dplyr::all_of(set_names_))
        out_[[1]][['estimate__']] <- outx[, 1]
        out_[[1]][['se__']]       <- outx[, 2]
        out_[[1]][['lower__']]    <- outx[, 3]
        out_[[1]][['upper__']]    <- outx[, 4]
        . <- out_
      } # else if(is.na(calling.args_ce$re_formula)) {
    } # if(!eval(full.args$deriv_model)) {
    

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
    # 6.03.2025
    suppressMessages({. <- plot(., plot = FALSE)[[1]] + ept(labels_ggfunx_str)})
    .
  }


#' @rdname plot_conditional_effects
#' @export
plot_conditional_effects <- function(model, ...) {
  UseMethod("plot_conditional_effects")
}





