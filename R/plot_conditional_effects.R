

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
           method = NULL,
           estimation_method = 'fitted',
           transform = NULL,
           transform_draws = NULL,
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
           ipts = NULL,
           deriv = 0,
           summary = FALSE,
           model_deriv = NULL,
           idata_method = NULL,
           verbose = FALSE,
           dummy_to_factor = NULL, 
           expose_function = FALSE,
           usesavedfuns = NULL,
           clearenvfuns = NULL,
           funlist = NULL,
           xvar = NULL,
           idvar = NULL,
           itransform = NULL,
           newdata_fixed = NULL,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- envir
    }
    
    if(!is.null(transform) & !is.null(transform_draws)) {
      stop("Please specify either transform or transform_draws, not both")
    }
    
    # 20.03.2025
    assign_function_to_environment(transform_draws, 'transform_draws', 
                                   envir = NULL)
    model$model_info[['transform_draws']] <- transform_draws
    
    # 20.03.2025
    # Depending on dpar 'mu' or 'sigma', subset model_info
    # This only when set_sigma_manual used to model a b c 
    # Not when a function such as splines::ns etc used in sigma_formula
    
    model <- getmodel_info(model = model, dpar = dpar, resp = resp)
    
    
    
    
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    # set_names_  <- c('Estimate', probtitles)
    
    

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
    
    
    rlang_trace_back <- rlang::trace_back()
    check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      # nothing
    } else {
      rlang_trace_back.bgmfit_i <- min(which(check_trace_back.bgmfit == TRUE))
      rlang_trace_back.bgmfit <- rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
      rlang_call_name <- rlang::call_name(rlang_trace_back.bgmfit)
      xcall <- rlang_call_name
    }
    
    
    
    check_if_package_installed(model, xcall = NULL)
    
    if(is.null(ndraws)) {
      ndraws <- brms::ndraws(model)
    }
    
    if(is.null(model_deriv)) {
      model_deriv <- TRUE
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
    if(is.null(levels_id) & is.null(idvar)) {
      idvar <- model$model_info[[groupvar_]]
      if (!is.null(model$model_info[[hierarchical_]])) {
        idvar <- model$model_info[[hierarchical_]]
      }
    } else if (!is.null(levels_id)) {
      idvar <- levels_id
    } else if (!is.null(idvar)) {
      idvar <- idvar
    }
    xvar  <- xvar
    idvar <- idvar
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
    full.args$model_deriv <- model_deriv
    
    
    
    
    full.args <- 
      sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                 arguments = full.args, 
                                 check_formalArgs = plot_conditional_effects.bgmfit,
                                 check_formalArgs_exceptions = c('object'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    if(!is.null(full.args[['transform_draws']])) {
      full.args[['transform']] <- full.args[['transform_draws']]
      if(verbose) message("'transform' set based on 'transform_draws'")
    }
    
    
    # Interpolation points
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)

    
    if(!is.null(model$xcall)) {
      arguments <- get_args_(as.list(match.call())[-1], model$xcall)
      newdata <- newdata
    } else {
      newdata <- CustomDoCall(get.newdata, full.args)
    }
    full.args$newdata <- newdata
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") model_deriv<- FALSE
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
    
    o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
    
    post_processing_checks_args[['all']]      <- TRUE
    oall <- CustomDoCall(post_processing_checks, post_processing_checks_args)
    post_processing_checks_args[['all']]      <- FALSE
    
   
    
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
        model_deriv <- FALSE
        call_slopes <- TRUE
        post_processing_checks_args[['deriv']]    <- 0
        o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    full.args$model_deriv <- model_deriv
    
    # 20.03.2025
    if(!is.null(model$model_info[['sigma_fun_mode']])) {
      sigma_fun_mode <- model$model_info[['sigma_fun_mode']]
      if(dpar == "sigma") {
        if(deriv > 0) {
          if(sigma_fun_mode == "inline") {
            check_fun    <- TRUE
            available_d1 <- FALSE
            model_deriv  <- FALSE
            call_slopes  <- TRUE
          }
        }
      }
    }
    
    
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall, 
                      usesavedfuns = usesavedfuns, 
                      deriv = post_processing_checks_args[['deriv']], 
                      envir = envir, 
                      model_deriv = model_deriv, 
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
    
    # 6.03.2025
    # see slopes will be mandatory
    check_fun <- FALSE
    if(deriv > 0) {
      available_d1 <- o[['available_d1']]
      if(!available_d1) {
        model_deriv <- FALSE
        call_slopes <- TRUE
        post_processing_checks_args[['deriv']]    <- 0
        o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    # post_processing_checks_args[['deriv']]    <- deriv
    
    
    
    
    
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
    
    # 6.03.2025
    calling.args <-
      sanitize_CustomDoCall_args(what = "CustomDoCallx",
                                 arguments = calling.args,
                                 check_formalArgs = plot_conditional_effects.bgmfit,
                                 check_formalArgs_exceptions = c('object', 'model'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    if(!is.null(calling.args[['transform_draws']])) {
      calling.args[['transform']] <- calling.args[['transform_draws']]
      if(verbose) message("'transform' set based on 'transform_draws'")
    }
    
    
    # Interpolation points
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    calling.args$ipts <- ipts <- check_ipts(ipts = calling.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
     
    
    if(check_fun) {
      if(!available_d1) {
        if(is.null(calling.args$re_formula)) {
          if(!spaghetti) stop("Please set 'spaghetti = TRUE'")
          # stop("For 'deriv = 1', use 'marginal_draws(..., deriv = 1)'",
          #      "\n ",
          #      " instead of 'plot_conditional_effects()'")
        }
      }
    }
    
    
    

    calling.args_ce <- calling.args_cefd <- calling.args
    calling.args_ce$newdata <- NULL
    calling.args_ce$x       <- calling.args_ce$object
    calling.args_ce$object  <- NULL
    if(!eval(full.args$model_deriv)) {
      if(is.null(calling.args_ce$re_formula)) {
        
        calling.args_ce_effects <- calling.args_ce
        calling.args_ce_effects[['effects']] <- paste0(xvar, ":", idvar)
        out_effects    <- CustomDoCall(brms::conditional_effects, calling.args_ce_effects)
        datace_effects <- out_effects[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
        
        out_    <- CustomDoCall(brms::conditional_effects, calling.args_ce)
        datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
        datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
        datace <- attrstrip(datace, keep = c('row.names', 'names', 'class'))
        
        datace <- attrstrip(datace_effects, keep = c('row.names', 'names', 'class'))

        # datace <- datace %>% dplyr::relocate(id, age)
        calling.args_cefd$newdata <- datace
        calling.args_cefd$model <- model
        calling.args_cefd$object <- NULL
        calling.args_cefd$xcall_str <- paste(deparse(xcallz), collapse = "")
        if (estimation_method == 'fitted') {
          outx <- CustomDoCall(fitted_draws, calling.args_cefd)
        } else if (estimation_method == 'predict') {
          outx <- CustomDoCall(predict_draws, calling.args_cefd)
        }
        # outx <-  CustomDoCall(fitted_draws, calling.args_cefd)
        # outx <- apply(outx, 2, mean)
        # outx0 <<- outx
        if(!summary) {
          outx <- brms::posterior_summary(outx , probs = probs, robust = robust)
          outx <- outx %>% data.frame()
        }
       
        if(spaghetti) {
          out_ <- out_effects
          out_[[1]][['estimate__']] <- outx[, 1] 
          out_[[1]][['se__']]       <- outx[, 2]
          out_[[1]][['lower__']]    <- outx[, 3]
          out_[[1]][['upper__']]    <- outx[, 4]
          . <- out_ <- plot(out_, plot = F)[[1]] + ggplot2::theme(legend.position = 'none')
        } else if(!spaghetti) {
          outx <-
          datace %>% dplyr::select(dplyr::all_of(c(xvar))) %>% dplyr::bind_cols(., outx) %>%
            dplyr::group_by_at(c(xvar)) %>%
            dplyr::summarise(dplyr::across(dplyr::where(is.numeric),
                                           ~ mean(.x, na.rm = TRUE))) %>%
            dplyr::ungroup() %>% dplyr::select(-dplyr::all_of(xvar)) %>%
            as.matrix()
          # outx <-
          #   datace %>% dplyr::select(dplyr::all_of(c(xvar))) %>% 
          #   dplyr::bind_cols(., outx) %>%
          #   dplyr::distinct(!! as.name('age'), .keep_all = T)
          
          out_[[1]][['estimate__']] <- outx[, 1] 
          out_[[1]][['se__']]       <- outx[, 2]
          out_[[1]][['lower__']]    <- outx[, 3]
          out_[[1]][['upper__']]    <- outx[, 4]
          # out_ <- plot(out_xx, plot = F)[[1]] + ggplot2::theme(legend.position = 'none')
          . <- out_
        } # if(spaghetti) {
         
        
      } else if(is.na(calling.args_ce$re_formula)) {
        out_    <- CustomDoCall(brms::conditional_effects, calling.args_ce)
        datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
        datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
        calling.args_cefd$newdata <- datace
        calling.args_cefd$model <- model
        calling.args_cefd$object <- NULL
        calling.args_cefd$xcall_str <- paste(deparse(xcallz), collapse = "")
        if (estimation_method == 'fitted') {
          outx <- CustomDoCall(fitted_draws, calling.args_cefd)
        } else if (estimation_method == 'predict') {
          outx <- CustomDoCall(predict_draws, calling.args_cefd)
        }
        # outx <-  CustomDoCall(fitted_draws, calling.args_cefd)
        if(!summary) {
          outx <- brms::posterior_summary(outx , probs = probs, robust = robust)
          outx <- outx %>% data.frame()
        }
        outx <- outx %>% dplyr::select(dplyr::all_of(set_names_))
        out_[[1]][['estimate__']] <- outx[, 1]
        out_[[1]][['se__']]       <- outx[, 2]
        out_[[1]][['lower__']]    <- outx[, 3]
        out_[[1]][['upper__']]    <- outx[, 4]
        . <- out_
      } # else if(is.na(calling.args_ce$re_formula)) {
    } # if(!eval(full.args$model_deriv)) {
    

    if(eval(full.args$model_deriv)) {
      calling.args_ce <- calling.args
      calling.args_ce$newdata <- NULL
      calling.args_ce$x <- calling.args_ce$object
      calling.args_ce$object <- NULL
      .   <- CustomDoCall(brms::conditional_effects, calling.args_ce)
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





