



#' @title Predicted values from the posterior predictive distribution
#'
#' @description The \strong{predict_draws()} function is a wrapper around the
#'   [brms::predict.brmsfit()] function, which obtains predicted values (and
#'   their summary) from the posterior distribution. See
#'   [brms::predict.brmsfit()] for details. An alternative approach is to
#'   [marginal_draws()] function which is based on the \pkg{marginaleffects}.
#'
#' @details The \strong{predict_draws()} function computes the expected values
#'   from the posterior distribution. The [brms::predict.brmsfit()] function
#'   from the \pkg{brms} package can be used to obtain predicted (distance)
#'   values when the outcome (e.g., height) is untransformed. However, when the
#'   outcome is log or square root transformed, the [brms::predict.brmsfit()]
#'   function will return the expected curve on the log or square root scale. In
#'   contrast, the \strong{predict_draws()} function returns the expected values
#'   on the original scale. Furthermore, \strong{predict_draws()} also computes
#'   the first derivative (velocity), again on the original scale, after making
#'   the necessary back-transformation. Aside from these differences, both
#'   functions ([brms::predict.brmsfit()] and \strong{predict_draws()}) work
#'   similarly. In other words, the user can specify all the options available
#'   in [brms::predict.brmsfit()].
#'
#' @param ... Additional arguments passed to the [brms::predict.brmsfit()]
#'   function. Please see [brms::predict.brmsfit()] for details on the various
#'   options available.
#'
#' @return An array of predicted response values. See [brms::predict.brmsfit()]
#'   for details.
#'
#' @inherit growthparameters.bgmfit params
#' @inherit fitted_draws.bgmfit params
#' @inherit brms::predict.brmsfit params
#'
#' @rdname predict_draws
#' @export
#'
#' @seealso [brms::predict.brmsfit()]
#'
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model
#'
#' # To avoid mode estimation, which takes time, the Bayesian SITAR model is fit 
#' # to the 'berkeley_exdata' and saved as an example fit ('berkeley_exfit').
#' # See the 'bsitar' function for details on 'berkeley_exdata' and 
#' # berkeley_exfit'.
#'
#' # Check and confirm whether the model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#'
#' model <- berkeley_exfit
#' 
#' # Population average distance curve
#' predict_draws(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' predict_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' predict_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' predict_draws(model, deriv = 1, re_formula = NULL)
#'  }
#'  
predict_draws.bgmfit <-
  function(model,
           newdata = NULL,
           resp = NULL,
           dpar = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           re_formula = NA,
           allow_new_levels = FALSE,
           sample_new_levels = "uncertainty",
           incl_autocor = TRUE,
           numeric_cov_at = NULL,
           levels_id = NULL,
           avg_reffects = NULL,
           aux_variables = NULL,
           ipts = NULL,
           deriv = 0,
           model_deriv = TRUE,
           summary = TRUE,
           robust = FALSE,
           transform_draws = NULL,
           probs = c(0.025, 0.975),
           xrange = NULL,
           xrange_search = NULL,
           parms_eval = FALSE,
           parms_method = 'getPeak',
           idata_method = NULL,
           verbose = FALSE,
           fullframe = NULL,
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
    
    # 20.03.2025
    assign_function_to_environment(transform_draws, 'transform_draws', 
                                   envir = NULL)
    model$model_info[['transform_draws']] <- transform_draws
    
    # 20.03.2025
    # Depending on dpar 'mu' or 'sigma', subset model_info
    # This only when set_sigma_manual used to model a b c 
    # Not when a function such as splines::ns etc used in sigma_formula
    
    if(is.null(dpar)) {
      dpar <- "mu"
    }
    
    model <- getmodel_info(model = model, dpar = dpar, resp = resp)
    
    
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
    
    if(is.null(model_deriv)) {
      model_deriv <- TRUE
    }
    
    if (is.null(idata_method)) {
      idata_method <- 'm2'
    }
    
    rlang_trace_back <- rlang::trace_back()
    
    check_trace_back.bgmfit <- grepl("plot_conditional_effects", 
                                     rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      plot_conditional_effects_calling <- FALSE
    } else {
      plot_conditional_effects_calling <- TRUE 
    }
    
    
    check_trace_back.bgmfit <- grepl("growthparameters", 
                                     rlang_trace_back[[1]])
    if(all(!check_trace_back.bgmfit)) {
      growthparameters_calling <- FALSE
    } else {
      growthparameters_calling <- TRUE 
    }
    
    if(any(grepl("modelbased_growthparameters", rlang_trace_back[[1]]))) {
      growthparameters_calling <- FALSE 
    }
    
    if(any(grepl("plot_curves", rlang_trace_back[[1]]))) {
      growthparameters_calling <- FALSE 
    }
    
    indirectcall <- FALSE
    if(!plot_conditional_effects_calling) {
      if(!is.null(model$xcall)) {
        arguments <- get_args_(as.list(match.call())[-1], model$xcall)
        full.args <- evaluate_call_args(cargs = arguments, 
                                        fargs = NULL, 
                                        dargs = NULL, 
                                        verbose = verbose)
        full.args$object <- full.args$model
        newdata          <- full.args$newdata
        indirectcall <- TRUE
      } else {
        full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                        fargs = formals(), 
                                        dargs = list(...), 
                                        verbose = verbose)
        full.args$model <- model
      }
    }
    
    xcall_str <- NULL
    if(plot_conditional_effects_calling) {
      full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                      fargs = formals(), 
                                      dargs = list(...), 
                                      verbose = verbose)
      # 6.03.2025
      xcall_str           <- full.args$xcall_str
      full.args$xcall_str <- NULL
    }
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") model_deriv<- FALSE
    }
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
    
    # 6.03.2025
    if(is.null(xcall_str)) {
      setxcall_   <- match.call()
    } else {
      setxcall_ <- xcall_str
    }
    
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
        o <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    
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
    
    
    # The deriv = 0/1 should also reflect in  setupfuns() 
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall,
                      usesavedfuns = usesavedfuns,
                      deriv = post_processing_checks_args[['deriv']],
                      envir = envir,
                      model_deriv = model_deriv,
                      ...)
    
    if(is.null(test)) {
      return(invisible(NULL))
    }
    
    
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
    
    misc <- c("verbose", "usesavedfuns", "clearenvfuns", 
              "envir", "fullframe")
    
    if(!indirectcall) {
      calling.args <- post_processing_args_sanitize(model = model,
                                                    xcall = match.call(),
                                                    resp = resp,
                                                    envir = envir,
                                                    deriv = deriv, 
                                                    dots = list(...),
                                                    misc = misc,
                                                    verbose = verbose)
    } else if(indirectcall) {
      calling.args <- full.args
    }
    
    calling.args$object <- full.args$model
    if(is.null(calling.args$newdata)) {
      if(!is.null(newdata)) calling.args$newdata <- newdata
    }
    
    
    if(growthparameters_calling) {
      if(calling.args$re_formula_opt == "V") {
        calling.args$re_formula <- NULL
      } else if(calling.args$re_formula_opt == "v") {
        calling.args$re_formula <- NA
      }
    }
    
    # 6.03.2025
    if(is.null(newdata) & !indirectcall ) {
      calling.args_newdata         <- calling.args
      calling.args_newdata$model   <- calling.args_newdata$object
      calling.args_newdata$newdata <- model$model_info$bgmfit.data 
      calling.args_newdata$dpar    <- dpar
      newdata <- CustomDoCall(get.newdata, calling.args_newdata)
      # if(growthparameters_calling) {
      #   newdata <- CustomDoCall(get.newdata, calling.args_newdata)
      # } else {
      #   newdata <- CustomDoCall(get.newdata, calling.args_newdata)
      # }
      rm('calling.args_newdata')
      calling.args$newdata <- newdata
    }
    
    
    # 6.03.2025
    if(!indirectcall & !growthparameters_calling) {
      if(is.null(newdata)) {
        calling.args_newdata         <- calling.args
        calling.args_newdata$dpar    <- dpar
        newdata <- CustomDoCall(get.newdata, calling.args_newdata)
        rm('calling.args_newdata')
        calling.args$newdata <- newdata
      }
    }
    
    # 6.03.2025
    if(is.null(full.args$newdata)) {
      full.args$newdata <- calling.args$newdata
    }
    
    # set up mesage
    if(check_fun) {
      if(!available_d1) {
        message_for_model_deriv_FALSE <- ""
        message_for_model_deriv_FALSE <- 
          paste0(message_for_model_deriv_FALSE, "\n",
                 "calculating deriv by differentiation of distance curve")
        
        if(is.null(ipts)) {
          message_for_model_deriv_FALSE <- 
            paste0(message_for_model_deriv_FALSE, "\n",
                   "It is strongly recommended not to set'ipts = NULL'", "\n",
                   "Typically, ipts > 100 is needed to get smooth deriv curve")
        }
        if(any(grepl("fitted_draws", rlang_trace_back[[1]])) |
           any(grepl("predict_draws", rlang_trace_back[[1]]))) {
          message_for_model_deriv_FALSE <- 
            paste0(message_for_model_deriv_FALSE, "\n",
                   "A better approach would be use 'marginal_draws()' ",
                   "instead")
        } # if(grepl("fitted_draws", rlang_trace_back[[1]]) |
      } # if(!available_d1) {
    } # if(check_fun) {
    
    
    calling.args <- 
      sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                 arguments = calling.args, 
                                 check_formalArgs = NULL,
                                 check_formalArgs_exceptions = c('object'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    
    # Interpolation points
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    calling.args$ipts <- ipts <- check_ipts(ipts = calling.args$ipts, 
                                            nipts = NULL, 
                                            check_fun  = check_fun, 
                                            available_d1 = available_d1, 
                                            xcall = NULL, verbose = verbose)
    
    
    if(!check_fun) {
      . <- CustomDoCall(predict, calling.args)
    }
    
    
    # TODO - check - Transformation applied to draw if deriv = 0, 
    # if(!check_fun) {
    #   if(deriv == 0) {
    #     if(!is.null(calling.args$model$model_info$transform_draws)) {
    #       . <- calling.args$model$model_info$transform_draws(.)
    #     } else if(!is.null(calling.args$object$model_info$transform_draws)) {
    #       . <- calling.args$object$model_info$transform_draws(.)
    #     } else {
    #       if(verbose) message("transform_draws function is NULL, check it")
    #     }
    #   } # if(deriv = 0) {
    # } # if(!check_fun) {
    
    
    # TODO - check - Transformation applied to draw if deriv = 0, 
    if(!check_fun) {
      if(deriv == 0) {
        if(!is.null(calling.args$model$
                    model_info$transform_draws)) {
          . <- calling.args$model$
            model_info$transform_draws(.)
        } else if(!is.null(calling.args$object$
                           model_info$transform_draws)) {
          . <- calling.args$object$
            model_info$transform_draws(.)
        } else {
          if(verbose) message("transform_draws function is NULL, check it")
        }
      } # if(deriv = 0) {
    } # if(!check_fun) {
    
    
    if(!plot_conditional_effects_calling) {
      if(check_fun) {
        if(deriv > 0) {
          if(available_d1) {
            . <- CustomDoCall(predict, calling.args)
          }
          if(!available_d1) {
            if(verbose) {
              message(message_for_model_deriv_FALSE)
            }
            calling.args_mapderivqr_args <- calling.args
            calling.args_mapderivqr_args[['summary']] <- FALSE
            y0 <- CustomDoCall(predict, calling.args_mapderivqr_args)
            # y0 <- calling.args$model$model_info$transform_draws(y0)
            if(!is.null(calling.args_mapderivqr_args$model$
                        model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$model$
                model_info$transform_draws(y0)
            } else if(!is.null(calling.args_mapderivqr_args$object$
                               model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$object$
                model_info$transform_draws(y0)
            } else {
              stop("transform_draws function is NULL, check it")
            }
            mapderivqr_args <- list()
            mapderivqr_args[['y0']] <- y0
            mapderivqr_args[['model']] <- calling.args[['object']]
            mapderivqr_args[['newdata']] <- calling.args[['newdata']]
            mapderivqr_args[['xvar']] <- calling.args[['xvar']]
            mapderivqr_args[['idvar']] <- calling.args[['idvar']]
            mapderivqr_args[['levels_id']] <- calling.args[['levels_id']]
            mapderivqr_args[['deriv']] <- calling.args[['deriv']]
            mapderivqr_args[['resp']] <- calling.args[['resp']]
            mapderivqr_args[['probs']] <- calling.args[['probs']]
            mapderivqr_args[['summary']] <- calling.args[['summary']]
            mapderivqr_args[['robust']] <- calling.args[['robust']]
            . <- CustomDoCall(mapderivqr, mapderivqr_args)
          }
        } # if(deriv > 0) {
      } # if(check_fun) {
    } # if(!plot_conditional_effects_calling) {
    
    
    full.args <-
      sanitize_CustomDoCall_args(what = "CustomDoCall",
                                 arguments = full.args,
                                 # check_formalArgs = plot_conditional_effects.bgmfit,
                                 check_formalArgs_exceptions = c('object'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    full.args$object <- full.args$model
    
    # Interpolation points
    if(!exists('check_fun')) check_fun <- FALSE
    if(!exists('available_d1')) available_d1 <- FALSE
    full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                         nipts = NULL, 
                                         check_fun  = check_fun, 
                                         available_d1 = available_d1, 
                                         xcall = NULL, verbose = verbose)
    
    
    if(plot_conditional_effects_calling) {
      if(check_fun) {
        if(deriv > 0) {
          if(available_d1) {
            . <- CustomDoCall(predict, full.args)
          }
          if(!available_d1) {
            if(verbose) {
              message(message_for_model_deriv_FALSE)
            }
            calling.args_mapderivqr_args <- full.args
            calling.args_mapderivqr_args[['summary']] <- FALSE
            y0 <- CustomDoCall(predict, calling.args_mapderivqr_args)
            # y0 <- calling.args$model$model_info$transform_draws(y0)
            if(!is.null(calling.args_mapderivqr_args$model$
                        model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$model$
                model_info$transform_draws(y0)
            } else if(!is.null(calling.args_mapderivqr_args$object$
                               model_info$transform_draws)) {
              y0 <- calling.args_mapderivqr_args$object$
                model_info$transform_draws(y0)
            } else {
              stop("transform_draws function is NULL, check it")
            }
            mapderivqr_args <- list()
            mapderivqr_args[['y0']] <- y0
            mapderivqr_args[['model']] <- calling.args[['object']]
            mapderivqr_args[['newdata']] <- calling.args[['newdata']]
            mapderivqr_args[['xvar']] <- calling.args[['xvar']]
            mapderivqr_args[['idvar']] <- calling.args[['idvar']]
            mapderivqr_args[['levels_id']] <- calling.args[['levels_id']]
            mapderivqr_args[['deriv']] <- calling.args[['deriv']]
            mapderivqr_args[['resp']] <- calling.args[['resp']]
            mapderivqr_args[['probs']] <- calling.args[['probs']]
            mapderivqr_args[['summary']] <- calling.args[['summary']]
            mapderivqr_args[['robust']] <- calling.args[['robust']]
            . <- CustomDoCall(mapderivqr, mapderivqr_args)
          }
        } # if(deriv > 0) {
      } # if(check_fun) {
    } # if(plot_conditional_effects_calling) {
    
    
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
    
    
    # fullframe
    full.args$idata_method <- idata_method
    full.args$fullframe <- eval(full.args$fullframe)
    
    # 6.03.2025 - extract info for if(is.na(model$model_info$univariate_by$by)){
    if(!is.null(eval(full.args$fullframe))) {
      if(eval(full.args$fullframe)) {
        if (is.null(resp)) {
          resp_rev_ <- resp
        } else if (!is.null(resp)) {
          resp_rev_ <- paste0("_", resp)
        }
        xvar_ <- paste0('xvar', resp_rev_)
        yvar_ <- paste0('yvar', resp_rev_)
        groupvar_ <- paste0('groupvar', resp_rev_)
        if(is.null(xvar)) xvar <- model$model_info[[xvar_]]
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
        # if (is.null(levels_id)) {
        #   IDvar <- model$model_info[[groupvar_]]
        #   if (!is.null(model$model_info[[hierarchical_]])) {
        #     IDvar <- model$model_info[[hierarchical_]]
        #   }
        # } else if (!is.null(levels_id)) {
        #   IDvar <- levels_id
        # }
        xvar  <- xvar
        idvar <- idvar
        if(length(idvar) > 1) idvar <- idvar[1]
        yvar  <- 'yvar'
      } # if(eval(full.args$fullframe)) {
    } # if(!is.null(eval(full.args$fullframe))) {
    
    
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
      if (!is.na(model$model_info$univariate_by$by)) {
        if(full.args$idata_method == 'm1') setfullframe <- FALSE
        if(full.args$idata_method == 'm2') setfullframe <- TRUE
      } else {
        setfullframe <- FALSE
      }
    }
    
    if (!is.na(model$model_info$univariate_by$by)) {
      if(is.null(full.args$fullframe)) {
        full.args$fullframe <- fullframe <- FALSE
      }
      if(full.args$fullframe & full.args$idata_method == 'm1') {
        setfullframe <- FALSE
      }
      if(full.args$fullframe & full.args$idata_method == 'm2') {
        setfullframe <- TRUE
      }
      if(!full.args$fullframe) {
        setfullframe <- FALSE
      }
      if(setfullframe) {
        uvarby <- model$model_info$univariate_by$by
        uvarbyresp <- paste0(uvarby, resp)
        uvarbynewdata <- eval(full.args$newdata) %>% 
          dplyr::filter(!!dplyr::sym(uvarbyresp) == 1)
        . <- cbind(., uvarbynewdata)
        
        # 6.03.2025
        # prepare_data2
        itransform_set <- get_itransform_call(itransform)
        if(any(itransform_set != "")) {
          . <- prepare_transformations(data = ., model = model, 
                                       itransform = itransform_set) 
        } # if(any(itransform_set != "")) {
      } # if(setfullframe) {
    } # if (!is.na(model$model_info$univariate_by$by)) {
    # 6.03.2025
    if (is.na(model$model_info$univariate_by$by)) {
      # Fot plot_curves() and growthparameters(), . must be be combined
      if(!is.null(eval(full.args$fullframe))) {
        if(eval(full.args$fullframe)) {
          cbindtonewdata <- eval(full.args$newdata)
          . <- cbind(cbindtonewdata, .)
        }
      }
      
      # prepare_data2
      itransform_set <- get_itransform_call(itransform)
      cbindtonewdata <- eval(full.args$newdata)
      if(any(itransform_set != "")) {
        . <- prepare_transformations(data = ., model = model, 
                                     itransform = itransform_set) 
      } # if(any(itransform_set != "")) {
    } # if (is.na(model$model_info$univariate_by$by)) {
    . 
  } # end predict_draws



#' @rdname predict_draws
#' @export
predict_draws <- function(model, ...) {
  UseMethod("predict_draws")
}


