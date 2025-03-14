


#' @title Fitted (Expected) Values from the Posterior Draws
#' 
#' @description The \strong{fitted_draws()} function is a wrapper around the
#'   [brms::fitted.brmsfit()] function, which allows users to obtain fitted
#'   values (and their summaries) from the posterior draws. For more details,
#'   refer to the documentation for [brms::fitted.brmsfit()]. An alternative
#'   approach is to [marginal_draws()] function which is based on the 
#'   \pkg{marginaleffects}.
#' 
#' @details The \strong{fitted_draws()} function computes the fitted values 
#'   from the posterior draws. While the [brms::fitted.brmsfit()] function 
#'   from the \pkg{brms} package can be used to obtain fitted (distance) 
#'   values when the outcome (e.g., height) is untransformed, it returns 
#'   fitted values on the log or square root scale if the outcome is 
#'   transformed. In contrast, \strong{fitted_draws()} returns fitted values 
#'   on the original scale. Additionally, \strong{fitted_draws()} computes 
#'   the first derivative (velocity) on the original scale, after applying 
#'   the necessary back-transformation. Apart from these differences, both 
#'   functions—[brms::fitted.brmsfit()] and [fitted_draws()]—operate in the 
#'   same manner, allowing users to specify all options available in 
#'   [brms::fitted.brmsfit()].
#'   
#' @param deriv An integer indicating whether to estimate the distance curve 
#'   or its derivative (velocity curve). The default \code{deriv = 0} is for 
#'   the distance curve, while \code{deriv = 1} is for the velocity curve.
#' 
#' @inherit growthparameters.bgmfit params
#' @inherit brms::fitted.brmsfit params
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit()] 
#'   function. For details on available options, please refer to 
#'   \code{brms::fitted.brmsfit()}.
#' 
#' @return An array of predicted mean response values when \code{summarise =
#'   FALSE}, or a \code{data.frame} when \code{summarise = TRUE}. For further
#'   details, refer to [brms::fitted.brmsfit].
#' 
#' @rdname fitted_draws
#' @export
#' 
#' @seealso [brms::fitted.brmsfit()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid time-consuming model estimation, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See the 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check if the model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Population average distance curve
#' fitted_draws(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' fitted_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' fitted_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' fitted_draws(model, deriv = 1, re_formula = NULL)
#' }
#' 
fitted_draws.bgmfit <-
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
           ipts = 10,
           deriv = 0,
           deriv_model = TRUE,
           summary = TRUE,
           robust = FALSE,
           transform = NULL,
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
           itransform = NULL,
           newdata_fixed = FALSE,
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
    
    # 6.03.2025
    
    
    
    check_CustomDoCall <- gsub_space(paste(deparse(sys.calls()), collapse = ""))
    eval_CustomDoCall <- FALSE
    if(grepl("CustomDoCall\\(", check_CustomDoCall)) {
      eval_CustomDoCall <- TRUE
      check_CustomDoCall <- regmatches(check_CustomDoCall, gregexpr("(?<=\\().*?(?=\\))", check_CustomDoCall, perl=T))[[1]]
      check_CustomDoCall <- strsplit(check_CustomDoCall[1], "\\(")[[1]][1]
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
    
   
    plot_conditional_effects_calling <- FALSE
    syscalls1 <- sys.calls()[[1]]
    syscallsall <- paste(deparse(syscalls1), collapse = "\n")
    for (xc in 1:length(syscallsall)) {
      if(any(grepl('plot_conditional_effects', syscallsall[[xc]]))) {
        plot_conditional_effects_calling <- TRUE
      }
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
        # 6.03.2025 - even when indirectcall = FALSE, need to create new data
        # But that must be based on the final calling.args
        # newdata <- do.call(get.newdata, full.args)
      }
      # 6.03.2025 - even when indirectcall = FALSE, need to create new data
      # But that must be based on the final calling.args
       # full.args$newdata <- newdata
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
      if(model$model_info$decomp == "QR") deriv_model<- FALSE
    }
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
    
    
    # 6.03.2025
    # xcall_str used for if(plot_conditional_effects_calling) {... 
    if(is.null(xcall_str)) {
      # setxcallarg <- paste(deparse(match.call()), collapse = "")
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
    
    o    <- do.call(post_processing_checks, post_processing_checks_args)

    post_processing_checks_args[['all']]      <- TRUE
    oall <- do.call(post_processing_checks, post_processing_checks_args)
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
        deriv_model <- FALSE
        call_slopes <- TRUE
        post_processing_checks_args[['deriv']]    <- 0
        o    <- do.call(post_processing_checks, post_processing_checks_args)
      }
      check_fun <- TRUE
    }
    

    # Unlike marginal_... functions where assign() works, for brms::fitted...
    # the environment is assigned to d0/d1 functions via setupfuns()
    # The deriv = 0/1 should also reflect in  setupfuns() 
   
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall,
                      usesavedfuns = usesavedfuns,
                      deriv = post_processing_checks_args[['deriv']],
                      envir = envir,
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

    

    growthparameters_calling <- FALSE
    syscalls1 <- sys.calls()[[1]]
    syscallsall <- paste(deparse(syscalls1), collapse = "\n")
    for (xc in 1:length(syscallsall)) {
      if(any(grepl('growthparameters', syscallsall[[xc]]))) {
        growthparameters_calling <- TRUE
        # 6.03.2025 - this "_growthparameters" to avoind false flag
        # modelbased_growthparameters
         if(any(grepl('_growthparameters', syscallsall[[xc]]))) {
          growthparameters_calling <- FALSE
         }
      }
    }
    if(growthparameters_calling) {
      if(calling.args$re_formula_opt == "V") {
        calling.args$re_formula <- NULL
      } else if(calling.args$re_formula_opt == "v") {
        calling.args$re_formula <- NA
      }
    }
    
    # 6.03.2025
    if(is.null(newdata)) {
      calling.args_newdata         <- calling.args
      calling.args_newdata$model   <- calling.args_newdata$object
      calling.args_newdata$newdata <- model$model_info$bgmfit.data 
      newdata <- do.call(get.newdata, calling.args_newdata)
      rm('calling.args_newdata')
      calling.args$newdata <- newdata
    }
    # 6.03.2025
    if(is.null(full.args$newdata)) {
      full.args$newdata <- calling.args$newdata
    }
    

    # set up mesage
    if(check_fun) {
      if(!available_d1) {
        message_for_deriv_model_FALSE <- ""
        message_for_deriv_model_FALSE <- 
          paste0(message_for_deriv_model_FALSE, "\n",
                 "calculating deriv by differentiation of distance curve")
        
        if(is.null(ipts)) {
          message_for_deriv_model_FALSE <- 
            paste0(message_for_deriv_model_FALSE, "\n",
                   "It is strongly recommended not to set'ipts = NULL'", "\n",
                   "Typically, ipts > 100 is needed to get smooth deriv curve")
        }
        
        for (xc in 1:length(syscallsall)) {
          syscallsall_names_ <- sub(" *\\(.*", "", syscallsall[[xc]])
          if(grepl('fitted_draws', syscallsall_names_) |
             grepl('predict_draws', syscallsall_names_)) {
            message_for_deriv_model_FALSE <- 
              paste0(message_for_deriv_model_FALSE, "\n",
                     "A better approach would be use 'marginal_draws()' ",
                     "instead of '", syscallsall_names_, "()'")
          }
        }
      }
    } # if(check_fun) {
    
    
   
    
    # 6.03.2025
    if(eval_CustomDoCall) {
        # 6.03.2025 - remove missing GOOD
      ownargs <- formalArgs(fitted_draws.bgmfit)
      for (i in setdiff(names(calling.args), ownargs)) {
        calling.args[[i]] <- NULL
      }
      for (i in names(calling.args)) {
        if(is.symbol(calling.args[[i]])) {
          if(deparse(calling.args[[i]]) == "") {
            calling.args[[i]] <- NULL
          }
        }
      }
      for (i in names(calling.args)) {
        calling.args[[i]] <- eval(calling.args[[i]], 
                                  envir = parent.frame())
        
      }
      calling.args$object <- calling.args$model
    } # if(eval_CustomDoCall) {
  
    

    if(!check_fun) {
      . <- do.call(fitted, calling.args)
    }
    if(check_fun) {
      if(deriv > 0) {
        if(available_d1) {
          . <- do.call(fitted, calling.args)
        } 
        if(!available_d1) {
          if(verbose) {
           message(message_for_deriv_model_FALSE)
          }
          calling.args_mapderivqr_args <- calling.args
          calling.args_mapderivqr_args[['summary']] <- FALSE
          y0 <- do.call(fitted, calling.args_mapderivqr_args)
          mapderivqr_args <- list()
          mapderivqr_args[['y0']] <- y0
          mapderivqr_args[['model']] <- calling.args[['object']]
          mapderivqr_args[['newdata']] <- calling.args[['newdata']]
          mapderivqr_args[['deriv']] <- calling.args[['deriv']]
          mapderivqr_args[['resp']] <- calling.args[['resp']]
          mapderivqr_args[['probs']] <- calling.args[['probs']]
          mapderivqr_args[['summary']] <- calling.args[['summary']]
          mapderivqr_args[['robust']] <- calling.args[['robust']]
          . <- do.call(mapderivqr, mapderivqr_args)
        }
      } # if(deriv > 0) {
    } # if(check_fun) {

    
    
    # if(!is.null((eval(full.args$deriv_model)))) {
    #   full.args$deriv_model <- FALSE
    # }

    if(!is.null((eval(full.args$deriv)))) {
      if(eval(full.args$deriv) > 0) { 
        # if(!is.null((eval(full.args$deriv_model)))) { # new layer if 6.03.2025
          if(!eval(full.args$deriv_model)) {
            full.args$. <- .
            . <- do.call(mapderivqr, full.args)
          } else {
            . <- .
          # }
        }
      }
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
        # why this if(setfullframe)? we are already in if(setfullframe)
        # if(setfullframe) . <- cbind(., uvarbynewdata)
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
  } # end fitted_draws


#' @rdname fitted_draws
#' @export
fitted_draws <- function(model, ...) {
  UseMethod("fitted_draws")
}


