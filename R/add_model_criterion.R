

#' @title Add Model Fit Criteria to Model
#' 
#' @description The \strong{add_model_criterion()} function is a wrapper around
#'   [brms::add_criterion()] that allows adding fit criteria to a model. Note
#'   that arguments such as \code{compare} and \code{pointwise} are relevant
#'   only for [brms::add_loo], while \code{summary}, \code{robust}, and
#'   \code{probs} are ignored except for the [brms::bayes_R2()].
#' 
#' @param model An object of class \code{bgmfit} representing the model to which
#'   the fit criteria will be added.
#' 
#' @inheritParams growthparameters.bgmfit
#' @inherit brms::add_criterion.brmsfit params description return
#' @inheritParams brms::bayes_R2.brmsfit
#' @inheritParams brms::waic.brmsfit
#' @inheritParams fitted_draws.bgmfit
#'  
#' @return An object of class \code{bgmfit} with the specified fit criteria added.
#' 
#' @rdname add_model_criterion
#' @export
#' 
#' @seealso [brms::add_loo], [brms::add_ic()], [brms::add_waic()],
#'   [brms::bayes_R2()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid model estimation which can take time, the Bayesian SITAR model fit
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' model <- berkeley_exfit
#' 
#' # Add model fit criteria (e.g., WAIC)
#' model <- add_model_criterion(model, criterion = c("waic"))
#' }
#' 
add_model_criterion.bgmfit <-
  function(model,
           criterion = c("loo", "waic"),
           ndraws = NULL,
           draw_ids = NULL,
           compare = TRUE,
           pointwise = FALSE,
           model_names = NULL,
           summary = TRUE,
           robust = FALSE,
           probs = c(0.025, 0.975),
           newdata = NULL,
           resp = NULL,
           cores = 1,
           model_deriv = NULL,
           verbose = FALSE,
           expose_function = FALSE,
           usesavedfuns = NULL,
           clearenvfuns = NULL,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      if(!is.null(model$model_info$exefuns[[1]])) {
        envir <- environment(model$model_info$exefuns[[1]])
      } else {
        envir <- envir
      }
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
    
    if(is.null(model_deriv)) {
      model_deriv <- TRUE
    }
    
    full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                           fargs = formals(), 
                                           dargs = list(...), 
                                           verbose = verbose)
    
    full.args$model <- model
    full.args$deriv <- deriv <- 0
    
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
    
    
    if(!is.null(model$xcall)) {
      arguments <- get_args_(as.list(match.call())[-1], model$xcall)
      newdata <- newdata
    } else {
      newdata <- CustomDoCall(get.newdata, full.args)
    }
    
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
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall, 
                      usesavedfuns = usesavedfuns, 
                      deriv = deriv, envir = envir, 
                      model_deriv = model_deriv, 
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
    if(!isTRUE(
      check_pkg_version_exists('brms', minversion = '2.20.17', 
                               prompt = FALSE,
                               stop = FALSE,
                               verbose = FALSE))) {
      if(is.null(check_if_functions_exists(model, o, model$xcall,
                                           usesavedfuns = usesavedfuns))) {
        return(invisible(NULL))
      }
    }
    
    
    misc <- c("verbose", "usesavedfuns", "clearenvfuns", 
              "envir", "fullframe", "dummy_to_factor")
    calling.args <- post_processing_args_sanitize(model = model,
                                          xcall = match.call(),
                                          resp = resp,
                                          envir = envir,
                                          deriv = deriv, 
                                          dots = list(...),
                                          misc = misc,
                                          verbose = verbose)
   
    
    calling.args$object    <- full.args$model
    calling.args$x         <- full.args$model
    calling.args$object    <- full.args$model <- NULL
    calling.args$criterion <- criterion

    calling.args <- 
      sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                 arguments = calling.args, 
                                 check_formalArgs = add_model_criterion.bgmfit,
                                 check_formalArgs_exceptions = c('x'),
                                 check_trace_back = NULL,
                                 envir = parent.frame())
    
    suppressWarnings({
      . <- CustomDoCall(brms::add_criterion, calling.args)
    })
    
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


#' @rdname add_model_criterion
#' @export
add_model_criterion <- function(model, ...) {
  UseMethod("add_model_criterion")
}


