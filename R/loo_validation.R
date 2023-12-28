

#' Perform leave-one-out (loo) cross-validation
#' 
#' @description The \strong{loo_validation()} is a wrapper around the
#'   [brms::loo()] function to perform approximate leave-one-out
#'   cross-validation based on the posterior likelihood. See [brms::loo()] for
#'   details.
#' 
#' @inherit brms::loo details 
#' 
#' @inheritParams growthparameters.bgmfit
#' 
#' @inheritParams plot_ppc.bgmfit
#' 
#' @param ... Additional arguments passed to the [brms::loo()] function. 
#' Please see \code{brms::loo} for details on various options available.
#' 
#' @return If only one model object is provided, then an object of class
#'   \code{loo} is returned. If multiple objects are provided, an object of
#'   class \code{loolist}.
#' 
#' @export loo_validation.bgmfit
#' @export
#' 
#' @seealso [brms::loo()] 
#' 
#' @inherit berkeley author
#'
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid fitting the model which takes time, the model  
#' # fit has already been saved as 'berkeley_mfit.rda' file.
#' # See examples section of the main function for details on the model fit.
#' 
#' model <- berkeley_mfit
#' 
#' \donttest{
#' loo_validation(model, cores = 1)
#' }
#' 
#' 
loo_validation.bgmfit <-
  function(model,
           resp = NULL,
           cores = 1,
           deriv = NULL,
           verbose = FALSE,
           usesavedfuns = FALSE,
           clearenvfuns = NULL,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- parent.frame()
    }
    
    if(is.null(deriv)) {
      deriv <- 0
    }
    
    full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                    fargs = formals(), 
                                    dargs = list(...), 
                                    verbose = verbose)
    
    
    if(!is.null(model$xcall)) {
      arguments <- get_args_(as.list(match.call())[-1], model$xcall)
      newdata <- newdata
    } else {
      newdata <- do.call(get.newdata, full.args)
    }
    
    
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
    
    # setupfuns.args <- full.args
    # setupfuns.args$o <- o
    # setupfuns.args$oall <- oall
    # setupfuns.args$deriv <- deriv
    # test <- do.call(setupfuns, setupfuns.args)
   
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall,
                      usesavedfuns = usesavedfuns,
                      deriv = deriv, envir = envir,
                      deriv_model = deriv_model,
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
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
    
    
    # . <- do.call(brms::loo, calling.args)
    
    
    . <- brms::loo(model, resp = resp, cores = cores ,...)
    
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
      if(eval(full.args$usesavedfuns)) {
        setcleanup <- TRUE 
      } else {
        setcleanup <- FALSE
      }
    }
    
    # Cleanup environment if requested
    if(setcleanup) {
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
    } # if(setcleanup) {
    
    .
  }


#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

