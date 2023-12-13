

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
#' @param deriv Must be \code{NULL}. 
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
           deriv = 0,
           usesavedfuns = FALSE,
           clearenvfuns = FALSE,
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- parent.frame()
    }
    
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv,
                             all = FALSE)
    
    oall <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv,
                             all = TRUE)
    
    getfunx1always <- model$model_info[['exefuns']][[o[[1]]]]
    
    if(deriv == 0) {
      getfunx <- model$model_info[['exefuns']][[o[[2]]]]
    } else if(deriv > 0) {
      stop("For loo_validation, the 'deriv' argument must be set as 0")
    }
    
    
    if(usesavedfuns) {
      setcleanup <- TRUE
      tempgenv <- parent.frame()
      oalli_c <- c()
      oalli_c <- c(oalli_c, paste0(o[[1]], "0"))
      for (oalli in names(oall)) {
        if(!grepl(o[[1]], oalli)) {
          oalli_c <- c(oalli_c, oalli)
        }
      }
      for (oalli in oalli_c) {
        assign(oalli, oall[[oalli]], envir = tempgenv)
      }
      assign(o[[1]], getfunx, envir = tempgenv)
    } else if(!usesavedfuns) {
      setcleanup <- FALSE
      if(!check_if_functions_exists(model, o, model$xcall)) {
        return(invisible(NULL))
      } else {
        setcleanup <- TRUE
        tempgenv <- parent.frame()
        if(exists(o[[1]], envir = tempgenv)) {
          assign(o[[1]], getfunx, envir = tempgenv)
        } else {
          assign(o[[1]], getfunx, envir = environment(getfunx))
        }
      }
    }
    
    
    
    . <- brms::loo(model, resp = resp, cores = cores ,...)
    
    assign(o[[1]], getfunx1always, environment(getfunx1always))
    
    if(!is.null(clearenvfuns)) {
      if(!is.logical(clearenvfuns)) {
        stop('clearenvfuns must be NULL or a logical')
      } else {
        setcleanup <- clearenvfuns
      }
    }
    
    if(setcleanup) {
      tempgenv <- parent.frame()
      for (oalli in names(oall)) {
        if(exists(oalli, envir = tempgenv )) {
          remove(list=oalli, envir = tempgenv)
        }
      }
    }
    
    .
  }


#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

