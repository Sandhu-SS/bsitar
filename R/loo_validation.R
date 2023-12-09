

#' Perform leave-one-out (loo) cross-validation
#' 
#' @description The \strong{loo_validation} function is a wrapper around the
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
                             deriv = deriv)
    
    getfunx1always <- model$model_info[['exefuns']][[o[[1]]]]
    
    if(deriv == 0) {
      getfunx <- model$model_info[['exefuns']][[o[[2]]]]
    } else if(deriv > 0) {
      stop("For loo_validation, the 'deriv' argument must be set as 0")
    }
    
    setcleanup <- FALSE
    if(!check_if_functions_exists(model, o, model$xcall)) {
      return(invisible(NULL))
    } else {
      setcleanup <- TRUE
      assign(o[[1]], getfunx, envir = environment(getfunx))
    }
    
    . <- brms::loo(model, resp = resp, cores = cores ,...)
    assign(o[[1]], getfunx1always, environment(getfunx1always))
    .
  }


#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

