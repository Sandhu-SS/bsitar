

#' Perform posterior predictive distribution check
#' 
#' @details The \strong{plot_ppc} function is a wrapper around 
#' the [brms::pp_check()]. 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @inherit plot_conditional_effects.bgmfit params
#' 
#' @param ... Additional arguments passed to the [brms::pp_check.brmsfit()] 
#' function. Please see [brms::pp_check.brmsfit()] for details.
#' 
#' @inherit brms::pp_check.brmsfit description 
#' 
#' @return A ggplot object that can be further customized using the
#' ggplot2 package.
#' 
#' @export
#' 
#' @inherit bgm author
#'
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid fitting the model which takes time, the model  
#' # fit has already been saved as 'berkeley_mfit.rda' file.
#' # See examples section of the bgm function for details on the model fit.
#' 
#' model <- berkeley_mfit
#' 
#' plot_ppc(model)
#' 
plot_ppc.bgmfit <-
  function(model,
           resp = NULL,
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
                             deriv = deriv,
                             all = FALSE)

    
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

    . <- brms::pp_check(model, resp = resp, ...)
    
    assign(o[[1]], getfunx1always, environment(getfunx1always))
    .
  }


#' @rdname plot_ppc.bgmfit
#' @export
plot_ppc <- function(model, ...) {
  UseMethod("plot_ppc")
}


