

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
                             deriv = 0,
                             all = FALSE)

    
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    
    
    if(!check_if_functions_exists(model, o)) return(invisible(NULL))

    . <- brms::pp_check(model, resp = resp, ...)
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname plot_ppc.bgmfit
#' @export
plot_ppc <- function(model, ...) {
  UseMethod("plot_ppc")
}


