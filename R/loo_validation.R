

#' Perform leave-one-out (loo) cross-validation
#' 
#' @details The \strong{loo_validation} function is a wrapper around 
#' the [brms::loo()] and works exactly the same way as [brms::loo()]. 
#' 
#' @inherit brms::loo description 
#' 
#' @inheritParams growthparameters.bgmfit
#' 
#' @inheritParams pp_check_bgm.bgmfit
#' 
#' @param deriv Must be \code{NULL}. 
#' 
#' @param ... Additional arguments passed to the [brms::loo()] 
#' function. Please see \code{brms::loo} for details on 
#' various options available.
#' 
#' @return If only one model object is provided, then an object of class
#'   \code{loo} is returned. If multiple objects are provided, an object of
#'   class \code{loolist}.
#' 
#' @export loo_validation.bgmfit
#' @export
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' # To avoid running the model which takes some time, model fit to the
#' # \code{berkeley_mdata} has already been saved as berkeley_mfit.rda object.
#' # Please see \code{bgm} examples.
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
           envir = globalenv(),
           deriv = NULL,
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv)
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    . <- brms::loo(model, resp = resp, cores = cores ,...)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

