

#' Leave-one-out (loo) cross-validation for \code{bgmfit} object
#' 
#' @details The \code{loo_bgm} function is a wrapper around 
#' the [brms::loo] and works exactly the same way as 
#' [brms::loo]. 
#' 
#' @inherit pp_check_bgm.bgmfit params
#' 
#' @inherit brms::loo description 
#' 
#' @inheritParams growthparameters.bgmfit
#' 
#' @param deriv Must be \code{NULL}. 
#' 
#' @param ... Additional arguments passed to the [brms::loo] 
#' function. Please see \code{brms::loo} for details on 
#' various options available.
#' 
#' @return If just one object is provided, an object of class \code{loo}. If 
#' multiple objects are provided, an object of class \code{loolist}.
#' 
#' @export loo_bgm.bgmfit
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
#' loo_bgm(model, cores = 1)
#' }
#' 
#' 
loo_bgm.bgmfit <-
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


#' @rdname loo_bgm.bgmfit
#' @export
loo_bgm <- function(model, ...) {
  UseMethod("loo_bgm")
}

