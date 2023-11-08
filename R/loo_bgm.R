

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
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @export
#'
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' # data <- berkeley
#' # berkeley_fit <- bgm(x = age, y = height, id = id, data = data, df = 4,
#' #                     chains = 2, iter = 1000, thin = 10)
#' 
#' # To avoid running the model which takes some time, the fitted model has 
#' # already been saved as berkeley_fit.rda object. The model is fitted using 2 
#' # chain  with 1000  iteration per chain (to save time) and setting thin as 1 
#' # (to save memory also).
#' 
#' model <- berkeley_mfit
#' 
#' \donttest{
#' loo_bgm(model, cores = 1)
#' }
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

