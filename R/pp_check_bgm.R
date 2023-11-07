

#' Posterior predictive distribution check for \code{bgmfit} object
#' 
#' @details The \code{pp_check_bgm} function is a wrapper around 
#' the [brms::loo] and works exactly the same way as 
#' [brms::pp_check.brmsfit]. 
#' 
#' @param model An object of class \code{bgmfit}.
#' function.
#' @param resp Response variable (default \code{NULL}) specified as a string 
#' character required during the post-processing of multivariate and 
#' univariate-by-subgroup model (see \code{bgmfit} function for details).
#' @param envir The calling environment. Deafault set to \code{parent.frame()}.
#' @param ... Additional arguments passed to the [brms::pp_check.brmsfit] 
#' function. Please see \code{brms::pp_check.brmsfit} for details on 
#' various options available.
#' 
#' @inherit brms::pp_check.brmsfit description 
#' 
#' @return A ggplot object that can be further customized using the
#' ggplot2 package.
#' 
#' @export
#'
#' @examples
#' #
#' # Fit Bayesian SITAR model 
#' # berkeley_fit <- bgm(x = age, y = height, id = id, data = data, df = 4,
#' #                     chains = 2, iter = 1000, thin = 10)
#' #
#' # To avoid running the model which takes some time, the fitted model has 
#' # already been saved as berkeley_fit.rda object. The model is fitted using 2 
#' # chain  with 1000  iteration per chain (to save time) and setting thin as 1 
#' # (to save memory also).
#' # 
#' model <- berkeley_fit
#' #
#' pp_check_bgm(model)
#' 

pp_check_bgm.bgmfit <-
  function(model,
           resp = NULL,
           envir = parent.frame(),
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = 0)
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    . <- brms::pp_check(model, resp = resp, ...)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname pp_check_bgm.bgmfit
#' @export
pp_check_bgm <- function(model, ...) {
  UseMethod("pp_check_bgm")
}


