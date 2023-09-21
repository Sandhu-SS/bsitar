

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
#' 
#' \dontrun{
#' pp_check_bgm(model)
#' }

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


