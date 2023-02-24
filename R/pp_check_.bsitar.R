

#' Posterior predictive distribution check for \code{bsitar} object
#' 
#' @details The \code{pp_check_} function is a wrapper around 
#' the [brms::loo] and works exactly the same way as 
#' [brms::pp_check.brmsfit]. 
#' 
#' @param model An object of class \code{bsitar}.
#' function.
#' @param resp Response variable (default \code{NULL}) specified as a string 
#' character required during the post-processing of multivariate and 
#' univariate-by-subgroup model (see \code{bsitar} function for details).
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
#' pp_check_(model)
#' }

pp_check_.bsitar <-
  function(model,
           resp = NULL,
           envir = .GlobalEnv,
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = deriv)
    assign(o[[1]], model$Spl_funs[[o[[2]]]], envir = envir)
    . <- pp_check(model, resp = resp, ...)
    assign(o[[1]], model$Spl_funs[[o[[1]]]], envir = envir)
    .
  }


#' @rdname pp_check_.bsitar
#' @export
pp_check_ <- function(model, ...) {
  UseMethod("pp_check_")
}


