

#' Leave-one-out (loo) cross-validation for \code{bsitar} object
#' 
#' @details The \code{loo_} function is a wrapper around 
#' the [brms::loo] and works exactly the same way as 
#' [bsitar::loo]. 
#' 
#' @inherit pp_check_.bsitar params
#' 
#' @inherit brms::loo description return
#' 
#' @param ... Additional arguments passed to the [brms::loo] 
#' function. Please see \code{brms::loo} for details on 
#' various options available.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' loo_(model)
#' }

loo_.bsitar <-
  function(model,
           resp = NULL,
           envir = parent.frame(),
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = deriv)
    assign(o[[1]], model$Spl_funs[[o[[2]]]], envir = envir)
    . <- loo(model, resp = resp, ...)
    assign(o[[1]], model$Spl_funs[[o[[1]]]], envir = envir)
    .
  }


#' @rdname loo_.bsitar
#' @export
loo_ <- function(model, ...) {
  UseMethod("loo_")
}

