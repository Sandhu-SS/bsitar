

#' Fitted (expected) values from the posterior predictive distribution
#' 
#' @details The \code{fitted_} function is a wrapper around 
#' the [brms::fitted.brmsfit]. The [brms::fitted.brmsfit]  
#' function from the \code{brms} package can used to plot the fitted (distance) 
#' curve for an *bsitar* model when outcome (e.g., height) is untransformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::fitted.brmsfit] will return the fitted curve on the log or 
#' square root scale whereas the [bsitar::fitted_] will 
#' return the fitted curve on the original scale. Furthermore, the 
#' fitted_ also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::fitted.brmsfit] and 
#' [bsitar::fitted_]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::fitted.brmsfit]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*fitted* to 
#' *fitted_*). 
#' 
#' @inherit conditional_effects_.bsitar params
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit] 
#' function. Please see \code{brms::fitted.brmsfit} for details on 
#' various options available.
#' 
#' @inherit brms::fitted.brmsfit description return
#' 
#' @export fitted_.bsitar
#' @export
#'
#' @examples
#' \dontrun{
#' # The following examples show the use of *fitted_* to estimate  
#' # population average and individual-specific distance and velocity 
#' # curves.
#' # Population average distance curve
#' fitted_(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' fitted_(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' fitted_(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' fitted_(model, deriv = 1, re_formula = NULL)
#'  
#' }

fitted_.bsitar <-
  function(model,
           resp = NULL,
           deriv = 0,
           envir = parent.frame(),
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = deriv)
    assign(o[[1]], model$Spl_funs[[o[[2]]]], envir = envir)
    . <- fitted(model, resp = resp, ...)
    assign(o[[1]], model$Spl_funs[[o[[1]]]], envir = envir)
    .
  }


#' @rdname fitted_.bsitar
#' @export
fitted_ <- function(model, ...) {
  UseMethod("fitted_")
}





