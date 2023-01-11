

#' Visualize the conditional effects of predictor
#' 
#' @details The \code{conditional_effects_} function is a wrapper around 
#' the [brms::conditional_effects]. The [brms::conditional_effects]  
#' function from the \code{brms} package can used to plot the fitted (distance) 
#' curve for an *bsitar* model when outcome (e.g., height) is untransformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::conditional_effects] will return the fitted curve on the log or 
#' square root scale whereas the [bsitar::conditional_effects_] will 
#' return the fitted curve on the original scale. Furthermore, the 
#' conditional_effects_ also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::conditional_effects] and 
#' [bsitar::conditional_effects_]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::conditional_effects]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*conditional_effects* to 
#' *conditional_effects_*). 
#'
#' @param model An object of class \code{bsitar}.
#' function.
#' @param resp Response variable (default \code{NULL}) specified as a string 
#' character required during the post-processing of multivariate and 
#' univariate-by-subgroup model (see \code{bsitar} function for details).
#' @param deriv An integer value to specify whether to estimate distance curve
#'  (i.e., model estimated curve(s)) or velocity curve (first derivative of the 
#'  model estimated curve(s)). A 0 value (default) is for distance curve and 
#'  1 for the velocity curve.
#' @param envir The calling environment. Deafault set to \code{parent.frame()}.
#' @param ... Additional arguments passed to the [brms::conditional_effects()] 
#' function. Please see [brms::conditional_effects()] for details.
#'
#' @inherit brms::conditional_effects description return
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # The examples below show the use of *conditional_effects_* to estimate  
#' # the population average and individual-specific distance and velocity 
#' # curves.
#' # Population average distance curve
#' conditional_effects_(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' conditional_effects_(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' conditional_effects_(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' conditional_effects_(model, deriv = 1, re_formula = NULL)
#'  
#' }

conditional_effects_.bsitar <-
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
    . <- conditional_effects(model, resp = resp, ...)
    assign(o[[1]], model$Spl_funs[[o[[1]]]], envir = envir)
    .
  }


#' @rdname conditional_effects_.bsitar
#' @export
conditional_effects_ <- function(model, ...) {
  UseMethod("conditional_effects_")
}





