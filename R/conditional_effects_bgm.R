

#' Visualize the conditional effects of predictor
#' 
#' @details The \code{conditional_effects_bgm} function is a wrapper around 
#' the [brms::conditional_effects()]. The [brms::conditional_effects()]  
#' function from the \code{brms} package can used to plot the fitted (distance) 
#' curve for an *bgmfit* model when response (e.g., height) is not transformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::conditional_effects] will return the fitted curve on the log or 
#' square root scale whereas the [bsitar::conditional_effects_bgm()] will 
#' return the fitted curve on the original scale. Furthermore, the 
#' conditional_effects_bgm also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::conditional_effects] and 
#' [bsitar::conditional_effects_bgm()]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::conditional_effects]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*conditional_effects* to 
#' *conditional_effects_bgm*). 
#'
#' @param model An object of class \code{bgmfit}.
#' function.
#' @param resp Response variable (default \code{NULL}) specified as a string 
#' character required during the post-processing of multivariate and 
#' univariate-by-subgroup model (see \code{bsitar::bgm()} for details).
#' @param deriv An integer value to specify whether to estimate distance curve
#'  (i.e., model estimated curve(s)) or velocity curve (first derivative of the 
#'  model estimated curve(s)). A 0 value (default) is for distance curve and 
#'  1 for the velocity curve.
#' @param envir The calling environment. Deafault set to \code{parent.frame()}.
#' @param ... Additional arguments passed to the [brms::conditional_effects()] 
#' function. Please see [brms::conditional_effects()] for details.
#'
#' @inherit brms::conditional_effects description 
#' 
#' @return An object of class 'brms_conditional_effects' which is a named list
#'   with one data.frame per effect containing all information required to
#'   generate conditional effects plots. See brms::conditional_effects for details.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # The examples below show the use of *conditional_effects_bgm* to estimate  
#' # the population average and individual-specific distance and velocity 
#' # curves.
#' # Population average distance curve
#' conditional_effects_bgm(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' conditional_effects_bgm(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' conditional_effects_bgm(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' conditional_effects_bgm(model, deriv = 1, re_formula = NULL)
#'  
#' }

conditional_effects_bgm.bgmfit <-
  function(model,
           resp = NULL,
           deriv = 0,
           envir = parent.frame(),
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = deriv,
                             envir = envir)
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    . <- brms::conditional_effects(model, resp = resp, ...)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname conditional_effects_bgm.bgmfit
#' @export
conditional_effects_bgm <- function(model, ...) {
  UseMethod("conditional_effects_bgm")
}





