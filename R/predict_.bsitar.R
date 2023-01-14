

#' Predicted values (draws) from the posterior predictive distribution
#' 
#' @details The \code{predict_} function is a wrapper around 
#' the [brms::predict.brmsfit]. The [brms::predict.brmsfit]  
#' function from the \code{brms} package can used to plot the predict (distance) 
#' curve for an *bsitar* model when outcome (e.g., height) is untransformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::predict.brmsfit] will return the predict curve on the log or 
#' square root scale whereas the [bsitar::predict_] will 
#' return the predict curve on the original scale. Furthermore, the 
#' predict_ also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::predict.brmsfit] and 
#' [bsitar::predict_]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::predict.brmsfit]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*predict* to 
#' *predict_*). 
#' 
#' @inherit conditional_effects_.bsitar params
#' 
#' @param ... Additional arguments passed to the [brms::predict.brmsfit] 
#' function. Please see \code{brms::predict.brmsfit} for details on 
#' various options available.
#' 
#' @inherit brms::predict.brmsfit description return
#' 
#' @export predict_.bsitar
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # The examples below show the use of *predict_* to estimate  
#' # population average and individual-specific distance and velocity 
#' # curves for the the predict model
#' # Population average distance curve
#' predict_(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' predict_(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' predict_(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' predict_(model, deriv = 1, re_formula = NULL)
#'  
#' }

predict_.bsitar <-
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
    . <- predict(model, resp = resp, ...)
    assign(o[[1]], model$Spl_funs[[o[[1]]]], envir = envir)
    .
  }


#' @rdname predict_.bsitar
#' @export
predict_ <- function(model, ...) {
  UseMethod("predict_")
}


