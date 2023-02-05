

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
#' @inherit gparameters.bsitar params
#' 
#' @inherit conditional_effects_.bsitar params
#' 
#' @inherit brms::predict.brmsfit description params
#' 
#' @param ... Additional arguments passed to the [brms::predict.brmsfit] 
#' function. Please see \code{brms::predict.brmsfit} for details on 
#' various options available.
#' 
#' @return An array of predicted response values. See [brms::predict.brmsfit] 
#' for details.
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
           newdata = NULL,
           resp = NULL,
           ndraws = NULL,
           re_formula = NULL,
           numeric_cov_at = NULL,
           levels_id = NULL,
           ipts = NULL,
           deriv = 0,
           summary = TRUE,
           robust = FALSE,
           probs = c(0.025, 0.975),
           envir = parent.frame(),
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = deriv)
    
    xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
    
    if(xcall == "gparameters" | xcall == "plot_bsitar") {
      arguments <- get_args_(as.list(match.call())[-1], xcall)
      newdata <- newdata
    } else {
      newdata <- get.newdata(model, 
                             newdata = newdata, 
                             resp = resp, 
                             numeric_cov_at = numeric_cov_at,
                             levels_id = levels_id,
                             ipts = ipts)
    }
    
    assign(o[[1]], model$Spl_funs[[o[[2]]]], envir = envir)
    
    . <- predict(model,
                newdata = newdata,
                resp = resp,
                ndraws = ndraws,
                re_formula = re_formula,
                numeric_cov_at = numeric_cov_at,
                levels_id = levels_id,
                ipts = ipts,
                deriv = 0,
                summary = TRUE,
                robust = FALSE,
                probs = c(0.025, 0.975), 
                ...)
    assign(o[[1]], model$Spl_funs[[o[[1]]]], envir = envir)
    .
  }


#' @rdname predict_.bsitar
#' @export
predict_ <- function(model, ...) {
  UseMethod("predict_")
}


