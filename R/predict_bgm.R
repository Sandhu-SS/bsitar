

#' Predicted values (draws) from the posterior predictive distribution
#' 
#' @details The \code{predict_bgm} function is a wrapper around 
#' the [brms::predict.brmsfit]. The [brms::predict.brmsfit]  
#' function from the \code{brms} package can used to plot the predict (distance) 
#' curve for an *bsitar* model when outcome (e.g., height) is untransformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::predict.brmsfit] will return the predict curve on the log or 
#' square root scale whereas the [bsitar::predict_bgm()] will 
#' return the predict curve on the original scale. Furthermore, the 
#' predict_bgm also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::predict.brmsfit] and 
#' [bsitar::predict_bgm()]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::predict.brmsfit]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*predict* to 
#' *predict_bgm*). 
#' 
#' @inherit growthparameters.bgmfit params
#' 
#' @inherit conditional_effects_bgm.bgmfit params
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
#' @export predict_bgm.bgmfit
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # The examples below show the use of *predict_bgm* to estimate  
#' # population average and individual-specific distance and velocity 
#' # curves for the the predict model
#' # Population average distance curve
#' predict_bgm(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' predict_bgm(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' predict_bgm(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' predict_bgm(model, deriv = 1, re_formula = NULL)
#'  
#' }

predict_bgm.bgmfit <-
  function(model,
           newdata = NULL,
           resp = NULL,
           ndraws = NULL,
           re_formula = NULL,
           numeric_cov_at = NULL,
           levels_id = NULL,
           ipts = NULL,
           deriv = 0,
           deriv_model = TRUE,
           summary = TRUE,
           robust = FALSE,
           probs = c(0.025, 0.975),
           xrange = NULL,
           parms_eval = FALSE,
           parms_method = 'getPeak',
           idata_method = 'm1',
           envir = parent.frame(),
           ...) {
    
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             deriv = deriv, 
                             envir = envir)
    
    if(!is.null(model$xcall)) {
      arguments <- get_args_(as.list(match.call())[-1], model$xcall)
      newdata <- newdata
    } else {
      newdata <- get.newdata(model, 
                             newdata = newdata, 
                             resp = resp, 
                             numeric_cov_at = numeric_cov_at,
                             levels_id = levels_id,
                             ipts = ipts,
                             xrange = xrange,
                             idata_method = idata_method)
    }
    
    if(deriv == 0) {
      assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
    }
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") deriv_model<- FALSE
    }
    
    if(!deriv_model) {
      if(deriv == 1 | deriv == 2) {
        summary <- FALSE
        assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
      }
    }
    
    if(deriv_model) {
      if(deriv == 1 | deriv == 2) {
        assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
      }
    }
    
    . <- predict(model,
                newdata = newdata,
                resp = resp,
                ndraws = ndraws,
                re_formula = re_formula,
                numeric_cov_at = numeric_cov_at,
                levels_id = levels_id,
                ipts = ipts,
                xrange = xrange,
                deriv = deriv,
                summary = summary,
                robust = robust,
                probs = probs,
                ...)
    
    if(!deriv_model) {
      if(deriv == 1 | deriv == 2) {
        . <- mapderivqr(model, ., newdata, deriv)
      }
    } 
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname predict_bgm.bgmfit
#' @export
predict_bgm <- function(model, ...) {
  UseMethod("predict_bgm")
}


