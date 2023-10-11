

#' Fitted (expected) values from the posterior predictive distribution
#' 
#' @details The \code{fitted_bgm} function is a wrapper around 
#' the [brms::fitted.brmsfit]. The [brms::fitted.brmsfit]  
#' function from the \code{brms} package can used to plot the fitted (distance) 
#' curve for an *bsitar* model when outcome (e.g., height) is untransformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::fitted.brmsfit] will return the fitted curve on the log or 
#' square root scale whereas the [bsitar::fitted_bgm()] will 
#' return the fitted curve on the original scale. Furthermore, the 
#' fitted_bgm also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::fitted.brmsfit] and 
#' [bsitar::fitted_bgm()]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::fitted.brmsfit]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*fitted* to 
#' *fitted_bgm*). 
#' 
#' @inherit growthparameters.bgmfit params
#' 
#' @inherit conditional_effects_bgm.bgmfit params
#' 
#' @inherit brms::fitted.brmsfit description params
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit] 
#' function. Please see \code{brms::fitted.brmsfit} for details on 
#' various options available.
#' 
#' @return An array of predicted mean response values. See [brms::fitted.brmsfit] 
#' for details.
#' 
#' @export fitted_bgm.bgmfit
#' @export
#'
#' @examples
#' \dontrun{
#' # The following examples show the use of *fitted_bgm* to estimate  
#' # population average and individual-specific distance and velocity 
#' # curves.
#' # Population average distance curve
#' fitted_bgm(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' fitted_bgm(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' fitted_bgm(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' fitted_bgm(model, deriv = 1, re_formula = NULL)
#'  
#' }

fitted_bgm.bgmfit <-
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
   
    
    . <- fitted(model,
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
        . <- mapderivqr(model, ., newdata = newdata, resp = resp, 
                        deriv = deriv, probs = probs, robust = robust)
      }
    } 
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname fitted_bgm.bgmfit
#' @export
fitted_bgm <- function(model, ...) {
  UseMethod("fitted_bgm")
}


