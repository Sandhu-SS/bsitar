

#' Predicted values (draws) from the posterior predictive distribution
#' 
#' @description The \strong{predict_draws} function is a wrapper around the
#'   [brms::predict.brmsfit()] function to obtain predicted values (and their
#'   summary) from the posterior distribution. See [brms::predict.brmsfit()] for
#'   details.
#' 
#' @details The \strong{fitted_draws} function computed the fitted values from
#'   the posterior distribution. The [brms::predict.brmsfit()] function
#'   from the \pkg{brms} package can used to get the predicted (distance) values
#'   when outcome (e.g., height) is untransformed. However, when the outcome is
#'   log or square root transformed, the [brms::predict.brmsfit()] function will
#'   return the fitted curve on the log or square root scale whereas the
#'   \strong{predict_draws} function returns the fitted values on the original
#'   scale. Furthermore, the \strong{predict_draws} also compute the first
#'   derivative of (velocity) that too on the original scale after making
#'   required back-transformation. Except for these differences, both these
#'   functions (i.e., [brms::predict.brmsfit()] and [predict_draws()]) work in
#'   the same manner. In other words, user can specify all the options available
#'   in the [brms::predict.brmsfit()].
#' 
#' 
#' @inherit growthparameters.bgmfit params
#' 
#' @inherit plot_conditional_effects.bgmfit params
#' 
#' @inherit brms::predict.brmsfit params
#' 
#' @param ... Additional arguments passed to the [brms::predict.brmsfit()]
#'   function. Please see [brms::predict.brmsfit()] for details on various
#'   options available.
#' 
#' @return An array of predicted response values. See [brms::predict.brmsfit()] 
#' for details.
#' 
#' @export predict_draws.bgmfit
#' @export
#' 
#' @seealso [brms::predict.brmsfit()] 
#' 
#' @inherit bgm author
#'
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid fitting the model which takes time, the model  
#' # fit has already been saved as 'berkeley_mfit.rda' file.
#' # See examples section of the bgm function for details on the model fit.
#' 
#' model <- berkeley_mfit
#' 
#' # Population average distance curve
#' predict_draws(model, deriv = 0, re_formula = NA)
#' 
#' \donttest{
#' # Individual-specific distance curves
#' predict_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' predict_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' predict_draws(model, deriv = 1, re_formula = NULL)
#'  }
#' 
predict_draws.bgmfit <-
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
           envir = globalenv(),
           ...) {
    
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv)
    
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
        . <- mapderivqr(model, ., newdata = newdata, resp = resp, 
                        deriv = deriv, probs = probs, robust = robust)
      }
    } 
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname predict_draws.bgmfit
#' @export
predict_draws <- function(model, ...) {
  UseMethod("predict_draws")
}


