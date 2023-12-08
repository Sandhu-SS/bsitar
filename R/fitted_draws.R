

#' Fitted (expected) values from the posterior predictive distribution
#' 
#' @description The \strong{fitted_draws} function is a wrapper around the
#' [brms::fitted.brmsfit()] function to obtain fitted values (and their summary) 
#' from the posterior draws. See [brms::fitted.brmsfit()] for details.
#'
#' @details The \strong{fitted_draws} function computed the fitted values from
#'   the posterior draws. The [brms::fitted.brmsfit()] function from the
#'   \pkg{brms} package can used to get the fitted (distance) values when
#'   outcome (e.g., height) is untransformed. However, when the outcome is log
#'   or square root transformed, the [brms::fitted.brmsfit()] function will
#'   return the fitted curve on the log or square root scale whereas the
#'   \strong{fitted_draws} function returns the fitted values on the original
#'   scale. Furthermore, the \strong{fitted_draws} also compute the first
#'   derivative of (velocity) that too on the original scale after making
#'   required back-transformation. Except for these differences, both these
#'   functions (i.e., [brms::fitted.brmsfit()] and [fitted_draws()]) work in the
#'   same manner. In other words, user can specify all the options available in
#'   the [brms::fitted.brmsfit()].
#' 
#' @inherit growthparameters.bgmfit params
#' 
#' @inherit plot_conditional_effects.bgmfit params
#' 
#' @inherit brms::fitted.brmsfit params
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit] 
#' function. Please see \code{brms::fitted.brmsfit} for details on 
#' various options available.
#' 
#' @return An array of predicted mean response values. See [brms::fitted.brmsfit] 
#' for details.
#' 
#' @export fitted_draws.bgmfit
#' @export
#' 
#' @seealso [brms::fitted.brmsfit()]
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
#' fitted_draws.bgmfit(model, deriv = 0, re_formula = NA)
#' 
#' \donttest{
#' # Individual-specific distance curves
#' fitted_draws(model, deriv = 0, re_formula = NULL)
#' 
#' # Population average velocity curve
#' fitted_draws(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' fitted_draws(model, deriv = 1, re_formula = NULL)
#' }
#' 
fitted_draws.bgmfit <-
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
           envir = NULL,
           ...) {
    
    if(is.null(envir)) {
      envir <- parent.frame()
    }
    
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") deriv_model<- FALSE
    }
    
 
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv,
                             all = FALSE)
    
  
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
    
    
    getfunx1always <- model$model_info[['exefuns']][[o[[1]]]]
    
    if(deriv == 0) {
      getfunx <- model$model_info[['exefuns']][[o[[2]]]]
    } else if(deriv > 0) {
      if(deriv_model) {
        getfunx <- model$model_info[['exefuns']][[o[[2]]]]
      } else if(!deriv_model) {
        summary <- FALSE
        getfunx <- model$model_info[['exefuns']][[o[[1]]]]
      }
    }
    
    setcleanup <- FALSE
    if(!check_if_functions_exists(model, o, model$xcall)) {
      #. <- NULL
       return(invisible(NULL))
      #return(.)
    } else {
      setcleanup <- TRUE
      assign(o[[1]], getfunx, envir = environment(getfunx))
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
    
    assign(o[[1]], getfunx1always, environment(getfunx1always))
    .
  }




#' @rdname fitted_draws.bgmfit
#' @export
fitted_draws <- function(model, ...) {
  UseMethod("fitted_draws")
}



