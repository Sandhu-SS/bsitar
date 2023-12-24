

#' Fitted (expected) values from the posterior predictive distribution
#' 
#' @description The \strong{fitted_draws()} is a wrapper around the
#'   [brms::fitted.brmsfit()] function to obtain fitted values (and their
#'   summary) from the posterior draws. See [brms::fitted.brmsfit()] for
#'   details.
#'
#' @details The \strong{fitted_draws()} computes the fitted values from the
#'   posterior draws. The [brms::fitted.brmsfit()] function from the \pkg{brms}
#'   package can used to get the fitted (distance) values when outcome (e.g.,
#'   height) is untransformed. However, when the outcome is log or square root
#'   transformed, the [brms::fitted.brmsfit()] function will return the fitted
#'   curve on the log or square root scale whereas the \strong{fitted_draws()}
#'   function returns the fitted values on the original scale. Furthermore, the
#'   \strong{fitted_draws()} also compute the first derivative of (velocity)
#'   that too on the original scale after making required back-transformation.
#'   Except for these differences, both these functions (i.e.,
#'   [brms::fitted.brmsfit()] and [fitted_draws()]) work in the same manner. In
#'   other words, user can specify all the options available in the
#'   [brms::fitted.brmsfit()].
#' 
#' @inherit growthparameters.bgmfit params
#' @inherit plot_conditional_effects.bgmfit params
#' @inherit brms::fitted.brmsfit params
#' 
#' @param ... Additional arguments passed to the [brms::fitted.brmsfit()] 
#' function. Please see \code{brms::fitted.brmsfit()} for details on 
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
#' @inherit berkeley author
#'
#' @examples
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid fitting the model which takes time, the model  
#' # fit has already been saved as 'berkeley_mfit.rda' file.
#' # See examples section of the main function for details on the model fit.
#' 
#' model <- berkeley_mfit
#' 
#' # Population average distance curve
#' fitted_draws(model, deriv = 0, re_formula = NA)
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
           verbose = FALSE,
           usesavedfuns = FALSE,
           clearenvfuns = NULL,
           envir = NULL,
           ...) {
   
  
    if(is.null(envir)) {
      envir <- parent.frame()
    }
    
    # For consistency across post processing functions
    charg <- ls()
    chcall <- match.call()
    if(!checkifargmiss(charg, chcall, 'deriv'))          deriv          <- NULL
    if(!checkifargmiss(charg, chcall, 'numeric_cov_at')) numeric_cov_at <- NULL
    if(!checkifargmiss(charg, chcall, 'levels_id'))      levels_id      <- NULL
    if(!checkifargmiss(charg, chcall, 'ipts'))           ipts           <- NULL
    if(!checkifargmiss(charg, chcall, 'idata_method'))   idata_method   <- NULL
    if(!checkifargmiss(charg, chcall, 'xrange'))         xrange         <- NULL
    if(!checkifargmiss(charg, chcall, 'probs'))          probs          <- NULL
    if(!checkifargmiss(charg, chcall, 'robust'))         robust         <- NULL
    if(!checkifargmiss(charg, chcall, 'newdata'))        newdata        <- NULL
    if(!checkifargmiss(charg, chcall, 'deriv_model'))    deriv_model    <- FALSE
    
    
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
    
    
    if(!is.null(model$model_info$decomp)) {
      if(model$model_info$decomp == "QR") deriv_model<- FALSE
    }
    
    expose_method_set <- model$model_info[['expose_method']]
    
    model$model_info[['expose_method']] <- 'NA' # Over ride method 'R'
    
    o <- post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv, 
                             all = FALSE)
    
    oall <- post_processing_checks(model = model,
                                xcall = match.call(),
                                resp = resp,
                                envir = envir,
                                deriv = deriv, 
                                all = TRUE)
    
    
    # if(!usesavedfuns) {
    #   if(is.null(check_if_functions_exists(model, o, model$xcall))) {
    #     return(invisible(NULL))
    #   }
    # }
    
    
    # if(usesavedfuns) {
    #   if(is.null(check_if_functions_exists(model, o, model$xcall))) {
    #     envir <- envir
    #   } else {
    #     envir <- getEnv(o[[1]], geteval = TRUE)
    #   }
    #   oall <- model$model_info[['exefuns']]
    #   oalli_c <- names(oall)
    #   for (oalli in oalli_c) {
    #     assign(oalli, oall[[oalli]], envir = envir)
    #   }
    # }
    

    # if(deriv == 0) {
    #   assignfun <- paste0(model$model_info[['namesexefuns']], deriv)
    #   assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
    # } else if(deriv > 0) {
    #   if(deriv_model) {
    #     assignfun <- paste0(model$model_info[['namesexefuns']], deriv)
    #   } else if(!deriv_model) {
    #     assignfun <- paste0(model$model_info[['namesexefuns']], '0')
    #   }
    #   assign(o[[1]], model$model_info[['exefuns']][[assignfun]], envir = envir)
    # }
    
    
    test <- setupfuns(model = model, resp = resp,
                      o = o, oall = oall, 
                      usesavedfuns = usesavedfuns, 
                      deriv = deriv, envir = envir, 
                      deriv_model = deriv_model, 
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
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
    
    if(!is.null(deriv)) {
      if(deriv > 0) { 
        if(!deriv_model) {
          . <- mapderivqr(model, ., newdata = newdata, resp = resp, 
                          deriv = deriv, probs = probs, robust = robust)
        } else {
          . <- .
        }
      }
    }
    

    
    # Restore function(s)
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    
    if(!is.null(clearenvfuns)) {
      if(!is.logical(clearenvfuns)) {
        stop('clearenvfuns must be NULL or a logical')
      } else {
        setcleanup <- clearenvfuns
      }
    }
    
    if(is.null(clearenvfuns)) {
      if(usesavedfuns) setcleanup <- TRUE else setcleanup <- FALSE
    }
      
    # Cleanup environment if requested
    if(setcleanup) {
      tempgenv <- envir
      for (oalli in names(oall)) {
        if(exists(oalli, envir = tempgenv )) {
          remove(list=oalli, envir = tempgenv)
        }
      }
      tempgenv <- test
      for (oalli in names(oall)) {
        if(exists(oalli, envir = tempgenv )) {
          remove(list=oalli, envir = tempgenv)
        }
      }
      
    } # if(setcleanup) {
    
    .
  }


#' @rdname fitted_draws.bgmfit
#' @export
fitted_draws <- function(model, ...) {
  UseMethod("fitted_draws")
}


