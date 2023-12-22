

#' Perform leave-one-out (loo) cross-validation
#' 
#' @description The \strong{loo_validation()} is a wrapper around the
#'   [brms::loo()] function to perform approximate leave-one-out
#'   cross-validation based on the posterior likelihood. See [brms::loo()] for
#'   details.
#' 
#' @inherit brms::loo details 
#' 
#' @inheritParams growthparameters.bgmfit
#' 
#' @inheritParams plot_ppc.bgmfit
#' 
#' @param ... Additional arguments passed to the [brms::loo()] function. 
#' Please see \code{brms::loo} for details on various options available.
#' 
#' @return If only one model object is provided, then an object of class
#'   \code{loo} is returned. If multiple objects are provided, an object of
#'   class \code{loolist}.
#' 
#' @export loo_validation.bgmfit
#' @export
#' 
#' @seealso [brms::loo()] 
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
#' \donttest{
#' loo_validation(model, cores = 1)
#' }
#' 
#' 
loo_validation.bgmfit <-
  function(model,
           resp = NULL,
           cores = 1,
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
    
    # setNULL <- as.character(expression(c(deriv, 
    #                                 numeric_cov_at, 
    #                                 levels_id,
    #                                 ipts,
    #                                 idata_method,
    #                                 xrange,
    #                                 probs,
    #                                 robust,
    #                                 )))[-1]
    # 
    # setFALSE <- as.character(expression(c(deriv_model
    # )))[-1]
    # 
    # for (setNULLi in setNULL) {
    #   if(!checkifargmiss(charg, chcall, setNULLi)) assign(setNULLi, NULL)
    # }
    # 
    # for (setFALSEi in setFALSE) {
    #   if(!checkifargmiss(charg, chcall, setFALSE)) assign(setFALSE, FALSE)
    # }
    
    
    
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
      if(model$model_info$decomp == "QR") deriv_model <- FALSE
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
    
    test <- setupfuns(model = model, o = o, oall = oall, 
                      usesavedfuns = usesavedfuns, 
                      deriv = deriv, 
                      envir = envir, 
                      deriv_model = deriv_model, 
                      ...)
    
    if(is.null(test)) return(invisible(NULL))
    
    
    . <- brms::loo(model, resp = resp, cores = cores ,...)
    
    
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


#' @rdname loo_validation.bgmfit
#' @export
loo_validation <- function(model, ...) {
  UseMethod("loo_validation")
}

