

#' Perform posterior predictive distribution checks
#' 
#' @details The \strong{plot_ppc()} is a wrapper around the [brms::pp_check()].
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @inherit growthparameters.bgmfit params
#' @inherit brms::pp_check.brmsfit description 
#' @inherit fitted_draws.bgmfit params
#' 
#' @param ... Additional arguments passed to the [brms::pp_check.brmsfit()] 
#' function. Please see [brms::pp_check.brmsfit()] for details.
#' 
#' @return A ggplot object that can be further customized using the
#' ggplot2 package.
#' 
#' @export
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
#' plot_ppc(model, ndraws = 100)
#' 
plot_ppc.bgmfit <-
  function(model,
           resp = NULL,
           deriv = 0,
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
    

    . <- brms::pp_check(model, resp = resp, ...)
    
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


#' @rdname plot_ppc.bgmfit
#' @export
plot_ppc <- function(model, ...) {
  UseMethod("plot_ppc")
}


