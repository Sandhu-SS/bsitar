

#' @title Estimate model based individual growth parameters
#' 
#' @description The \strong{modelbased_growthparameters()} function estimates
#'   individual growth parameters by mapping population average estimate of age
#'   of interest (such as age at peak growth velocity or age at take off) on
#'   individual velocity curves defined by individual level random effects.
#'   After the individual level age parameter is found, the corresponding
#'   distance and velocity is also estimated.
#'   
#' @details Since SITAR is a shape-invariant model, each individual curve has a  
#' peak velocity point that can be mapped by knowing the population average age 
#' at peak velocity. This hold true even when a individual lacks measurements at 
#' the expected turning point.
#' 
#' @inheritParams growthparameters.bgmfit
#' @inheritParams marginal_growthparameters.bgmfit
#' @inheritParams marginal_comparisons.bgmfit
#' @inheritParams marginaleffects::predictions
#' @inheritParams marginaleffects::plot_predictions
#' @inheritParams brms::fitted.brmsfit
#' 
#' @param ... Additional arguments passed to the function. 
#' 
#' @return A named list of 3 comprising individual level estimate of
#'   \strong{age}, \strong{distance} and \strong{velocity}. Each of the list is
#'   a data frame with one row per individual and six columns.\cr
#'   
#'   \strong{age} \cr \item{id}{subject identifier} \item{Estimate}{subject's
#'   age corresponding to \code{x}.} \item{Est.Error}{SD of Estimate} \item{Q2.5
#'   }{Lower CI} \item{Q97.5}{Upper CI}
#'   \item{missing}{logical flags where TRUE means subject's specified age lies
#'   outside their measurement range}
#'   
#'   \strong{distance} \cr \item{id}{subject identifier}
#'   \item{Estimate}{distance corresponding to subject's age}
#'   \item{Est.Error}{SD of Estimate} \item{Q2.5 }{Lower CI} \item{Q97.5}{Upper
#'   CI}
#'   \item{missing}{logical flags where TRUE means subject's specified age lies
#'   outside their measurement range}
#'   
#'   \strong{velocity} \cr \item{id}{subject identifier}
#'   \item{Estimate}{velocity corresponding to subject's age}
#'   \item{Est.Error}{SD of Estimate} \item{Q2.5 }{Lower CI} \item{Q97.5}{Upper
#'   CI}
#'   \item{missing}{logical flags where TRUE means subject's specified age lies
#'   outside their measurement range}
#' 
#' @rdname modelbased_growthparameters
#' @export
#' 
#' @seealso [marginaleffects::predictions()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' 
#' \donttest{
#' 
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' modelbased_growthparameters(model, ndraws = 10)
#' 
#' }
#' 
modelbased_growthparameters.bgmfit <-
  function(model,
           resp = NULL,
           dpar = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           newdata = NULL,
           re_formula = NA,
           parameter = NULL,
           acg_velocity = 0.10,
           digits = 2,
           seed = 123,
           future = FALSE,
           future_session = 'multisession',
           future_splits = NULL,
           future_method = 'future',
           future_re_expose = NULL,
           envir = NULL, 
           ...) {
    
    
   
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- parent.frame()
    }
    
    
    # Depending on dpar 'mu' or 'sigma', subset model_info
    model <- getmodel_info(model = model, dpar = dpar)
    
    
    ndraws_org <- ndraws
    ndraws_exe <- FALSE
    if(!is.null(ndraws)) {
      ndraws_exe <- TRUE
      ndraws <- ndraws
    }
    
    if(is.null(ndraws)) {
      ndraws <- brms::ndraws(model)
    }
    
    out <- modelbased_growthparameters_call(model,
                    resp = resp,
                    dpar = dpar,
                    ndraws = ndraws,
                    draw_ids = draw_ids,
                    newdata = newdata,
                    re_formula = re_formula,
                    parameter = parameter,
                    acg_velocity = acg_velocity,
                    digits = digits,
                    seed = seed,
                    future = future,
                    future_session = future_session,
                    future_splits = future_splits,
                    future_method = future_method,
                    future_re_expose = future_re_expose,
                    envir = envir, 
                    ...) 
    
    
    return(out) 
    
  }


#' @rdname modelbased_growthparameters
#' @export
modelbased_growthparameters <- function(model, ...) {
  UseMethod("modelbased_growthparameters")
}


