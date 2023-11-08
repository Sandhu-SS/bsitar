

#' Visualize the conditional effects of predictor
#' 
#' @details The \code{conditional_effects_bgm} function is a wrapper around 
#' the [brms::conditional_effects()]. The [brms::conditional_effects()]  
#' function from the \code{brms} package can used to plot the fitted (distance) 
#' curve for an *bgmfit* model when response (e.g., height) is not transformed.
#' However, when the outcome is log or square root transformed, the 
#' [brms::conditional_effects] will return the fitted curve on the log or 
#' square root scale whereas the [bsitar::conditional_effects_bgm()] will 
#' return the fitted curve on the original scale. Furthermore, the 
#' conditional_effects_bgm also displays the velocity curve on the original scale
#' after making required back-transformation. Apart from these differences, 
#' both these functions ([brms::conditional_effects] and 
#' [bsitar::conditional_effects_bgm()]) work in the same manner. In other words, 
#' user can specify all the arguments which are available in the 
#' [brms::conditional_effects]. Because of this, the name is kept same except 
#' for adding an underscore at the end of the name (*conditional_effects* to 
#' *conditional_effects_bgm*). 
#'
#' @param model An object of class \code{bgmfit}. function.
#' 
#' @param resp Response variable (default \code{NULL}) specified as a string
#'   character required during the post-processing of multivariate and
#'   univariate-by-subgroup model (see \code{bsitar::bgm()} for details).
#'   
#' @param deriv An integer to specify whether to estimate distance curve or
#'   derivatives (velocity and acceleration curves). Default \code{deriv = 0} is
#'   for the distance curve whereas \code{deriv = 1} for velocity curve and
#'   \code{deriv = 2} for the acceleration curve.
#'   
#' @param deriv_model A logical (default \code{TRUE}) to indicate whether to
#'   estimate model based derivatives or from the differentiation of the
#'   distance curve. When model is fit with \code{decomp = 'QR'}, the only
#'   approach available to estimate derivatives is the  differentiation of the
#'   distance curve.
#'   
#' @param envir The calling environment. Deafault set to \code{globalenv()}.
#' 
#' @param ... Additional arguments passed to the [brms::conditional_effects()]
#'   function. Please see [brms::conditional_effects()] for details.
#'
#' @inherit brms::conditional_effects description 
#' 
#' @return An object of class 'brms_conditional_effects' which is a named list
#'   with one data.frame per effect containing all information required to
#'   generate conditional effects plots. See brms::conditional_effects for
#'   details.
#'
#' @export
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @examples
#' 
#' # The examples below show the use of *conditional_effects_bgm* to plot  
#' # the population average and individual-specific distance and velocity 
#' # curves.
#' 
#' # Fit Bayesian SITAR model 
#' # To avoid running the model which takes some time, model fit to the
#' # \code{berkeley_mdata} has already been saved as berkeley_mfit.rda object.
#' # Please see \code{bgm} examples.
#' 
#' model <- berkeley_mfit
#' 
#' # Population average distance curve
#' conditional_effects_bgm(model, deriv = 0, re_formula = NA)
#' 
#' # Individual-specific distance curves
#' conditional_effects_bgm(model, deriv = 0, re_formula = NULL)
#' 
#' \donttest{
#' # Population average velocity curve
#' conditional_effects_bgm(model, deriv = 1, re_formula = NA)
#' 
#' # Individual-specific velocity curves
#' conditional_effects_bgm(model, deriv = 1, re_formula = NULL)
#' }
#' 
conditional_effects_bgm.bgmfit <-
  function(model,
           resp = NULL,
           deriv = 0,
           deriv_model = TRUE,
           envir = globalenv(),
           ...) {
    o <-
      post_processing_checks(model = model,
                             xcall = match.call(),
                             resp = resp,
                             envir = envir,
                             deriv = deriv)
    
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
    
    
    if(!deriv_model) {
      xvar  <- model$model_info$xvar
      idvar <- model$model_info$groupvar
      if(length(idvar) > 1) idvar <- idvar[1]
      yvar  <- 'yvar'
      out_    <- brms::conditional_effects(model, resp = resp, ...)
      datace <- out_[[1]] %>% dplyr::select(dplyr::all_of(names(model$data)))
      datace[[idvar]] <- unique(levels(model$data[[idvar]]))[1]
      outx <- fitted_bgm(model, resp = resp, newdata = datace, deriv = deriv, ...)
      out_[[1]][['estimate__']] <- outx[, 1]
      out_[[1]][['se__']] <- outx[, 2]
      out_[[1]][['lower__']] <- outx[, 3]
      out_[[1]][['upper__']] <- outx[, 4]
      # plot(out_, plot = FALSE)[[1]]
      . <- out_
    }
    
    if(deriv_model) . <- brms::conditional_effects(model, resp = resp, ...)
    
    assign(o[[1]], model$model_info[['exefuns']][[o[[1]]]], envir = envir)
    .
  }


#' @rdname conditional_effects_bgm.bgmfit
#' @export
conditional_effects_bgm <- function(model, ...) {
  UseMethod("conditional_effects_bgm")
}





