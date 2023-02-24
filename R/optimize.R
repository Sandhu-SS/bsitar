

#' optimize \pkg{bsitar} models
#'
#' @param model An object of class \code{bsitar}.
#' @param newdata newdata Optional \code{data.frame} to update the model with
#'   new data. Note that data-dependent default priors will not be updated
#'   automatically.
#' @param recompile A logical to indicate whether the Stan model should be
#'   recompiled. When \code{NULL} (the default), \code{update} tries to figure
#'   out internally, if recompilation is necessary. Setting it to \code{FALSE}
#'   will cause all Stan code changing arguments to be ignored.
#'
#' @param optimize A list specifying the transformations of predictor (typically
#'   \code{age}) and the outcome.
#'
#' @param add_fit_criteria An optional (default \code{NULL}) indicator to add
#'   fit criteria to the model fit. options are \code{loo} and \code{waic}.
#'   Please see [brms::add_criterion()] for details.
#'
#' @param add_fit_bayes_R An optional (default \code{NULL}) to add Bayesian R
#'   square.
#'
#' @param ... Other arguments passed to \code{\link{brms}}.
#'
#' @return An updated object of class \code{brmsfit, bsiatr}, that contains the
#'   posterior draws and other useful information about the model.
#'
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @export update.bsitar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(heights)
#' data_males <- heights %>% filter(sex == 'Male')
#' fit_males <- bsitar(x=age, y=height, id=id, data=heights, df=4)
#' fit_males2 <- update(df=5)
#' }
#' 
optimize.bsitar <- function(model, newdata = NULL, recompile = T, 
                            optimize = list(x = c('NULL', 'log',  'sqrt'),
                                            y = c('NULL', 'log',  'sqrt')),
                            add_fit_criteria = NULL,
                            add_fit_bayes_R = NULL,
                            ...) {
  
  # optimize = list(x = c('NULL', 'log',  'sqrt'),
  #                 y = c('NULL', 'log',  'sqrt'))
  
  
  scale_set_comb <- 
    with(expand.grid(optimize[[1]], optimize[[2]]), paste0(Var1,'_',Var2))
  
  
  brms_fitfun <- function(.x, model) {
    xsplit <- strsplit(.x, "_")[[1]]
    xfun <- xsplit[1]
    yfun <- xsplit[2]
    if(xfun == 'NULL') xfun <- NULL else xfun <- xfun
    if(yfun == 'NULL') yfun <- NULL else yfun <- yfun
    
    #  cat(paste0('Degree of freedom ', .x), '\n')
    if(is.null(xfun)) xfun_print <- deparse(xfun) else xfun_print <- xfun
    if(is.null(yfun)) yfun_print <- deparse(yfun) else yfun_print <- yfun
    cat("\n")
    cat(paste0("Transformations: ", "x = ", xfun_print, "; y = ", yfun_print), "\n")
    
    fit <- update(model, xfun = xfun, yfun = yfun)
    
    ##########
    # Very important to set cores = 1 on windows
    if(!is.null(add_fit_criteria)) {
      message("Adding criterion....\n")
      cat(" ", add_fit_criteria)
      fit <- add_criterion(fit, add_fit_criteria, cores = 1)
      cat(" \n ")
    }
    
    if(!is.null(add_fit_bayes_R) & !is.null(add_fit_criteria)) {
      message("Adding R square...\n")
      cat(" ", add_fit_bayes_R)
      fit$criteria$bayes_R <- bayes_R2(fit)
      cat(" \n ")
    }
    
    if(!is.null(add_fit_bayes_R) & is.null(add_fit_criteria)) {
      message("Adding R square...\n")
      cat(" ", add_fit_bayes_R)
      fit$bayes_R <- bayes_R2(fit)
      cat(" \n ")
    }
    
  }
  
  # brmsmodels <-
  #   scale_set_comb %>% 
  #   purrr::set_names() %>% 
  #   purrr::map(~ brms_fitfun(.x, model))
  
  brmsmodels <- lapply(scale_set_comb, function(.x) brms_fitfun(.x, model))

  
  brmsmodels
  
} # end


# ccc <- update_bsitar(fitx, cores = 2, x = age, a_init_beta = 2, y = y, xfun = NULL, df = 5, data = mmm)


#' @rdname update.bsitar
#' @export
optimize <- function(model, ...) {
  UseMethod("optimize")
}
