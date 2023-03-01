

#'Optimize \pkg{bsitar} model
#'
#'@param model An object of class \code{bsitar}.
#'@param newdata An optional \code{data.frame} to be used when optimising the
#'  model. If \code{NULL} (default), the same same data used for original model
#'  fit is used. Note that data-dependent default priors will not be updated
#'  automatically.
#'@param optimize_df A vector specifying the updation of degree of freedom
#'  \code{df}. For \code{univariate-by-sungroup} and \code{multivariate} model
#'  (see [bsitar::bsitar] for details on these arguments), \code{optimize_df}
#'  can be a single integer (e.g., \code{optimize_df = 4}) or a list as a string
#'  (e.g., \code{optimize_df = 'list(4,5)'}). If \code{NULL}, then \code{df} is
#'  taken from the original model.
#'@param optimize_x A vector specifying the transformations of predictor
#'  (typically \code{age}) variable (via \code{xvar}). The option are 'NULL',
#'  'log' and  'sqrt' or their combinations. The default \code{optimize_x =
#'  c('NULL', 'log',  'sqrt')} is to explore all possible combinations of
#'  'NULL', 'log' and  'sqrt'.
#'
#'@param optimize_y A vector specifying the transformations of the outcome
#'  variable (via \code{yvar}). The options are 'NULL', 'log' and  'sqrt' or
#'  their combinations. The default \code{optimize_y = c('NULL', 'log',
#'  'sqrt')} is to explore all possible combinations of 'NULL', 'log' and
#'  'sqrt'.
#'
#'
#'@param exclude_default_funs A logical (deafult \code{FALSE}) to indicate
#'  whether transformations (\code{xvar} and \code{yvar}) used in the original
#'  model fit should be excluded. If \code{TRUE}, then the \code{xvar} and
#'  \code{yvar} transformations spevified for the original model fit are
#'  excluded from the \code{optimize_x} and \code{optimize_y}. From example, if
#'  original model is fit with \code{xvar = log} and \code{yvar = NULL}, then
#'  the \code{optimize_x} is translated into \code{optimize_x = c('NULL',
#'  'sqrt')} and  \code{optimize_y} as \code{optimize_y = c('log', 'sqrt')}.
#'
#'
#'@param add_fit_criteria An optional (default \code{NULL}) indicator to add fit
#'  criteria to the model fit. options are \code{loo} and \code{waic}. Please
#'  see [brms::add_criterion()] for details.
#'
#'@param add_fit_bayes_R An optional (default \code{NULL}) to add Bayesian R
#'  square.
#'
#'@param ... Other arguments passed to \code{\link{update_bsitar}}.
#'
#'@return An updated object of class \code{brmsfit, bsiatr}, that contains the
#'  posterior draws and other useful information about the model.
#'
#'@author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#'@export optimize_bsitar.bsitar
#'
#'@export
#'
#' @examples
#' \dontrun{
#' data(heights)
#' data_males <- heights %>% filter(sex == 'Male')
#' fit_males <- bsitar(x=age, y=height, id=id, data=heights, df=4)
#' fit_males2 <- optimize_bsitar(fit_males)
#' }
#' 
optimize_bsitar.bsitar <- function(model, newdata = NULL, 
                                   optimize_df = NULL,
                                   optimize_x = c('NULL', 'log',  'sqrt'),
                                   optimize_y = c('NULL', 'log',  'sqrt'),
                                   exclude_default_funs = FALSE,
                                   add_fit_criteria = NULL,
                                   add_fit_bayes_R = NULL,
                                   ...) {
  
  
  call_o <- match.call()
  call_o_args <- as.list(call_o)[-1]
  
  args_o <- as.list(model$model_info$call.full.bsitar)[-1]
  args_o_dots_ <- list(...)
  if (length(args_o_dots_) > 0) {
    for (i in names(args_o_dots_)) {
      args_o[[i]] <- args_o_dots_[[i]]
    }
  }
  
  # if(is.null(newdata)) {
  #   data <- eval.parent(model$model_info$call.bsitar$data)
  # }

  if(exclude_default_funs) {
    optimize_x <- optimize_x[!optimize_x %in% model$model_info$xfuns]
    optimize_y <- optimize_y[!optimize_y %in% model$model_info$xfuns]
  }
  
  if(is.null(optimize_df)) optimize_df <- "NULL"
  if(is.null(optimize_x)) optimize_x <- "NULL"
  if(is.null(optimize_y)) optimize_y <- "NULL"
  
  optimize_df_x_y <- 
    expand.grid(optimize_df, optimize_x, optimize_y)
  
  colnames(optimize_df_x_y) <- c("df", "xfun", "yfun")

  optimize_fun <- function(.x, model ) {
    # xsplit <- strsplit(.x, "_")[[1]]
    exe_row <- optimize_df_x_y[.x, ]
    df <- levels(droplevels(exe_row$df)) 
    xfun <- levels(droplevels(exe_row$xfun)) # xsplit[1]
    yfun <- levels(droplevels(exe_row$xfun)) # xsplit[2]
    if(df == 'NULL') df <- paste0("list(", paste(model$model_info$dfs, collapse = ","), ")") else df <- df
    if(xfun == 'NULL') xfun <- NULL else xfun <- xfun
    if(yfun == 'NULL') yfun <- NULL else yfun <- yfun

    if(is.null(xfun)) xfun_print <- deparse(xfun) else xfun_print <- xfun
    if(is.null(yfun)) yfun_print <- deparse(yfun) else yfun_print <- yfun
    cat("\n")
    cat(paste0("df = ", df, "; xfun = ", xfun_print, "; yfun = ", yfun_print), "\n")

    optimization_info <- paste0("df = ", df, "; xfun = ", xfun_print, "; yfun = ", yfun_print)    

    
    args_o$model <- model
    args_o$df <- eval(parse(text = df))
    args_o$xfun <- xfun
    args_o$yfun <- yfun
    
    if(!is.null(newdata)) {
      args_o$data <- call_o_args$newdata
    }
    
    fit <- do.call(update_bsitar, args_o)
    fit$model_info$optimization_info <- optimization_info
   
    ##########
    # Very important to set cores = 1 on windows
    if(!is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_criteria, collapse = ", ")
      message("Adding", " ", what_, " ", "...")
      cat("\n")
      suppressWarnings(fit <- add_criterion(fit, add_fit_criteria, cores = 1))
    }
    
    if(!is.null(add_fit_bayes_R) & !is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_bayes_R, collapse = ", ")
      message("Adding", " ", what_, " ", "...")
      cat("\n")
      suppressWarnings(fit$criteria$bayes_R <- bayes_R2(fit))
    }
    
    if(!is.null(add_fit_bayes_R) & is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_bayes_R, collapse = ", ")
      message("Adding", " ", what_, " ", "...")
      cat("\n")
      suppressWarnings(fit$bayes_R <- bayes_R2(fit))
    }
    return(fit)
  }
  lapply(1:nrow(optimize_df_x_y), function(.x) optimize_fun(.x, model))
} 




#' @rdname optimize_bsitar.bsitar
#' @export
optimize_bsitar <- function(model, ...) {
  UseMethod("optimize_bsitar")
}
