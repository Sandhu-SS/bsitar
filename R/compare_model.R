

#' Model comparison with the loo package
#' 
#' @description 
#' A wrapper around [brms::loo_compare()] to compare models using fit criteria.
#' 
#' @details
#' All models should have pre computed fit criterion See [add_model_criterion()]
#' for more help.
#' 
#'
#' @param model An \code{R} object of class \code{bsitar}
#' 
#' @param check_criterion A logical (default \code{FALSE}) indicating whether to
#'   check and add fit criterion if not already added to model objects
#'   
#' @param add_criterion_args A named list to set up the arguments for the
#'   [add_model_criterion()] function. Ignored when \code{check_criterion =
#'   FALSE} or when each model objects has precomputed fit criteria.
#'   
#' @param verbose A logical (default \code{FALSE}) to print some useful
#'   information.
#'   
#' @inheritParams brms::loo_compare
#' @inherit brms::loo_compare return
#'
#'
#' @rdname compare_model
#' @export
#' 
#' @seealso [brms::loo_compare()]
#' 
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid model estimation which can take time, the Bayesian SITAR model fit
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' model <- getNsObject(berkeley_exfit)
#' 
#' # Add model fit criteria (e.g., WAIC). For illustration purposes, we compare 
#' # model with itself
#' model_1 <- add_model_criterion(model, criterion = c("waic"))
#' model_2 <- add_model_criterion(model, criterion = c("waic"))
#' 
#' # compare models model_1 and model_2
#' compare_model_12 <- compare_model(model_1, model_2, criterion = c("waic"))
#' }
#' 
compare_model.bgmfit <- function(model, 
                                 ..., 
                                 criterion = "loo", 
                                 model_names = NULL,
                                 check_criterion = TRUE,
                                 add_criterion_args = list(),
                                 verbose = FALSE) {
  
  if(!is.list(add_criterion_args)) {
    stop("Argument 'add_criterion_args' must be a named list")
  }
  
  is_brmsfit <- function(x) inherits(x, "brmsfit")
  
  flatten_models <- function(x) {
    if (is.null(x)) return(list())
    if (is_brmsfit(x)) return(list(x))
    if (is.list(x)) return(unlist(lapply(x, flatten_models), recursive = FALSE))
    stop("All inputs must be brmsfit objects or lists of brmsfit objects.")
  }
  
  flatten_exprs <- function(expr, value) {
    if (is.null(value)) return(character())
    if (is_brmsfit(value)) return(deparse(expr, nlines = 1))
    
    if (is.list(value)) {
      out <- vector("list", length(value))
      nm <- names(value)
      
      for (i in seq_along(value)) {
        child_expr <- if (!is.null(nm) && nzchar(nm[i])) {
          as.name(nm[i])
        } else if (is.call(expr) && identical(expr[[1]], as.name("list"))) {
          expr[[i + 1]]
        } else {
          as.call(list(as.name("[["), expr, i))
        }
        out[[i]] <- flatten_exprs(child_expr, value[[i]])
      }
      return(unlist(out, use.names = FALSE))
    }
    stop("All inputs must be brmsfit objects or lists of brmsfit objects.")
  }
  
  has_criterion <- function(fit, criterion) {
    !is.null(fit$criteria) && !is.null(fit$criteria[[criterion]])
  }
  
  exprs <- as.list(substitute(list(model, ...)))[-1]
  vals  <- c(list(model), list(...))
  
  models <- unlist(lapply(vals, flatten_models), recursive = FALSE)
  
  if (is.null(model_names)) {
    model_names <- unlist(Map(flatten_exprs, exprs, vals), use.names = FALSE)
  }
  
  if (length(models) < 2) {
    stop("Need at least two models for loo_compare().")
  }
  
  
  if(!check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion(fit, criterion)) {
        stop("No precomputed criterion availabel for one or more models.",
             " Either add criterion before hand using 'add_model_criterion()'",
             " or else set check_criterion = TRUE. Note that arguments to",
             " 'add_model_criterion()' function can be set by using",
             " 'add_criterion_args' which must be a named list")
      }
    })
  } # if(!check_criterion) {
  
  
  if(check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion(fit, criterion)) {
        suppressWarnings({
        fit <- do.call(
          brms::add_criterion,
          c(list(x = fit, criterion = criterion), add_criterion_args))
        })
      }
      fit
    })
  } # if(check_criterion) {
  
  
  if (length(model_names) != length(models)) {
    # model_names <- paste0("model", seq_along(models))
    nnames <- length(model_names)
    model_names_all <- paste0("model", seq_along(models)-length(model_names))
    model_names_all[1:nnames] <- model_names
    model_names <- model_names_all
    message2c("The number of model names is not same as the number of models.  
              The remaining models are named sequentially as model1,...")
  }
  
  do.call(brms::loo_compare, c(models, 
                               list(criterion = criterion, 
                                    model_names = model_names)))

}



#' @rdname compare_model
#' @export
compare_model <- function(model, ...) {
  UseMethod("compare_model")
}


#' An alias of 'compare_model()'
#' @rdname compare_model
#' @export
#' 
compare_models <- function(model, ...) {
  UseMethod("compare_model")
}



