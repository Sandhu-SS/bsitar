

#' Update \pkg{bsitar} model
#'
#' @param model An object of class \code{bsitar}.
#' @param newdata An optional \code{data.frame} to be used when updating the
#'   model. If \code{NULL} (default), the same same data used for original model
#'   fit is used. Note that data-dependent default priors will not be updated
#'   automatically.
#' @param recompile A logical to indicate whether the Stan model should be
#'   recompiled. When \code{NULL} (the default), \code{update_bsitar} tries to
#'   figure out internally, if recompilation is necessary. Setting it to
#'   \code{FALSE} will cause all Stan code changing arguments to be ignored.
#' @param ... Other arguments passed to \code{\link{brms}}.
#'
#' @return An updated object of class \code{brmsfit, bsiatr}, that contains the
#'   posterior draws and other useful information about the model.
#'
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @export update_bsitar.bsitar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(heights)
#' data_males <- heights %>% filter(sex == 'Male')
#' fit_males <- bsitar(x=age, y=height, id=id, data=heights, df=4)
#' fit_males2 <- update_bsitar(df=5)
#' }
#' 
update_bsitar.bsitar <- function(model, newdata = NULL, recompile = NULL, ...) {
  
  
  mcall_ <- match.call()
  call_ <- model$model_info$call.full.bsitar
  if(!is.null(newdata)) {
    call_$data <- mcall_$newdata
  }
  
  mcall_$newdata <- NULL
  
  
  
  if(!identical(names(mcall_)[-c(1:2)], character(0))) {
    for (i in names(names(mcall_)[-c(1:2)])) {
      if(!i %in% names(call_)) {
        stop("Argument ", i, " is not a valid arguments",
             " \n ",
             " Please see 'bsitar' function ")
      }
    }
  }
  
  
  
  
  eval_identical <- function(x, y) {
    identical(eval(x), eval(y))
  }
  
  if(!is.null(mcall_$backend)) {
    if(eval_identical(mcall_$backend, call_$backend)) {
      same_backend <- TRUE
    } else {
      same_backend <- FALSE
    }
  } else if(is.null(mcall_$backend)) {
    same_backend <- TRUE
  }
  
  
  for (i in names(mcall_)[-c(1:2)]) {
    call_[[i]] <- mcall_[[i]]
  }
  
  
  
  iter <- call_$iter
  silent <- call_$silent
  
  as_one_logical <- is_equal <- needs_recompilation <- NULL
  as_one_logical <- utils::getFromNamespace("as_one_logical", "brms")
  is_equal <- utils::getFromNamespace("is_equal", "brms")
  needs_recompilation <- utils::getFromNamespace("needs_recompilation", "brms")
  
  
  if (is.null(recompile)) {
    call_stancode <- call_
    call_stancode$get_stancode <- TRUE
    new_stancode <- eval.parent(call_stancode)
    new_stancode <- sub("^[^\n]+[[:digit:]]\\.[^\n]+\n", "", new_stancode)
    old_stancode <- model$bmodel
    old_stancode <- sub("^[^\n]+[[:digit:]]\\.[^\n]+\n", "", old_stancode)
    
    # Not using but good for later use
    
    # call_stancode$get_stancode <- FALSE
    # call_stancode$get_standata <- TRUE
    # new_standata <- eval.parent(call_stancode)
    # old_standata <- standata(model)
    # if(new_standata$prior_only == old_standata$prior_only) {
    #   sample_prior_arg_same <- TRUE
    # } else {
    #   sample_prior_arg_same <- FALSE
    # }
    
    # !sample_prior_arg_same || 
    
    recompile <- needs_recompilation(model) || !same_backend || 
      !is_equal(new_stancode, old_stancode)
    if (recompile && silent < 2) {
      message("Update requires model recompilation")
    }
  }
  
  recompile <- as_one_logical(recompile)
  
  # This controls the sample prior argument
  attr(model$prior, "sample_prior") <- call_$sample_prior 
  
  if (recompile) {
    call_$fit <- NA
    model <- eval.parent(call_) 
  } else {
    call_$fit <- model
    model <- eval.parent(call_)
  }
  # if(call_$expose_function) model <- expose_bsitar_functions(model, model$bmodel)
  return(model)
}



#' @rdname update_bsitar.bsitar
#' @export
update_bsitar <- function(model, ...) {
  UseMethod("update_bsitar")
}

