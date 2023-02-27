


#' Update \pkg{bsitar} models
#'
#' @param model An object of class \code{bsitar}.
#' @param data An optional \code{data.frame} to be used when updating the model.
#'   If \code{NULL} (default), the same same data used for original model fit is
#'   used. Note that data-dependent default priors will not be updated
#'   automatically.
#' @param recompile A logical to indicate whether the Stan model should be
#'   recompiled. When \code{NULL} (the default), \code{update} tries to figure
#'   out internally, if recompilation is necessary. Setting it to \code{FALSE}
#'   will cause all Stan code changing arguments to be ignored.
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
update.bsitar <- function(model, data = NULL, recompile = NULL, ...) {
  
  mcall <- match.call()
  mcall_ <- mcall[3:length(mcall)]
  
  call_ <- model$model_info$call.full.bsitar #[-1]
  
  for (i in names(mcall_[-1])) {
    if(!i %in% names(call_)) {
      stop("Argument ", i, " is not a valid arguments",
           " \n ",
           " Please see 'bsitar' function ")
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
  
  
  for (i in names(mcall_)) {
    call_[[i]] <- mcall_[[i]]
  }
  
  iter <- call_$iter
  silent <- call_$silent
  
  
  as_one_logical <- NULL
  backend_choices <- backend_choices <- first_not_null <- as_one_logical
  get_drop_unused_levels <- is_equal <- is_normalized <- as_one_logical
  needs_recompilation <- validate_silent <- as_one_logical
  
  as_one_logical <- utils::getFromNamespace("as_one_logical", "brms")
  backend_choices <- utils::getFromNamespace("backend_choices", "brms")
  first_not_null <- utils::getFromNamespace("first_not_null", "brms")
  get_drop_unused_levels <- utils::getFromNamespace("get_drop_unused_levels", "brms")
  is_equal <- utils::getFromNamespace("is_equal", "brms")
  is_normalized <- utils::getFromNamespace("is_normalized", "brms")
  needs_recompilation <- utils::getFromNamespace("needs_recompilation", "brms")
  validate_silent <- utils::getFromNamespace("validate_silent", "brms")
  
  
  validate_threads <- NULL
  algorithm_choices <- validate_sample_prior <- validate_threads
  validate_save_pars <- validate_threads
  
  
  if (is.null(recompile)) {
    call_stancode <- call_
    call_stancode$get_stancode <- TRUE
    
    new_stancode <- eval.parent(call_stancode)
    new_stancode <- sub("^[^\n]+\n", "", new_stancode)
    
    old_stancode <- model$bmodel
    old_stancode <- sub("^[^\n]+\n", "", old_stancode)
    
    recompile <- needs_recompilation(model) || !same_backend || 
      !is_equal(new_stancode, old_stancode)
    
    if (recompile && silent < 2) {
      message("The desired updates require recompiling the model")
    }
  }
  
  
  recompile <- as_one_logical(recompile)
  
  
  if (recompile) {
    call_$fit <- NA
    model <- eval.parent(call_) 
  } else {
    call_$fit <- model
    model <- eval.parent(call_)
  }
  
  return(model)
}



#' @rdname update.bsitar
#' @export
update <- function(model, ...) {
  UseMethod("update")
}


