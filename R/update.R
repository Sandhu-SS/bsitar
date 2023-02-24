


#' Update \pkg{bsitar} models
#'
#' @param model An object of class \code{bsitar}.
#' @param newdata newdata Optional \code{data.frame} to update the model with
#'   new data. Note that data-dependent default priors will not be updated
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
update.bsitar <- function(model, newdata = NULL, recompile = T, ...) {
  
  mcall_ <- match.call()
  
  call_ <- model$model_info$call.bsitar#[-1]
  
  arguments <- as.list(mcall_)[-1]
  
  arguments <- arguments[-1]

  for (i in names(arguments)) {
    call_[[i]] <- arguments[[i]]
  }
  
  args_all_ <- formalArgs(bsitar)
  args_all_ <- args_all_[!args_all_ %in% "..."]
  
  for (i in names(call_[-1])) {
    if(!i %in% args_all_) {
      # stop("Argument ", i, " is not a valid arguments",
      #      "\ n",
      #      " valid arguments are: ", args_all_)
      stop("Argument ", i, " is not a valid arguments",
           " \n ",
           " Please see 'bsitar' function ")
    }
  }
  
  # call_ <- call_
  
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
  
  
  ############################
  
  # do not use 'is.null' to allow updating arguments to NULL
  if (!"data2" %in% names(arguments)) {
    arguments$data2 <- model$data2
  }
  if (!"stanvars" %in% names(arguments)) {
    arguments$stanvars <- model$stanvars
  }
  if (!"algorithm" %in% names(arguments)) {
    arguments$algorithm <- model$algorithm
  }
  if (!"backend" %in% names(arguments)) {
    arguments$backend <- model$backend
  }
  if (!"threads" %in% names(arguments)) {
    arguments$threads <- model$threads
  }
  if (!"save_pars" %in% names(arguments)) {
    arguments$save_pars <- model$save_pars
  }
 
  if (!"drop_unused_levels" %in% names(arguments)) {
    #brms:::
    arguments$drop_unused_levels <- get_drop_unused_levels(model$data)
  }
  if (!"normalize" %in% names(arguments)) {
    #brms:::
    arguments$normalize <- is_normalized(model$model)
  }
  
  # update arguments controlling the sampling process
  if (is.null(arguments$iter)) {
    # only keep old 'warmup' if also keeping old 'iter'
    arguments$warmup <- first_not_null(arguments$warmup, model$fit@sim$warmup)
  }
  #brms:::
  arguments$iter <- first_not_null(arguments$iter, model$fit@sim$iter)
  arguments$chains <- first_not_null(arguments$chains, model$fit@sim$chains)
  arguments$thin <- first_not_null(arguments$thin, model$fit@sim$thin)
  arguments$backend <- match.arg(arguments$backend, backend_choices())
  same_backend <- is_equal(arguments$backend, model$backend)
  if (same_backend) {
    # reusing control arguments in other backends may cause errors #1259
    control <- attr(model$fit@sim$samples[[1]], "args")$control
    control <- control[setdiff(names(control), names(arguments$control))]
    arguments$control[names(control)] <- control
    # reuse backend arguments originally passed to brm #1373
    names_old_stan_args <- setdiff(names(model$stan_args), names(arguments))
    arguments[names_old_stan_args] <- model$stan_args[names_old_stan_args]
  }
  
  
  ##########
  arguments$formula <- model$formula
  if(is.null(newdata)) {
    arguments$data <- model$data
  } else {
    arguments$data <- newdata
  }
  
  arguments$prior <- model$prior
  
  arguments$formula <- model$formula
  
  
  if ("silent" %in% names(arguments)) {
    #brms:::
    arguments$silent <- validate_silent(arguments$silent)
  } else {
    arguments$silent <- model$stan_args$silent %||% 1L
  }
  silent <- arguments$silent
  
  ##########
  
  if (is.null(recompile)) {
    # only recompile if new and old stan code do not match
    new_stancode <- suppressMessages(do_call(make_stancode, arguments))
    # stan code may differ just because of the version number (#288)
    new_stancode <- sub("^[^\n]+\n", "", new_stancode)
    old_stancode <- stancode(model, version = FALSE)
    #brms:::
    recompile <- needs_recompilation(model) || !same_backend ||
      !is_equal(new_stancode, old_stancode)
    if (recompile && silent < 2) {
      message("The desired updates require recompiling the model")
    }
  }
  #brms:::
  recompile <- as_one_logical(recompile)
  if (recompile) {
    # recompliation is necessary
    arguments$fit <- NA
    # if (!testmode) {
      model <- do_call(brm, arguments)
    # }
  } else {
    # refit the model without compiling it again
    # if (!is.null(arguments$formula)) {
    #   model$formula <- arguments$formula
    #   arguments$formula <- NULL
    # }
    # bterms <- brmsterms(model$formula)
    # model$data2 <- validate_data2(arguments$data2, bterms = bterms)
    # model$data <- validate_data(
    #   arguments$data, bterms = bterms, data2 = model$data2,
    #   knots = arguments$knots, drop_unused_levels = arguments$drop_unused_levels
    # )
    # model$prior <- .validate_prior(
    #   arguments$prior, bterms = bterms, data = model$data,
    #   sample_prior = arguments$sample_prior
    # )
    # model$family <- get_element(model$formula, "family")
    # model$autocor <- get_element(model$formula, "autocor")
    # model$ranef <- tidy_ranef(bterms, data = model$data)
    # model$stanvars <- validate_stanvars(arguments$stanvars)
    model$threads <- validate_threads(arguments$threads)
    if ("sample_prior" %in% names(arguments)) {
      arguments$sample_prior <- validate_sample_prior(arguments$sample_prior)
      attr(model$prior, "sample_prior") <- arguments$sample_prior
    }
    model$save_pars <- validate_save_pars(
      save_pars = arguments$save_pars,
      save_ranef = arguments$save_ranef,
      save_mevars = arguments$save_mevars,
      save_all_pars = arguments$save_all_pars
    )
    algorithm <- match.arg(arguments$algorithm, algorithm_choices())
    arguments$algorithm <- model$algorithm <- algorithm
    # can only avoid recompilation when using the old backend
    arguments$backend <- model$backend
    # if (!testmode) {
      arguments$fit <- model
      model <- do_call(brm, arguments)
    # }
  }
  
} # end


# ccc <- update_bsitar(fitx, cores = 2, x = age, a_init_beta = 2, y = y, xfun = NULL, df = 5, data = mmm)


#' @rdname update.bsitar
#' @export
update <- function(model, ...) {
  UseMethod("update")
}
