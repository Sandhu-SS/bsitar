


#' Update \pkg{bsitar} model
#'
#' @param model An object of class \code{bsitar}.
#' 
#' @details This is an adapted version of the \code{update} function of 
#' \code{brms} package.
#'
#' @param newdata An optional \code{data.frame} to be used when updating the
#'   model. If \code{NULL} (default), the same same data used for original model
#'   fit is used. Note that data-dependent default priors will not be updated
#'   automatically.
#'
#' @param recompile A logical to indicate whether the Stan model should be
#'   recompiled. When \code{NULL} (the default), \code{update_bsitar} tries to
#'   figure out internally, if recompilation is necessary. Setting it to
#'   \code{FALSE} will cause all Stan code changing arguments to be ignored.
#'
#' @param ... Other arguments passed to \code{\link{brms}}.
#'
#' @return An updated model of class \code{brmsfit, bsiatr}, that contains the
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
update_bsitar.bsitar <-
  function(model,
           newdata = NULL,
           recompile = NULL,
           ...) {
    formula. <- NULL
    args <- formalArgs(bsitar)
    args <- args[!args == "..."]
    
    call_ <- model$model_info$call.full.bsitar[-1] %>% as.list()
    
    call_$data <- NULL
    mcall_ <- list(...)
    
    if (length(mcall_) != 0) {
      for (i in names(mcall_)) {
        if (!i %in% args) {
          stop("Argument ",
               i,
               " is not a valid arguments",
               " \n ",
               " Please see 'bsitar' function ")
        } else {
          call_[[i]] <- mcall_[[i]]
        }
      }
    }
    
    
    dot_and_call_intersect <-
      intersect(names(list(...)), names(call_))
    
    exclude_args_names <- c(model$model_info[['brms_arguments_list']])
    
    exclude_args_names <-
      c(exclude_args_names, dot_and_call_intersect)
    
    if ("init" %in% dot_and_call_intersect)
      new_init_arg <- TRUE
    else
      new_init_arg <- FALSE
    
    
    for (ix in  exclude_args_names) {
      call_[[ix]] <- NULL
    }
    
    dots <- list(...)
    
    dots$data <- NULL
    
    # dots$select_model <- model$model_info[['select_model']]
    # dots$multivariate <- model$model_info[['multivariate']]
    
    
    as_one_logical         <-
      is_equal <- needs_recompilation <- substitute_name <- NULL
    as_one_logical         <-
      utils::getFromNamespace("as_one_logical", "brms")
    is_equal               <-
      utils::getFromNamespace("is_equal", "brms")
    needs_recompilation    <-
      utils::getFromNamespace("needs_recompilation", "brms")
    substitute_name        <-
      utils::getFromNamespace("substitute_name", "brms")
    get_drop_unused_levels <-
      utils::getFromNamespace("get_drop_unused_levels", "brms")
    validate_data          <-
      utils::getFromNamespace("validate_data", "brms")
    get_data_name          <-
      utils::getFromNamespace("get_data_name", "brms")
    validate_formula       <-
      utils::getFromNamespace("validate_formula", "brms")
    get_knots              <-
      utils::getFromNamespace("get_knots", "brms")
    is_normalized          <-
      utils::getFromNamespace("is_normalized", "brms")
    first_not_null         <-
      utils::getFromNamespace("first_not_null", "brms")
    backend_choices        <-
      utils::getFromNamespace("backend_choices", "brms")
    validate_data2         <-
      utils::getFromNamespace("validate_data2", "brms")
    .validate_prior        <-
      utils::getFromNamespace(".validate_prior", "brms")
    get_element            <-
      utils::getFromNamespace("get_element", "brms")
    tidy_ranef             <-
      utils::getFromNamespace("tidy_ranef", "brms")
    validate_stanvars      <-
      utils::getFromNamespace("validate_stanvars", "brms")
    validate_threads       <-
      utils::getFromNamespace("validate_threads", "brms")
    validate_sample_prior  <-
      utils::getFromNamespace("validate_sample_prior", "brms")
    validate_save_pars     <-
      utils::getFromNamespace("validate_save_pars", "brms")
    standata_basis         <-
      utils::getFromNamespace("standata_basis", "brms")
    algorithm_choices      <-
      utils::getFromNamespace("algorithm_choices", "brms")
    get_nl                 <-
      utils::getFromNamespace("get_nl", "brms")
    get_arg                <-
      utils::getFromNamespace("get_arg", "brms")
    is_nonlinear           <-
      utils::getFromNamespace("is_nonlinear", "brms")
    subset2                <-
      utils::getFromNamespace("subset2", "brms")
    rcols_prior            <-
      utils::getFromNamespace("rcols_prior", "brms")
    stop2                  <- utils::getFromNamespace("stop2", "brms")
    
    validate_silent        <- utils::getFromNamespace("validate_silent", "brms")
    
    
    
    
    
    testmode <- isTRUE(dots[["testmode"]])
    dots$testmode <- NULL
    if ("silent" %in% names(dots)) {
      dots$silent <- validate_silent(dots$silent)
    } else {
      dots$silent <- model$stan_args$silent %||% 1L
    }
    silent <- dots$silent
    model <- brms::restructure(model)
    if (isTRUE(model$version$brms < "2.0.0")) {
      warning2("Updating models fitted with older versions of brms may fail.")
    }
    model$file <- NULL
    
    if ("data" %in% names(dots)) {
      # otherwise the data name cannot be found by substitute
      stop2("Please use argument 'newdata' to update the data.")
    }
    if (!is.null(newdata)) {
      dots$data <- newdata
      data_name <- substitute_name(newdata)
    } else {
      dots$data <- model$data
      data_name <- get_data_name(model$data)
    }
    
    if (missing(formula.) || is.null(formula.)) {
      dots$formula <- model$formula
      if (!is.null(dots[["family"]])) {
        dots$formula <- bf(dots$formula, family = dots$family)
      }
      if (!is.null(dots[["autocor"]])) {
        dots$formula <- bf(dots$formula, autocor = dots$autocor)
      }
    } else {
      # TODO: restructure updating of the formula
      if (is.mvbrmsformula(formula.) ||
          is.mvbrmsformula(model$formula)) {
        stop2("Updating formulas of multivariate models is not yet possible.")
      }
      if (is.brmsformula(formula.)) {
        nl <- get_nl(formula.)
      } else {
        formula. <- as.formula(formula.)
        nl <- get_nl(formula(model))
      }
      family <- get_arg("family", formula., dots, model)
      autocor <- get_arg("autocor", formula., dots, model)
      dots$formula <-
        bf(formula.,
           family = family,
           autocor = autocor,
           nl = nl)
      if (is_nonlinear(model)) {
        # if (length(setdiff(all.vars(dots$formula$formula), ".")) == 0L) {
        #   dots$formula <- update(model$formula, dots$formula, mode = "keep")
        # } else {
        #   dots$formula <- update(model$formula,dots$formula, mode = "replace")
        #   if (silent < 2) {
        #     message("Argument 'formula.' will completely replace the ",
        #             "original formula in non-linear models.")
        #   }
        # }
      } else {
        mvars <- all.vars(dots$formula$formula)
        mvars <- setdiff(mvars, c(names(model$data), "."))
        if (length(mvars) && is.null(newdata)) {
          stop2(
            "New variables found: ",
            collapse_comma(mvars),
            "\nPlease supply your data again via argument 'newdata'."
          )
        }
        dots$formula <- update(formula(model), dots$formula)
      }
    }
    # update response categories and ordinal thresholds
    dots$formula <- validate_formula(dots$formula, data = dots$data)
    
    if (is.null(dots$prior)) {
      dots$prior <- model$prior
    } else {
      if (!is.brmsprior(dots$prior)) {
        stop2("Argument 'prior' needs to be a 'brmsprior' model.")
      }
      # update existing priors manually and keep only user-specified ones
      # default priors are recomputed base on newdata if provided
      
      
      # old_user_prior <- subset2(model$prior, source = "user")
      # dots$prior <- rbind(dots$prior, old_user_prior)
      # dupl_priors <- duplicated(dots$prior[, rcols_prior()])
      # dots$prior <- dots$prior[!dupl_priors, ]
    }
    # make sure potentially updated priors pass 'validate_prior'
    attr(dots$prior, "allow_invalid_prior") <- TRUE
    if (!"sample_prior" %in% names(dots)) {
      dots$sample_prior <- attr(model$prior, "sample_prior")
      if (is.null(dots$sample_prior)) {
        has_prior_pars <- any(grepl("^prior_", variables(model)))
        dots$sample_prior <- if (has_prior_pars)
          "yes"
        else
          "no"
      }
    }
    # do not use 'is.null' to allow updating arguments to NULL
    if (!"data2" %in% names(dots)) {
      dots$data2 <- model$data2
    }
    if (!"stanvars" %in% names(dots)) {
      dots$stanvars <- model$stanvars
    }
    if (!"algorithm" %in% names(dots)) {
      dots$algorithm <- model$algorithm
    }
    if (!"backend" %in% names(dots)) {
      dots$backend <- model$backend
    }
    if (!"threads" %in% names(dots)) {
      dots$threads <- model$threads
    }
    if (!"save_pars" %in% names(dots)) {
      dots$save_pars <- model$save_pars
    }
    if (!"knots" %in% names(dots)) {
      dots$knots <- get_knots(model$data)
    }
    if (!"drop_unused_levels" %in% names(dots)) {
      dots$drop_unused_levels <- get_drop_unused_levels(model$data)
    }
    if (!"normalize" %in% names(dots)) {
      dots$normalize <- is_normalized(model$model)
    }
    
    # update arguments controlling the sampling process
    if (is.null(dots$iter)) {
      # only keep old 'warmup' if also keeping old 'iter'
      dots$warmup <- first_not_null(dots$warmup, model$fit@sim$warmup)
    }
    dots$iter <- first_not_null(dots$iter, model$fit@sim$iter)
    dots$chains <- first_not_null(dots$chains, model$fit@sim$chains)
    dots$thin <- first_not_null(dots$thin, model$fit@sim$thin)
    dots$backend <- match.arg(dots$backend, backend_choices())
    same_backend <- is_equal(dots$backend, model$backend)
    if (same_backend) {
      # reusing control arguments in other backends may cause errors #1259
      control <- attr(model$fit@sim$samples[[1]], "args")$control
      control <- control[setdiff(names(control), names(dots$control))]
      dots$control[names(control)] <- control
      # reuse backend arguments originally passed to brm #1373
      names_old_stan_args <-
        setdiff(names(model$stan_args), names(dots))
      dots[names_old_stan_args] <-
        model$stan_args[names_old_stan_args]
    }
    
    if (is.null(recompile)) {
      # only recompile if new and old stan code do not match
      
      # new_stancode <- suppressMessages(do_call(make_stancode, dots))
      
      dots_for_scode              <- dots
      dots_for_scode$prior        <- NULL
      dots_for_scode$stanvars     <- NULL
      dots_for_scode$formula      <- NULL
      if (!new_init_arg)
        dots_for_scode$init         <- NULL
      dots_for_scode              <- c(dots_for_scode, call_)
      dots_for_scode$get_stancode <- TRUE
      new_stancode <-
        suppressMessages(do.call(bsitar, dots_for_scode))
      
      
      # stan code may differ just because of the version number (#288)
      new_stancode <- sub("^[^\n]+\n", "", new_stancode)
      old_stancode <- brms::stancode(model, version = FALSE)
      recompile <- needs_recompilation(model) || !same_backend ||
        !is_equal(new_stancode, old_stancode)
      if (recompile && silent < 2) {
        message("The desired updates require recompiling")
      }
    }
    recompile <- as_one_logical(recompile)
    if (recompile) {
      # recompliation is necessary
      dots$fit <- NA
      if (!testmode) {
        dots_for_recompile          <- dots
        dots_for_recompile$prior    <- NULL
        dots_for_recompile$stanvars <- NULL
        dots_for_recompile$formula  <- NULL
        if (!new_init_arg)
          dots_for_recompile$init     <- NULL
        dots_for_recompile          <- c(dots_for_recompile, call_)
        model <- do.call(bsitar, dots_for_recompile)
      }
    } else {
      # refit without compiling it again
      if (!is.null(dots$formula)) {
        model$formula <- dots$formula
        dots$formula <- NULL
      }
      bterms <- brms::brmsterms(model$formula)
      model$data2 <- validate_data2(dots$data2, bterms = bterms)
      model$data <- validate_data(
        dots$data,
        bterms = bterms,
        data2 = model$data2,
        knots = dots$knots,
        drop_unused_levels = dots$drop_unused_levels
      )
      model$prior <- .validate_prior(
        dots$prior,
        bterms = bterms,
        data = model$data,
        sample_prior = dots$sample_prior
      )
      model$family <- get_element(model$formula, "family")
      model$autocor <- get_element(model$formula, "autocor")
      model$ranef <- tidy_ranef(bterms, data = model$data)
      model$stanvars <- validate_stanvars(dots$stanvars)
      model$threads <- validate_threads(dots$threads)
      if ("sample_prior" %in% names(dots)) {
        dots$sample_prior <- validate_sample_prior(dots$sample_prior)
        attr(model$prior, "sample_prior") <- dots$sample_prior
      }
      model$save_pars <- validate_save_pars(
        save_pars = dots$save_pars,
        save_ranef = dots$save_ranef,
        save_mevars = dots$save_mevars,
        save_all_pars = dots$save_all_pars
      )
      model$basis <- standata_basis(bterms, data = model$data)
      algorithm <- match.arg(dots$algorithm, algorithm_choices())
      dots$algorithm <- model$algorithm <- algorithm
      # can only avoid recompilation when using the old backend
      dots$backend <- model$backend
      if (!testmode) {
        dots$fit <- model
        dots$file_refit <- NULL
        dots_for_norecompile          <- dots
        dots_for_norecompile$prior    <- NULL
        dots_for_norecompile$stanvars <- NULL
        dots_for_norecompile$formula  <- NULL
        if (!new_init_arg | new_init_arg) {
          dots_for_norecompile$init     <- NULL
          dots_for_norecompile          <-
            c(dots_for_norecompile, call_)
          model <- do.call(bsitar, dots_for_norecompile)
        } # if(!new_init_arg) {
        if (new_init_arg) {
          # set_brms_args            <- dots
          # set_brms_args$prior      <- model$prior
          # set_brms_args$stanvars   <- model$stanvars
          # set_brms_args$formula    <- model$formula
          # # set_brms_args$model <- NULL
          # # set_brms_args$newdata <- NULL
          # set_brms_args$data <- newdata
          # set_brms_args$init <- random
          # set_brms_args          <- c(set_brms_args)
          # model <- do.call(brm, set_brms_args)
        } # if(new_init_arg) {
      } # if (!testmode) {
    }
    attr(model$data, "data_name") <- data_name
    model
  }

#' @rdname update_bsitar.bsitar
#' @export
update_bsitar <- function(model, ...) {
  UseMethod("update_bsitar")
}
