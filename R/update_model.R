


#' @title Update model
#'
#' @description The \strong{update_model()} function is a wrapper around the
#'   \code{update()} function from the \pkg{brms} package, which refits the model
#'   based on the user-specified updated arguments.
#' 
#' @details This function is an adapted version of the \strong{update()} function 
#'   from the \pkg{brms} package.
#' 
#' @param model An object of class \code{bgmfit}.
#'
#' @param newdata An optional \code{data.frame} to be used when updating the
#'   model. If \code{NULL} (default), the data used in the original model fit is
#'   reused. Note that data-dependent default priors are not automatically updated.
#'
#' @param recompile A logical value indicating whether the Stan model should be
#'   recompiled. When \code{NULL} (default), \strong{update_model()} tries to
#'   internally determine whether recompilation is required. Setting
#'   \code{recompile} to \code{FALSE} will ignore any changes in the Stan code.
#'   
#' @param check_newargs A logical value (default \code{FALSE}) indicating whether
#'   to check if the arguments in the original \code{model} fit and the 
#'   \code{update_model} are identical. When \code{check_newargs = TRUE} and the 
#'   arguments are identical, it indicates that an update is unnecessary. In this 
#'   case, the original \code{model} object is returned, along with a message if 
#'   \code{verbose = TRUE}.
#' 
#' @inherit growthparameters.bgmfit params
#'
#' @param ... Other arguments passed to \code{[brms::brm()]}.
#'
#' @return An updated object of class \code{brmsfit}.
#'   
#' @rdname update_model
#' @export
#'
#' @inherit berkeley author
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
#' # the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' # Check and confirm whether model fit object 'berkeley_exfit' exists
#'  berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' # Update model
#' # Note that in case all arguments supplied to the update_model() call are 
#' # same as the original model fit (checked via check_newargs = TRUE), then  
#' # original model object is returned.   
#' # To explicitly get this information whether model is being updated or not, 
#' # user can set verbose = TRUE. The verbose = TRUE also useful in getting the
#' # information regarding what all arguments have been changed as compared to
#' # the original model.
#' 
#' model2 <- update_model(model, df = 5, check_newargs = TRUE, verbose = TRUE)
#' 
#' }
#'
update_model.bgmfit <-
  function(model,
           newdata = NULL,
           recompile = NULL,
           expose_function = FALSE,
           verbose = FALSE,
           check_newargs = FALSE,
           envir = NULL,
           ...) {
    # changes on 21.07.2024 for brms version ‘2.21.6’
    # deletions/insertions -> tidy_ranef standata_basis 
    # insertion -> bframe
    # edited -> model$prior
    if(is.null(envir)) {
      envir <- model$model_info$envir
    } else {
      envir <- envir
    }
    

    if(check_newargs) {
      call_o <- match.call()
      call_o_args <- as.list(call_o)[-1]
      
      args_o <- as.list(model$model_info$call.full.bgmfit)[-1]
      
      args_o_dots_ <- list(...)
      if (length(args_o_dots_) > 0) {
        for (i in names(args_o_dots_)) {
          args_o[[i]] <- args_o_dots_[[i]]
        }
      }
      
      # This to evaluate T/F to TRUE/FALSE
      for (i in names(args_o)) {
        if (is.symbol(args_o[[i]])) {
          if (args_o[[i]] == "T")
            args_o[[i]] <- eval(args_o[[i]])
          if (args_o[[i]] == "F")
            args_o[[i]] <- eval(args_o[[i]])
        }
      }
      
      args_o_new <- args_o_dots_
      args_o_new[['expose_function']] <- expose_function
      
      calling    <- model$model_info$call.full.bgmfit
      calling[['verbose']] <- NULL
      
      args_o_org <- calling
      args_o_new$data <- NULL
      args_o_org$data <- NULL
      
      # This to evaluate T/F to TRUE/FALSE
      for (i in names(args_o_org)) {
        if (is.symbol(args_o_org[[i]])) {
          if (args_o_org[[i]] == "T")
            args_o_org[[i]] <- eval(args_o_org[[i]])
          if (args_o_org[[i]] == "F")
            args_o_org[[i]] <- eval(args_o_org[[i]])
        }
      }
      
      all_same_args_c <- all_same_args <- c()
      # args_o_org_updated <- list()
      for (args_oi in names(args_o_new)) {
        all_same_args_c <- c(all_same_args_c, identical(args_o_org[[args_oi]],
                                                        args_o_new[[args_oi]]) 
        )
      }
      
      all_same_args_c_diffs <- args_o_new[!all_same_args_c]
      if(length(all_same_args_c_diffs) > 0) {
        all_same_args <- FALSE 
      } else {
        all_same_args <- TRUE
      }
      
      if(all_same_args) {
        if(verbose) {
          cat("\n")
          message("Arguemnets supplied for 'update_model()' call are same as ",
              "the original model fit.", 
              "\n ",
              "Therefore, returning the original model object")
          cat("\n")
        }
      }
      return(model)
    } # if(check_newargs) {
    
    
    
    
    check_if_package_installed(model, xcall = NULL)
    
    formula. <- NULL
    args <- formalArgs(bsitar)
    args <- args[!args == "..."]
    
    call_ <- model$model_info$call.full.bgmfit[-1] %>% as.list()
    
    call_$data <- NULL
    mcall_ <- list(...)
    
    if (length(mcall_) != 0) {
      for (i in names(mcall_)) {
        if (!i %in% args) {
          stop("Argument ",
               i,
               " is not a valid arguments",
               " \n ",
               " Please see the main calling function ")
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
    
    
    
    
    as_one_logical <- is_equal <- NULL
    needs_recompilation <- substitute_name <- NULL
    
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
    # tidy_ranef             <-
    #   utils::getFromNamespace("tidy_ranef", "brms")
    getframe_re      <-
      utils::getFromNamespace("frame_re", "brms")
    validate_stanvars      <-
      utils::getFromNamespace("validate_stanvars", "brms")
    validate_threads       <-
      utils::getFromNamespace("validate_threads", "brms")
    validate_sample_prior  <-
      utils::getFromNamespace("validate_sample_prior", "brms")
    validate_save_pars     <-
      utils::getFromNamespace("validate_save_pars", "brms")
    # standata_basis         <-
    #   utils::getFromNamespace("standata_basis", "brms")
    getframe_basis      <-
      utils::getFromNamespace("frame_basis", "brms")
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
    
    getbrmsframe        <- utils::getFromNamespace("brmsframe", "brms")
    
    
    testmode <- isTRUE(dots[["testmode"]])
    dots$testmode <- NULL
    if ("silent" %in% names(dots)) {
      dots$silent <- validate_silent(dots$silent)
    } else {
      dots$silent <- model$stan_args$silent %||% 1L
    }
    silent <- dots$silent
    model <- brms::restructure(model)
    # if (isTRUE(model$version$brms < "2.0.0")) {
    #   warning2("Updating models fitted with older versions of brms may fail.")
    # }
    model$file <- NULL
    
    if ("data" %in% names(dots)) {
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
    dots$formula <- validate_formula(dots$formula, data = dots$data)
    
    if (is.null(dots$prior)) {
      dots$prior <- model$prior
    } else {
      if (!is.brmsprior(dots$prior)) {
        stop2("Argument 'prior' needs to be a 'brmsprior' model.")
      }
     
    }
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
    
    if (is.null(dots$iter)) {
      dots$warmup <- first_not_null(dots$warmup, model$fit@sim$warmup)
    }
    dots$iter <- first_not_null(dots$iter, model$fit@sim$iter)
    dots$chains <- first_not_null(dots$chains, model$fit@sim$chains)
    dots$thin <- first_not_null(dots$thin, model$fit@sim$thin)
    dots$backend <- match.arg(dots$backend, backend_choices())
    same_backend <- is_equal(dots$backend, model$backend)
    if (same_backend) {
      control <- attr(model$fit@sim$samples[[1]], "args")$control
      control <- control[setdiff(names(control), names(dots$control))]
      dots$control[names(control)] <- control
      names_old_stan_args <-
        setdiff(names(model$stan_args), names(dots))
      dots[names_old_stan_args] <-
        model$stan_args[names_old_stan_args]
    }
    
    
    
    if (is.null(recompile)) {
      dots_for_scode              <- dots
      dots_for_scode$prior        <- NULL
      dots_for_scode$stanvars     <- NULL
      dots_for_scode$formula      <- NULL
      if (!new_init_arg)
        dots_for_scode$init         <- NULL
      dots_for_scode              <- c(dots_for_scode, call_)
      dots_for_scode$get_stancode <- TRUE
      new_stancode <- suppressMessages(do.call(bsitar, dots_for_scode))
      # new_stancode <- suppressMessages(CustomDoCall(bsitar, dots_for_scode))
      
      
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
        # model <- CustomDoCall(bsitar, dots_for_recompile)
      }
    } else {
      if (!is.null(dots$formula)) {
        model$formula <- dots$formula
        dots$formula <- NULL
      }
      bterms <- brms::brmsterms(model$formula)
      bframe <- getbrmsframe(bterms, data = model$data)
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
        # bterms = bterms,
        bframe = bframe,
        # data = model$data,
        sample_prior = dots$sample_prior
      )
      model$family <- get_element(model$formula, "family")
      model$autocor <- get_element(model$formula, "autocor")
      # model$ranef <- tidy_ranef(bterms, data = model$data)
      model$ranef <- getframe_re(bterms, data = model$data)
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
      # model$basis <- standata_basis(bterms, data = model$data)
      model$basis <- getframe_basis(bframe, data = model$data)
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
          # model <- CustomDoCall(bsitar, dots_for_norecompile)
        } # if(!new_init_arg) {
        if (new_init_arg) {
          # TODO
        } # if(new_init_arg) {
      } # if (!testmode) {
    }
    if(expose_function) model <- expose_model_functions(model, envir = envir)
    attr(model$data, "data_name") <- data_name
    model
  }



#' @rdname update_model
#' @export
update_model <- function(model, ...) {
  UseMethod("update_model")
}
