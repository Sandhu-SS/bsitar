

#' @title Return Model-Fit Criteria for a Bayesian SITAR Model
#'
#' @description
#' \code{get_model_criterion()} is a wrapper around [add_model_criterion()] that
#' computes and returns model-fit criteria. See [add_model_criterion()] for
#' details and available arguments.
#' 
#' @param reformat Logical indicating whether to apply [base::round()] to the
#'   numeric variables in the output \code{data.frame}. The default is
#'   \code{NULL}, which is treated as \code{TRUE}. If \code{TRUE}, then numeric
#'   variables are rounded using \code{digits}.
#' 
#' @param add_attr A logical indicating whether complex list elements of
#'   the criteria should be added as attributes to the returned
#'   \code{data.frame}. Defaults to \code{FALSE}.
#'
#' @param ... Additional arguments passed to [add_model_criterion()].
#' 
#' @inheritParams compare_models.bgmfit
#' @inheritParams add_model_criterion.bgmfit
#' @inheritParams growthparameters.bgmfit
#' @inheritParams brms::bayes_R2.brmsfit
#' @inheritParams brms::add_criterion.brmsfit
#' @inheritParams brms::waic.brmsfit
#' @inheritParams fitted_draws.bgmfit
#' 
#' @return A \code{data.frame}. If \code{add_attr = TRUE}, an additional
#'   attribute named \code{"attr_object"} is attached to the returned data
#'   frame. This attribute contains the complex, nested list components from
#'   \code{add_model_criterion()} that cannot be represented as regular data
#'   frame columns.
#' 
#' @rdname get_model_criterion
#' @export
#' 
#' @seealso [brms::add_loo], [brms::add_ic()], [brms::add_waic()],
#'   [brms::bayes_R2()]
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
#' # For illustration purposes, we make a copy of model with itself
#' # In the example below, get_model_criterion() should indicate no difference 
#' # between model_1 and model_2 as both these models are exactly identical
#' model_1 <- model
#' model_2 <- model
#' 
#' # Add model fit criteria (e.g., WAIC). 
#' out_1 <- get_model_criterion(model_1, criterion = c("waic"))
#' out_2 <- get_model_criterion(model_2, criterion = c("loo"))
#' 
#' # compare models model_1 and model_2
#' out_12 <- get_model_criterion(model_1, model_2, criterion = c("waic"))
#' 
#' # compare models model_1 and model_2 - model names: "mods[[1L]]" "mods[[2L]]"
#' mods <- list(model_1, model_2)
#' out_12 <- get_model_criterion(mods)
#' 
#' # compare models model_1 and model_2 - model names extracted as such
#' # Note that list() could be supplied as a string "list()" also.
#' out_12 <- get_model_criterion(list(model_1, model_2))
#' out_12 <- get_model_criterion("list(model_1, model_2)")
#' 
#' }
#' 
get_model_criterion.bgmfit <- function(model,
                                       ...,
                                       criterion = "loo",
                                       ndraws = NULL,
                                       draw_ids = NULL,
                                       pointwise = FALSE,
                                       model_name = NULL,
                                       summary = TRUE,
                                       robust = FALSE,
                                       probs = c(0.025, 0.975),
                                       newdata = NULL,
                                       resp = NULL,
                                       cores = 1,
                                       expose_function = FALSE, 
                                       verbose = FALSE,
                                       reformat = NULL,
                                       digits = 2,
                                       add_attr = FALSE) {

  only_object <- NULL
  if (is.character(model) && length(model) == 1 && 
      grepl("^\\s*list\\s*\\(", model)) {
    only_object <- FALSE
  } 
  if(is.list(model)) {
    only_object <- FALSE
  } 
  if(is.bgmfit(model)) {
    only_object <- TRUE
  }
  if(length(c(list(model), list(...))) > 1) {
    only_object <- FALSE
  }
  
  if (is.character(model) && length(model) == 1 && 
      grepl("^\\s*list\\s*\\(", model)) {
    model_list_str <-  model
    model <- eval(parse(text = model), envir = parent.frame())
  } else {
    model_list_str <-  deparse(substitute(model))
  }
  gsub_namespace <- FALSE
  model_list_names <- NULL
  if(grepl("list\\(", model_list_str)) {
    model_list_str <- gsub("list", "", model_list_str)
    model_list_str <- strsplit(model_list_str, ",")[[1]]
    # replace :: with _ , not when :::
    model_list_str <- gsub("[^A-Za-z0-9_:]", "", model_list_str)
    if(gsub_namespace) {
      tmp <- gsub(":::", "@@@COLONPAIR@@@", model_list_str)
      tmp <- gsub("::", "_", tmp)
      model_list_str <- gsub("@@@COLONPAIR@@@", ":::", tmp)
    }
    model_list_names <- model_list_str
  }

  if (is.null(model_name)) {
    if(!is.null(model_list_names)) model_names <- model_list_names
    if( is.null(model_list_names)) model_names <- NULL
  } else {
    model_names <- model_name
  }
  
  add_args <- as.list(match.call(expand.dots = FALSE))
  defaults_it <- base::as.list(base::formals(get_model_criterion))
  for (i in names(defaults_it)) {
    if(is.null(add_args[[i]])) add_args[[i]] <- defaults_it[[i]]
  }
  
  defaults <- base::as.list(base::formals(get_model_criterion.bgmfit))
  defaults[['model']] <- NULL
  build_args <- utils::modifyList(defaults, add_args)

  check_criterion <- TRUE
  add_criterion_args                      <- build_args
  add_criterion_args[["model_name"]]      <- NULL
  add_criterion_args[["check_criterion"]] <- FALSE
  add_criterion_args[["compare"]]         <- FALSE
  add_criterion_args[["return_criteria"]] <- TRUE
  add_criterion_args[["return_model"]]    <- FALSE
  add_criterion_args[["add_attr"]]         <- NULL
  
  add_criterion_args <- move_to_front_list(add_criterion_args,c("model", "..."))
  
  add_criterion_args[['model']] <- NULL
  add_criterion_args[['...']] <- NULL
  add_criterion_args[['compare']] <- FALSE
  
  if (!is.list(add_criterion_args)) {
    stop("Argument 'add_criterion_args' must be a named list")
  }
  
  exprs <- as.list(substitute(list(model, ...)))[-1]
  vals <- c(list(model), list(...))
  models <- unlist(lapply(vals, flatten_models), recursive = FALSE)
  if (is.null(model_names)) {
    model_names <- unlist(Map(flatten_exprs, exprs, vals), 
                          use.names = FALSE)
  }
  
  only_one_model <- FALSE
  if (length(models) < 2) {
    only_one_model <- TRUE
  }
  
  only_one_criterion <- FALSE
  if (length(criterion) < 2) {
    only_one_criterion <- TRUE
  }
  
  if (!check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion_multiple(fit, criterion)) {
        stop("No precomputed criterion availabel for one or more models.", 
             " Either add criterion before hand using 'add_model_criterion()'", 
             " or else set check_criterion = TRUE. Note that arguments to", 
             " 'add_model_criterion()' function can be set by using", 
             " 'add_criterion_args' which must be a named list")
      }
    })
  }
  
  if (check_criterion) {
    models <- lapply(models, function(fit) {
      if (!has_criterion_multiple(fit, criterion)) {
        fit <- do.call(expose_model_functions, c(list(model = fit, 
                                                              expose = expose_function) ))
        suppressWarnings({
          fit <- do.call(add_model_criterion, c(list(model = fit), add_criterion_args))
        })
      }
      fit
    })
  }
  
  if (length(model_names) != length(models)) {
    nnames <- length(model_names)
    model_names_all <- paste0("model", seq_along(models) - 
                                length(model_names))
    model_names_all[1:nnames] <- model_names
    model_names <- model_names_all
    message2c("The number of model names is not same as the number of models.\n              The remaining models are named sequentially as model1,...")
  }

  out <- nested_to_df(models, model_names = model_names, add_attr = T)
  
  if(is.null(reformat)) {
    reformat <- TRUE
  }
  
  if(reformat) {
    out <- out %>% 
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), 
                                  ~ round(.x, digits = digits)))
  }
  
  return(out)
}




#' @rdname get_model_criterion
#' @export
get_model_criterion <- function(model, ...) {
  UseMethod("get_model_criterion")
}


#' @rdname get_model_criterion
#' @export
get_model_criterion.list <- function(model, ...) {
  get_model_criterion.bgmfit(model, ...)
}


#' @rdname get_model_criterion
#' @export
get_model_criterion.character <- function(model, ...) {
  get_model_criterion.bgmfit(model, ...)
}


#' @rdname get_model_criterion
#' @export
get_model_criterion.default <- function(model, ...) {
  if(!inherits(model, 'bgmfit') & 
     !inherits(model, 'list') &
     !inherits(model, 'character'))
  stop(
    "`model` must be an object of class 'bgmfit', a list, or a string",
    call. = FALSE
  )
}



