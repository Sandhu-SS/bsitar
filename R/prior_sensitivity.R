

#' Prior sensitivity analysis workflow for \pkg{bsitar} models
#'
#' Run a \pkg{priorsense} power-scaling sensitivity analysis for a fitted
#' \code{bsitar} model that returns a \code{brmsfit} object, or for a
#' \code{brmsfit} directly. Optionally, append derived quantities to a
#' posterior draws object in the same style as the \pkg{priorsense}
#' airquality workflow and immediately return a requested diagnostic plot.
#'
#' This function is intended as a package-level convenience wrapper around
#' [priorsense::powerscale_sensitivity()] and the associated plotting tools.
#' Since \pkg{priorsense} already supports \code{brmsfit} objects directly,
#' the default implementation uses the fitted \pkg{brms} object itself,
#' which is the most stable route for routine package use.
#'
#' If \code{add_draws = TRUE}, the function also constructs an augmented
#' draws object using posterior draws, log-likelihood, optional log-prior,
#' joint log-likelihood, optional Bayesian R-squared, and optional posterior
#' expected predictions for representative covariate values. This follows the
#' general workflow demonstrated in the \pkg{priorsense} vignette, where
#' derived quantities are bound to the draws object before sensitivity
#' analysis.
#'
#' @param model A fitted model object from \code{bsitar}.
#' 
#' @param variable Optional character vector of parameter or quantity names
#'   to assess. Passed to [priorsense::powerscale_sensitivity()] and, when
#'   plotting, to the corresponding \pkg{priorsense} plotting function.
#'   
#' @param component Character vector indicating which component(s) to assess.
#'   Typically one or both of \code{"prior"} and \code{"likelihood"}.
#'   Passed to [priorsense::powerscale_sensitivity()].
#'   
#' @param prior_selection Passed to [priorsense::powerscale_sensitivity()]
#'   for selecting prior terms to power-scale.
#'   
#' @param likelihood_selection Passed to
#'   [priorsense::powerscale_sensitivity()] for selecting likelihood terms
#'   to power-scale.
#'   
#' @param resample Logical; passed to [priorsense::powerscale_sensitivity()].
#' 
#' @param moment_match Logical; if \code{TRUE}, passed to
#'   [priorsense::powerscale_sensitivity()] to request moment matching when
#'   supported. This can improve unstable importance-sampling estimates in some
#'   cases, at additional computational cost.
#'   
#' @param add_draws Logical; if \code{TRUE}, build and store an augmented draws
#'   object with additional derived quantities. If \code{FALSE} (default), rely
#'   on the native \code{brmsfit} support in \pkg{priorsense}.
#'   
#' @param include_criterion Logical; if \code{TRUE}, append Bayesian R-squared
#'   draws using [brms::bayes_R2()].
#'   
#' @param include_jointlik Logical; if \code{TRUE} and \code{add_draws = TRUE},
#'   append a \code{joint log likelihood} quantity using
#'   [priorsense::predictions_as_draws()] with the customized
#'   \code{jointlik_fun} function.
#'   
#' @param include_predictor A character string or a named list to append a
#'   \code{expected values at predictor values} quantity using
#'   [priorsense::predictions_as_draws()] with the customized
#'   [brms::posterior_epred()].
#'    
#' @param name_predictor A character string or \code{NULL} (default). if
#'   \code{name_predictor = NULL}, then \code{name_predictor} set internally as
#'   \code{xvar}.
#' 
#' @param newdata Optional \code{data.frame} of representative predictor values
#'   used to compute posterior expected predictions with
#'   [brms::posterior_epred()] when \code{add_draws = TRUE}. Ignored.
#'
#' @param prediction_names Optional character vector of names for prediction
#'   quantities added from \code{newdata}. If \code{NULL}, names are created as
#'   \code{"pred_1"}, \code{"pred_2"}, and so on.
#'   
#' @param plot One of \code{NULL}, \code{FALSE}, \code{"dens"}, \code{"ecdf"},
#'   \code{"quantities"}, or \code{"both"}. If \code{NULL} or \code{FALSE}, no
#'   plot is produced. If a string is supplied, the corresponding
#'   \pkg{priorsense} plot is created and returned in the output object.
#'   
#' @param plot_variable Optional character string naming the variable to plot.
#'   If \code{NULL}, the first element of \code{variable} is used when
#'   available.
#'   
#' @param plot_print Optional logical indicating whether to print plot object.
#'
#' @param plot_return Optional logical indicating whether to print plot object.
#'
#' @param facet Optional value passed to \code{facet_rows} in
#'   [priorsense::powerscale_plot_dens()] and
#'   [priorsense::powerscale_plot_ecdf()].
#'   
#' @param return_ps Optional logical indicating whether to return a simple
#'   sensitivity object from the [priorsense::powerscale_sensitivity()] or a
#'   structured list with class attributes. Note that \code{return_ps} should be
#'   \code{FALSE} (default) in case user later want to call
#'   [prior_sensitivity_conflict()].
#'   
#' @param ... Additional arguments passed to
#'   [priorsense::powerscale_sensitivity()].
#'   
#' @inheritParams add_model_criterion
#' @inheritParams fitted_draws
#' @inheritParams priorsense::powerscale_sensitivity
#'
#' @details
#' The wrapper supports two closely related workflows.
#'
#' \strong{Default workflow (\code{add_draws = FALSE}):}
#' \pkg{priorsense} is run directly on the extracted \code{brmsfit}, which
#' is the most robust and package-friendly path because \pkg{brms} already
#' implements [brms::create_priorsense_data.brmsfit()] for this purpose.
#'
#' \strong{Extended workflow (\code{add_draws = TRUE}):}
#' a posterior draws object is created with [posterior::as_draws_df()] and
#' optionally extended with:
#' \itemize{
#'   \item pointwise log-likelihood draws via
#'         [priorsense::log_lik_draws()],
#'   \item log-prior draws via [priorsense::log_prior_draws()] if available,
#'   \item Bayesian R-squared via [brms::bayes_R2()],
#'   \item a \code{jointlik} quantity via \code{jointlik_fun},
#'   \item posterior expected predictions via [brms::posterior_epred()].
#' }
#'
#' This extended draws object can be useful when sensitivity is to be
#' explored not only for parameters but also for model-level summaries and
#' predictive quantities, consistent with the \pkg{priorsense} vignette
#' workflow.
#'
#' In the \pkg{priorsense} framework, sensitivity to both prior and
#' likelihood perturbations can indicate possible prior-data conflict,
#' whereas strong prior sensitivity with relatively weak likelihood
#' sensitivity can suggest that the likelihood is only weakly informative
#' for that quantity.
#'
#' @return
#' An object of class \code{"prior_sensitivity"}, a list with components:
#' \describe{
#'   \item{fit}{The extracted \code{bsitar} object.}
#'   \item{post_draws}{The augmented draws object if
#'     \code{add_draws = TRUE}, otherwise \code{NULL}.}
#'   \item{sensitivity}{The object returned by
#'     [priorsense::powerscale_sensitivity()].}
#'   \item{plot}{A plot object, or list of plot objects when
#'     \code{plot = "both"}, otherwise \code{NULL}.}
#'   \item{plot_type}{The requested plot type, if any.}
#'   \item{plot_variable}{The plotted variable, if any.}
#'   \item{newdata}{The supplied \code{newdata}.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @seealso
#' [priorsense::powerscale_sensitivity()],
#' [priorsense::powerscale_plot_dens()],
#' [priorsense::powerscale_plot_ecdf()],
#' [priorsense::powerscale_plot_quantities()],
#' [priorsense::predictions_as_draws()],
#' [prior_sensitivity_conflict()]
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
#' # Basic workflow: use brmsfit support directly
#' ps <- prior_sensitivity(
#'   model = model,
#'   variable = c("b_a_Intercept", "sigma")
#' )
#'
#' # Directly request a density plot
#' ps2 <- prior_sensitivity(
#'   model = model,
#'   variable = c("sigma"),
#'   plot = "dens"
#' )
#' ps2$plot
#'
#' # Extended workflow with derived predictive quantities
#' ps3 <- prior_sensitivity(
#'   model = model,
#'   variable = c("sigma"),
#'   criterion = c('bayes_R2', "jointlik"),
#'   include_predictor = c("min", "max"),
#'   plot = "quantities",
#'   plot_variable = c("sigma", "jointlik", "agemin", "agemax")
#' )
#' }
#' 
#' @rdname prior_sensitivity
#' @export
#' 
prior_sensitivity.bgmfit <- function(
    model,
    variable = NULL,
    lower_alpha = 0.99,
    upper_alpha = 1.01,
    div_measure = "cjs_dist",
    measure_args = list(),
    component = c("prior", "likelihood"),
    sensitivity_threshold = 0.05,
    moment_match = FALSE,
    k_threshold = 0.5,
    resample = FALSE,
    transform = NULL,
    prediction = NULL,
    prior_selection = NULL,
    likelihood_selection = NULL,
    num_args = NULL,
    newdata = NULL,
    add_draws = TRUE,
    include_criterion = FALSE,
    include_jointlik = TRUE,
    include_predictor = NULL,
    name_predictor = NULL,
    prediction_names = NULL,
    plot = NULL,
    plot_variable = NULL,
    plot_print = FALSE,
    plot_return = FALSE,
    facet = NULL,
    return_ps = FALSE,
    criterion = "bayes_R2",
    return_model = FALSE,
    return_criteria = TRUE,
    pointwise = FALSE,
    model_name = NULL,
    overwrite = FALSE,
    force_save = FALSE,
    file = NULL,
    resp = NULL,
    dpar = NULL,
    ndraws = NULL,
    draw_ids = NULL,
    re_formula = NULL,
    allow_new_levels = FALSE,
    sample_new_levels = "uncertainty",
    incl_autocor = TRUE,
    numeric_cov_at = NULL,
    levels_id = NULL,
    avg_reffects = NULL,
    aux_variables = NULL,
    grid_add = NULL,
    summary = FALSE,
    robust = FALSE,
    transform_draws = NULL,
    scale = c("response", "linear"),
    probs = c(0.025, 0.975),
    xrange = NULL,
    xrange_search = NULL,
    idata_method = NULL,
    verbose = FALSE,
    fullframe = NULL,
    dummy_to_factor = NULL, 
    expose_function = FALSE,
    usesavedfuns = NULL,
    clearenvfuns = NULL,
    funlist = NULL,
    xvar = NULL,
    idvar = NULL,
    itransform = NULL,
    envir = NULL,
    ...) {
  
  component <- match.arg(component, several.ok = TRUE)
  insight::check_if_installed('priorsense', prompt = FALSE)
  
  if(is.null(name_predictor)) {
    if(is.null(xvar)) {
      name_predictor <- get_basic_info(model, dpar = dpar, resp = resp)
    } else {
      name_predictor <- xvar
    }
  }
  
  if (!is.null(plot) && !identical(plot, FALSE)) {
    plot <- match.arg(plot, c("dens", "ecdf", "quantities", "both"))
  } else {
    plot <- NULL
  }
  
  if (is.null(plot_variable) && !is.null(variable) && length(variable) >= 1L) {
    plot_variable <- variable[[1]]
  }
  
  if(is.null(envir)) {
    if(!is.null(model$model_info$exefuns[[1]])) {
      envir <- environment(model$model_info$exefuns[[1]])
    } else {
      envir <- envir
    }
  }
  
  if(is.null(ndraws)) {
    ndraws <- brms::ndraws(model)
  }
  
  fitted_draws_args <-list()
  add_model_criterion_args <-list()
  
  prior_sensitivity_args <- get_args_(as.list(match.call())[-1], 
                                      "prior_sensitivity.bgmfit")
  
  prior_sensitivity_args[['...']] <- NULL
  prior_sensitivity_args <- c(prior_sensitivity_args, list(...))
  
  fitted_draws_args_names <- methods::formalArgs(fitted_draws.bgmfit)
  fitted_draws_args_names <- 
    fitted_draws_args_names[!fitted_draws_args_names == "..."]
  for (argsi in fitted_draws_args_names) {
    fitted_draws_args[[argsi]] <- prior_sensitivity_args[[argsi]]
  }
  
  add_model_criterion_args_names <- 
    methods::formalArgs(add_model_criterion.bgmfit)
  add_model_criterion_args_names <- 
    add_model_criterion_args_names[!add_model_criterion_args_names == "..."]
  for (argsi in add_model_criterion_args_names) {
    add_model_criterion_args[[argsi]] <- prior_sensitivity_args[[argsi]]
  }
  
  fitted_draws_args[['ipts']] <- FALSE
  fitted_draws_args[['deriv']] <- 0
  fitted_draws_args[['model_deriv']] <- TRUE
  fitted_draws_args[['difx']] <- NULL
  fitted_draws_args[['newdata_fixed']] <- TRUE
  
  add_model_criterion_args[['ipts']] <- FALSE
  add_model_criterion_args[['deriv']] <- 0
  add_model_criterion_args[['model_deriv']] <- TRUE
  add_model_criterion_args[['difx']] <- NULL
  add_model_criterion_args[['newdata_fixed']] <- TRUE
  
  fitted_draws_args <- 
    sanitize_CustomDoCall_args(what = "CustomDoCall", 
                               arguments = fitted_draws_args, 
                               check_formalArgs = fitted_draws.bgmfit,
                               check_formalArgs_exceptions = c('x'),
                               check_trace_back = NULL,
                               envir = parent.frame())
  
  add_model_criterion_args <- 
    sanitize_CustomDoCall_args(what = "CustomDoCall", 
                               arguments = add_model_criterion_args, 
                               check_formalArgs = add_model_criterion.bgmfit,
                               check_formalArgs_exceptions = c('x'),
                               check_trace_back = NULL,
                               envir = parent.frame())
  
  posterior_draw_args <- fitted_draws_args
  posterior_draw_args[['x']] <- posterior_draw_args[['model']]
  posterior_draw_args[['model']] <- NULL
  
  posterior_draw_args[['x']] <- 
    expose_model_functions(posterior_draw_args[['x']],
                           expose = expose_function)
  
  get_post_draws <-  CustomDoCall(posterior::as_draws_df, 
                                  posterior_draw_args)
  
  get_log_lik_draws <-  CustomDoCall(priorsense::log_lik_draws, 
                                     posterior_draw_args)
  
  post_draws <- get_post_draws %>% posterior::bind_draws(get_log_lik_draws)
  
  if(!is.null(add_model_criterion_args[['criterion']])) {
    allowed_criterion <- c("bayes_R2", "loo_R2", "marglik", "jointlik")
    for (izx in add_model_criterion_args[['criterion']]) {
      if(!izx %in% allowed_criterion) {
        paste0(izx, " is not allowed. The allowed criterion are",
               collapse_comma(allowed_criterion))
      }
    }
    if("jointlik" %in% allowed_criterion) {
      add_model_criterion_args[['criterion']] <- 
        add_model_criterion_args[['criterion']][
          add_model_criterion_args[['criterion']] != "jointlik"
        ]
      if(is.null(include_jointlik)) include_jointlik <- TRUE
    }
  }
  
  if(is.null(include_criterion)) {
    if(!is_emptyx(add_model_criterion_args[['criterion']])) {
      include_criterion <- TRUE
    }
  }
  
  
  get_criterion_draws_fun <- function(add_model_criterion_args, 
                                      predict_fn = NULL,
                                      ...) {
    add_model_criterion_args[['model']] <- NULL
    add_model_criterion_name <- add_model_criterion_args[['criterion']]
    
    if(is.null(predict_fn)) {
      predict_fn <- 'add_model_criterion'
    }
    
    if(is.function(predict_fn)) {
      predict_fn <- predict_fn
    } else if(predict_fn == 'add_model_criterion') {
      prediction_names <- add_model_criterion_name
      predict_fn <- function(x, ...) {
        out <- do.call(
          add_model_criterion,
          c(list(model = x), add_model_criterion_args, list(...))
        )
        return(out[[add_model_criterion_name]])
      }
    } else if(predict_fn == 'jointlik') {
      prediction_names <- 'jointlik'
      jointlik_fun <- function(x) {
        as.matrix(rowSums(log_lik(x)))
      }
      predict_fn <- jointlik_fun
    }
    
    predictions_as_draws_args <- list(
      x = model,
      predict_fn = predict_fn,
      prediction_names =prediction_names)
    
    get_criterion_draws <- do.call(
      priorsense::predictions_as_draws,
      predictions_as_draws_args) 
  }
  
  if (isTRUE(include_criterion)) {
    for (ctri in add_model_criterion_args[['criterion']]) {
      add_model_criterion_args[['criterion']] <- ctri
      add_model_criterion_args[['clearenvfuns']] <- FALSE
      criterion_draws_df <- get_criterion_draws_fun(add_model_criterion_args,
                                                    predict_fn = 
                                                      'add_model_criterion',
                                                    ...)
      post_draws <- post_draws %>% posterior::bind_draws(criterion_draws_df)
    }
  }
  
  if (isTRUE(include_jointlik)) {
    add_model_criterion_args[['model']] <- 
      add_model_criterion_args[['clearenvfuns']] <- FALSE
    jointlik_draws_df <- get_criterion_draws_fun(add_model_criterion_args,
                                                 predict_fn = 'jointlik',
                                                 ...)
    post_draws <- post_draws %>% posterior::bind_draws(jointlik_draws_df)
  }
  
  call_include_predictor <- FALSE
  if (!is.null(include_predictor)) {
    call_include_predictor <- TRUE
    include_predictor_df <- fitted_draws_args[['model']][['data']]
    include_predictor_list <- list()
    if(is.character(include_predictor)) {
      allowed_include_predictor <- c('minmax', 'range', "min", 'max', 
                                     "mean", 'median')
      for (izx in include_predictor) {
        if(!izx %in% allowed_include_predictor) {
          paste0(izx, " is not allowed. The allowed include_predictor are",
                 collapse_comma(allowed_include_predictor))
        }
      }
      include_predictor_vals <- c()
      include_predictor_vals_names <- c()
      if('minmax' %in% include_predictor | 'range' %in% include_predictor) {
        include_predictor_vals <- 
          c(include_predictor_vals,
            range(include_predictor_df[[name_predictor]]))
        include_predictor_vals_names <- c(include_predictor_vals_names,
                                          "min", 'max')
        
      }
      if('min' %in% include_predictor) {
        include_predictor_vals <- c(include_predictor_vals,
                                    min(include_predictor_df[[name_predictor]]))
        include_predictor_vals_names <- c(include_predictor_vals_names,
                                          'min')
      }
      if('max' %in% include_predictor) {
        include_predictor_vals <- c(include_predictor_vals,
                                    max(include_predictor_df[[name_predictor]]))
        include_predictor_vals_names <- c(include_predictor_vals_names,
                                          'max')
      }
      if('mean' %in% include_predictor) {
        include_predictor_vals <- 
          c(include_predictor_vals,
            mean(include_predictor_df[[name_predictor]]))
        include_predictor_vals_names <- c(include_predictor_vals_names,
                                          'mean')
      }
      if('median' %in% include_predictor) {
        include_predictor_vals <- 
          c(include_predictor_vals,
            median(include_predictor_df[[name_predictor]]))
        include_predictor_vals_names <- c(include_predictor_vals_names,
                                          'median')
      }
      include_predictor_vals_names <- paste0(name_predictor, "", 
                                             include_predictor_vals_names)
      include_predictor_list[[name_predictor]] <- include_predictor_vals
    } else if(is.list(include_predictor)) {
      err_msg <- "'include_predictor' must be named 
      list with each elemnet namaed"
      if(is.null(names(include_predictor))) stop2c(err_msg)
      if(length(names(include_predictor)) != length(include_predictor)) {
        stop2c(err_msg)
      }
      include_predictor_vals <- unlist(include_predictor)
      if(!is.numeric(include_predictor_vals)) 
        stop("include_predictor must be numeric")
      include_predictor_vals_names <- paste0(name_predictor, "", 
                                             include_predictor_vals)
      include_predictor_list[[name_predictor]] <- include_predictor
    }
  }
  
  if(call_include_predictor) {
    get_fitted_draws_fun <- function(fitted_draws_args, 
                                     predict_fn = NULL,
                                     include_predictor_df = NULL,
                                     name_predictor = NULL,
                                     prediction_names = NULL,
                                     prediction_vals = NULL,
                                     ...) {
      fitted_draws_args[['summary']] <- TRUE
      fitted_draws_args[['newdata']] <- include_predictor_df
      fitted_draws_args[['newdata']][[name_predictor]] <- prediction_vals
      fitted_draws_args[['model']] <- 
        expose_model_functions(fitted_draws_args[['model']], expose = F)
      
      if(is.null(predict_fn)) {
        predict_fn <- brms::posterior_epred
        fitted_draws_args[['model']] <- NULL
      } else if(predict_fn == 'fitted_draws') {
        predict_fn <- function(x, ...) {
          fitted_draws_args[['model']] <- NULL
          out <- do.call(
            fitted_draws,
            c(list(model = x), fitted_draws_args, list(...))
          )[,1] %>% as.matrix()
          return(out)
        }
      }
      
      predictions_as_draws_args <- list(
        x = model,
        predict_fn = predict_fn,
        prediction_names =prediction_names)
      
      get_criterion_draws <- do.call(
        priorsense::predictions_as_draws,
        predictions_as_draws_args) 
    }
    
    for (ix in 1:length(include_predictor_vals_names)) {
      fitted_draws_df <-  get_fitted_draws_fun (
        fitted_draws_args, 
        predict_fn = 'fitted_draws',
        include_predictor_df = include_predictor_df,
        name_predictor = name_predictor,
        prediction_names = include_predictor_vals_names[ix],
        prediction_vals = include_predictor_vals[ix]
      )
      
      post_draws <- post_draws %>% posterior::bind_draws(fitted_draws_df)
    }
  }
  
  sens <- priorsense::powerscale_sensitivity(
    x = post_draws,
    variable = variable,
    lower_alpha = lower_alpha,
    upper_alpha = upper_alpha,
    div_measure = div_measure,
    measure_args = measure_args,
    component = component,
    sensitivity_threshold = sensitivity_threshold,
    moment_match = moment_match,
    k_threshold = k_threshold,
    resample = resample,
    transform = transform,
    prediction = prediction,
    prior_selection = prior_selection,
    likelihood_selection = likelihood_selection,
    num_args = num_args,
    ...
  )
  
  plot_obj <- NULL
  if (!is.null(plot)) {
    if (is.null(plot_variable)) {
      stop2c("A plotting variable must be supplied via 
             `plot_variable` or `variable`.")
    }
    
    if(is.null(facet)) facet <- 'variable'
    
    if(call_include_predictor | include_jointlik | include_criterion) {
      plot_input <- post_draws
    } else {
      plot_input <- model
    }
    
    if (identical(plot, "dens")) {
      plot_obj <- priorsense::powerscale_plot_dens(
        plot_input,
        variable = plot_variable,
        facet_rows = facet
      )
    }
    
    if (identical(plot, "ecdf")) {
      plot_obj <- priorsense::powerscale_plot_ecdf(
        plot_input,
        variable = plot_variable,
        facet_rows = facet
      )
    }
    
    if (identical(plot, "quantities")) {
      plot_obj <- priorsense::powerscale_plot_quantities(
        plot_input,
        variable = plot_variable
      )
    }
    
    if (identical(plot, "both")) {
      plot_obj <- list(
        dens = priorsense::powerscale_plot_dens(
          plot_input,
          variable = plot_variable,
          facet_rows = facet
        ),
        ecdf = priorsense::powerscale_plot_ecdf(
          plot_input,
          variable = plot_variable,
          facet_rows = facet
        )
      )
    }
  }
  
  
  if(!is.null(plot_obj)) {
    if(plot_print) print(plot_obj)
    if(plot_return) print(plot_obj)
  }
  
  if(return_ps) return(sens)
  
  out <- 
    structure(
      list(
        model = model,
        post_draws = post_draws,
        sensitivity = sens,
        plot = plot_obj,
        plot_type = plot,
        plot_variable = plot_variable,
        call = match.call()
      ),
      class = c("prior_sensitivity", "bgmfit")
    )
  
  return(out)
}





#' @rdname prior_sensitivity
#' @export
prior_sensitivity <- function(model, ...) {
  UseMethod("prior_sensitivity")
}


