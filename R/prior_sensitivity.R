

#' Prior sensitivity analysis workflow for \pkg{bsitar} models
#'
#' Run a \pkg{priorsense} power-scaling sensitivity analysis for a fitted
#' \code{bsitar} model. Optionally, append derived quantities to a posterior
#' draws object in the same style as the recommended \pkg{priorsense} workflow.
#' Also support returning diagnostic plot(s).
#'
#' This function is intended as a package-level convenience wrapper around
#' [priorsense::powerscale_sensitivity()] and the associated plotting tools.
#'
#' If requested, the function also constructs an augmented draws object using
#' fit criterion \code{see include_criterion}, joint log-likelihood \code{see
#' include_criterion}, or posterior expected predictions for representative
#' covariate values \code{see include_predictor}. This follows the general
#' workflow recommended in the \pkg{priorsense} package.
#'
#' @param model A fitted model object from \code{bsitar}.
#'
#' @param variable Optional character vector of parameter or quantity names to
#'   assess. Passed to [priorsense::powerscale_sensitivity()] and, when
#'   plotting, to the corresponding \pkg{priorsense} plotting function.
#' 
#' @param include_criterion Logical; if \code{TRUE}, appends one or more of the
#'   following fit criteria using [add_model_criterion()]: \code{"bayes_R2"},
#'   \code{"loo_R2"}, \code{"jointlik"}, or \code{"margliklik"}.
#'
#' @param include_jointlik Logical; if \code{TRUE}, appends a \code{joint log-
#'   likelihood} quantity using [priorsense::predictions_as_draws()] with the
#'   custom \code{jointlik_fun} function. If \code{include_criterion} contains
#'   \code{"jointlik"}, \code{include_jointlik} is automatically set to
#'   \code{TRUE}.
#'
#' @param include_predictor A character string or named list specifying
#'   predictor-based quantities to append as \code{expected values at predictor
#'   values}, using [priorsense::predictions_as_draws()] with a custom call to
#'   [brms::posterior_epred()].
#'
#' @param name_predictor A character string or \code{NULL} (default). If
#'   \code{NULL}, \code{name_predictor} is set internally to \code{xvar}. Note
#'   that \code{prediction_names} are generated internally as
#'   \code{paste0(name_predictor, include_predictor)}. These
#'   \code{prediction_names} can then be included as \code{plot_variable} or to
#'   get potential conflicts by running the [prior_conflict()].
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
#' @param return_table Optional logical indicating whether to return the raw
#'   sensitivity object from [priorsense::powerscale_sensitivity()] or a
#'   structured list with class attributes. Set \code{return_table = FALSE} (the
#'   default) if you plan to call [prior_conflict()] later.
#'   
#' @param return_conflict Optional logical indicating whether to return the raw
#'   sensitivity results \code{return_conflict = FALSE} from
#'   [priorsense::powerscale_sensitivity()] or the \code{conflict}
#'   \code{return_conflict = TRUE} as evaluated using the [prior_conflict()].
#'   The \code{return_table} behavior is same whether \code{return_conflict =
#'   TRUE} or \code{return_conflict = FALSE}.
#'   
#' @param digits Optional integer to pass on to [base::round()] applied to the
#'   numeric variables returned from the [priorsense::powerscale_sensitivity()].
#'   When \code{digits = NULL} (default), \code{digits} is ignored.
#' 
#' @param ... Additional arguments passed to
#'   [priorsense::powerscale_sensitivity()].
#'   
#' @inheritParams add_model_criterion
#' @inheritParams fitted_draws
#' @inheritParams priorsense::powerscale_sensitivity
#' @inheritParams prior_conflict
#'
#' @details
#' The wrapper supports two closely related workflows.
#'
#' \strong{Extended workflow:} A posterior draws object is created with
#' [posterior::as_draws_df()] and can optionally be extended with:
#' \itemize{
#'   \item point wise log-likelihood draws via [priorsense::log_lik_draws()]
#'   (currently ignored),
#'   \item log-prior draws via [priorsense::log_prior_draws()]
#'   (currently ignored),
#'   \item fit criteria via [add_model_criterion()]. See
#'   \code{include_criterion} for details.
#'   \item a \code{jointlik} quantity computed with \code{jointlik_fun},
#'   \item posterior expected predictions via [brms::posterior_epred()].
#' }
#'
#' This extended draws object can be useful when sensitivity is to be explored
#' not only for parameters but also for model-level summaries and predictive
#' quantities, consistent with the \pkg{priorsense} workflow.
#'
#' In the \pkg{priorsense} framework, sensitivity to both prior and
#' likelihood perturbations can indicate possible prior-data conflict,
#' whereas strong prior sensitivity with relatively weak likelihood
#' sensitivity can suggest that the likelihood is only weakly informative
#' for that quantity.
#'
#' @return 
#' If \code{return_table = TRUE}, a table of sensitivity values for each specified
#' variable is returned. If \code{return_table = FALSE}, a
#' \code{"prior_sensitivity"} object is returned as a list with the following
#' components:
#' \describe{
#'   \item{model}{The extracted \code{bsitar} object.}
#'   \item{post_draws}{The augmented draws object if any of 
#'   \code{include_criterion}, \code{include_jointlik}, or 
#'   \code{include_predictor} are specified, otherwise \code{NULL}.}
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
#' [prior_conflict()]
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
#' # In `variable` and `plot_variable`, the custom variables `bayes_R2`,
#' # `jointlik`, `agemin`, and `agemax` are requested via the `criterion`
#' # and `include_predictor` arguments. Here, the predictor of interest is
#' # `age`, so `min` and `max` are prefixed with `age`, giving `agemin`
#' # and `agemax`. The predictor `age` is automatically inferred from the
#' # `model` object, but it can also be specified explicitly using `xvar`.
#' 
#' ps3 <- prior_sensitivity(
#'   model = model,
#'   variable = c("sigma", "bayes_R2", "jointlik", "agemin", "agemax"),
#'   criterion = c('bayes_R2', "jointlik"),
#'   include_predictor = c("min", "max"),
#'   plot = "ecdf",
#'   plot_variable = c("sigma", "bayes_R2", "jointlik", "agemin", "agemax")
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
    include_criterion = NULL,
    include_jointlik = NULL,
    include_predictor = NULL,
    name_predictor = NULL,
    plot = NULL,
    plot_variable = NULL,
    plot_print = FALSE,
    plot_return = FALSE,
    facet = NULL,
    empty = "-",
    print = FALSE,
    return_table = FALSE,
    return_file = NULL,
    flex_table = FALSE,
    return_conflict = FALSE,
    digits = NULL,
    path = NULL,
    title = NULL,
    align = "center",
    sheet_name = "table",
    criterion = "bayes_R2",
    pointwise = FALSE,
    model_name = NULL,
    overwrite = FALSE,
    force_save = FALSE,
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
  add_model_criterion_args[['file']] <- NULL
  add_model_criterion_args[['newdata_fixed']] <- TRUE
  add_model_criterion_args[['return_model']] <- FALSE
  add_model_criterion_args[['return_criteria']] <- TRUE
  
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
  

  
  suppressWarnings({
    suppressMessages({
      
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
      
    })
  })
  
  

  if(!is.null(digits)) {
    sens <- sens %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), 
                                  ~ round(.x, digits)))
  }
  
  if(!return_conflict) {
    if(print) print(sens)
    if(return_table) {
      if(!flex_table) {
        return(sens)
      } else if( flex_table) {
        sens_flex <- flextable::as_flextable(sens, max_row = nrow(sens))
        if (!is.null(title)) {
          sens_flex <- flextable::set_caption(sens_flex, caption = title)
        }
        return(sens_flex)
      }
    }
  }
  
  
  
  
  
  out_list <- 
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
  
  
  
  
  if(return_conflict) {
    out_list <- prior_conflict(out_list,  print = FALSE)
  }
  
  if(!return_table) {
    return(out_list)
  }
 
  out_flex <- list_to_flextable(out_list, align = align, empty = empty)
  
  if (!is.null(title)) {
    out_flex <- flextable::set_caption(out_flex, caption = title)
  }
  
  if(print) print(out_flex$body$dataset)
  
  out <- export_flextable(ft = out_flex,
                          return_file = return_file,
                          path = path,
                          title = title,
                          align = align,
                          sheet_name = sheet_name)
  
  if(is.null(return_file)) {
    #
  } else {
    return_table <- FALSE
  }
  
  
  if(return_table) {
    if(!flex_table) {
      return(out$body$dataset)
    } else if( flex_table) {
      if (!is.null(title)) {
        out <- flextable::set_caption(out, caption = title)
      }
      return(out)
    }
  }
  
  if(!return_table) {
    export_flextable(ft = out_flex,
                     return_file = return_file,
                     path = path,
                     title = title,
                     align = align,
                     sheet_name = sheet_name)
  }
  
  
  
  
}





#' @rdname prior_sensitivity
#' @export
prior_sensitivity <- function(model, ...) {
  UseMethod("prior_sensitivity")
}


