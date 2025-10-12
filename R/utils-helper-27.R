

##############################################################################

# Greedy knot selection algorithm for restricted cubic spline regression

# https://pmc.ncbi.nlm.nih.gov/articles/PMC10910934/

# knutar

# https://github.com/jo-inge-arnes/knutar-experiments

##############################################################################

# knutar - example

# https://github.com/jo-inge-arnes/knutar-experiments/blob/main/mcycle-example.R

##############################################################################


# keywords internal
# noRd
# export

# https://github.com/jo-inge-arnes/knutar/blob/main/R/choose_removal.R

##############################################################################
# choose_model
# importFrom mfp fp
##############################################################################

#' Function assessing different types of regression models for a range of
#' degrees of freedom (knots in the case of splines), returning the one yielding
#' the best results according to an information criterion. The function uses
#' fractional polynomials, restricted cubic splines, and restricted cubic
#' splines with non-uniform knot placements.
#' 
#' @param dataset The data frame
#' 
#' @param dependent The dependent variable in the formula
#' 
#' @param independents The independent variables in the formula
#' 
#' @param icr_fn The information criterion function for comparing different
#'   models with different degress for freedom or knots (default BIC)
#'   
#' @param cost_fn The criterion used to choose which knots to remove, passed to
#'   the function choose_removal. Defaults to BIC.
#'   
#' @param fp_alpha The relax factor for multivariate fractional polynomials.
#'   Ignored
#'   
#' @param max_nsknots The max number of inner knots for restricted cubic splines
#'   (default 4)
#'   
#' @param max_fp_df The max degrees of freedom for fractional polynomials
#'   (default 4, i.e., which is the same as a second-degree/two-terms FP).
#'   Ignored
#'   
#' @param verbose Verbose output, default FALSE
#' 
#' @param boundary_knots The boundary knot placements for restricted cubic
#'   splines or NA if not specified.
#'   
#' @return A list with named elements, such as 'model', 'type', 'score'. The
#'   function returns a list with named elements and sublists, see examples for
#'   full overview of the returned values.
#' 
#' @keywords internal
#' @noRd
#' @examples
#' my_model <- choose_model(d, y, x)$model
#' result <- choose_model(d, y, x, icr_fn = BIC, verbose = FALSE)
#'
#' ret <- choose_model(d, y, x)
#'
#' ret$labels[[ret$type]] # Human readable name of chosen model type
#' ret[[ret$type]]        # Gives more values for the chosen model if available
#'
#' ret$model      # The chosen model
#' ret$score      # The chosen model's score
#' ret$type       # The type of model chosen as a string
#' ret$score_name # The type of score used as a string
#' ret$score_fn   # The function used for scores
#'
#' ret$labels     # Description string for the types of models
#'
#' ret$mfp        # Multivariate fractional polynomial
#' ret$mfp$model  # The model
#' ret$mfp$score  # The score
#'
#' ret$ns         # Restricted cubic splines from regular sequence of quantiles
#' ret$ns$model   # The model
#' ret$ns$score   # The score
#' ret$ns$knot_cnt_arg      # The number of inner knots as input argument
#' ret$ns$knot_cnt_distinct # The number of distinct placements in the result
#' ret$ns$knot_placements   # Knots and boundary knots
#' ret$ns$knot_placements$knots           # The knot placements as a list
#' ret$ns$knot_placements$Boundary.knots  # The boundary knots as a list
#'
#' ret$ns_nu         # Restricted cubic splines w/non-uniform placements
#' ret$ns_nu$model   # The model
#' ret$ns_nu$score   # The score
#' ret$ns_nu$knot_cnt_distinct # The number of distinct placements in the result
#' ret$ns_nu$knot_placements   # Knots and boundary knots
#' ret$ns_nu$knot_placements$knots          # The knot placements as a list
#' ret$ns_nu$knot_placements$Boundary.knots # The boundary knots as a list
#' 
choose_model <- function(dataset,
                         dependent,
                         independents,
                         ...,
                         icr_fn = stats::BIC,
                         cost_fn = stats::BIC,
                         fp_alpha = NA,
                         max_nsknots = 4,
                         max_fp_df = 4,
                         verbose = TRUE,
                         boundary_knots = NA) {
  
  dependent <- rlang::enquo(dependent)
  independents <- rlang::enquo(independents)
  
  score_type <- NULL
  if (missing(icr_fn)) {
    icr_fn <- stats::BIC
    score_type <- "BIC"
  } else {
    score_type <- deparse(substitute(icr_fn))
  }
  if (missing(cost_fn)) cost_fn <- stats::BIC
  # if (missing(fp_alpha)) fp_alpha <- NA
  if (missing(max_nsknots)) max_nsknots <- 4
  # if (missing(max_fp_df)) max_fp_df <- 4
  if (missing(verbose)) verbose <- TRUE
  if (missing(boundary_knots)) boundary_knots <- NA
  
  ret_desc <- list(
    # "mfp" = "Multivariate fractional polynomials",
    "ns_nu" = "Restricted cubic splines with freely placed knots",
    "ns" = "Restricted cubic splines with knots placed at quantiles")
  
  ret <- list(labels = ret_desc, score_fn = icr_fn, score_name = score_type)
  
  only_positive_independents <- 
    all(as.vector(model.matrix(independents, dataset)[, 2]) > 0)
  
  
  # if (only_positive_independents) {
  #   independents_str <- sub("~", "", deparse(independents))
  #   # Multivariate fractional polynomials (move to separate func)
  #   fp_formula <- stats::formula(paste0(
  #     rlang::as_name(dependent),
  #     " ~ fp(",
  #     independents_str,
  #     ", df = ",
  #     max_fp_df,
  #     ")"
  #   ))
  #   mfp_res <- NA
  #   if (is.na(fp_alpha)) {
  #     mfp_res <-
  #       mfp::mfp(fp_formula, data = dataset, verbose = verbose)
  #   } else {
  #     mfp_res <-
  #       mfp::mfp(fp_formula,
  #                alpha = fp_alpha, data = dataset, verbose = verbose)
  #   }
  #   mfp_mod <- eval(summary(mfp_res)$call)
  #   mfp_score <- icr_fn(mfp_mod)
  # 
  #   ret <- append(ret, list(mfp = list(model = mfp_mod, score = mfp_score)))
  # 
  #   if (verbose) {
  #     R.utils::printf("-----------------------------------------------------\n")
  #     R.utils::printf("%s\n%s: %f\n\n",
  #                     ret_desc[["mfp"]], score_type, mfp_score)
  #   }
  # } else {
  #   if (verbose) {
  #     R.utils::printf("-----------------------------------------------------\n")
  #     R.utils::printf(paste("Skipped fractional polynomials because some",
  #                           "independent values were not larger than zero.\n", sep = " "))
  #   }
  # 
  #   mfp_mod <- NA
  #   mfp_score <- Inf
  # }


  
  suppressWarnings({
    # Restricted cubic splines with knots distanced by regular sequence
    # of quantiles between the boundary knots
    knotcnt_suggestion <-
      suggest_knotcount(dataset, !!dependent, !!independents, max_nsknots,
                        icr_fn = icr_fn, boundary_knots = boundary_knots)
    ns_mod <- model_by_count(dataset, !!dependent, !!independents,
                             knotcnt_suggestion$nknots, boundary_knots = boundary_knots)
    ns_score <- icr_fn(ns_mod)
    
    extracted_knots <- extract_knots(ns_mod)
    ret <-
      append(ret, list(ns =
                         list(model = ns_mod,
                              score = ns_score,
                              knot_cnt_arg = knotcnt_suggestion$nknots,
                              knot_cnt_distinct = length(extracted_knots$knots),
                              knot_placements = extracted_knots)))
    
    if (verbose) {
      R.utils::printf("%s\n%s: %f\n", ret_desc[["ns"]], score_type, ns_score)
      R.utils::printf("Suggested knot count: %d\n", knotcnt_suggestion$nknots)
      print_knots(extract_knots(ns_mod))
      R.utils::printf("\n")
    }
    
    # Restricted cubic splines with freely placed knots between the boundaries
    knutar_res <- choose_splines(dataset, !!dependent, !!independents,
                                 max_nsknots, icr_fn = icr_fn, cost_fn = cost_fn,
                                 boundary_knots = boundary_knots)
    
    ret <-
      append(ret, list(ns_nu =
                         list(model = knutar_res$model,
                              score = knutar_res$score,
                              knot_cnt_distinct = length(knutar_res$knots$knots),
                              knot_placements = knutar_res$knots)))
    
    if (verbose) {
      R.utils::printf("%s\n%s: %f\n",
                      ret_desc[["ns_nu"]], score_type, knutar_res$score)
      print_knots(knutar_res$knots)
      R.utils::printf("\n")
    }
    
    # if ((mfp_score <= ns_score) && (mfp_score <= knutar_res$score)) {
    #   ret <- append(ret, list(model = mfp_mod, type = "mfp", score = mfp_score))
    # } else if (ns_score <= knutar_res$score) {
    #   ret <- append(ret, list(model = ns_mod, type = "ns", score = ns_score))
    # } else {
    #   ret <- append(ret, list(model = knutar_res$model, type = "ns_nu",
    #                           score = knutar_res$score))
    # }
    
    if (ns_score <= knutar_res$score) {
      ret <- append(ret, list(model = ns_mod, type = "ns", score = ns_score))
    } else {
      ret <- append(ret, list(model = knutar_res$model, type = "ns_nu",
                              score = knutar_res$score))
    }
    
    if (verbose) {
      R.utils::printf("Chosen model type:\n%s\n", ret_desc[[ret$type]])
    }
  })
  
  return(ret)
}


##############################################################################
# choose_removal 
##############################################################################

#' Finds the best knot to remove from the given model
#'
#' @details
#' The function searches through the knots of the given model to find the inner
#' knot yielding the best resulting score if one knot has to be removed. The
#' selection criterion function can be seen as a cost function, where lower
#' scores are better.
#' 
#' @param dataset The data frame
#' 
#' @param dependent The dependent variable in the formula
#' 
#' @param independents The independent variables in the formula
#' 
#' @param knots The knot placements
#' 
#' @param boundary_knots The boundary knot placements
#' 
#' @param cost_fn The function for the selection criterion score (BIC default)
#' 
#' @return A named list with 'index' of the chosen knot, 'model', and 'score'
#' 
#' @keywords internal
#' @noRd
#' @examples
#' removal_index <- choose_removal(my_data, y, x, knots, boundary_knots)$index
#' removal_index <- choose_removal(my_data, y, x, knots, boundary_knots,
#'  cost_fn = function(model) { return(-summary(model)$fstatistic[[1]]) })$index
#'  
choose_removal <- function(dataset,
                           dependent,
                           independents,
                           knots,
                           boundary_knots,
                           cost_fn = stats::BIC) {
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  
  if (missing(cost_fn)) cost_fn <- stats::BIC
  
  model_scores <- lapply(seq_along(knots), function(i) {
    mod <- model_by_knots(dataset, !!dependent, !!independents,
                          knots = knots[-i], boundary_knots = boundary_knots)
    mod_score <- cost_fn(mod)
    return(list(model = mod, score = mod_score))
  })
  
  scores <- unlist(lapply(model_scores, "[[", "score"))
  index <- which.min(scores)
  min_score <- scores[[index]]
  
  return(list(model = model_scores[[index]][["model"]],
              score = min_score, index = index))
}


##############################################################################
# choose_splines 
##############################################################################

#' Chooses the best from a set of restricted cubic spline (RCS) models with
#' an inner knot count lower than or equal to a specified maximum number
#'
#' @details
#' The maximum number of inner knots is given as an input argument.
#' 
#' @param dataset The data frame
#' 
#' @param dependent The dependent variable in the formula
#' 
#' @param independents The independent variables in the formula
#' 
#' @param max_nknots The maximum number of inner knots wanted
#' 
#' @param icr_fn The information criterion function comparing models with
#'   different knot counts (BIC default)
#'   
#' @param cost_fn The criterion used to choose which inner knot to remove, used
#'   by the function 'choose_removal'. Default is BIC.
#'   
#' @param initial_nknots The initial high number inner of knots for the
#'   algorithm (default is the value from the 'suggest_knotcount'-function)
#'   
#' @param diff_better How much lower must the score be to make a higher knot
#'   model be deemed a better model than an alternative lower knot model?
#'   
#' @param all_models If TRUE, the function will include all intermediate models
#'   in the results as 'all_models'. Default is FALSE.
#'   
#' @param boundary_knots The boundary knot placements or NA if not specified
#' 
#' @return The chosen 'model', 'score', 'knots', and 'all_models'
#' 
#' @keywords internal
#' @noRd
#' @examples
#' my_model <- choose_splines(d, y, x, 7)
#' my_model <- choose_splines(d, y, x, 7, BIC)
#' 
choose_splines <- function(dataset,
                           dependent,
                           independents,
                           max_nknots = 10,
                           ...,
                           icr_fn = stats::BIC,
                           cost_fn = stats::BIC,
                           initial_nknots = -1,
                           diff_better = 0,
                           all_models = FALSE,
                           boundary_knots = NA) {
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  
  
  if (missing(max_nknots)) max_nknots <- 10
  if (missing(icr_fn)) icr_fn <- stats::BIC
  if (missing(cost_fn)) cost_fn <- stats::BIC
  if (missing(initial_nknots)) initial_nknots <- -1
  if (missing(diff_better)) diff_better <- 0
  if (missing(all_models)) all_models <- FALSE
  if (missing(boundary_knots)) boundary_knots <- NA
  
  if (initial_nknots == -1) {
    initial_nknots <-
      suggest_knotcount(dataset, !!dependent, !!independents,
                        boundary_knots = boundary_knots)$nknots
  }
  
  upper_model <- suggest_splines(dataset, !!dependent, !!independents,
                                 max_nknots,
                                 initial_knots = initial_nknots,
                                 cost_fn = cost_fn,
                                 boundary_knots = boundary_knots)
  
  cur_model <- upper_model
  best_model <- cur_model
  cur_score <- icr_fn(best_model)
  best_score <- cur_score
  best_knots <- extract_knots(best_model)
  cur_nknots <- length(best_knots$knots)
  intermediate_models <- list()
  
  while (cur_nknots > 0) {
    these_knots <- extract_knots(cur_model)
    chosen <- choose_removal(dataset, 
                             !!dependent, 
                             !!independents,
                             these_knots$knots,
                             these_knots$Boundary.knots, 
                             cost_fn)
    
    cur_score <- icr_fn(chosen$model)
    if (cur_score <= (best_score + diff_better)) {
      best_model <- chosen$model
      best_score <- cur_score
      best_knots <- extract_knots(best_model)
    }
    cur_model <- chosen$model
    cur_nknots <- length(extract_knots(cur_model)$knots)
    if (all_models) {
      intermediate_models <- append(intermediate_models, list(cur_model))
    }
  }
  
  return(
    list(model = best_model, score = best_score, knots = best_knots,
         all_models = intermediate_models))
}

##############################################################################
# extract_knots 
##############################################################################

#' Extracts the distinct knot placements of natural spline regression models
#'
#' @param ns_model The restricted cubic spline regression model
#' 
#' @return A list with named elements 'knots' (inner) and 'Boundary.knots'
#' 
#' @keywords internal
#' @noRd
#' @examples
#' extract_knots(my_model)
#' 
extract_knots <- function(ns_model) {
  knots <- attr(ns_model$model[[2]], "knots")
  knots <- unique(knots)
  bknots <- attr(ns_model$model[[2]], "Boundary.knots")
  return(list(knots = knots, Boundary.knots = bknots))
}


##############################################################################
# generate_data 
##############################################################################

#' Generates synthetic data
#'
#' Se the examples section for information about the returned column names
#'
#' @param n The number of units in the sample
#' 
#' @param x_accr The accuracy (number of decimals) of independent variable
#' 
#' @param y_accr The accuracy (number of decimals) of dependent variable
#' 
#' @param f_x_dist Function for the distribution of the independent variable
#' 
#' @param f_signal The true function for the population means. Takes the x
#'   values as input.
#'   
#' @param f_noise The noise function, i.e., variance distribution(s) around the
#'   population means for y, which can be heteroscedastic.
#'   
#' @return A dataframe with generated data
#' 
#' @keywords internal
#' @noRd
#' @examples
#' The returned data frame has the following columns:
#' ID               # IDs
#' Independent      # The measured X values (rounded according to accuracy)
#' Dependent        # The measured Y values (rounded according to accuracy)
#' IndependentRaw   # The raw X values (not rounded)
#' DependentRaw     # The raw Y values (not rounded)
#' SignalRaw        # The raw signal part of the population mean curve for Y
#' Noise            # The raw noise (variance) for the DependentRaw values
#' SignalMeasured   # The signal rounded according to accuracy
#' 
generate_data <- function(n, x_accr, y_accr, f_x_dist,
                          f_signal, f_noise) {
  ids <- 1:n
  xs_raw <- f_x_dist(n)
  ys_signal <- f_signal(xs_raw)
  ys_noise <- f_noise(xs_raw)
  ys_raw <- ys_signal + ys_noise
  
  xs_measured <- xs_raw
  if (!missing(x_accr) && !is.null(x_accr) && x_accr > -1) {
    xs_measured <- round(xs_raw, x_accr)
  }
  
  ys_measured <- ys_raw
  if (!missing(y_accr) && !is.null(y_accr) && y_accr > -1) {
    ys_measured <- round(ys_raw, y_accr)
  }
  
  xs_measured_signal <- f_signal(xs_measured)
  
  return(data.frame(
    ID = ids,
    Independent = xs_measured,
    Dependent = ys_measured,
    IndependentRaw = xs_raw,
    DependentRaw = ys_raw,
    SignalRaw = ys_signal,
    Noise = ys_noise,
    SignalMeasured = xs_measured_signal
  ))
}


##############################################################################
# plot_model 
##############################################################################

#' Plots the (glm) model with data, curve, confidence bands, etc.
#'
#' @param dataset The data frame
#' 
#' @param mod The model
#' 
#' @keywords internal
#' @noRd
#' @examples
#' plot_model(d, y, x, my_mod)
#' 
plot_model <- function(dataset, dependent, independent, model) {
  fit_link <- NULL;
  se_link  <- NULL;
  lwr  <- NULL;
  upr  <- NULL;
  pred  <- NULL;
  
  # The confidence bands has to be in the scale of the link function,
  # and we need to get it's inverse for plotting
  ilink <- stats::family(model)$linkinv
  
  d <- dataset
  d <- d %>%
    dplyr::mutate(pred = stats::predict(model)) %>%
    dplyr::bind_cols(stats::setNames(dplyr::as_tibble(
      stats::predict(model, newdata = d, se.fit = TRUE)[1:2]),
      c("fit_link", "se_link"))) %>%
    dplyr::mutate(fit_resp = ilink(fit_link),
                  upr = ilink(fit_link + (2 * se_link)),
                  lwr = ilink(fit_link - (2 * se_link)))
  
  fig <- ggplot2::ggplot(d,
                         ggplot2::aes(x = {{ independent }}, y = {{ dependent }})) +
    ggplot2::geom_point() +
    ggplot2::geom_ribbon(data = d,
                         ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = pred),
                       color = "blue", size = 1)
  
  return(fig)
}



##############################################################################
# plot_rstandard 
##############################################################################

#' Plots the standardized residuals for the model
#'
#' @param dataset The data frame
#' 
#' @param mod The model
#' 
#' @keywords internal
#' @noRd
#' @examples
#' plot_rstandard(d, my_mod)
#' 
plot_rstandard <- function(dataset, model) {
  ggplot2::ggplot(dataset, ggplot2::aes(x = stats::predict(model, dataset),
                                        y = stats::rstandard(model))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x)
}



##############################################################################
# print_knots 
##############################################################################

#' Print knot placements
#'
#' @param knot_placements A list with two named lists, 'knots' and
#'   'Boundary.knots'
#'   
#' @keywords internal
#' @noRd
#' @examples
#' print_knots(knot_placements)
#' 
print_knots <- function(knot_placements) {
  R.utils::printf("Inner knots count: %d\n", length(knot_placements$knots))
  knots_str <- paste0("[", paste0(knot_placements$knots, collapse = ", "), "]")
  boundary_str <-
    paste0("[", paste0(knot_placements$Boundary.knots, collapse = ", "), "]")
  R.utils::printf("Inner knots: %s\nBoundary knots: %s\n", knots_str,
                  boundary_str)
}


##############################################################################
# restricted_cubic_splines
##############################################################################

#' Utility function for restricted cubic splines that wraps the 'glm' and 'ns'
#' functions
#'
#' @details
#' Creates a restricted cubic splines regression model given the wanted number
#' of inner knots
#' 
#' @param dataset The data frame
#' 
#' @param dependent The dependent variable in the formula
#' 
#' @param independents The independent variables in the formula
#' 
#' @param nknots The requested number of inner knots, excluding the boundary
#' knots
#' 
#' @param boundary_knots The boundary knot placements
#' 
#' @return The regression model
#' 
#' @keywords internal
#' @noRd
#' @examples
#' my_model <- create_model(my_data, y, x, 7)
model_by_count <- function(dataset, 
                           dependent, 
                           independents, 
                           nknots,
                           boundary_knots = NA,
                           stype = 'ns'
                           ) {
  
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  
  independents_str <- sub("~", "", deparse(independents))
  
  if (missing(boundary_knots)) boundary_knots <- NA
  
  if (length(boundary_knots) != 2) {
    model_formula <- stats::as.formula(
      paste0(rlang::as_name(dependent), 
             " ~ splines::ns(", 
             independents_str,
             ", df = ", 
             nknots + 1, 
             ")"))
  } else {
    boundary_knots_str <- paste0(
      "c(", paste0(boundary_knots, collapse = ", "), ")")
    formula_str <- paste0(
      rlang::as_name(dependent), 
      " ~ splines::ns(", 
      independents_str,
      ", df = ", 
      nknots + 1, 
      ", Boundary.knots = ", boundary_knots_str, 
      ")")
    model_formula <- stats::as.formula(formula_str)
  }
  
  ns_model <- stats::glm(model_formula, data = dataset)
  
  return(ns_model)
}

#' Utility function for restricted cubic splines that wraps the 'glm' and 'ns'
#' functions
#'
#' Creates a restricted cubic spline regression model from the given placements
#' of the inner knots and boundary knots
#' @param dataset The data frame
#' @param dependent The dependent variable in the formula
#' @param independents The independent variables in the formula
#' @param knots The inner knot placements
#' @param boundary_knots The boundary knot placements
#' @return The regression model
#' @keywords internal
#' @noRd
#' @examples
#' my_model <- model_by_knots(my_data, y, x, c(0.1, 0.2), c(0.0, 0.3))
model_by_knots <- function(dataset,
                           dependent,
                           independents,
                           knots,
                           boundary_knots,
                           stype = 'ns'
                           ) {
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  
  independents_str <- sub("~", "", deparse(independents))
  
  knots_str <- paste0(
    "c(", paste0(knots, collapse = ", "), ")")
  boundary_knots_str <- paste0(
    "c(", paste0(boundary_knots, collapse = ", "), ")")
  formula_str <- paste0(
    rlang::as_name(dependent), 
    " ~ splines::ns(", 
    independents_str,
    ", knots = ", 
    knots_str, 
    ", Boundary.knots = ", 
    boundary_knots_str, 
    ")")
  model_formula <- stats::as.formula(formula_str)
  
  ns_model <- stats::glm(model_formula, data = dataset)
  
  return(ns_model)
}



##############################################################################
# suggest_knotcount 
##############################################################################

#' Finds the number of restricted cubic spline inner knots that gives the
#' lowest score for an information criterion and a data set.
#'
#' @param dataset The data frame
#' @param dependent The dependent variable in the formula
#' @param independents The independent variable(s) in the formula
#' @param maximum_knots The highest inner knot count to assess.
#' Defaults the lower of 50 and n / 2
#' @param icr_fn The information criterion function. Defaults to BIC
#' @param all_scores If TRUE, all scores are returned in a list 'all_scores'
#' @param boundary_knots The boundary knot placements
#' @return A list with named elements 'nknots', 'score', and 'all_scores'
#' @keywords internal
#' @noRd
#' @examples
#' suggest_knotcount(d, nwsize, age_dec)
suggest_knotcount <- function(dataset,
                              dependent,
                              independents,
                              max_nknots = -1,
                              ...,
                              icr_fn = stats::BIC,
                              all_scores = FALSE,
                              boundary_knots = NA,
                              stype = 'ns'
                              ) {
  dependent <- rlang::enquo(dependent)
  independents <- rlang::enquo(independents)
  
  if (missing(max_nknots) || max_nknots == -1) {
    max_nknots <- min(50, nrow(dataset) %/% 2)
  }
  
  if (missing(icr_fn)) icr_fn <- stats::BIC
  if (missing(all_scores)) all_scores <- FALSE
  if (missing(boundary_knots)) boundary_knots <- NA
  
  min_icr <- Inf
  min_ndf <- Inf
  
  n_knots <- list()
  scores <- list()
  
  independents_str <- sub("~", "", deparse(independents))
  
  if (length(boundary_knots) != 2) {
    boundary_knots_str <- ""
  } else {
    boundary_knots_str <- paste0(
      ", Boundary.knots = c(", paste0(boundary_knots, collapse = ", "), ")"
    )
  }
  
  consecutive_non_convergance <- 0
  
  for (i in 1:(max_nknots + 1)) {
    # Because of a change between R 4.2.0 and 4.3.0, quantiles coinciding with
    # boundary knots must now be removed.
    #
    # This caused a bit of rewriting...
    #
    # See: https://stat.ethz.ch/pipermail/r-announce/2023/000691.html
    # "bs() and ns() in the (typical) case of automatic knot construction, when
    # some of the supposedly inner knots coincide
    # with boundary knots, now moves them inside (with a warning),
    # building on PR#18442 by Ben Bolker."
    
    subset_data <- dataset %>%
      dplyr::filter(
        !!independents >= boundary_knots[[1]],
        !!independents <= boundary_knots[[2]]
      )
    
    n <- i - 1
    knots <- c()
    
    if (n > 0) {
      quantiles <- unique(quantile(subset_data[[rlang::as_label(independents)]],
                                   probs = c(0, seq(1 / n, 1, by = 1 / n))
      ))
      
      if (quantiles[[1]] == boundary_knots[[1]]) {
        quantiles <- quantiles[-1]
      }
      
      if (!length(quantiles) == 0 &&
          quantiles[[length(quantiles)]] == boundary_knots[[2]]) {
        quantiles <- quantiles[-length(quantiles)]
      }
      
      knots <- quantiles
    }
    
    model_formula_str <- NULL
    
    if (length(knots) > 0) {
      model_formula_str <- paste0(
        rlang::as_label(dependent),
        " ~ splines::ns(",
        independents_str,
        ", knots = c(",
        paste(knots, collapse = ", "),
        ")",
        boundary_knots_str,
        ")"
      )
    } else {
      model_formula_str <- paste0(
        rlang::as_label(dependent),
        " ~ splines::ns(",
        independents_str,
        ", df = 1",
        boundary_knots_str,
        ")"
      )
    }
    
    model_formula <- formula(model_formula_str)
    
    mod_spline <- NULL
    
    try(mod_spline <- stats::glm(model_formula, data = dataset))
    
    if (!is.null(mod_spline) && mod_spline$converged) {
      consecutive_non_convergance <- 0
    } else {
      consecutive_non_convergance <- consecutive_non_convergance + 1
    }
    
    if (consecutive_non_convergance == 0) {
      icr_score <- icr_fn(mod_spline)
      
      if (all_scores) {
        scores <- append(scores, icr_score)
        n_knots <- append(n_knots, i - 1)
      }
      
      if (icr_score < min_icr) {
        min_icr <- icr_score
        min_ndf <- i
      }
    } else if (consecutive_non_convergance >= 3) {
      warning(paste(
        "Models failed to converge three consecutive times,",
        "will not assess any higher knot counts."
      ))
      break
    }
  }
  
  return(list(
    nknots = min_ndf - 1, score = min_icr,
    all_scores = list(scores = scores, n_knots = n_knots)
  ))
}


##############################################################################
# suggest_splines 
##############################################################################

#' Suggests a restricted cubic splines regression model with inner knot
#' placements that can be non-uniform with respect to quantiles and widths
#'
#' The target number of inner knots for the model is given as an input argument.
#' The algorithm starts with a model with a high number of knots and
#' systematically removes inner knots until the target number of inner knots is
#' reached.
#' The initial number of inner knots, before starting to remove knots, can be
#' given as an argument as well, but  defaults to the suggested number
#' obtained from the function 'suggest_knotcount'.
#' @param dataset The data frame
#' @param dependent The dependent variable in the formula
#' @param independents The independent variables in the formula
#' @param target_nknots The target and maximum number of inner knots for the
#' model
#' @param initial_nknots The number of inner knots initially, defaults to the
#' result from the function 'suggest_knotcount'
#' @param cost_fn The function for the selection criterion score (BIC default)
#' used to compare which inner knot should be removed, passed to choose_removal
#' @param all_knots If TRUE, then knots for all intermediate models will be
#' included in return value. Default is FALSE.
#' @param boundary_knots The boundary knot placements or NA if not specified
#' @return Returns the suggested natural splines model, or if the 'all_knots'
#' argument was TRUE, then a list with named elements 'model', 'all_knots', and
#' 'Boundary.knots' is returned.
#' @keywords internal
#' @noRd
#' @examples
#' my_model <- suggest_splines(d, y, x, 4)
#' my_model <- suggest_splines(d, y, x, 4, initial_nknots = 100, cost_fn = BIC)
suggest_splines <- function(dataset,
                            dependent,
                            independents,
                            target_nknots,
                            ...,
                            initial_nknots = -1,
                            cost_fn = stats::BIC,
                            all_knots = FALSE,
                            boundary_knots = NA) {
  independents <- rlang::enquo(independents)
  dependent <- rlang::enquo(dependent)
  if (missing(boundary_knots)) boundary_knots <- NA
  
  if (initial_nknots == -1) {
    initial_nknots <-
      suggest_knotcount(dataset, !!dependent, !!independents,
                        boundary_knots = boundary_knots)$nknots
  }
  
  if (missing(cost_fn)) cost_fn <- stats::BIC
  if (missing(all_knots)) all_knots <- FALSE
  
  # Find the initial model with a high number of knots, and get the distinct
  # knot placements
  ns_model <-
    model_by_count(dataset, !!dependent, !!independents, initial_nknots,
                   boundary_knots)
  knots <- extract_knots(ns_model)
  
  intermediate_knots <- list()
  
  # Initialize the variables that will hold the final knot placements
  final_knots <- knots$knots
  boundary_knots <- knots$Boundary.knots
  
  if (all_knots) {
    intermediate_knots <- append(intermediate_knots, list(final_knots))
  }
  
  # As long as there are more inner knots left than the target number, remove
  # inner knots one by one, by always removing the knot that gives the best
  # resulting model
  if (length(knots$knots) > target_nknots) {
    for (i in 1:(length(knots$knots) - target_nknots)) {
      rm_index <- choose_removal(dataset, 
                                 !!dependent, 
                                 !!independents,
                                 knots = final_knots, 
                                 boundary_knots = boundary_knots, 
                                 cost_fn)$index
      
      final_knots <- final_knots[-rm_index]
      
      if (all_knots) {
        intermediate_knots <- append(intermediate_knots, list(final_knots))
      }
    }
  }
  
  final_mod <- model_by_knots(dataset, !!dependent, !!independents,
                              knots = final_knots, boundary_knots = boundary_knots)
  
  if (all_knots) {
    return(list(model = final_mod, all_knots = intermediate_knots,
                Boundary.knots = boundary_knots))
  } else {
    return(final_mod)
  }
}



##############################################################################
# choose_model
##############################################################################

# https://github.com/jo-inge-arnes/knutar-experiments

create_figure <- function(d, mod, title, x, y) {
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x <- sub("~", "", deparse(x))
  y <- sub("~", "", deparse(y))
  fig <- ggplot2::ggplot()
  fig <- fig + ggplot2::theme_bw()
  fig <- fig + ggplot2::geom_point(data = d, 
                                   ggplot2::aes(.data[[x]], .data[[y]]), 
                                   shape = 1,
                                   color = "gray50")
 
  fig <- fig + ggplot2::ggtitle(title)
  fig <- fig + ggplot2::xlab(scales::parse_format()("'Predictor'~X"))
  fig <- fig + ggplot2::ylab(scales::parse_format()("'Response'~Y"))
  
  knots <- extract_knots(mod)
 
  fig <- fig + ggplot2::geom_vline(xintercept = knots$knots, linetype = "dashed")
  fig <- fig + ggplot2::geom_vline(xintercept = knots$Boundary.knots, linetype = "solid")
  
  fig <- fig + ggplot2::geom_line(ggplot2::aes(x = mod$data[[x]], y = mod$fitted.values),
                         linetype = "solid", color = "black", size = 0.5)
  
  # fig <- fig + force_panelsizes(rows = unit(3.5, "in"), cols = unit(7, "in"))
  
  return(fig)
}




##############################################################################
##############################################################################

# https://github.com/davidaarmstrong/damisc


##############################################################################
##############################################################################




##############################################################################
# NKnots 
##############################################################################

# https://github.com/davidaarmstrong/damisc/blob/master/R/DAMisc_functions.R

#' AIC and BIC selection of number of spline knots
#' 
#' Calculates AIC and BIC for the selection of knots in a spline over values
#' (potentially including polynomials) up to a user-defined maximum.
#' 
#' 
#' @param form A formula detailing the model for which smoothing is to be
#' evaluated.
#' @param var A character string identifying the variable for which smoothing
#' is to be evaluated.
#' @param data Data frame providing values of all variables in \code{form}.
#' @param degree Degree of polynomial in B-spline basis functions.
#' @param min.knots Minimum number of internal B-spline knots to be evaluated.
#' @param max.knots Maximum number of internal B-spline knots to be evaluated.
#' @param includePoly Include linear and polynomial models up to, and including
#' \code{degree}-th order polynomials.
#' @param plot Logical indicating whether a plot should be returned.
#' @param criterion Statistical criterion to minimize in order to find the best
#' number of knots - AIC, BIC or Cross-validation. Ignored here
#' @param cvk Number of groups for cross-validation. cross-validation ignored
#'   here.
#' @param cviter Number of iterations of cross-validation to average over.10
#'   Ignored here.
#' is the default but in real-world applications, this should be somewhere
#' around 200.
#' @returns A plot, if \code{plot=TRUE}, otherwise a data frame with the degrees
#' of freedom and corresponding fit measure.
#' @author Dave Armstrong
#' 
#' @keywords internal
#' @noRd
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' NKnots(prestige ~ education + type, var="income", data=na.omit(Prestige), plot=FALSE)
#' 
NKnots <- function(form, var, data, degree=3, min.knots=1,
                   max.knots=10, includePoly = FALSE, plot=FALSE, 
                   criterion=c("AIC", "BIC", "CV"),
                   cvk=10, cviter=10, stype = 'nsk'){
  crit <- match.arg(criterion)
  k <- seq(min.knots, max.knots, by=1)
  forms <- vector("list", ifelse(includePoly, length(k)+3, length(k)))
  df_poly <- NULL
  m <- 1
  if(includePoly){
    df_poly <- 1:3
    forms[[1]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], " + ", var, sep=""))
    forms[[2]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], "+ poly(", var,  ", 2)", sep=""))
    forms[[3]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], "+ poly(", var,  ", 3)", sep=""))
    m <- 4
  }
  df_spline <- NULL
  for(i in 1:length(k)){
    df_spline <- c(df_spline, degree+i)
    forms[[m]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], "+ splines2::bsp(", var, ", df=", degree+k[i],
                                  ", Boundary.knots=c(", min(data[[var]], na.rm=TRUE),", ", max(data[[var]], na.rm=TRUE), "))", sep=""))
    m <- m+1
  }
  if(crit %in% c("AIC", "BIC")){
    mods <- lapply(forms, function(x)lm(x, data=data))
    stats <- sapply(mods, function(x)do.call(crit, list(object=x)))
  }
  if(crit == "CV"){
    stop("criterion 'CV' not supported")
  }
  # if(crit == "CV"){
  #   tmp.stats <- NULL
  #   for(j in 1:cviter){
  #     mods <- list()
  #     for(i in 1:length(forms)){
  #       tmp <- stats::glm(forms[[i]], data=data, family=gaussian)
  #       tmpdat <- data[rownames(stats::model.frame(tmp)), ]
  #       mods[[i]] <- boot::cv.glm(tmpdat, tmp, K=cvk)
  #     }
  #     tmp.stats <- rbind(tmp.stats, sapply(mods, function(x)x$delta[1]))
  #   }
  #   stats <- colMeans(tmp.stats)
  # }
  if(plot){
    k <- k+3
    if(includePoly){k <- c(1:3, k)}
    plot(k, stats, type="o", pch=16, col="black", xlab="# Degrees of Freedom", ylab = crit)
    graphics::points(k[which.min(stats)], min(stats), pch=16, col="red")
  }else{
    return(data.frame(df = c(df_poly, df_spline), stat=stats))
  }
}



##############################################################################
# NKnotsTest 
##############################################################################

#' Test of functional form assumption using B-splines
#' 
#' Estimate hypothesis test of lower- and higher-order non-linear relationships
#' against an assumed target relationship.
#' 
#' 
#' @param form A formula detailing the model for which smoothing is to be
#' evaluated.
#' @param var A character string identifying the variable for which smoothing
#' is to be evaluated.
#' @param data Data frame providing values of all variables in \code{form}.
#' @param targetdf The assumed degrees of freedom against which the tests will
#' be conducted.
#' @param degree Degree of polynomial in B-spline basis functions.
#' @param min.knots Minimum number of internal B-spline knots to be evaluated.
#' @param max.knots Maximum number of internal B-spline knots to be evaluated.
#' @param adjust Method by which p-values will be adjusted (see
#' \code{\link{p.adjust}})
#' @returns A matrix with the following columns: \item{F}{F statistics of test
#' of candidate models against target model} \item{DF1}{Numerator DF from
#' F-test} \item{DF2}{Denominator DF from F-test} \item{p(F)}{p-value from the
#' F-test} \item{Clarke}{Test statistic from the Clarke test}
#' \item{Pr(Better)}{The Clarke statistic divided by the number of
#' observations} \item{p(Clarke)}{p-value from the Clarke test.  (T) means that
#' the significant p-value is in favor of the Target model and (C) means the
#' significant p-value is in favor of the candidate (alternative) model.}
#' \item{Delta_AIC}{AIC(candidate model) - AIC(target model)}
#' \item{Delta_AICc}{AICc(candidate model) - AICc(target model)}
#' \item{Delta_BIC}{BIC(candidate model) - BIC(target model)}
#' @author Dave Armstrong
#' 
#' @keywords internal
#' @noRd
#' 
#' @examples
#' 
#' data(Prestige, package="carData")
#' NKnotsTest(prestige ~ education + type, var="income", data=na.omit(Prestige), targetdf=3)
#' 
NKnotsTest <- function(form, var, data, targetdf = 1, degree=3, min.knots=1,
                       max.knots=10, adjust="none"){
  k <- seq(min.knots, max.knots, by=1)
  forms <- vector("list", length(k)+3)
  m <- 1
  forms[[1]]<- as.formula(paste(as.character(form)[2], "~",
                                as.character(form)[3], " + ", var, sep=""))
  forms[[2]]<- as.formula(paste(as.character(form)[2], "~",
                                as.character(form)[3], "+ poly(", var,  ", 2)", sep=""))
  forms[[3]]<- as.formula(paste(as.character(form)[2], "~",
                                as.character(form)[3], "+ poly(", var,  ", 3)", sep=""))
  m <- 4
  for(i in 1:length(k)){
    forms[[m]]<- as.formula(paste(as.character(form)[2], "~",
                                  as.character(form)[3], "+ bs(", var, ", df=", degree+k[i], ")", sep=""))
    m <- m+1
  }
  mods <- lapply(forms, function(x)lm(x, data=data))
  mods.df <- c(1:3, k+3)
  
  target.mod <- mods[[which(mods.df == targetdf)]]
  cand.mods <- mods
  cand.mods[[which(mods.df == targetdf)]] <- NULL
  tests <- lapply(cand.mods, function(x)as.matrix(stats::anova(target.mod, x)))
  # tests2 <- lapply(cand.mods, function(x)clarke_test(target.mod, x))
  num.df <- sapply(tests, function(x)abs(diff(x[,1])))
  denom.df <- sapply(tests, function(x)min(x[,1]))
  Fstats <- sapply(tests, function(x)x[2,5])
  pval <- stats::p.adjust(sapply(tests, function(x)x[2,6]), method=adjust)
  # cstats <- sapply(tests2, function(x)x$stat)
  # cprobs <- sapply(tests2, function(x)x$stat/x$nobs)
  # cminstat <- sapply(tests2, function(x)min(x$stat, x$nobs - x$stat))
  # cbetter <- cp <- 2 * stats::pbinom(cminstat, sapply(tests2, function(x)x$nobs), 0.5)
  # pref <- sapply(tests2, function(x)ifelse(x$stat > x$nobs - x$stat,  "(T)", "(C)"))
  # pref <- ifelse(cp > .05, "", pref)
  cstats <- NULL
  cprobs <- NULL
  pref   <- NULL
  cp     <- NULL
  
  delta.aic <- sapply(cand.mods, stats::AIC) - stats::AIC(target.mod)
  delta.bic <- sapply(cand.mods, stats::BIC) - stats::BIC(target.mod)
  # delta.aicc <- sapply(cand.mods, AICcmodavg::AICc) - AICcmodavg::AICc(target.mod)
  delta.aicc <- NULL
  res <- cbind(Fstats, num.df, denom.df, pval, cstats, cprobs, cp, delta.aic, delta.aicc, delta.bic)
  sigchar <- ifelse(res[,4] < .05, "*", " ")
  sigchar2 <- ifelse(res[,7] < .05, "*", " ")
  strres <- NULL
  digs <- c(3,0,0,3, 0, 3, 3, 3, 3, 3)
  for(i in 1:10){
    tmp <- sprintf(paste("%.", digs[i], "f", sep=""), res[,i])
    if(i == 1){
      tmp <- paste(tmp, sigchar, sep="")
    }
    if(i == 5){
      tmp <- paste(tmp, sigchar2, sep="")
    }
    if(i == 7){
      tmp <- paste(tmp, pref,  sep=" ")
    }
    
    strres <- cbind(strres,tmp )
  }
  colnames(strres) <- c("F", "DF1", "DF2", "p(F)", "Clarke", "Pr(Better)", "p(Clarke)", "Delta_AIC", "Delta_AICc", "Delta_BIC")
  rownames(strres) <- paste("DF=", targetdf, " vs. DF=", mods.df[-targetdf], sep="")
  if(targetdf > 1){
    below <- strres[1:(targetdf-1), , drop=F]
    above <- strres[targetdf:nrow(strres),, drop=F]
    strres <- rbind(below, rep("", 10), above)
    rownames(strres)[targetdf] <- "   Target"
  }
  print(strres, quote=FALSE)
}


##############################################################################
# 
##############################################################################




##############################################################################
# 
##############################################################################




##############################################################################
# 
##############################################################################





##############################################################################
# 
##############################################################################





##############################################################################
# 
##############################################################################
