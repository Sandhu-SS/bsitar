

#'Optimize \pkg{bsitar} model
#'
#'@param model An object of class \code{bsitar}.
#'
#'@param newdata An optional \code{data.frame} to be used when optimizing the
#'  model. If \code{NULL} (default), the same same data used for original model
#'  fit is used. Note that data-dependent default priors will not be updated
#'  automatically.
#'
#'@param optimize_df A vector specifying the degree of freedom (\code{df}) 
#'  update. If \code{NULL} (default), the \code{df} is taken from the
#'  original model. For \code{univariate-by-sungroup} and \code{multivariate}
#'  models (see [bsitar::bsitar()] for details on these arguments),
#'  \code{optimize_df} can be a single integer (e.g., \code{optimize_df = 4}) or
#'  a list (e.g., \code{optimize_df = list(4,5)}). For optimization over
#'  different \code{df}, say for example df 4 and 5 for univariate model, the
#'  corresponding code is \code{optimize_df = list(4,5)}. For a multivariate
#'  model fit to two outcomes with different \code{df}, the optimization over
#'  \code{df} 4 and 5 for the first sub model and 5 and 6 for the second
#'  sub model, the corresponding \code{optimize_df} code is \code{optimize_df =
#'  list(list(4,5), list(5,6))} i.e, a list of lists.
#'
#'@param optimize_x A vector specifying the transformations of predictor
#'  (typically \code{age}) variable (via \code{xvar}). The option are 
#'  \code{NULL}, \code{log}, \code{sqrt} and their combinations. Note that user 
#'  need not to enclose these options in a single or double quotes as they are 
#'  take care of internally. The default setting is to explore all possible 
#'  combination i.e., \code{optimize_x = list(NULL, log,  sqrt)}. Similar to
#'  the \code{optimize_df}, user can specify different \code{optimize_x} for
#'  \code{univariate-by-sungroup} and \code{multivariate} sub models.
#'
#'@param optimize_y A vector specifying the transformations for the response
#'  variable (via \code{yvar}). The approach and options available for 
#'  \code{optimize_y} are identical to the \code{optimize_x} (see above).
#'
#'@param exclude_default_funs A logical to indicate whether transformations for
#'  (\code{xvar} and \code{yvar}) used in the original model fit should be
#'  excluded. If \code{TRUE} (default), the \code{xvar} and \code{yvar}
#'  transformations specified for the original model fit are excluded from the
#'  \code{optimize_x} and \code{optimize_y}. From example, if original model is
#'  fit with \code{xvar = log} and \code{yvar = NULL}, then \code{optimize_x} is
#'  translated into \code{optimize_x = list(NULL, sqrt)} and \code{optimize_y}
#'  as \code{optimize_y = list(log, sqrt)}.
#'
#'
#'@param add_fit_criteria An optional (default \code{NULL}) indicator to add fit
#'  criteria to the model fit. options are \code{loo} and \code{waic}. Please
#'  see [brms::add_criterion()] for details.
#'
#'@param add_fit_bayes_R An optional (default \code{NULL}) to add Bayesian R
#'  square.
#'
#'@param byresp A logical (default \code{FALSE}) to indicate if response wise
#'  fit criteria to be calculated. This argument is evaluated only for the
#'  \code{multivariate} model for which options are available for joint
#'  calculation of pointwise log likelihood or response specific. For,
#'  \code{univariate-by-subgroup} model, the only option available is to
#'  calculate separate pointwise log likelihood for each sub-model.
#'
#'@param digits An integer to set the number of decimal places.
#'
#'@param ... Other arguments passed to \code{\link{update_bsitar}}.
#'
#'@return A list containing the optimized models of class \code{brmsfit,
#'  bsiatr}, and the the combined summary statistics if \code{add_fit_criteria}
#'  and/or \code{add_fit_bayes_R} are specified.
#'
#'@author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#'@importFrom loo pareto_k_table
#'
#'@export optimize_bsitar.bsitar
#'
#'@export
#'
#' @examples
#' \dontrun{
#' data(heights)
#' data_males <- heights %>% filter(sex == 'Male')
#' fit_males <- bsitar(x=age, y=height, id=id, data=heights, df=4)
#' fit_males2 <- optimize_bsitar(fit_males)
#' }
#'
optimize_bsitar.bsitar <- function(model,
                                   newdata = NULL,
                                   optimize_df = NULL,
                                   optimize_x = list(NULL, log,  sqrt),
                                   optimize_y = list(NULL, log,  sqrt),
                                   exclude_default_funs = TRUE,
                                   add_fit_criteria = NULL,
                                   add_fit_bayes_R = NULL,
                                   byresp = FALSE,
                                   digits = 2,
                                   ...) {
  call_o <- match.call()
  call_o_args <- as.list(call_o)[-1]
  
  args_o <- as.list(model$model_info$call.full.bsitar)[-1]
  args_o_dots_ <- list(...)
  if (length(args_o_dots_) > 0) {
    for (i in names(args_o_dots_)) {
      args_o[[i]] <- args_o_dots_[[i]]
    }
  }
  
  # Global R cmd check
  outcome <- xfun <- yfun <- NULL
  
  
  # This to evaluate T/F to TRUE/FALSE
  for (i in names(args_o)) {
    if (is.symbol(args_o[[i]])) {
      if (args_o[[i]] == "T")
        args_o[[i]] <- eval(args_o[[i]])
      if (args_o[[i]] == "F")
        args_o[[i]] <- eval(args_o[[i]])
    }
  }
  
  # print(sort(names(args_o)))
  # print(args_o$expose_function)
  # xxx <<- args_o$expose_function
  # stop()
  
  for (add_fit_criteriai in add_fit_criteria) {
    if (!add_fit_criteriai %in% c("loo", "waic")) {
      stop("only loo and waic criteria are supported")
    }
  }
  
  for (bayes_Ri in add_fit_bayes_R) {
    if (!bayes_Ri %in% c("bayes_R2")) {
      stop("only bayes_R2 as R square measure is supported")
    }
  }
 
  if (!args_o$expose_function) {
    if (!is.null(add_fit_criteria) | !is.null(add_fit_bayes_R)) {
      stop(
        "Argument expose_function must be set to TRUE when ",
        "\n ",
        " adding fit criteria and bayes_R2"
      )
    }
  }
  
  get_args_opt <- function(xo) {
    get_within_fist_last_paranthesese <- function(x__) {
      x__ <- sub('\\(', '[', x__)
      x__ <- sub("\\)([^)]*)$", "]\\1", x__)
      x__ <-
        gsub("[\\[\\]]", "", regmatches(x__, gregexpr("\\[.*?\\]", x__))[[1]])
      x__ <- gsub("\\[|\\]", "", x__)
      x__
    }
    gsub_comma_within_paranthesese <-
      function(x__, replace_comma_by) {
        tt <-
          gsub("[\\(\\)]", "", regmatches(x__, gregexpr("\\(.*?\\)", x__))[[1]])
        tt2 <- gsub(",", replace_comma_by, tt, fixed = T)
        j <- 0
        for (i in tt) {
          j <- j + 1
          x__ <- gsub(tt[j], tt2[j], x__, fixed = T)
        }
        x__
      }
    xxo <- gsub("[[:space:]]", "", xo)
    
    numeric_dx <-
      is.numeric(eval(parse(text = gsub('\"', "", xxo))))
    if (xxo != "NULL" & xxo != "\"NULL\"" & !numeric_dx) {
      xxo <- get_within_fist_last_paranthesese(xxo)
      xxo <- gsub_comma_within_paranthesese(xxo, "_comma_")
      xxo <- strsplit(xxo, ",")[[1]]
      xxo <- gsub("_comma_" , ",", xxo)
      xxo <- gsub('\"', "", xxo)
    } else {
      xxo <- xxo
      xxo <- gsub('\"', "", xxo)
    }
    xxo
  }
  
  optimize_df <- get_args_opt(deparse(substitute(optimize_df)))
  optimize_x  <- get_args_opt(deparse(substitute(optimize_x)))
  optimize_y  <- get_args_opt(deparse(substitute(optimize_y)))
  
  if (exclude_default_funs) {
    optimize_x <- optimize_x[!optimize_x %in% model$model_info$xfuns]
    optimize_y <-
      optimize_y[!optimize_y %in% model$model_info$xfuns]
    if (identical(optimize_x, character(0)))
      optimize_x <- "NULL"
    if (identical(optimize_y, character(0)))
      optimize_y <- "NULL"
  }
  
  optimize_df_x_y <-
    expand.grid(optimize_df, optimize_x, optimize_y)
  
  colnames(optimize_df_x_y) <- c("df", "xfun", "yfun")
  
  add_summary_waic <- NULL
  Count <- Est.Error <- Inference <- Min..n_eff <- where <- NULL
  Min.n_eff <- Percent <- Proportion <- Range <- SE <- NULL
  
  
  
  
  combine_summaries <- function(model_list, summary_obj) {
    ic = 0
    list_c <- list()
    for (model_listi in 1:length(model_list)) {
      if (!is.null(model_list[[model_listi]][[summary_obj]])) {
        ic <- ic + 1
        list_c[[ic]] <- model_list[[model_listi]][[summary_obj]]
      }
      summary_of_obj <-
        list_c %>% do.call(rbind, .) %>% data.frame()
    }
    if (nrow(summary_of_obj) < 1)
      summary_of_obj <- NULL
    summary_of_obj
  }
  
  
  
  
  
  # resp = NULL is only placeholder that too only for multivariate
  # if NULL, then combined log likelihood used for multivariate model
  # if anything else e.g., resp = 'NULL' or anything, '
  # then separate likelihood for responses
  
  add_citeria_fun <- function(fit,
                              add_fit_criteria = NULL,
                              add_fit_bayes_R = NULL,
                              resp = NULL,
                              digits = 2,
                              df,
                              xfun_print,
                              yfun_print) {
    if (!is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_criteria, collapse = ", ")
      message(" Adding", " ", what_, " ", "...")
      cat("\n")
      if (is.na(fit$model_info$univariate_by) |
          !fit$model_info$multivariate) {
        if (!fit$model_info$multivariate) {
          suppressWarnings(fit <- add_criterion(fit,
                                                add_fit_criteria, cores = 1))
        }
        if (fit$model_info$multivariate) {
          if (is.null(resp)) {
            suppressWarnings(fit <- add_criterion(fit,
                                                  add_fit_criteria, cores = 1))
          }
          if (!is.null(resp)) {
            for (aci in fit$model_info$ys) {
              suppressWarnings(fit <- add_criterion(
                fit,
                add_fit_criteria,
                resp = aci,
                cores = 1
              ))
              aci_names <- paste0(names(fit$criteria), aci)
              names(fit$criteria) <- aci_names
            }
            aci_names <- c()
            for (aci in fit$model_info$ys) {
              aci_names <- c(aci_names, paste0(add_fit_criteria, aci))
            }
            names(fit$criteria) <- aci_names
          }
        }
      }
      
      if (!is.na(fit$model_info$univariate_by)) {
        for (aci in fit$model_info$ys) {
          suppressWarnings(fit <- add_criterion(
            fit,
            add_fit_criteria,
            resp = aci,
            cores = 1
          ))
          aci_names <- paste0(names(fit$criteria), aci)
          names(fit$criteria) <- aci_names
        }
        aci_names <- c()
        for (aci in fit$model_info$ys) {
          aci_names <- c(aci_names, paste0(add_fit_criteria, aci))
        }
        names(fit$criteria) <- aci_names
      }
    } # if (!is.null(add_fit_criteria))
    
    
    if (!is.null(add_fit_bayes_R)) {
      what_ <- paste(add_fit_bayes_R, collapse = ", ")
      message(" Adding", " ", what_, " ", "...")
      cat("\n")
      if (is.na(fit$model_info$univariate_by)) {
        if (!fit$model_info$multivariate) {
          aci_names <- paste0(add_fit_bayes_R, '')
          suppressWarnings(fit$criteria[[aci_names]] <-
                             bayes_R2(fit, cores = 1))
          fit$criteria[[aci_names]] <-
            fit$criteria[[aci_names]] %>%
            data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
            dplyr::relocate(Parameter)
          rownames(fit$criteria[[aci_names]]) <- NULL
        }
        if (fit$model_info$multivariate) {
          if (is.null(resp)) {
            aci_names <- paste0(add_fit_bayes_R, '')
            suppressWarnings(fit$criteria[[aci_names]] <-
                               bayes_R2(fit,
                                        cores = 1))
            fit$criteria[[aci_names]] <-
              fit$criteria[[aci_names]] %>%
              data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
              dplyr::relocate(Parameter)
            rownames(fit$criteria[[aci_names]]) <- NULL
          }
          if (!is.null(resp)) {
            for (aci in fit$model_info$ys) {
              aci_names <- paste0(add_fit_bayes_R, aci)
              suppressWarnings(fit$criteria[[aci_names]] <-
                                 bayes_R2(fit,
                                          resp = aci,
                                          cores = 1))
              fit$criteria[[aci_names]] <-
                fit$criteria[[aci_names]] %>%
                data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
                dplyr::relocate(Parameter)
              rownames(fit$criteria[[aci_names]]) <- NULL
            }
          }
        }
      }
      
      
      
      
      if (!is.na(fit$model_info$univariate_by)) {
        for (aci in fit$model_info$ys) {
          aci_names <- paste0(add_fit_bayes_R, aci)
          suppressWarnings(fit$criteria[[aci_names]] <-
                             bayes_R2(fit,
                                      resp = aci,
                                      cores = 1))
          fit$criteria[[aci_names]] <-
            fit$criteria[[aci_names]] %>%
            data.frame() %>% dplyr::mutate(Parameter = rownames(.)) %>%
            dplyr::relocate(Parameter)
          rownames(fit$criteria[[aci_names]]) <- NULL
        }
      }
      
      # xx <- fit$criteria$bayes_R %>% data.frame()
      # names(xx) <- sub('^bayes_R.', '', names(xx))
      # xx$Parameter <- row.names(xx)
      # row.names(xx) <- NULL
      # xx <- xx %>% dplyr::relocate(Parameter)
      # fit$criteria$bayes_R <- xx
    } # if (!is.null(add_fit_bayes_R)) {
    
    
    
    ################
    add_summary_waic <- function(x, round_digits = 1) {
      summary_waic <- x
      summary_waic$pointwise <- NULL
      summary_waic <- summary_waic$estimates
      summary_waic <- summary_waic %>% data.frame()
      summary_waic <- summary_waic %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, round_digits))
      summary_waic$Parameter <- row.names(summary_waic)
      row.names(summary_waic) <- NULL
      summary_waic <-
        summary_waic %>% dplyr::relocate(Parameter, Estimate, SE)
      summary_waic
    }
    
    add_summary_bayes_R2 <- function(x, round_digits = 2) {
      summary_bayes_R <- x
      summary_bayes_R <- summary_bayes_R %>% data.frame()
      summary_bayes_R <- summary_bayes_R %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, round_digits))
      # summary_bayes_R$Parameter <- row.names(summary_bayes_R)
      row.names(summary_bayes_R) <- NULL
      summary_bayes_R$SE <- summary_bayes_R$Est.Error
      summary_bayes_R <-
        summary_bayes_R %>% dplyr::select(-c(Est.Error))
      summary_bayes_R <- summary_bayes_R %>%
        dplyr::relocate(Parameter, Estimate, SE)
      summary_bayes_R
    }
    
    
    
    add_summary_loo <- function(x, round_digits = 1) {
      summary_loo <- x
      summary_loo$pointwise <- NULL
      summary_loo <- summary_loo$estimates
      summary_loo <- summary_loo %>% data.frame()
      summary_loo <- summary_loo %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, round_digits))
      summary_loo$Parameter <- row.names(summary_loo)
      row.names(summary_loo) <- NULL
      summary_loo <-
        summary_loo %>% dplyr::relocate(Parameter, Estimate, SE)
      summary_loo
    }
    
    add_diagnostic_loo <- function(x, round_digits = 1) {
      summary_loo_diagnostic <- loo::pareto_k_table(x) %>% data.frame()
      row.names(summary_loo_diagnostic) <- NULL
      summary_loo_diagnostic$Range <- attr(loo::pareto_k_table(x),
                                           "dimnames")[[1]]
      summary_loo_diagnostic$Inference <-
        c('Good', "Ok", "Bad", "Very bad")
      summary_loo_diagnostic$Percent <-
        round(summary_loo_diagnostic$Proportion * 100, round_digits)
      summary_loo_diagnostic$Min.n_eff  <-
        summary_loo_diagnostic$Min..n_eff
      summary_loo_diagnostic$Min.n_eff <-
        round(summary_loo_diagnostic$Min.n_eff)
      summary_loo_diagnostic <- summary_loo_diagnostic %>%
        dplyr::select(-c(Proportion, Min..n_eff))
      summary_loo_diagnostic <- summary_loo_diagnostic %>%
        dplyr::relocate(Range, Inference, Count, Percent, Min.n_eff)
      summary_loo_diagnostic
    }
    
    
    if ('waic' %in% add_fit_criteria) {
      err. <- FALSE
      tryCatch(
        expr = {
          if (!is.na(fit$model_info$univariate_by)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0('waic', aci)
              list_c_[[aci]] <-
                add_summary_waic(fit$criteria[[getit_]], 
                                 round_digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_waic <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & !is.null(resp)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0('waic', aci)
              list_c_[[aci]] <-
                add_summary_waic(fit$criteria[[getit_]], 
                                 round_digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_waic <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & is.null(resp)) {
            getit_ <- paste0('waic', '')
            summary_waic <-
              add_summary_waic(fit$criteria[[getit_]], round_digits = digits)
          } else if (is.na(fit$model_info$univariate_by) &
                     !fit$model_info$multivariate) {
            getit_ <- paste0('waic', '')
            summary_waic <-
              add_summary_waic(fit$criteria[[getit_]], round_digits = digits)
          }
          summary_waic$df <- df
          summary_waic$xfun <- xfun_print
          summary_waic$yfun <- yfun_print
          summary_waic <-
            summary_waic %>% dplyr::relocate(df, xfun, yfun)
          rownames(summary_waic) <- NULL
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (err.) {
        summary_waic <- NULL
      } else {
        summary_waic <- summary_waic
      }
      fit$summary_waic <- summary_waic
    }
    
    
    
    
    
    if ('bayes_R2' %in% add_fit_bayes_R) {
      err. <- FALSE
      tryCatch(
        expr = {
          if (!is.na(fit$model_info$univariate_by)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0(add_fit_bayes_R, aci)
              list_c_[[aci]] <-
                add_summary_bayes_R2(fit$criteria[[getit_]], 
                                     round_digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_bayes_R2 <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & !is.null(resp)) {
            list_c_ <- list()
            for (aci in fit$model_info$ys) {
              getit_ <- paste0(add_fit_bayes_R, aci)
              list_c_[[aci]] <-
                add_summary_bayes_R2(fit$criteria[[getit_]], 
                                     round_digits = digits) %>%
                dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
            }
            summary_bayes_R2 <-
              list_c_ %>%  do.call(rbind, .) %>% data.frame()
          } else if (fit$model_info$multivariate & is.null(resp)) {
            getit_ <- paste0('bayes_R2', '')
            summary_bayes_R2 <-
              add_summary_bayes_R2(fit$criteria[[getit_]], 
                                   round_digits = digits)
          } else if (is.na(fit$model_info$univariate_by) &
                     !fit$model_info$multivariate) {
            getit_ <- paste0('bayes_R2', '')
            summary_bayes_R2 <-
              add_summary_bayes_R2(fit$criteria[[getit_]], 
                                   round_digits = digits)
          }
          summary_bayes_R2$df <- df
          summary_bayes_R2$xfun <- xfun_print
          summary_bayes_R2$yfun <- yfun_print
          summary_bayes_R2 <-
            summary_bayes_R2 %>% dplyr::relocate(df, xfun, yfun)
          rownames(summary_bayes_R2) <- NULL
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (err.) {
        summary_bayes_R2 <- NULL
      } else {
        summary_bayes_R2 <- summary_bayes_R2
      }
      fit$summary_bayes_R2 <-
        summary_bayes_R2 %>% dplyr::select(-Parameter)
    }
    
    
    
    if ('loo' %in% add_fit_criteria) {
      if ('loo' %in% add_fit_criteria) {
        err. <- FALSE
        tryCatch(
          expr = {
            if (!is.na(fit$model_info$univariate_by)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_summary_loo(fit$criteria[[getit_]], 
                                  round_digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              summary_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       !is.null(resp)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_summary_loo(fit$criteria[[getit_]], 
                                  round_digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              summary_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       is.null(resp)) {
              getit_ <- paste0('loo', '')
              
              summary_loo <-
                add_summary_loo(fit$criteria[[getit_]], round_digits = digits)
            } else if (is.na(fit$model_info$univariate_by) &
                       !fit$model_info$multivariate) {
              getit_ <- paste0('loo', '')
              summary_loo <-
                add_summary_loo(fit$criteria[[getit_]], round_digits = digits)
            }
            summary_loo$df <- df
            summary_loo$xfun <- xfun_print
            summary_loo$yfun <- yfun_print
            summary_loo <-
              summary_loo %>% dplyr::relocate(df, xfun, yfun)
            rownames(summary_loo) <- NULL
          },
          error = function(e) {
            err. <<- TRUE
          }
        )
        if (err.) {
          summary_loo <- NULL
        } else {
          summary_loo <- summary_loo
        }
        fit$summary_loo <- summary_loo
      }
      
      if ('loo' %in% add_fit_criteria) {
        err. <- FALSE
        tryCatch(
          expr = {
            if (!is.na(fit$model_info$univariate_by)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_diagnostic_loo(fit$criteria[[getit_]], 
                                     round_digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              diagnostic_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       !is.null(resp)) {
              list_c_ <- list()
              for (aci in fit$model_info$ys) {
                getit_ <- paste0('loo', aci)
                list_c_[[aci]] <-
                  add_diagnostic_loo(fit$criteria[[getit_]], 
                                     round_digits = digits) %>%
                  dplyr::mutate(outcome = aci) %>% dplyr::relocate(outcome)
              }
              diagnostic_loo <-
                list_c_ %>%  do.call(rbind, .) %>% data.frame()
            } else if (fit$model_info$multivariate &
                       is.null(resp)) {
              getit_ <- paste0('loo', '')
              diagnostic_loo <-
                add_diagnostic_loo(fit$criteria[[getit_]], 
                                   round_digits = digits)
            } else if (is.na(fit$model_info$univariate_by) &
                       !fit$model_info$multivariate) {
              diagnostic_loo <-
                add_diagnostic_loo(fit$criteria[[getit_]], 
                                   round_digits = digits)
            }
            diagnostic_loo$df <- df
            diagnostic_loo$xfun <- xfun_print
            diagnostic_loo$yfun <- yfun_print
            diagnostic_loo <-
              diagnostic_loo %>% dplyr::relocate(df, xfun, yfun)
            rownames(diagnostic_loo) <- NULL
          },
          error = function(e) {
            err. <<- TRUE
          }
        )
        if (err.) {
          diagnostic_loo <- NULL
        } else {
          diagnostic_loo <- diagnostic_loo
        }
        fit$diagnostic_loo <- diagnostic_loo
      }
    } # if('loo' %in% add_fit_criteria) {
    
    return(fit)
  } # add_citeria_fun
  
  
  
  
  
  
  
  
  
  
  
  
  optimize_fun <- function(.x, model) {
    message("\nOptimizing model no. ",
            .x,
            " (total ",
            nrow(optimize_df_x_y),
            " models)")
    # cat("\n")
    exe_row <- optimize_df_x_y[.x, ]
    df <- levels(droplevels(exe_row$df))
    xfun <- levels(droplevels(exe_row$xfun))
    yfun <- levels(droplevels(exe_row$yfun))
    if (df == 'NULL')
      df <-
      paste0("list(", paste(model$model_info$dfs, collapse = ","), ")")
    else
      df <- df
    if (xfun == 'NULL')
      xfun <- NULL
    else
      xfun <- xfun
    if (yfun == 'NULL')
      yfun <- NULL
    else
      yfun <- yfun
    
    if (is.null(xfun))
      xfun_print <- deparse(xfun)
    else
      xfun_print <- xfun
    if (is.null(yfun))
      yfun_print <- deparse(yfun)
    else
      yfun_print <- yfun
    
    cat("\n")
    cat(paste0("df = ", df, "; xfun = ", xfun_print, "; yfun = ", yfun_print),
        "\n")
    
    optimization_info <-
      paste0("df = ", df, "; xfun = ", xfun_print, "; yfun = ", yfun_print)
    
    
    args_o$model <- model
    args_o$df    <- eval(parse(text = df))
    args_o$xfun  <- xfun
    args_o$yfun  <- yfun
    
    if (!is.null(newdata)) {
      args_o$data <- call_o_args$newdata
    }
    
    # Somehow update_bsitar is not accepting symbol as arg e.g, x = age
    
    # fit <- do.call(update_bsitar, args_o)
    fit <- do.call(bsitar, args_o)
    
    
    fit$model_info$optimization_info <- optimization_info
    fit$model_info$optimize_df <- df
    fit$model_info$optimize_x <- xfun_print
    fit$model_info$optimize_y <- yfun_print
    
    # Add fit_criteria and bares_R to the fit
    # Add summary data frames for criteria and R square
    
    # setresp to anything so that even multivariate will be response wise
    # if desired, this behavious
    # if(length(fit$model_info$ys) == 1) setresp <- NULL
    # if(length(fit$model_info$ys) > 1) setresp <- 'TRUE'
    
    if (fit$model_info$multivariate) {
      if (byresp) {
        setresp <- 'TRUE'
      } else if (!byresp) {
        setresp <- NULL
      }
    } else if (!fit$model_info$multivariate) {
      setresp <- NULL
    }
    
    if (!is.null(add_fit_criteria)) {
      fit <- add_citeria_fun(
        fit,
        add_fit_criteria = add_fit_criteria,
        add_fit_bayes_R =  NULL,
        resp = setresp,
        digits = digits,
        df = df,
        xfun_print = xfun_print,
        yfun_print = yfun_print
      )
    }
    
    if (!is.null(add_fit_bayes_R)) {
      fit <- add_citeria_fun(
        fit,
        add_fit_criteria = NULL,
        add_fit_bayes_R =  add_fit_bayes_R,
        resp = setresp,
        digits = digits,
        df = df,
        xfun_print = xfun_print,
        yfun_print = yfun_print
      )
    }
    return(fit)
  }
  
  optimize_list <- lapply(1:nrow(optimize_df_x_y), function(.x)
    optimize_fun(.x, model))
  
  
  loo             <- combine_summaries(optimize_list, 'summary_loo')
  loo_diagnostic  <-
    combine_summaries(optimize_list, 'diagnostic_loo')
  waic            <-
    combine_summaries(optimize_list, 'summary_waic')
  bayes_R2        <-
    combine_summaries(optimize_list, 'summary_bayes_R2')
  
  list(
    models = optimize_list,
    loo = loo,
    loo_diagnostic = loo_diagnostic,
    waic = waic,
    bayes_R2 = bayes_R2
  )
  
}




#' @rdname optimize_bsitar.bsitar
#' @export
optimize_bsitar <- function(model, ...) {
  UseMethod("optimize_bsitar")
}
