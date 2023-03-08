

#'Optimize \pkg{bsitar} model
#'
#'@param model An object of class \code{bsitar}.
#'
#'@param newdata An optional \code{data.frame} to be used when optimizing the
#'  model. If \code{NULL} (default), the same same data used for original model
#'  fit is used. Note that data-dependent default priors will not be updated
#'  automatically.
#'  
#'@param optimize_df A vector specifying the updation of degree of freedom
#'  \code{df}. If \code{NULL} (default), the \code{df} is taken from the
#'  original model. For \code{univariate-by-sungroup} and \code{multivariate}
#'  models (see [bsitar::bsitar] for details on these arguments),
#'  \code{optimize_df} can be a single integer (e.g., \code{optimize_df = 4}) or
#'  a list (e.g., \code{optimize_df = list(4,5)}). For optimization over
#'  different \code{df}, say for example df 4 and 5 for univariate model, the
#'  corresponding code is \code{optimize_df = list(4,5)}. For a multivariate
#'  model fit to two outcomes with different \code{df}, the optimization over
#'  \code{df} 4 and 5 for the first submodel and 5 and 6 for the second
#'  submodel, the  corresponding \code{optimize_df} code is \code{optimize_df =
#'  kist(list(4,5), list(5,6))}.
#'  
#'@param optimize_x A vector specifying the transformations of predictor
#'  (typically \code{age}) variable (via \code{xvar}). The option are 'NULL',
#'  'log' and  'sqrt' or their combinations. Note that user need not to enclose
#'  these options in either single or double quotes as they are take care of
#'  internally. The default \code{optimize_x = list(NULL, log,  sqrt)} is to
#'  explore all possible combinations of 'NULL', 'log' and  'sqrt'. Similar to
#'  the \code{optimize_df}, user can specify different  \code{optimize_x} for
#'  \code{univariate-by-sungroup} and \code{multivariate} submodels.
#'
#'@param optimize_y A vector specifying the transformations of the outcome
#'  variable (via \code{yvar}). The option are 'NULL', 'log' and  'sqrt' or
#'  their combinations. Note that user need not to enclose these options in
#'  either single or double quotes as they are take care of internally. The
#'  default \code{optimize_y = list(NULL, log,  sqrt)} is to explore all
#'  possible combinations of 'NULL', 'log' and  'sqrt'. Similar to the
#'  \code{optimize_df}, user can specify different  \code{optimize_y} for
#'  \code{univariate-by-sungroup} and \code{multivariate} submodels.
#'
#'@param exclude_default_funs A logical to indicate
#'  whether transformations (\code{xvar} and \code{yvar}) used in the original
#'  model fit should be excluded. If \code{TRUE} (deafult), the the \code{xvar} 
#'  and \code{yvar} transformations specified for the original model fit are
#'  excluded from the \code{optimize_x} and \code{optimize_y}. From example, if
#'  original model is fit with \code{xvar = log} and \code{yvar = NULL}, then
#'  the \code{optimize_x} is translated into \code{optimize_x = list(NULL,
#'  sqrt)} and  \code{optimize_y} as \code{optimize_y = list(log, sqrt)}.
#'
#'
#'@param add_fit_criteria An optional (default \code{NULL}) indicator to add fit
#'  criteria to the model fit. options are \code{loo} and \code{waic}. Please
#'  see [brms::add_criterion()] for details.
#'
#'@param add_fit_bayes_R An optional (default \code{NULL}) to add Bayesian R
#'  square.
#'
#'@param ... Other arguments passed to \code{\link{update_bsitar}}.
#'
#'@return A list containing the optimized models of class \code{brmsfit, bsiatr},
#' and the the combined summary statistics if \code{add_fit_criteria} and/or
#' \code{add_fit_bayes_R} are specified. 
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
  
  # This to evaluate T/F to TRUE/FALSE
  for (i in names(args_o)) {
    if(is.symbol(args_o[[i]])) {
      if(args_o[[i]] == "T") args_o[[i]] <- eval(args_o[[i]])
      if(args_o[[i]] == "F") args_o[[i]] <- eval(args_o[[i]])
    }
  }
  
  # print(sort(names(args_o)))
  # print(args_o$expose_function)
  # xxx <<- args_o$expose_function
  # stop()
  
  for (add_fit_criteriai in add_fit_criteria) {
    if(!add_fit_criteriai %in% c("loo", "waic")) {
      stop("only loo and waic criteria are supported")
    }
  }
  
  for (bayes_Ri in add_fit_bayes_R) {
    if(!bayes_Ri %in% c("bayes_R")) {
      stop("only bayes_R as R square measure is supported")
    }
  }
  # print(args_o$expose_function)
  # stop()
  if(!args_o$expose_function) {
    if(!is.null(add_fit_criteria) | !is.null(add_fit_bayes_R)) {
      stop("Argument expose_function must be set to TRUE when ",
           "\n ",
           " adding fit criteria and bayes_R")
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
    
    numeric_dx <-  is.numeric(eval(parse(text = gsub('\"', "", xxo))))
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
  
  
  add_summary_waic <- function(x, digits = 1) {
    summary_waic <- x
    summary_waic$ pointwise <- NULL
    summary_waic <- summary_waic$ estimates
    summary_waic <- summary_waic %>% data.frame()
    summary_waic <- summary_waic %>% 
      dplyr::mutate(dplyr::across(where(is.numeric), round, digits))
    summary_waic$Parameter <- row.names(summary_waic)
    row.names(summary_waic) <- NULL
    summary_waic <- summary_waic %>% dplyr::relocate(Parameter, Estimate, SE)
    summary_waic
  }
  
  add_summary_loo <- function(x, digits = 1) {
    summary_loo <- x
    summary_loo$ pointwise <- NULL
    summary_loo <- summary_loo $ estimates
    summary_loo <- summary_loo %>% data.frame()
    summary_loo <- summary_loo %>% 
      dplyr::mutate(dplyr::across(where(is.numeric), round, digits))
    summary_loo$Parameter <- row.names(summary_loo)
    row.names(summary_loo) <- NULL
    summary_loo <- summary_loo %>% dplyr::relocate(Parameter, Estimate, SE)
    summary_loo
  }
  
  add_diagnostic_loo <- function(x, digits = 1) {
    summary_loo_diagnostic <- loo::pareto_k_table(x) %>% data.frame()
    row.names(summary_loo_diagnostic) <- NULL
    summary_loo_diagnostic$Range <- attr(loo::pareto_k_table(x), "dimnames")[[1]]
    summary_loo_diagnostic$Inference <-  c('Good', "Ok", "Bad", "Very bad")
    summary_loo_diagnostic$Percent <- round(summary_loo_diagnostic$Proportion*100, digits)
    summary_loo_diagnostic$Min.n_eff  <- summary_loo_diagnostic$Min..n_eff 
    summary_loo_diagnostic$Min.n_eff <- round(summary_loo_diagnostic$Min.n_eff)
    summary_loo_diagnostic <- summary_loo_diagnostic %>% 
      dplyr::select(-c(Proportion, Min..n_eff))
    summary_loo_diagnostic <- summary_loo_diagnostic %>% 
      dplyr::relocate(Range, Inference, Count, Percent, Min.n_eff)
    summary_loo_diagnostic
  }
  
  add_summary_bayes_R <- function(x, digits = 2) {
    summary_bayes_R <- x
    summary_bayes_R <- summary_bayes_R %>% data.frame()
    summary_bayes_R <- summary_bayes_R %>% 
      dplyr::mutate(dplyr::across(where(is.numeric), round, digits))
    summary_bayes_R$Parameter <- row.names(summary_bayes_R)
    row.names(summary_bayes_R) <- NULL
    summary_bayes_R$SE <- summary_bayes_R$Est.Error  
    summary_bayes_R <- summary_bayes_R %>% dplyr::select(-c(Est.Error))
    summary_bayes_R <- summary_bayes_R %>% 
      dplyr::relocate(Parameter, Estimate, SE)
    summary_bayes_R
  }
  
  combine_summaries <- function(model_list, summary_obj) {
    ic = 0
    list_c <- list()
    for (model_listi in 1:length(model_list)) {
      if(!is.null(model_list[[model_listi]][[summary_obj]])) {
        ic <- ic + 1
        list_c[[ic]] <- model_list[[model_listi]][[summary_obj]]
      }
      summary_of_obj <- list_c %>% do.call(rbind, .) %>% data.frame()
    }
    if(nrow(summary_of_obj) < 1) summary_of_obj <- NULL
    summary_of_obj
  }
  
  optimize_fun <- function(.x, model) {
    message("Working on model no. ", .x, " of ", nrow(optimize_df_x_y), " models")
    cat("\n")
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
    args_o$df <- eval(parse(text = df))
    args_o$xfun <- xfun
    args_o$yfun <- yfun
    
    if (!is.null(newdata)) {
      args_o$data <- call_o_args$newdata
    }
    
    # print(sort(names(args_o)))
    # print(args_o$expose_function)
    # xxx <<- args_o$sample_prior
    # stop()
    
    fit <- do.call(update_bsitar, args_o) 
    
    fit$model_info$optimization_info <- optimization_info
    fit$model_info$optimize_df <- df
    fit$model_info$optimize_x <- xfun_print
    fit$model_info$optimize_y <- yfun_print
    
    # Add fit_criteria and bares_R to the fit 
    
    # Very important to set cores = 1 on windows otherwise loo hangs
    if (!is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_criteria, collapse = ", ")
      message(" Adding", " ", what_, " ", "...")
      cat("\n")
      suppressWarnings(fit <-
                         add_criterion(fit, add_fit_criteria, cores = 1))
    }
    
    if (!is.null(add_fit_bayes_R) & !is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_bayes_R, collapse = ", ")
      message(" Adding", " ", what_, " ", "...")
      cat("\n")
      suppressWarnings(fit$criteria$bayes_R <- bayes_R2(fit))
    }
    
    if (!is.null(add_fit_bayes_R) & is.null(add_fit_criteria)) {
      what_ <- paste(add_fit_bayes_R, collapse = ", ")
      message(" Adding", " ", what_, " ", "...")
      cat("\n")
      suppressWarnings(fit$bayes_R <- bayes_R2(fit))
    }
    
    # Add summary data frames for criteria and R square
    
    if('loo' %in% add_fit_criteria) {
      summary_loo <- add_summary_loo(fit$criteria$loo, digits = 1)
      summary_loo$df <- df
      summary_loo$xfun <- xfun_print
      summary_loo$yfun <- yfun_print
      summary_loo <- summary_loo %>% dplyr::relocate(df, xfun, yfun)
      diagnostic_loo <- add_diagnostic_loo(fit$criteria$loo, digits = 1)
      diagnostic_loo$df <- df
      diagnostic_loo$xfun <- xfun_print
      diagnostic_loo$yfun <- yfun_print
      diagnostic_loo <- diagnostic_loo %>% dplyr::relocate(df, xfun, yfun)
      fit$summary_loo <- summary_loo
      fit$diagnostic_loo <- diagnostic_loo
    }
    if('waic' %in% add_fit_criteria) {
      summary_waic <- add_summary_waic(fit$criteria$waic, digits = 1)
      summary_waic$df <- df
      summary_waic$xfun <- xfun_print
      summary_waic$yfun <- yfun_print
      summary_waic <- summary_waic %>% dplyr::relocate(df, xfun, yfun)
      fit$summary_waic <- summary_waic
    }
    if('bayes_R' %in% add_fit_bayes_R) {
      summary_bayes_R <- add_summary_bayes_R(fit$criteria$bayes_R, digits = 1)
      summary_bayes_R$df <- df
      summary_bayes_R$xfun <- xfun_print
      summary_bayes_R$yfun <- yfun_print
      summary_bayes_R <- summary_bayes_R %>% dplyr::relocate(df, xfun, yfun)
      fit$summary_bayes_R <- summary_bayes_R
    }
    
    return(fit)
  }
  
  optimize_list <- lapply(1:nrow(optimize_df_x_y), function(.x)
    optimize_fun(.x, model))
  
  # summary_loo_all     <- combine_summaries(optimize_list, 'summary_loo')
  # diagnostic_loo_all  <- combine_summaries(optimize_list, 'diagnostic_loo')
  # summary_waic_all    <- combine_summaries(optimize_list, 'summary_waic')
  # summary_bayes_R_all <- combine_summaries(optimize_list, 'summary_bayes_R')
  
  loo_summary     <- combine_summaries(optimize_list, 'summary_loo')
  loo_diagnostic  <- combine_summaries(optimize_list, 'diagnostic_loo')
  waic_summary    <- combine_summaries(optimize_list, 'summary_waic')
  bayes_R_summary <- combine_summaries(optimize_list, 'summary_bayes_R')
  
  list(models = optimize_list, loo_summary = loo_summary,
       loo_diagnostic = loo_diagnostic, waic_summary = waic_summary,
       bayes_R_summary = bayes_R_summary)
}




#' @rdname optimize_bsitar.bsitar
#' @export
optimize_bsitar <- function(model, ...) {
  UseMethod("optimize_bsitar")
}

