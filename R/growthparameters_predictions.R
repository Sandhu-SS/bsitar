




#' Compute growth parameters for the \code{bgmfit} model
#' 
#'@description The \code{growthparameters_predictions} function provides
#'  estimates of growth parameters (such as peak velocity and the age at peak
#'  velocity) derived from the growth velocity curve. This function is a wrapper
#'  around the [marginaleffects::predictions()] and
#'  [marginaleffects::avg_predictions()]. A detailed explanation is available at
#'  <https://marginaleffects.com/articles/predictions.html> and
#'  <https://marginaleffects.com/>.
#'  
#' @details The \code{growthparameters_predictions} function estimates and
#'   returns the following growth paramaeters:
#' \itemize{
#'   \item pv  - peak velocity
#'   \item apv - age at peak velocity
#'   \item tv  - take off velocity
#'   \item atv - age at take off velocity
#'   \item cv  - cessation velocity
#'   \item acv - age at cessation velocity
#' }
#' 
#' The take off velocity is the minimum velocity before the peak velocity and it
#' indicates the beginning of the pubertal growth spurt. The cessation velocity
#' indicates the end of the active pubertal growth spurt and is calculated as
#' the percentage of the peak velocity (\code{pv}). The percentage of the peak
#' velocity used in the calculation of the the cessation velocity (\code{'cv'})
#' is controlled via the \code{acg_velocity} argument. The \code{acg_velocity}
#' takes a real value as an input (greater than 0 and less than 1). Typically, a
#' 10 percent of \code{pv} (i.e., \code{acg_velocity = 0.1}) is considered as a
#' good indicator of the cessation of the active pubertal growth spurt
#' \insertCite{Anna2022}{bsitar}.
#' 
#' 
#' @param model An object of class \code{bgmfit}.
#' 
#' @param parameter A single character string, or a vector of character strings 
#' indicating the parameter to be estimated. Options available are \code{ato}, 
#' \code{vto}, \code{vpv}, \code{apv}, \code{vcg}, and \code{acg}. For 
#' \code{parameter = NULL} (default), age at peak growth velocity is estimated.
#' 
#' @param acg_velocity A real number to set the percentage of peak growth growth
#'   velocity as the cessation velocity when estimating the \code{vcg} and
#'   \code{acg} growth parameters. The \code{acg_velocity} should be greater
#'   than \code{0} and less than \code{1}. The default \code{acg_velocity =
#'   0.10} indicates that a 10 per cent of the peak growth velocity will be used
#'   to get the cessation velocity and the corresponding age at the cessation
#'   velocity. For example if peak growth velocity estimate is \code{10
#'   mm/year}, then cessation velocity is \code{1 mm/year}.
#' 
#' @param digits An integer (default \code{2}) to set the decimal places for the
#'   estimated growth parameters. The \code{digits} is passed on to the
#'   [base::round()] function.
#'   
#' @param average A logical to indicate whether to internally call the
#'    [marginaleffects::predictions()] or the
#'    [marginaleffects::avg_predictions()] function. If \code{FALSE} (default),
#'    [marginaleffects::predictions()] is called otherwise
#'    [marginaleffects::avg_predictions()] when \code{average = TRUE}.
#' 
#'@param at_age An integer to indicate the age at which size, (distance) and
#'  velocity to be calculated. For the \code{'logistic3'} model, the possible
#'  \code{at_age} are \code{1} (first maximum growth velocity), \code{2} (at
#'  childhood growth  spurt), and \code{3} (pubertal growth spurt). For
#'  \code{at_age = NULL} (default), parameters are calculated at the pubertal
#'  growth spurt, i.e., \code{at_age = 3}.
#'
#'@param tsv An integer to indicate whether to get the timing, size, (distance)
#'  at the given timing or velocity at give timing. If \code{tsv = NULL}
#'  (default), timing is returned otherwise distance if \code{tsv = 0}, or
#'  velocity if \code{tsv = 1}.
#' 
#' @param reformat A logical (default \code{TRUE}) to reformat returned output
#' as a data.frame with colnames renamed as follows: \code{estimate} as 
#' \code{Estimate}, \code{conf.low} as \code{Q2.5}, and \code{conf.high} as 
#' \code{Q97.5} (assuming that \code{conf_int = 0.95}). Also, following 
#' columns are dropped from the data frame: \code{term}, \code{contrast},
#' \code{tmp_idx}, \code{predicted_lo}, \code{predicted_hi}, \code{predicted}.
#' 
#' @inheritParams  growthparameters_comparison.bgmfit
#' @inheritParams  growthparameters.bgmfit
#' @inheritParams  marginaleffects::predictions
#' @inheritParams  marginaleffects::avg_predictions
#'
#' @return A data frame objects with estimates and CIs for computed parameter(s)
#' 
#' @export growthparameters_predictions.bgmfit
#' 
#' @export
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' \dontrun{
#' growthparameters_predictions(model, parameter = 'apv')
#' }
#' 
growthparameters_predictions.bgmfit <- function(model,
                                                resp = NULL,
                                                ndraws = NULL,
                                                draw_ids = NULL,
                                                newdata = NULL,
                                                re_formula = NA,
                                                at_age = NULL,
                                                tsv = NULL,
                                                parameter = NULL,
                                                xrange = 1,
                                                acg_velocity = 0.10,
                                                digits = 2,
                                                levels_id = NULL,
                                                ipts = NULL,
                                                seed = 123,
                                                future = FALSE,
                                                future_session = 'multisession',
                                                cores = NULL,
                                                average = FALSE, 
                                                variables = NULL,
                                                type = NULL,
                                                by = FALSE,
                                                conf_level = 0.95,
                                                transform = NULL,
                                                wts = NULL,
                                                hypothesis = NULL,
                                                equivalence = NULL,
                                                reformat = TRUE,
                                                envir = parent.frame(),
                                                ...) {
  
  if(system.file(package='tidyr') == "") {
    stop("Please install 'tidyr' package before calling the function",
         "\n ",
         "'growthparameters_predictions'")
  }
  
  if (is.null(ndraws))
    ndraws  <- brms::ndraws(model)
  else
    ndraws <- ndraws
  
  
  if(is.null(at_age)) {
    if(model$model_info$select_model == "logistic3") {
      at_age <- 3
    } else if(model$model_info$select_model == "logistic2") {
      
    }
  }
  
  
  
  if(!is.null(tsv)) {
    if(tsv == 0) {
      deriv <- 0 
    } else if(tsv == 1) {
      deriv <- 1
    } else {
      stop("'tsv' should be between 0 and 1")
    }
  } 
  
  if(is.null(tsv)) {
    deriv <- 0 # placeholder
  }

  o <- post_processing_checks(model = model,
                              xcall = match.call(),
                              deriv = deriv,
                              resp  = resp, 
                              envir = envir)
  
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  
  get_xcall <- function(xcall, scall) {
    scall <- scall[[length(scall)]]
    if(any(grepl("growthparameters_predictions", scall, fixed = T)) |
       any(grepl("growthparameters_predictions.bgmfit", scall, fixed = T))) {
      xcall <- "growthparameters_predictions"
    } else {
      xcall <- xcall
    } 
  }
  
  if(!is.null(model$xcall)) {
    if(model$xcall == "growthparameters_predictions") {
      xcall <- "growthparameters_predictions"
    }
  } else {
    scall <- sys.calls()
    xcall <- get_xcall(xcall, scall)
  }
  
  
  model$xcall <- xcall
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  
  
  # This arguments$model <- model required when using pipe %>% to use gparameter
 # arguments$model <- model
  
  # https://michaelbarrowman.co.uk/post/getting-a-variable-name-in-a-pipeline/
  if(deparse(arguments$model) == ".") {
    first_call <- sys.calls()[[1]] # get the first entry on the call stack
    lhs <- first_call[[2]] # get the second element of this entry
    mymodel <- lhs # rlang::as_name(lhs)
    arguments$model <- mymodel
  }
  
  
  
  conf <- conf_level
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  # set_names_  <- c('Estimate', 'Est.Error', probtitles)
  set_names_  <- c('Estimate', probtitles)
  
  
  
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- setincores <-  get.cores_[['max.cores']]
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if (future_session == 'multisession') {
      future::plan('multisession', workers = setincores)
    } else if (future_session == 'multicore') {
      future::plan('multicore', workers = setincores)
    }
  }
  
  
  if (is.null(newdata)) {
    newdata <- model$data
  }
  
  # if (is.null(newdata)) {
  #   newdata <- get.newdata(
  #     model,
  #     newdata = newdata,
  #     resp = resp,
  #     numeric_cov_at = NULL,
  #     aux_variables = NULL,
  #     levels_id = NULL,
  #     ipts = ipts,
  #     xrange = xrange
  #   )
  # }
  
  
  
  arguments$newdata <- newdata
  
  
  arguments[["..."]] <- NULL
  
  
  ##############################################
  # Initiate non formalArgs()
  ##############################################
  term <- NULL;
  contrast <- NULL;
  tmp_idx <- NULL;
  predicted_lo <- NULL;
  predicted_hi <- NULL;
  predicted <- NULL;
  conf.high <- NULL;
  conf.low <- NULL;
  estimate <- NULL;
  `:=` <- NULL;
  `.` <- NULL;
  
  
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  
  if (is.null(levels_id)) {
    IDvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      IDvar <- model$model_info[[hierarchical_]]
    }
  } else if (!is.null(levels_id)) {
    IDvar <- levels_id
  }
  
  xfun_ <- paste0('xfun', resp_rev_)
  yfun_ <- paste0('yfun', resp_rev_)
  xfun <- model$model_info[[xfun_]]
  yfun <- model$model_info[[yfun_]]
  
  
  cov_ <- paste0('cov', resp_rev_)
  cov  <- model$model_info[[cov_]]
  
  cov_sigma_ <- paste0('cov_sigma', resp_rev_)
  cov_sigma  <- model$model_info[[cov_sigma_]]
  
  
  uvarby <- model$model_info$univariate_by
  
  
  
  #####################
  
  
  
  predictions_arguments <- arguments
  exclude_args <- as.character(quote(
    c(
      xrange,
      acg_velocity,
      digits,
      levels_id,
      ipts,
      seed,
      future,
      future_session
    )
  ))[-1]
  
  for (exclude_argsi in exclude_args) {
    predictions_arguments[[exclude_argsi]] <- NULL
  }
  
  
  # if (!is.null(variables)) {
  #   if (!is.list(variables)) {
  #     stop('variables must be a named list')
  #   } else if (is.list(variables)) {
  #     set_variables <- variables
  #   }
  # }
  # 
  # 
  # if (is.null(variables)) {
  #   set_variables <- list(eps)
  #   names(set_variables) <- xvar
  # }
  
  
  
  if(is.null(by)) {
    if(is.null(cov)) {
      set_group <- FALSE
    } else if(!is.null(cov)) {
      set_group <- cov
      if (!set_group %in% cov) {
        stop('by must be one of the ', cov)
      }
    }
  } else if(!is.null(by)) {
    if (!isFALSE(by)) {
      set_group <- by
    } else if (isFALSE(by)) {
      set_group <- FALSE
    }
  }
  
  
  
  if (acg_velocity >= 1 | acg_velocity <= 0) {
    stop("The acg_velocity should be set between 0.01 and 0.99")
  }
  
  
  
  newdata_org <- newdata
  newdata     <- marginaleffects::datagrid(model = model)
  
  get_coefs_brms      <-
    utils::getFromNamespace("coef.brmsfit", "brms")
  
  assign(o[[1]], model$model_info[['exefuns']][[o[[2]]]], envir = envir)
  
  if(model$model_info$select_model == "logistic3") {
    if(is.null(re_formula)) {
      xvar_names_       <- c(xvar, IDvar)
      all_names_        <- names(newdata)
      all_names_as_such <- setdiff(all_names_, xvar_names_)
      coef_   <- get_coefs_brms(model)
      coef_   <- coef_[[IDvar]]
      coef_d1 <- coef_[ , , 'c_Intercept']
      coef_d2 <- coef_[ , , 'f_Intercept']
      coef_d3 <- coef_[ , , 'i_Intercept']
      setid   <- levels(newdata_org[[IDvar]])
      if(at_age == 1) {
        setage  <- coef_d1
      } else if(at_age == 2) {
        setage  <- coef_d2
      } else if(at_age == 3) {
        setage  <- coef_d3
      }
      setage_expnad <- setage[,1]
      newdata <- newdata %>% 
        tidyr::expand(tidyr::nesting(!!xvar := setage_expnad, 
                              !!IDvar := setid, 
                              dplyr::select(
                                newdata, 
                                dplyr::all_of(!!all_names_as_such))))
      
      predictions_arguments$newdata    <- newdata
      if(tsv == 0 | tsv == 1) {
        suppressWarnings({
          if(!average) {
            out <- do.call(marginaleffects::predictions, 
                           predictions_arguments)
          } else if(average) {
            out <- do.call(marginaleffects::avg_predictions, 
                           predictions_arguments)
          }
        })
      } else {
        out <- setage
      }
      # out <- out %>% dplyr::mutate(!!IDvar := setid) 
    } else if(is.na(re_formula)) { # if(is.null(re_formula)) {
      xvar_names_       <- c(xvar)
      all_names_        <- names(newdata_org)
      all_names_as_such <- setdiff(all_names_, xvar_names_)
      fixed_  <- brms::fixef(model)
      fixed_2 <- fixed_[,1]
      if(at_age == 1) {
        setage  <- fixed_2[grepl("^c",names(fixed_2))]
      } else if(at_age == 2) {
        setage  <- fixed_2[grepl("^f",names(fixed_2))]
      } else if(at_age == 3) {
        setage  <- fixed_2[grepl("^i",names(fixed_2))]
      }
      # setage  <- fixed_2[grepl("^c|^f|^i",names(fixed_2))]
      setage_expnad <- setage
      newdata <- newdata %>% 
        tidyr::expand(tidyr::nesting(!!xvar := setage_expnad, 
                              dplyr::select(
                                newdata, 
                                dplyr::all_of(!!all_names_as_such))))
      predictions_arguments$newdata    <- newdata
      if(!is.null(tsv)) {
        suppressWarnings({
          if(!average) {
            out <- do.call(marginaleffects::predictions, 
                           predictions_arguments)
          } else if(average) {
            out <- do.call(marginaleffects::avg_predictions, 
                           predictions_arguments)
          }
        })
      } else if(is.null(tsv)) {
        setagenames <- names(setage)
        out <- brms::posterior_summary(model, 
                                       variable = paste0('b_', setagenames))
        row.names(out) <- setagenames
      }
      
    } # if(is.na(re_formula)) {
  } # if(model$model_info$select_model == "logistic3") {
  
  
  out_sf <- out %>% data.frame() # %>% dplyr::select(!dplyr::all_of(!!all_names_))
  
  out_sf <- out_sf %>% 
    dplyr::rename(!!as.symbol('Parameter') := parameter) %>% 
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ round(., digits = digits))) %>% 
    data.frame()
  
  
  if (reformat) {
    out_sf <- out_sf %>% 
      dplyr::rename(!!as.symbol(set_names_[1]) := estimate) %>% 
      dplyr::rename(!!as.symbol(set_names_[2]) := conf.low) %>% 
      dplyr::rename(!!as.symbol(set_names_[3]) := conf.high) 
    data.frame()
    
    remove_cols_ <- c('term', 'contrast', 'tmp_idx', 'predicted_lo', 
                      'predicted_hi', 'predicted')
    
    out_sf <- out_sf[,!names(out_sf) %in% remove_cols_]
  }
  
  
  return(out_sf)
}




#' @rdname growthparameters_predictions.bgmfit
#' @export
growthparameters_predictions <- function(model, ...) {
  UseMethod("growthparameters_predictions")
}



# growthparameters_predictions(berkeley_fit, tsv = 0, re_formula= NA, variables = 'sex',by = 'sex')



# get_contrasts <- marginaleffects:::get_contrasts
# for (i in 1:1) {
#  # print(i)
# c <- growthparameters_predictions(female_1444, draw_ids = i,
# by = NULL, parameter = c( 'vpv'))
#   #print(c)
# }


