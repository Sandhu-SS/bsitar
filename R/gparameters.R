


#' Growth parameter estimation for \code{bsitar} model
#'
#' @description The \code{gparameters} computes the growth parameter estimates
#'   and the uncertainty (standard errors, SE and credible intervals, CI) for
#'   population average and individual-specific parameters such as the peak
#'   growth velocity and the age at peak growth velocity.
#'
#' @details The \code{gparameters} functions first call the appropriate function
#'   (fitted or predict) to estimate velocity curve for each draw from the
#'   posterior distribution and then estimates various growth parameters such as
#'   age at peak growth velocity (APGV), peak growth velocity (PGV), age at
#'   takeoff growth velocity (ATGV), takeoff growth velocity (TGV), age at
#'   cessation of growth velocity (ACGV), and the cessation growth velocity
#'   (CGV). The growth parameters APGV and PGV are estimated by calling the
#'   [sitar::getPeak] function whereas the ATGV and TGV are estimated by using
#'   the [sitar::getTakeoff] function. The [sitar::getTrough] function is used
#'   to estimates ACGV and CGV parameters. These parameters obtained for each
#'   posterior draw are then summarized appropriately to get the estimate and
#'   the uncertainty (SE and CI) around these estimates. Please note that it is
#'   not always possible to estimate cessation and takeoff growth parameters
#'   when there are no distinct pre-peak or post-peak troughs.
#'
#'
#' @param model An object of class \code{bsitar}.
#' @param resp Optional names of response variables. If specified, predictions
#'   are performed only for the specified response variables.
#' @param ndraws Positive integer indicating the number of posterior draws to be
#'   used in estimation. If \code{NULL} (default), all draws are used.
#' @param newdata An optional data.frame for which to evaluate predictions. If
#'   \code{NULL} (default), the original data of the model is used.
#' @param summary A logical (default \code{TRUE}) to indicate whether only the
#'   Estimate should be returned or Estimate along with SE and CI should be
#'   computed. Setting this option to \code{FALSE} will reduce the computation
#'   time but no SE or CI estimates will be available.
#' @param robust If \code{FALSE} (the default) the mean is used as the measure
#'   of central tendency and the standard deviation as the measure of
#'   variability. If \code{TRUE}, the median and the median absolute deviation
#'   (MAD) are applied instead. Ignored if summary is \code{FALSE}
#' @param re_formula Option to indicate whether or not to include the
#'   individual/group-level effects in the estimation. When \code{NA} (default),
#'   the individual-level effects are excluded and therefore population average
#'   growth parameters are computed. When \code{NULL}, individual-level effects
#'   are included in the computation and hence the growth parameters estimates
#'   returned are individual-specific. In both sitations, (i.e,, \code{NA} or
#'   \code{NULL}), continuous and factor covariate(s) are appropriately
#'   included. When the continuous covariates by default are set to their means
#'   (see \code{numeric_cov_at} for details), the population average as well as
#'   individual-specific growth parameter estimates are returned for each level
#'   of the factor cocariate.
#' @param peak Optional logical specifying whether or not to calculate the age
#'   at peak velocity from the velocity curve. If \code{TRUE} (default), age at
#'   peak velocity (APGV) and the peak velocity (PGV) are calculated along with
#'   the uncertainty (SE and CI) around the APGV and PGV parameters. See
#'   @details for further information.
#' @param takeoff Optional logical (default \code{FALSE}) specifying whether or
#'   not to calculate the age at takeoff velocity from the velocity curve. If
#'   \code{TRUE}, age at takeoff velocity (ATGV) and the takeoff growth velocity
#'   (TGV) are  calculated along with the the uncertainty (SE and CI) around the
#'   ATGV and TGV parameters. See @details for further information.
#' @param trough Optional logical (default \code{FALSE}) specifying whether or
#'   not to calculate the age at takeoff velocity from the velocity curve. If
#'   \code{TRUE}, age at cessation of growth velocity (ACGV) and the cessation
#'   growth velocity (CGV) are  calculated along with the the uncertainty (SE
#'   and CI) around the ACGV and CGV parameters. See @details for further
#'   information.
#' @param estimation_method A character string to specify the estimation method
#'   used to calculate velocity from the posterior draws. The \code{'fitted'}
#'   method internally calls the [bsitar::fitted_.bsitar] function whereas the
#'   option \code{'predict'} calls the [bsitar::predict_.bsitar] function. See
#'   [brms::fitted.brmsfit] and [brms::predict.brmsfit] for derails on fitted
#'   versus predict estimations as they are the backend functions called by the
#'   [bsitar::fitted_.bsitar] and [bsitar::predict_.bsitar] functions.
#' @param numeric_cov_at An optional argument to specify the value of continuous
#'   covariate(s) before calling the [bsitar::fitted_.bsitar] and
#'   [bsitar::predict_.bsitar] functions. The default \code{NULL} option set the
#'   continuous covariate(s) at mean. Alternatively, a named list can be
#'   supplied to manualy set the values. For example, \code{numeric_cov_at =
#'   list(xx = 2)} will set the continuous covariate varibale 'xx' at 2. The
#'   argument \code{numeric_cov_at} is ignored when no continuous covariate is
#'   included in the model.
#' @param levels_id An optional argument to specify the ids for hierarchical
#'   model (default \code{NULL}. It is used only when a with 3 or higher level
#'   structures are fit. For two level model, the the id for second level is
#'   automatically inferred from the fitted model. Even for 3 or higher level
#'   model, ids are inferred from the fitted model but under the assumption that
#'   hierarchy is specified from lower to upper levels i.e, id, study assuming
#'   that id is nested within the studies. However, it is not gaunated that
#'   these ids are sorted correctly. Therefore, it is better to set it manually.
#' @param conf A numeric value (default \code{0.95}) to compute CI. Internally,
#'   this is translated into a paired probability values as \code{c((1 - conf) /
#'   2, 1 - (1 - conf) / 2)}. For \code{conf = 0.95}, this will compute 95% CI
#'   with CI varibales named as Q.2.5 and Q.97.5.
#' @param ipts A numeric value to interpolate the predictor value to get a
#'   smooth velocity curve. The \code{NULL} (default) will return original
#'   values whereas a positive real value (e.g., \code{ipts = 10}) will
#'   interpolate the predictor. It is important to note that these
#'   interpolations do not alter the range of predictor when calculating
#'   population average and the individual specific velocity curves.
#' @param seed An integer (default \code{123}) that is passed to the estimation
#'   method.
#' @param future A logical (default \code{FALSE}) to specify whether or not to
#'   perform parallel computations. If set to \code{TRUE}, the
#'   [future.apply::future_sapply] function is used for summarizing the draws.
#' @param future_session A character string to set the session type when
#'   \code{future = TRUE}. The \code{'multisession'} (default) options sets the
#'   multisession whereas the \code{'multicore'} set the multicore session. Note
#'   that multicore session are not supported on Windows systems. For more
#'   details, see [future.apply] and [future].
#' @param cores Number of cores to be used when running the parallel
#'   computations by setting the option \code{future = TRUE}. On non-Windows,
#'   systems this argument can be set globally via the mc.cores option. For the
#'   default \code{NULL} option, the number of cores are set automatically by
#'   calling the [future::availableCores()]. The number of cores used are the
#'   maximum number of cores avaialble minus one, i.e.,
#'   \code{future::availableCores() - 1}.
#' @param ... Further arguments passed to \code{brms::fitted} and
#'   \code{brms::predict} functions that control several aspects of data
#'   validation and prediction. See [brms::fitted.brmsfit] and
#'   [brms::predict.brmsfit] for details.
#'
#' @return a data frame with five columns when \code{summary = TRUE}, and two
#'   columns when \code{summary = False} (assuming \code{re_formual = NULL}).
#'   The first two columns common to both these approaches \code{summary =
#'   TRUE/False} are 'Parameter' and 'Estimate' which indicates the name of the
#'   growth parameter (e.g., APGV, PGV etc), and its value. When \code{summary =
#'   TRUE}, three additional columns are added, the 'Est.Error' and a pair of
#'   columns showing the CI. The CI columns are named as Q with appropriate
#'   suffix indicating the percentiles used to construct these intervals (such
#'   as Q.2.5 and Q.97.5 where 2.5 and 97.5 are the 0.025 and 0.975 percentiles
#'   used to compute by the 95% CI by calling the quantile function. When
#'   \code{re_formual = NULL}, an additional column is added that shows the
#'   individuals/groups corresponding to the estimates.
#'
#' @importFrom rlang .data
#'
#' @export gparameters.bsitar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Population average APGV and PGV
#' gparameters(model, re_formula = NA)
#'
#' #' # Population average APGV, PGV, ATGV, TGV
#' gparameters(model, re_formula = NA, peak = TRUE, takeoff = TRUE)
#'
#' # Individual-specific APGV, PGV, ATGV, TGV
#' gparameters(model, re_formula = NULL, peak = TRUE, takeoff = TRUE)
#'
#' }
#' 
gparameters.bsitar <- function(model,
                               resp = NULL,
                               ndraws = NULL,
                               newdata = NULL,
                               summary = TRUE,
                               robust = FALSE,
                               re_formula = NA,
                               peak = TRUE,
                               takeoff = FALSE,
                               trough = FALSE,
                               estimation_method = 'fitted',
                               numeric_cov_at = NULL,
                               levels_id = NULL,
                               conf = 0.95,
                               ipts = NULL,
                               seed = 123,
                               future = FALSE,
                               future_session = 'multisession',
                               cores = NULL,
                               ...) {
  if (is.null(ndraws))
    ndraws  <- ndraws(model)
  else
    ndraws <- ndraws
  
  
  oo <- post_processing_checks(model = model,
                               xcall = match.call(),
                               resp = resp)
  
  
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  
  if(xcall == 'plot_bsitar') arguments$plot <- TRUE else arguments$plot <- FALSE
  
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', 'Est.Error', probtitles)
  
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
  
  
  call_posterior_summary <- function(dat) {
    if (!robust) {
      . <- posterior::summarise_draws(dat,
                                      ~ mean(.x, na.rm = T),
                                      ~ sd(.x, na.rm = T),
                                      ~ quantile(.x,
                                                 probs = probs,
                                                 na.rm = T),
                                      .cores = .cores_ps)[-1] %>%
        data.frame() %>% stats::setNames(set_names_)
    } else if (robust) {
      . <- posterior::summarise_draws(
        dat,
        ~ median(.x, na.rm = T),
        ~ mad(.x, na.rm = T),
        ~ quantile(.x, probs = probs, na.rm = T),
        .cores = .cores_ps
      )[-1] %>%
        data.frame() %>% stats::setNames(set_names_)
    }
    as.data.frame(.)
  }
  
  #
  summarise_gp <-
    function(.x,
             probs,
             future,
             cores,
             peak,
             takeoff,
             trough,
             summary,
             robust) {
      Xnames <- names(.x)[grepl("^P._D.", names(.x))]
      .x <- .x %>% data.frame()
      if (peak) {
        if (future)
          out_1 <-
            future.apply::future_sapply(Xnames, function(x)
              sitar::getPeak(.x[[xvar]], as.numeric(.x[[x]])))
        if (!future)
          out_1 <-
            sapply(Xnames, function(x)
              sitar::getPeak(.x[[xvar]], as.numeric(.x[[x]])))
        out_1 <- t(out_1)
        colnames(out_1) <- c("APGV", "PGV")
      } else if (!peak) {
        out_1 <- NULL
      }
      
      if (takeoff) {
        if (future)
          out_2 <-
            future.apply::future_sapply(Xnames, function(x)
              sitar::getTakeoff(.x[[xvar]], as.numeric(.x[[x]])))
        if (!future)
          out_2 <-
            sapply(Xnames, function(x)
              sitar::getTakeoff(.x[[xvar]], as.numeric(.x[[x]])))
        out_2 <- t(out_2)
        colnames(out_2) <- c("ATGV", "TGV")
      } else if (takeoff) {
        out_2 <- NULL
      }
      
      if (trough) {
        if (future)
          out_3 <-
            future.apply::future_sapply(Xnames, function(x)
              sitar::getTrough(.x[[xvar]], as.numeric(.x[[x]])))
        if (!future)
          out_3 <-
            sapply(Xnames, function(x)
              sitar::getTrough(.x[[xvar]], as.numeric(.x[[x]])))
        out_3 <- t(out_3)
        colnames(out_3) <- c("ACGV", "CGV")
      } else if (trough) {
        out_3 <- NULL
      }
      
      xframe <- out_1
      if (exists('out_2'))
        xframe <- cbind(xframe, out_2)
      if (exists('out_3'))
        xframe <- cbind(xframe, out_3)
      xframe <- xframe %>% as.matrix()
      pnames <- colnames(xframe)[!grepl("^P._D.", colnames(xframe))]
      
      if (!robust) {
        o_ <- posterior::summarise_draws(
          xframe,
          ~ mean(.x, na.rm = T),
          ~ sd(.x, na.rm = T),
          ~ quantile(.x, probs = probs, na.rm = T),
          .cores = .cores_ps
        )[-1] %>%
          data.frame() %>% stats::setNames(set_names_) %>%
          dplyr::mutate(Parameter = pnames) %>%
          dplyr::relocate(Parameter)
      } else if (robust) {
        o_ <- posterior::summarise_draws(
          xframe,
          ~ median(.x, na.rm = T),
          ~ mad(.x, na.rm = T),
          ~ quantile(.x, probs = probs, na.rm = T),
          .cores = .cores_ps
        )[-1] %>%
          data.frame() %>% stats::setNames(set_names_) %>%
          dplyr::mutate(Parameter = pnames) %>%
          dplyr::relocate(Parameter)
      }
      
      if (summary) {
        . <- o_[, 1:2]
      } else if (!summary) {
        . <-  o_
      }
      .
    }
  
  
  multiNewVar <- function(df, df2, varname){
    df %>% dplyr::mutate(., !!varname := df2[[varname]])
  }
  
  get_gparameters <-
    function(out_v_,
             newdata,
             groupby_str,
             summary,
             ...) {
      
      # if (!summary) {
      #   out__ <- out_v_ %>%
      #     data.frame() %>%
      #     stats::setNames(paste0('P._D.', names(.))) %>%
      #     dplyr::mutate(!!IDvar := newdata[[IDvar]]) %>%
      #     dplyr::mutate(!!xvar := newdata[[xvar]])
      # } else if (summary) {
      #   out__ <- out_v %>% data.frame() %>%
      #     dplyr::select(1) %>%  data.frame() %>%
      #     stats::setNames(paste0('P._D.', names(.))) %>%
      #     dplyr::mutate(!!IDvar := newdata[[IDvar]]) %>%
      #     dplyr::mutate(!!xvar := newdata[[xvar]])
      # }
      
      if (!summary) {
        out__ <- out_v_ %>%
          data.frame() %>%
          stats::setNames(paste0('P._D.', names(.))) %>%
          # dplyr::mutate(!!IDvar := newdata[[IDvar]]) %>%
          dplyr::mutate(!!xvar := newdata[[xvar]])
        for(i in IDvar) {
          out__ <- out__ %>% multiNewVar(df=., df2 = newdata, varname=i)
        } 
      } else if (summary) {
        out__ <- out_v %>% data.frame() %>%
          dplyr::select(1) %>%  data.frame() %>%
          stats::setNames(paste0('P._D.', names(.))) %>%
          # dplyr::mutate(!!IDvar := newdata[[IDvar]]) %>%
          dplyr::mutate(!!xvar := newdata[[xvar]])
        for(i in IDvar) {
          out__ <- out__ %>% multiNewVar(df=., df2 = newdata, varname=i)
        } 
      }
      
      
      if (!is.null(groupby_str)) {
        out__ <-
          cbind(out__, newdata %>% dplyr::select(all_of(groupby_str))) %>%
          data.frame()
        out__ <-
          out__ %>% dplyr::group_by(across(all_of(groupby_str)))
      } else if (is.null(groupby_str)) {
        out__ <- cbind(out__, newdata)
      }
      
      if (is.null(re_formula)) {
        if (!summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                summary = summary,
                robust = robust
              )
            ) %>% dplyr::ungroup()
        } else if (summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                summary = summary,
                robust = robust
              )
            ) %>% dplyr::ungroup()
        }
      } else if (!is.null(re_formula)) {
        if (!summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                summary = summary,
                robust = robust
              )
            ) %>% dplyr::ungroup()
        } else if (summary) {
          parameters <- out__ %>%
            dplyr::group_modify(
              ~ summarise_gp(
                .x,
                probs = probs,
                future = future,
                cores = setincores,
                peak = peak,
                takeoff = takeoff,
                trough = trough,
                summary = summary,
                robust = robust
              )
            ) %>% dplyr::ungroup()
        }
      }
      parameters
    }
  ###################################
  
  
  if (arguments$plot) {
    out_summary <- list()
    if(!is.null(arguments$...)) {
      arguments <- c(arguments, list(arguments$...))
    }
    
    for (argumentsi in names(arguments)) {
      if (length(arguments[[argumentsi]]) != 0) {
        if (argumentsi != "...") {
          arguments[[argumentsi]] <- eval(arguments[[argumentsi]])
        } 
      }
    }
    
    arguments[which(names(arguments) %in% "")] <- NULL
    
    list2env(arguments, envir = parent.env(new.env()))
    
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    
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
    
    
    if(is.null(newdata)) {
      newdata <- get.newdata(model, newdata = newdata, 
                             resp = resp, 
                             numeric_cov_at = numeric_cov_at,
                             levels_id = levels_id,
                             ipts = ipts)
    }
    
    
    list_c <- attr(newdata, 'list_c')
    
    for (list_ci in names(list_c)) {
      assign(list_ci, list_c[[list_ci]])
    }
    check__ <- c('xvar', 'yvar', 'IDvar', 'cov_vars', 'cov_factor_vars', 
                 'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 'uvarby', 'subindicatorsi')
    
    for (check___ in check__) {
      if(!exists(check___)) assign(check___, NULL)
    }
    
    
    # newdata <- i_data(model, newdata, resp = resp,
    #                            cov_factor_vars = cov_factor_vars, 
    #                            cov_numeric_vars = cov_numeric_vars,
    #                            ipts = ipts)
    
    newdata___ <- newdata
    
    if (grepl("d", opt, ignore.case = T)) {
      index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
      dist.. <- substr(opt, index_opt, index_opt)
    } else if (!grepl("d", opt, ignore.case = T)) {
      dist.. <- ""
    }
    
    if (grepl("v", opt, ignore.case = T)) {
      index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
      velc.. <- substr(opt, index_opt, index_opt)
    } else if (!grepl("v", opt, ignore.case = T)) {
      velc.. <- ""
    }
    
    if ((apv) & velc.. == "") {
      stop("You have set apv = TRUE but your opt argument",
           "\n ",
           "contains no 'v' or 'V' option")
    }
    
    if (dist.. != "") {
      if (grepl("^[[:upper:]]+$", dist..)) {
        groupby_str_d <- groupby_fistr
      } else  if (!grepl("^[[:upper:]]+$", dist..)) {
        groupby_str_d <- groupby_fstr
      }
      if (identical(groupby_str_d, character(0)))
        groupby_str_d <- NULL
      groupby_str_d <- groupby_str_d
    } else {
      groupby_str_d <- NULL
    }
    
    
    if (velc.. != "") {
      if (grepl("^[[:upper:]]+$", velc..)) {
        groupby_str_v <- groupby_fistr
      } else  if (!grepl("^[[:upper:]]+$", velc..)) {
        groupby_str_v <- groupby_fstr
      }
      if (identical(groupby_str_v, character(0)))
        groupby_str_v <- NULL
    } else {
      groupby_str_v <- NULL
    }
    
    if (dist.. != "") {
      newdata <- newdata___
      if (grepl("^[[:upper:]]+$", dist..)) {
        arguments$re_formula <- NULL
      } else if (!grepl("^[[:upper:]]+$", dist..)) {
        arguments$re_formula <- NA
        if (!is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(groupby_fstr, xvar)
        } else if (is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(xvar)
        }
        newdata <- newdata %>%
          dplyr::group_by(across(all_of(groupby_fstr_xvars))) %>%
          dplyr::slice(1) %>% dplyr::ungroup()
      }
      arguments$newdata <- newdata
      arguments$deriv <- 0
      
      arguments$envir <-arguments$envir_
      
      if (estimation_method == 'fitted') {
        out_d_ <- do.call(fitted_.bsitar, arguments)
      } else if (estimation_method == 'predict') {
        out_d_ <- do.call(predict_.bsitar, arguments)
      }
      # No need to summarise by calling call_posterior_summary
      # if (!summary) {
      #   out_d <- call_posterior_summary(out_d_)
      # } else if (summary) {
      #   out_d <- out_d_
      # }
      
      out_d <- out_d_
      
      if(!is.na(model$model_info$univariate_by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
     
      out_summary[['distance']] <-  cbind(newdata,
                                          out_d %>% data.frame() %>%
                                            dplyr::mutate(curve = 'distance')) %>%
        data.frame()
      
    } else if (dist.. == "") {
      out_summary[['distance']] <- NULL
    }
    
    
    if (velc.. != "") {
      newdata <- newdata___
      if (grepl("^[[:upper:]]+$", velc..)) {
        arguments$re_formula <- NULL
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        arguments$re_formula <- NA
        if (!is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(groupby_fstr, xvar)
        } else if (is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(xvar)
          
        }
        newdata <- newdata %>%
          dplyr::group_by(across(all_of(groupby_fstr_xvars))) %>%
          dplyr::slice(1) %>% dplyr::ungroup()
      }
      arguments$newdata <- newdata
      arguments$deriv <- 1
      
      arguments$envir <-arguments$envir_
      
      if (estimation_method == 'fitted') {
        out_v_ <- do.call(fitted_.bsitar, arguments)
      } else if (estimation_method == 'predict') {
        out_v_ <- do.call(predict_.bsitar, arguments)
      }
      # No need to summa call_posterior_summary
      # if (!summary) {
      #   out_v <- call_posterior_summary(out_v_)
      # } else if (summary) {
      #   out_v <- out_v_
      # }
      
      out_v <- out_v_
      
      if(!is.na(model$model_info$univariate_by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
      
      out_summary[['velocity']] <-
        cbind(newdata,
              out_v %>% data.frame() %>%
                dplyr::mutate(curve = 'velocity')) %>%
        data.frame()
    } else if (velc.. == "") {
      out_summary[['velocity']] <- NULL
    }
    
    
    if (apv) {
      if(!is.na(model$model_info$univariate_by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
      # out_summary[['parameters']] <-
      #   get_gparameters(t(out_v_), newdata, groupby_str_v, summary)
      out_summary[['parameters']] <-
        get_gparameters(out_v_, newdata, groupby_str_v, summary)
    }
    out_summary[['groupby_str_d']] <- groupby_str_d
    out_summary[['groupby_str_v']] <- groupby_str_v
    out_summary[['probtitles']] <- probtitles
    return(out_summary)
  }
  
  
  if (!arguments$plot) {
    
    newdata <- get.newdata(model, newdata = newdata, 
                           resp = resp, 
                           numeric_cov_at = numeric_cov_at,
                           levels_id = levels_id,
                           ipts = ipts)
    
    list_c <- attr(newdata, 'list_c')
    for (list_ci in names(list_c)) {
      assign(list_ci, list_c[[list_ci]])
    }
    check__ <- c('xvar', 'yvar', 'IDvar', 'cov_vars', 'cov_factor_vars', 
                 'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 'uvarby', 'subindicatorsi')
    
    for (check___ in check__) {
      if(!exists(check___)) assign(check___, NULL)
    }
    
    # newdata <- i_data(model, newdata, resp = resp, 
    #                            cov_factor_vars = cov_factor_vars, 
    #                            cov_numeric_vars = cov_numeric_vars,
    #                            ipts = ipts)
    
    if (is.null(re_formula)) {
      groupby_str <- groupby_fistr
    } else  if (!is.null(re_formula)) {
      groupby_str <- groupby_fstr
    }
    if (identical(groupby_str, character(0))) groupby_str <- NULL
    
    if (!is.null(re_formula)) {
      if (!is.null(groupby_fstr)) {
        groupby_fstr_xvars <- c(groupby_fstr, xvar)
      } else if (is.null(groupby_fstr)) {
        groupby_fstr_xvars <- c(xvar)
      }
      newdata <- newdata %>%
        dplyr::group_by(across(all_of(groupby_fstr_xvars))) %>%
        dplyr::slice(1) %>% dplyr::ungroup()
    }
    
    arguments$newdata <- newdata
    arguments$deriv <- 1
    
    arguments$envir <- parent.frame()
    
    if (estimation_method == 'fitted') {
      out_v_ <- do.call(fitted_.bsitar, arguments)
    } else if (estimation_method == 'predict') {
      out_v_ <- do.call(predict_.bsitar, arguments)
    }
    if (!summary) {
      out_v <- call_posterior_summary(out_v_)
    } else if (summary) {
      out_v <- out_v_
    }
    
    if(!is.na(model$model_info$univariate_by)) {
      newdata <- newdata %>%
        dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
    }
    
    parameters <-
      get_gparameters(out_v, newdata, groupby_str, summary)
    return(parameters)
  }
}

#' @rdname gparameters.bsitar
#' @export
gparameters <- function(model, ...) {
  UseMethod("gparameters")
}
