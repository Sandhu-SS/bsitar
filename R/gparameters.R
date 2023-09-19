


#' Growth parameter estimation for \code{bsitar} model
#'
#' @description The \code{gparameters} computes growth parameters and the 
#'   uncertainty (standard error, SE and the credible interval, CI) for
#'   population average and individual-specific parameters (see @details).
#'
#' @details The \code{gparameters} calls the appropriate function (fitted or
#'   predict) to estimate the first derivative (velocity curve) for each
#'   posterior draw (posterior distribution) and then computes growth parameters
#'   such as age at peak growth velocity (APGV), peak growth velocity (PGV), age
#'   at takeoff growth velocity (ATGV), takeoff growth velocity (TGV), age at
#'   cessation of growth velocity (ACGV), and the cessation growth velocity
#'   (CGV). The growth parameters APGV and PGV are estimated by calling the
#'   [sitar::getPeak()] function whereas the ATGV and TGV are estimated by using
#'   the [sitar::getTakeoff()] function. The [sitar::getTrough()] function is
#'   used to estimates ACGV and CGV parameters. The parameters obtained from
#'   each posterior draw are then summarized appropriately to get the estimate
#'   and the uncertainty (SE and CI) around these estimates. Please note that it
#'   is not always possible to estimate cessation and takeoff growth parameters
#'   when there are no distinct pre-peak or post-peak troughs.
#'
#'
#' @param model An object of class \code{bsitar}.
#' 
#' @param resp An optional character string to specify response variable when
#'   estimating growth parameter for the univariate-by-subgroup and multivariate
#'   models (see [bsitar::bsitar()] for details on univariate-by-subgroup and
#'   multivariate models). For univariate model, \code{resp = NULL} (default).
#'   
#' @param ndraws Positive integer indicating the number of posterior draws to be
#'   used in estimation. If \code{NULL} (default), all draws are used.
#'   
#' @param draw_ids An integer indicating the specif posterior draw(s) 
#' to be used. If \code{NULL} (default), all draws are used.
#' 
#' @param newdata An optional data frame to be used in predictions. If
#'   \code{NULL} (default), the original data from the fitted model is used.
#'   
#' @param summary A logical (default \code{TRUE}) indicating whether only the
#'   Estimate should be returned or Estimate along with SE and CI should be
#'   computed. Setting this option to \code{FALSE} will reduce the computation
#'   time but no SE or CI estimates will be available.
#'   
#' @param robust If \code{FALSE} (the default) the mean is used as the measure
#'   of central tendency and the standard deviation as the measure of
#'   variability. If \code{TRUE}, the median and the median absolute deviation
#'   (MAD) are applied instead. Ignored if summary is \code{FALSE}. 
#'   
#' @param re_formula Option to indicate whether or not to include the
#'   individual/group-level effects in the estimation. When \code{NA} (default),
#'   the individual-level effects are excluded and therefore population average
#'   growth parameters are computed. When \code{NULL}, individual-level effects
#'   are included in the computation and hence the growth parameters estimates
#'   returned are individual-specific. In both situations, (i.e,, \code{NA} or
#'   \code{NULL}), continuous and factor covariate(s) are appropriately included
#'   in the estimation. The continuous covariates by default are set to their
#'   means (see \code{numeric_cov_at} for details) whereas factor covariates are
#'   left unaltered thereby allowing estimation of (factor) covariate specific
#'   population average and individual-specific growth parameter at mean value
#'   of continous covaristes(s).
#'   
#' @param peak Optional (logical, default \code{TRUE}) to indicate whether or
#'   not to calculate the age at peak velocity (APGV) and the peak velocity
#'   (PGV). See @details for further information.
#'   
#' @param takeoff  Optional (logical, default \code{FALSE}) to indicate whether
#'   or not to calculate the age at takeoff velocity (ATGV) and the takeoff
#'   growth velocity (TGV). See @details for further information.
#'
#' @param trough Optional (logical, default \code{FALSE}) to indicate whether or
#'   not to calculate the age at cessation of growth velocity (ACGV) and the
#'   cessation of growth velocity (CGV). See @details for further information.
#'   See @details for further information.
#' 
#' @param acgv Optional (logical, default \code{FALSE}) to indicate whether or
#'   not to calculate the age at cessation of growth velocity from the velocity
#'   curve. If \code{TRUE}, age at cessation of growth velocity (ACGV) and the
#'   cessation growth velocity (CGV) are  calculated based on the percentage of
#'   the peak growth velocity as defined by the \code{acgv_velocity} argument
#'   (see below) which is typically set at 10 percent of the peak growth
#'   velocity. The ACGV and CGV are calculated along with the the uncertainty
#'   (SE and CI) around the ACGV and CGV parameters. 
#'   
#' @param acgv_velocity Specify the percentage of the peak growth velocity to be 
#'  used when estimating \code{acgv}. The default value is \code{0.10} i.e., 
#'  10 percent of the peak growth velocity.
#'   
#' @param estimation_method A character string to specify the estimation method
#'   when calculating the velocity from the posterior draws. The \code{fitted}
#'   method internally calls the [bsitar::fitted_.bsitar()] function whereas the
#'   option \code{predict} calls the [bsitar::predict_.bsitar()] function. See
#'   [brms::fitted.brmsfit()] and [brms::predict.brmsfit()] for derails on
#'   fitted versus predict estimations.
#'   
#' @param numeric_cov_at An optional argument to specify the value of continuous
#'   covariate(s). The default \code{NULL} option set the continuous
#'   covariate(s) at their mean. Alternatively, a named list can be
#'   supplied to manually set these values. For example, \code{numeric_cov_at =
#'   list(xx = 2)} will set the continuous covariate varibale 'xx' at 2. The
#'   argument \code{numeric_cov_at} is ignored when no continuous covariate is
#'   included in the model.
#'   
#' @param levels_id An optional argument to specify the ids for hierarchical
#'   model (default \code{NULL}. It is used only when model is fitted to the
#'   data with 3 or more levels of hierarchy. For a two level model, the id for
#'   second level is automatically inferred from the fitted model. Even for 3 or
#'   higher level model, ids are inferred from the fitted model but under the
#'   assumption that hierarchy is specified from lower to upper levels i.e, id,
#'   study assuming that id is nested within the studies. However, it is not
#'   gauranted that these ids are sorted correctly. Therefore, it is better to
#'   set them manually.
#'   
#' @param avg_reffects An optional argument (default \code{NULL} to calculate
#'   (marginal/average) curves and growth parameters (such as APGV and PGV). If
#'   specified, it must be a named list indicating the \code{over} and the fixed
#'   efects \code{feby} and random effects efects \code{reby} arguments e.g.,
#'   \code{avg_reffects = list(feby = 'study', reby = NULL, over = 'age'}. The
#'   \code{over} is typically age and is used to average over the random
#'   effects. The second argument is \code{by} that specifies the factor
#'   variable by which \code{over} is executed.
#'   
#'@param aux_variables An optional argument to specify the variables that can be
#'  passed to the \code{ipts} argument (see below). This is useful when fitting
#'  location scale models and the measurement error models.
#'   
#' @param ipts An integer to set the length of the predictor variable to get a
#'   smooth velocity curve. The \code{NULL} (default) will return original
#'   values whereas an integer (e.g., \code{ipts = 10}) will interpolate the
#'   predictor. It is important to note that these interpolations do not alter
#'   the range of predictor when calculating population average and the
#'   individual specific velocity curves.
#'   
#'@param conf A numeric value (default \code{0.95}) to compute CI. Internally,
#'   this is translated into a paired probability values as \code{c((1 - conf) /
#'   2, 1 - (1 - conf) / 2)}. For \code{conf = 0.95}, this will compute 95% CI
#'   with CI varibales named as Q.2.5 and Q.97.5.
#'   
#' @param xrange An integer to set the predictor range (i.e., age) when
#'   executing the interpolation via \code{ipts}. The default \code{NULL} sets
#'   the individual specific predictor range whereas code \code{xrange = 1} sets
#'   same range for all individuals within the higher order grouping variable
#'   (e.g., study). Code \code{xrange  = 2} sets the identical range dplyr::across the
#'   entire sample. Lastly, a paired numeric values can be supplied e.g.,
#'   \code{xrange = c(6, 20)} will set the range between 6 and 20.
#'   
#' @param seed An integer (default \code{123}) that is passed to the estimation
#'   method.
#'   
#' @param future A logical (default \code{FALSE}) to specify whether or not to
#'   perform parallel computations. If set to \code{TRUE}, the
#'   [future.apply::future_sapply()] function is used for summarizing the draws.
#'   
#' @param future_session A character string to set the session type when
#'   \code{future = TRUE}. The \code{multisession} (default) options sets the
#'   multisession whereas the \code{multicore} sets the multicore session. Note
#'   that multicore session is not supported on Windows systems. For more
#'   details, see [future.apply::future_sapply()].
#'   
#' @param cores Number of cores to be used when running the parallel
#'   computations by setting the option \code{future = TRUE}. On non-Windows,
#'   systems this argument can be set globally via the mc.cores option. For the
#'   default \code{NULL} option, the number of cores are set automatically by
#'   calling the [future::availableCores()]. The number of cores used are the
#'   maximum number of cores avaialble minus one, i.e.,
#'   \code{future::availableCores() - 1}.
#'  
#' @param parms_eval A logical (default \code{FALSE}) to specify whether or not 
#' to get growth parameters on the fly. 
#' 
#' @param parms_method A character to specify the method used to when 
#'  evaluating \code{parms_eval}. The default is \code{getPeak} which uses
#'  the [sitar::getPeak()] function from the \code{sitar} package. The 
#'  alternative option is \code{findpeaks} that uses the [pracma::findpeaks()] 
#'  function function from the \code{pracma} package.
#'  
#' @param envir Indicator to set the environment of function evaluation.  
#'  The default is \code{parent.frame}. 
#'   
#' @param ... Further arguments passed to \code{brms::fitted()} and
#'   \code{brms::predict()} functions that control several aspects of data
#'   validation and prediction. See [brms::fitted.brmsfit()] and
#'   [brms::predict.brmsfit()] for details.
#'
#' @return A data frame with five columns when \code{summary = TRUE}, and two
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
#' @inherit brms::prepare_predictions.brmsfit params 
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
                               draw_ids = NULL,
                               newdata = NULL,
                               summary = TRUE,
                               robust = FALSE,
                               re_formula = NA,
                               peak = TRUE,
                               takeoff = FALSE,
                               trough = FALSE,
                               acgv = FALSE,
                               acgv_velocity = 0.10,
                               estimation_method = 'fitted',
                               allow_new_levels = FALSE,
                               sample_new_levels = "uncertainty",
                               incl_autocor = TRUE,
                               numeric_cov_at = NULL,
                               levels_id = NULL,
                               avg_reffects = NULL,
                               aux_variables = NULL,
                               ipts = NULL,
                               conf = 0.95,
                               xrange = NULL,
                               seed = 123,
                               future = FALSE,
                               future_session = 'multisession',
                               cores = NULL,
                               parms_eval = FALSE,
                               parms_method = 'getPeak',
                               envir = NULL,
                               ...) {
  if (is.null(ndraws))
    ndraws  <- ndraws(model)
  else
    ndraws <- ndraws
  
  if(is.null(envir)) envir <- parent.frame()
  
  ###
  xvar <- NULL;
  acgv_asymptote <- NULL;
  apv <- NULL;
  Parameter <- NULL;
  IDvar <- NULL;
  groupby_fistr <- NULL;
  groupby_fstr <- NULL;
  subindicatorsi <- NULL;
  Estimate <- NULL;
  ':=' <- NULL;
  . <- NULL;
  ###
  
  oo <- post_processing_checks(model = model,
                               xcall = match.call(),
                               resp = resp,
                               envir = envir)
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  
  get_xcall <- function(xcall, scall) {
    scall <- scall[[length(scall)]]
    if(any(grepl("plot_bsitar", scall, fixed = T)) |
       any(grepl("plot_bsitar.bsitar", scall, fixed = T))) {
      xcall <- "plot_bsitar"
    } else if(any(grepl("gparameters", scall, fixed = T)) |
              any(grepl("gparameters.bsitar", scall, fixed = T))) {
      xcall <- "gparameters"
    } else {
      xcall <- xcall
    } 
  }
  
  if(!is.null(model$xcall)) {
    if(model$xcall == "plot_bsitar") {
      xcall <- "plot_bsitar"
    }
  } else {
    scall <- sys.calls()
    xcall <- get_xcall(xcall, scall)
  }
  
  
  model$xcall <- xcall
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  
  # This arguments$model <- model required when using pipe %>% to use gparameter
  arguments$model <- model
  
  if(xcall == 'plot_bsitar') arguments$plot <- TRUE else arguments$plot <- FALSE
  
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', 'Est.Error', probtitles)
  
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- setincores <-  get.cores_[['max.cores']] 
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if(is.null(cores)) stop("Please set the number of cores for 'future' by  
                            using the the 'cores' argument, e.g. cores = 4")
    if (arguments$future_session == 'multisession') {
      future::plan('multisession', workers = setincores)
    } else if (arguments$future_session == 'multicore') {
      future::plan('multicore', workers = setincores)
    }
  }
  
  
  
  #
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
             acgv,
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
        # out_1x <<- out_1
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
      } else if (!takeoff) {
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
      } else if (!trough) {
        out_3 <- NULL
      }
      
      #######
      
      if (acgv) {
        # acgv_asymptote not supported becuase this currect function used velocity
        if(!is.null(acgv_asymptote)) stop("Currently acgv_asymptote not 
                                           supported. Please use acgv_velocity")
        if(!is.null(acgv_asymptote) & !is.null(acgv_velocity) ) {
          stop("Specify either acgv_asymptote or acgv_velocity but not both")
        }
        if(is.null(acgv_asymptote) & is.null(acgv_velocity) ) {
          stop("Specify either acgv_asymptote or acgv_velocity")
        }
        
        set_get___fun <- function(x) {
          if (!peak) pkkk <- sitar::getPeak(.x[[xvar]], as.numeric(.x[[x]]))
          if ( peak | apv) pkkk <- out_1
          get__ <- as.numeric(.x[[x]])
          set_x_for_afo <- .x[[xvar]]
          tempbind <- cbind(get__, set_x_for_afo) %>% 
            data.frame() %>% 
            filter(set_x_for_afo > pkkk [1])
          get__ <- tempbind[,1]
          set_x_for_afo <-  tempbind[,2]
          if(length(get__) != 0) {
            get___pct <- max(get__) * acgv_velocity
            target.index <- which(abs(get__ - get___pct) == min(abs(get__ - get___pct)))
            get_asymptote_pct_x <- set_x_for_afo[target.index]
            if(length(get_asymptote_pct_x) > 1) get_asymptote_pct_x <- mean(get_asymptote_pct_x)
            out_3_temp <- c(get_asymptote_pct_x, get___pct)
          } else if(length(get__) == 0) {
            out_3_temp <- c(NA, NA)
          }
          # print(length(out_3_temp))
          # print(out_3_temp)
          out_3_temp
        }
        if (future)
          out_3 <-
            future.apply::future_sapply(Xnames, function(x)
              set_get___fun(x))
        if (!future)
          out_3 <-
            sapply(Xnames, function(x)
              set_get___fun(x))
        out_3 <- t(out_3)
        # out_3x <<- out_3
        # print("kkkkkkkkkkkkk")
       #  if(!summary) out_3 <- out_3 %>% do.call(rbind, .) 
        colnames(out_3) <- c("ACGV", "CGV")
      } else if (!acgv) {
        out_3 <- NULL
      }
      

      ######
      
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
          cbind(out__, newdata %>% dplyr::select(dplyr::all_of(groupby_str))) %>%
          data.frame()
        out__ <-
          out__ %>% dplyr::group_by(dplyr::across(dplyr::all_of(groupby_str)))
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
                acgv = acgv,
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
                acgv = acgv,
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
                acgv = acgv,
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
                acgv = acgv,
                summary = summary,
                robust = robust
              )
            ) %>% dplyr::ungroup()
        }
      }
      parameters
    }
  
  
  #
  get_avg_over <- function(raw_re, newdata, by, probs, robust) {
    raw_re_c <- c()
    getitEstimate <- getitarray <- NULL
    for (i in 1:(dim(raw_re)[1])) {
      getitEstimate <- raw_re[i,]
      raw_re_c[i] <- cbind(newdata, getitEstimate) %>% data.frame() %>% 
        dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
        dplyr::summarise(getitEstimate = mean(getitEstimate), .groups = 'drop') %>%
        dplyr::ungroup() %>%
        dplyr::select(getitEstimate) 
    }
    
    getitarray <- array(unlist(raw_re_c), 
                        dim=c(length(raw_re_c[[1]]), length(raw_re_c)  ))
    getitarray <- t(getitarray)
    # brms::posterior_summary(getitarray, probs = probs, robust = robust)
    getitarray
  }
  
  ###################################
  
  
  if (arguments$plot) {
    out_summary <- list()
    if(!is.null(arguments$...)) {
      arguments <- c(arguments, list(arguments$...))
    }
    arguments$model <- model
    
    
    
    for (argumentsi in names(arguments)) {
      if (length(arguments[[argumentsi]]) != 0) {
        if (argumentsi != "...") {
          arguments[[argumentsi]] <- eval(arguments[[argumentsi]])
        } 
      }
    }
    
    arguments[which(names(arguments) %in% "")] <- NULL
    
    list2env(arguments, envir = parent.env(new.env()))
    
    # arguments$model <- arguments$model_
    
    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    probtitles <- probs[order(probs)] * 100
    probtitles <- paste("Q", probtitles, sep = "")
    set_names_  <- c('Estimate', 'Est.Error', probtitles)
    
    get.cores_ <- get.cores(arguments$cores)
    arguments$cores <- setincores <-  get.cores_[['max.cores']] 
    .cores_ps <- get.cores_[['.cores_ps']]
    
    if (future) {
      if(is.null(cores)) stop("Please set the number of cores for 'future' by  
                            using the the 'cores' argument, e.g. cores = 4")
      if (arguments$future_session == 'multisession') {
        future::plan('multisession', workers = setincores)
      } else if (arguments$future_session == 'multicore') {
        future::plan('multicore', workers = setincores)
      }
    }
    
    
    newdata <- get.newdata(model, newdata = newdata, 
                           resp = resp, 
                           numeric_cov_at = numeric_cov_at,
                           aux_variables = aux_variables,
                           levels_id = levels_id,
                           ipts = ipts,
                           xrange = xrange)
    
    
    list_c <- attr(newdata, 'list_c')
    
    
    for (list_ci in names(list_c)) {
      assign(list_ci, list_c[[list_ci]])
    }
    check__ <- c('xvar', 'yvar', 'IDvar', 'cov_vars', 'cov_factor_vars', 
                 'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 'uvarby', 'subindicatorsi')
    
    for (check___ in check__) {
      if(!exists(check___)) assign(check___, NULL)
    }
    
    
    
    newdata___ <- newdata
    
    if(!is.null(avg_reffects)) {
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
        if (dist.. != "" & grepl("^[[:upper:]]+$", dist..) )
          stop("use option 'd' (and not 'D') with avg_reffects" )
      }
      if (grepl("v", opt, ignore.case = T) ) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (velc.. != "" & grepl("^[[:upper:]]+$", velc..) )
          stop("use option 'v' (and not 'V') with avg_reffects" )
      }
    }
    
    
    if(is.null(avg_reffects)) {
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
        # groupby_str_d <- c(avg_reffects_groupby_str_d, groupby_str_d)
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
        # groupby_str_v <- c(avg_reffects_groupby_str_v, groupby_str_v)
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
            dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars))) %>%
            dplyr::slice(1) %>% dplyr::ungroup()
        }
        arguments$newdata <- newdata
        arguments$deriv <- 0
        # don't let the ipts to pass again to the fitted_.bsitar
        arguments$ipts <- NULL 
        arguments$envir <- .GlobalEnv # arguments$envir_
        arguments$probs <- probs
        
        if (estimation_method == 'fitted') {
          out_d_ <- do.call(fitted_.bsitar, arguments)
        } else if (estimation_method == 'predict') {
          out_d_ <- do.call(predict_.bsitar, arguments)
        }
        probs
        if (!summary) {
          out_d <- call_posterior_summary((out_d_))
        } else if (summary) {
          out_d <- out_d_
        }
        
        # out_d <- out_d_
       
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
            dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars))) %>%
            dplyr::slice(1) %>% dplyr::ungroup()
        }
        arguments$newdata <- newdata
        arguments$deriv <- 1
        # don't let the ipts to pass again to the fitted_.bsitar
        arguments$ipts <- NULL 
        arguments$envir <- .GlobalEnv # arguments$envir_
        arguments$probs <- probs
        
        if (estimation_method == 'fitted') {
          out_v_ <- do.call(fitted_.bsitar, arguments)
        } else if (estimation_method == 'predict') {
          out_v_ <- do.call(predict_.bsitar, arguments)
        }
        out_v__apv_ <- out_v_
        if (!summary) {
          out_v <- call_posterior_summary((out_v_))
        } else if (summary) {
          out_v <- out_v_
        }
        
        # out_v <- out_v_
        
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
      
      
      if (apv | takeoff | trough | acgv) {
        if(!is.na(model$model_info$univariate_by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
        }
        # out_summary[['parameters']] <-
        #   get_gparameters(t(out_v_), newdata, groupby_str_v, summary)
        
        out_v__apv_ <- t(out_v__apv_)
        out_summary[['parameters']] <-
          get_gparameters(out_v__apv_, newdata, groupby_str_v, summary)
      }
      out_summary[['groupby_str_d']] <- groupby_str_d
      out_summary[['groupby_str_v']] <- groupby_str_v
      out_summary[['probtitles']] <- probtitles
    } # if(is.null(avg_reffects)) {
    
    
    
    
    if(!is.null(avg_reffects)) {
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
      
      # if (dist.. != "") {
      #   if (grepl("^[[:upper:]]+$", dist..)) {
      #     groupby_str_d <- groupby_fistr
      #   } else  if (!grepl("^[[:upper:]]+$", dist..)) {
      #     groupby_str_d <- groupby_fstr
      #   }
      #   if (identical(groupby_str_d, character(0)))
      #     groupby_str_d <- NULL
      #   groupby_str_d <- groupby_str_d
      # } else {
      #   groupby_str_d <- NULL
      # }
      
      if (dist.. != "") {
        groupby_str_d <- groupby_fstr
        if (identical(groupby_str_d, character(0)))
          groupby_str_d <- NULL
        groupby_str_d <- groupby_str_d
      } else {
        groupby_str_d <- NULL
      }
      
      # if (velc.. != "") {
      #   if (grepl("^[[:upper:]]+$", velc..)) {
      #     groupby_str_v <- groupby_fistr
      #   } else  if (!grepl("^[[:upper:]]+$", velc..)) {
      #     groupby_str_v <- groupby_fstr
      #   }
      #   if (identical(groupby_str_v, character(0)))
      #     groupby_str_v <- NULL
      # } else {
      #   groupby_str_v <- NULL
      # }
      
      if (velc.. != "") {
        groupby_str_v <- groupby_fstr
        if (identical(groupby_str_v, character(0)))
          groupby_str_v <- NULL
      } else {
        groupby_str_v <- NULL
      }
      
      summary_org <- arguments$summary
      arguments$summary <- FALSE
      arguments$re_formula <- NULL
      
      # avg_reffects_groupby_str_d <- avg_reffects[['feby']]
      # avg_reffects_groupby_str_v <- avg_reffects[['feby']]
      
      groupby_str_d <- c(groupby_str_d, avg_reffects[['feby']])
      groupby_str_v <- c(groupby_str_v, avg_reffects[['feby']])
      

      # Don't let below arguments$re_formula override the above 'd' and 'v'
      if (dist.. != "") {
        newdata <- newdata___
        if (grepl("^[[:upper:]]+$", dist..)) {
          # arguments$re_formula <- NULL
        } else if (!grepl("^[[:upper:]]+$", dist..)) {
          #  arguments$re_formula <- NA
          if (!is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(groupby_fstr, xvar)
          } else if (is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(xvar)
          }
          
          # newdata <- newdata %>%
          #   dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars))) %>%
          #   dplyr::slice(1) %>% dplyr::ungroup()
        }
        
        arguments$newdata <- newdata
        arguments$deriv <- 0
        # don't let the ipts to pass again to the fitted_.bsitar
        arguments$ipts <- NULL 
        arguments$envir <- .GlobalEnv # arguments$envir_
        arguments$probs <- probs
        
        
        if (estimation_method == 'fitted') {
          out_d_ <- do.call(fitted_.bsitar, arguments)
        } else if (estimation_method == 'predict') {
          out_d_ <- do.call(predict_.bsitar, arguments)
        }
        
        arguments$summary <- summary_org
        
        # moved here from below for avg_reffects to work with univariate_by
        if(!is.na(model$model_info$univariate_by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
        }
        
         selectby <- avg_reffects[['reby']]
        # selectby <- c(avg_reffects[['reby']], avg_reffects[['feby']])
        selectover <- avg_reffects[['over']]
        selectby_over <- c(selectby, selectover)
        out_d_ <- get_avg_over(out_d_, newdata = newdata, by = selectby_over,
                               probs = probs, robust = robust)
        
        out_d <- brms::posterior_summary(out_d_, probs = probs, robust = robust)
        # out_d <- out_d_
        
        # if(!is.na(model$model_info$univariate_by)) {
        #   newdata <- newdata %>%
        #     dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
        # }
        
        newdata <- newdata %>%
          dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), .keep_all = T) %>% 
          dplyr::arrange(!! as.name(selectby_over)) %>% 
          droplevels()
        
        # out_d <<- out_d_
        # newdata <<- newdata
        # cbind(newdata %>% arrange(age), out_d) %>% ggplot(., aes(x = age)) + geom_line(aes(y = Estimate))
        
        
        out_summary[['distance']] <-  cbind(newdata,
                                            out_d %>% data.frame() %>%
                                              dplyr::mutate(curve = 'distance')) %>%
          data.frame()
        
      } else if (dist.. == "") {
        out_summary[['distance']] <- NULL
      }
      
      
      summary_org <- arguments$summary
      arguments$summary <- FALSE
      arguments$re_formula <- NULL
      
      # Don't let below arguments$re_formula override the above 'd' and 'v'
      if (velc.. != "") {
        newdata <- newdata___
        if (grepl("^[[:upper:]]+$", velc..)) {
          # arguments$re_formula <- NULL
        } else if (!grepl("^[[:upper:]]+$", velc..)) {
          # arguments$re_formula <- NA
          if (!is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(groupby_fstr, xvar)
          } else if (is.null(groupby_fstr)) {
            groupby_fstr_xvars <- c(xvar)
            
          }
          # newdata <- newdata %>%
          #   dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars))) %>%
          #   dplyr::slice(1) %>% dplyr::ungroup()
        }
        arguments$newdata <- newdata
        arguments$deriv <- 1
        # don't let the ipts to pass again to the fitted_.bsitar
        arguments$ipts <- NULL 
        arguments$envir <- .GlobalEnv # arguments$envir_
        arguments$probs <- probs
        
        
        if (estimation_method == 'fitted') {
          out_v_ <- do.call(fitted_.bsitar, arguments)
        } else if (estimation_method == 'predict') {
          out_v_ <- do.call(predict_.bsitar, arguments)
        }
        
        arguments$summary <- summary_org
        
        # moved here from below for avg_reffects to work with univariate_by
        if(!is.na(model$model_info$univariate_by)) {
          newdata <- newdata %>%
            dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
        }
        
        selectby <- avg_reffects[['reby']]
        # selectby <- c(avg_reffects[['reby']], avg_reffects[['feby']])
        selectover <- avg_reffects[['over']]
        selectby_over <- c(selectby, selectover)
        out_v_ <- get_avg_over(out_v_, newdata = newdata, by = selectby_over,
                               probs = probs, robust = robust)
        
        out_v <- brms::posterior_summary(out_v_, probs = probs, robust = robust)
        
        # out_v <- out_v_
        
        # if(!is.na(model$model_info$univariate_by)) {
        #   newdata <- newdata %>%
        #     dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
        # }
        
        
        newdata <- newdata %>%
          dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), .keep_all = T) %>% 
          dplyr::arrange(!! as.name(selectby_over)) %>% 
          droplevels()
        
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
          # newdata <- newdata %>%
          #   dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
        }
        
        newdata <- newdata %>%
          dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), .keep_all = T) %>% 
          dplyr::arrange(!! as.name(selectby_over)) %>% 
          droplevels()
        
        out_v_ <- t(out_v_)
        out_summary[['parameters']] <-
          get_gparameters(out_v_, newdata, groupby_str_v, summary) # out_v_
      }
      out_summary[['groupby_str_d']] <- groupby_str_d
      out_summary[['groupby_str_v']] <- groupby_str_v
      out_summary[['probtitles']] <- probtitles
    } # if(!is.null(avg_reffects)) {
    
    return(out_summary)
  } # if (arguments$plot) {
  
  
  if (!arguments$plot) {
    newdata <- get.newdata(model, newdata = newdata, 
                           resp = resp, 
                           numeric_cov_at = numeric_cov_at,
                           aux_variables = aux_variables,
                           levels_id = levels_id,
                           ipts = ipts,
                           xrange = xrange)
    

    list_c <- attr(newdata, 'list_c')
    for (list_ci in names(list_c)) {
      assign(list_ci, list_c[[list_ci]])
    }
    check__ <- c('xvar', 'yvar', 'IDvar', 'cov_vars', 'cov_factor_vars', 
                 'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 
                 'uvarby', 'subindicatorsi')
    
    for (check___ in check__) {
      if(!exists(check___)) assign(check___, NULL)
    }
    
    
    if(is.null(arguments$re_formula)) {
      opt <- 'V'
    }
    if(!is.null(arguments$re_formula)) {
      opt <- 'v'
    }
    
    if(!is.null(avg_reffects)) {
      # if (grepl("d", opt, ignore.case = T)) {
      #   index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
      #   dist.. <- substr(opt, index_opt, index_opt)
      #   if (dist.. != "" & grepl("^[[:upper:]]+$", dist..) )
      #     stop("use option 'd' (and not 'D') with avg_reffects" )
      # }
      if (grepl("v", opt, ignore.case = T) ) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (velc.. != "" & grepl("^[[:upper:]]+$", velc..) )
          stop("use option 'v' (and not 'V') with avg_reffects" )
      }
    }
    
    
    if(is.null(avg_reffects)) {
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("v", opt, ignore.case = T)) {
        velc.. <- ""
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
          dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars))) %>%
          dplyr::slice(1) %>% dplyr::ungroup()
      }
      arguments$newdata <- newdata
      arguments$deriv <- 1
      # don't let the ipts to pass again to the fitted_.bsitar
      arguments$ipts <- NULL 
      arguments$envir <- .GlobalEnv # arguments$envir_
      arguments$probs <- probs
     
      if (estimation_method == 'fitted') {
        out_v_ <- do.call(fitted_.bsitar, arguments)
      } else if (estimation_method == 'predict') {
        out_v_ <- do.call(predict_.bsitar, arguments)
      }
      
      out_v__apv_ <- out_v_
      if (!summary) {
        out_v <- call_posterior_summary((out_v_))
      } else if (summary) {
        out_v <- out_v_
      }
      
      if(!is.na(model$model_info$univariate_by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% 
          droplevels()
      }
      
     
      
      # out_v <<- out_v
      # newdata <<- newdata
      # cbind(newdata , out_v) %>% ggplot(., aes(x = age)) + 
      #   geom_line(aes(y = Estimate, group = id))
      
      out_v__apv_ <- t(out_v__apv_)
      parameters <-
        get_gparameters(out_v__apv_, newdata, groupby_str_v, summary)
    } # if(is.null(avg_reffects)) {
    
    
    
    
    if(!is.null(avg_reffects)) {
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
      } else if (!grepl("v", opt, ignore.case = T)) {
        velc.. <- ""
      }
      
      summary_org <- arguments$summary
      arguments$summary <- FALSE
      arguments$re_formula <- NULL
      
      groupby_str_d <- avg_reffects[['feby']]
      groupby_str_v <- avg_reffects[['feby']]
      
      if (grepl("^[[:upper:]]+$", velc..)) {
        # arguments$re_formula <- NULL
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        # arguments$re_formula <- NA
        if (!is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(groupby_fstr, xvar)
        } else if (is.null(groupby_fstr)) {
          groupby_fstr_xvars <- c(xvar)
          
        }
        # newdata <- newdata %>%
        #   dplyr::group_by(dplyr::across(dplyr::all_of(groupby_fstr_xvars))) %>%
        #   dplyr::slice(1) %>% dplyr::ungroup()
      }
      arguments$newdata <- newdata
      arguments$deriv <- 1
      # don't let the ipts to pass again to the fitted_.bsitar
      arguments$ipts <- NULL 
      arguments$envir <- .GlobalEnv # arguments$envir_
      arguments$probs <- probs
      
      if (estimation_method == 'fitted') {
        out_v_ <- do.call(fitted_.bsitar, arguments)
      } else if (estimation_method == 'predict') {
        out_v_ <- do.call(predict_.bsitar, arguments)
      }
      
      arguments$summary <- summary_org
      # out_v_1 <<- out_v_
      
      selectby <- avg_reffects[['reby']]
      # selectby <- c(avg_reffects[['reby']], avg_reffects[['feby']])
      selectover <- avg_reffects[['over']]
      selectby_over <- c(selectby, selectover)
      out_v_ <- get_avg_over(out_v_, newdata = newdata, by = selectby_over,
                             probs = probs, robust = robust)
      
      
      # out_v_ <<- out_v_
      out_v <- out_v_
      
      if(!is.na(model$model_info$univariate_by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
      
      newdata <- newdata %>%
        dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), .keep_all = T) %>% 
        dplyr::arrange(!! as.name(selectby_over)) %>% 
        droplevels()
      
      if(!is.na(model$model_info$univariate_by)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
      
      newdata <- newdata %>%
        dplyr::distinct(., dplyr::across(dplyr::all_of(selectby_over)), .keep_all = T) %>% 
        dplyr::arrange(!! as.name(selectby_over)) %>% 
        droplevels()
     
      
      # out_v <<- out_v
      # newdata <<- newdata
      # cbind(newdata , out_v) %>% ggplot(., aes(x = age)) + 
      #   geom_line(aes(y = Estimate))
      out_v_ <- t(out_v_)
      parameters <-
        get_gparameters(out_v_, newdata, groupby_str_v, summary)
    } # if(!is.null(avg_reffects)) {
    
    return(parameters)
  } # if (!arguments$plot) {
  
} # end gparameters


#' @rdname gparameters.bsitar
#' @export
gparameters <- function(model, ...) {
  UseMethod("gparameters")
}




