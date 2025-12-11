

#' @title Estimate model-based growth parameters for the Bayesian SITAR model
#' 
#' @description The \strong{modelbased_growthparameters()} function estimates
#'   individual growth parameters by mapping population average estimate of age
#'   of interest (such as age at peak growth velocity or age at take off) on to
#'   the individual velocity curves defined by individual level random effects.
#'   Note that option \code{'dpar'} can not be used along with \code{'nlpar'} in 
#'   [brms::posterior_linpred()].
#' 
#' @details Since SITAR is a shape-invariant model, each individual curve has a  
#' peak velocity point that can be mapped by knowing the population average age 
#' at peak velocity. This hold true even when a individual lacks measurements at 
#' the expected turning point.
#' 
#' @inheritParams growthparameters.bgmfit
#' @inheritParams marginal_growthparameters.bgmfit
#' @inheritParams brms::prepare_predictions
#' @inheritParams brms::fitted.brmsfit
#' 
#' @param parameter A single character string or a character vector specifying
#'   the growth parameter(s) to be estimated. Options include age at peak growth
#'   velocity (\code{'apgv'}) and \code{'atgv'} (age at takeoff growth
#'   velocity). The corresponding distance and velocity at \code{'apgv'} and
#'   \code{'atgv'} are computed by default.
#' 
#' @param parameter_method An integer (default = \code{1}) that specifies the
#'   method used to compute individual-level growth parameters from fitted
#'   spline models.
#'
#'   Two methods are available:
#'
#'   \describe{
#'     \item{\code{1} (Model-based differentiation)}{
#'       This method estimates growth parameters by directly leveraging the structure
#'       of the fitted spline model. The spline curve is segmented into pieces of
#'       cubic polynomials. Each segment is then analytically differentiated to obtain
#'       first and second derivatives, which represent the growth velocity and
#'       acceleration, respectively. This approach is more faithful to the underlying
#'       model and is recommended when the goal is to derive precise, model-consistent
#'       growth characteristics.
#'     }
#'     \item{\code{2} (Plug-in adjustment using random effects)}{
#'       This method takes a two-step approach. First, it calculates the population-level
#'       average estimate of the age or time point of interest. Then, it adjusts this
#'       estimate for individual subjects by incorporating random effects from the model.
#'       This approach is computationally simpler and may be preferable when derivative
#'       estimation from the spline model is not feasible or when interpretation in terms
#'       of deviations from the population norm is of primary interest.
#'     }
#'   }
#'
#' @param re_formula Option to indicate whether or not to include
#'   individual/group-level effects in the estimation. When \code{NA} (default),
#'   individual-level effects are excluded, and population average growth
#'   parameters are computed. When \code{NULL}, individual-level effects are
#'   included in the computation, and the resulting growth parameters are
#'   individual-specific.
#' 
#' @param add_xtm A logical (default \code{FALSE}) to indicate whether to 
#'   compute \code{x} and \code{y} adjusted to the mean. Ignored if
#'   \code{parameter_method == 1}. Note that \code{add_xtm} does not affect the 
#'   estimation of parameters.
#' 
#' @param subset_by A logical or character string (default = \code{NULL}) that
#'   determines how to subset the data to retain a single unique row per
#'   \code{id}. This parameter is only used when \code{parameter_method == 2}
#'   and is ignored if \code{add_xtm = TRUE}.
#'
#'   When \code{subset_by = NULL}, the function automatically sets
#'   \code{subset_by} to the value of the \code{by} argument (i.e.,
#'   \code{subset_by = by}). To override this default behavior and skip
#'   subsetting altogether, set \code{subset_by = ""} (an empty string).
#'
#'   This option is useful when multiple rows per \code{id} are present and a
#'   reduction to one representative row per individual is needed for downstream
#'   calculations.
#'   
#'   For \code{parameter_method == 2} with \code{re_formula == NA}, the 
#'   \code{subset_by = "one-row"} will provide one row per parameter.
#' 
#' @param call_function A character string indicating the source of the function
#'   used for computation—either a native \code{R} implementation or a compiled
#'   function exposed from \code{Stan}. Valid options are \code{"R"} and
#'   \code{"Stan"}. The \code{call_function} is ignored when
#'   \code{parameter_method == 2}.
#'   
#'   Though \code{"Stan"} is much faster than \code{R}, The native \code{R}
#'   implementation (\code{call_function = "R"}, default) is recommended when the
#'   number of posterior draws is small. This is because the \code{Stan}
#'   function needs compilation which takes time. The future plan is to allow
#'   integration with \code{Stan}-based computational back ends which then be
#'   exposed along with other functions from \code{Stan} function block.
#' 
#' @param ... Additional arguments passed to the function. 
#' 
#' @return A data frame comprising growth parameter estimates for \strong{age},
#'   \strong{distance} and \strong{velocity}.
#' 
#' @rdname modelbased_growthparameters
#' @export
#' 
#' @seealso [marginaleffects::predictions()]
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
#' berkeley_exfit <- getNsObject(berkeley_exfit)
#' 
#' model <- berkeley_exfit
#' 
#' modelbased_growthparameters(model, ndraws = 10)
#' 
#' }
#' 
modelbased_growthparameters.bgmfit <-
  function(model,
           resp = NULL,
           dpar = NULL,
           ndraws = NULL,
           draw_ids = NULL,
           newdata = NULL,
           datagrid = NULL,
           re_formula = NA,
           newdata2 = NULL,
           allow_new_levels = FALSE,
           sample_new_levels = "gaussian",
           parameter = NULL,
           xrange = 1,
           acg_velocity = 0.10,
           digits = 2,
           numeric_cov_at = NULL,
           aux_variables = NULL,
           levels_id = NULL,
           avg_reffects = NULL,
           idata_method = NULL,
           ipts = NULL,
           seed = 123,
           future = FALSE,
           future_session = 'multisession',
           future_splits = NULL,
           future_method = 'future',
           future_re_expose = NULL,
           usedtplyr = FALSE,
           usecollapse = TRUE,
           parallel = FALSE,
           cores = NULL,
           fullframe = FALSE, 
           average = FALSE, 
           plot = FALSE, 
           showlegends = NULL, 
           variables = NULL,
           deriv = NULL,
           model_deriv = NULL,
           method = 'custom',
           marginals = NULL, 
           preparms = NULL,
           pdraws = FALSE, 
           pdrawso = FALSE,
           pdrawsp = FALSE, 
           pdrawsh = FALSE, 
           comparison = "difference",
           type = NULL,
           by = FALSE,
           bys = NULL,
           conf_level = 0.95,
           transform = NULL,
           transform_draws = NULL,
           cross = FALSE,
           wts = NULL,
           hypothesis = NULL,
           equivalence = NULL,
           eps = NULL,
           constrats_by = FALSE,
           constrats_at = FALSE,
           constrats_subset = FALSE,
           reformat = NULL,
           estimate_center = NULL,
           estimate_interval = NULL,
           dummy_to_factor = NULL, 
           verbose = FALSE,
           expose_function = FALSE,
           usesavedfuns = NULL,
           clearenvfuns = NULL,
           funlist = NULL,
           xvar = NULL,
           idvar = NULL,
           itransform = NULL,
           incl_autocor = TRUE,
           parameter_method = 1,
           subset_by = NULL,
           add_xtm = FALSE,
           call_function = "R",
           newdata_fixed = NULL,
           envir = NULL, 
           ...) {
          
    
    insight::check_if_installed('cheapr', prompt = FALSE, stop = FALSE)
    
    
      try(zz <- insight::check_if_installed(c("marginaleffects"), 
                                            minimum_version = 
                                              get_package_minversion(
                                                'marginaleffects'
                                              ), 
                                            prompt = FALSE,
                                            stop = FALSE))
      
      if(!isTRUE(zz)) {
        message("Please install the latest version of the 'marginaleffects' package",
                "\n ",
                "remotes::install_github('vincentarelbundock/marginaleffects')")
        return(invisible(NULL))
      }
      
      if(usecollapse) {
        usedtplyrcheck <- usedtplyr
        usedtplyr <- FALSE
        if(verbose & usedtplyrcheck) message("Setting usedtplyr = FALSE because ",
                                             "usecollapse = TRUE")
      } else {
        usedtplyr <- usedtplyr
      }
      
      if(usedtplyr) {
        try(zz <- insight::check_if_installed(c("dtplyr"), 
                                              minimum_version =
                                                get_package_minversion(
                                                  'dtplyr'
                                                ),
                                              prompt = FALSE,
                                              stop = FALSE))
        
        try(zz <- insight::check_if_installed(c("data.table"), 
                                              minimum_version =
                                                get_package_minversion(
                                                  'data.table'
                                                ),
                                              prompt = FALSE,
                                              stop = FALSE))
      }
      
      if(usecollapse) {
        try(zz <- insight::check_if_installed(c("collapse"), 
                                              minimum_version =
                                                get_package_minversion(
                                                  'collapse'
                                                ),
                                              prompt = FALSE,
                                              stop = FALSE))
        
        try(zz <- insight::check_if_installed(c("doParallel"), 
                                              minimum_version =
                                                get_package_minversion(
                                                  'doParallel'
                                                ),
                                              prompt = FALSE,
                                              stop = FALSE))
        
        try(zz <- insight::check_if_installed(c("foreach"), 
                                              minimum_version =
                                                get_package_minversion(
                                                  'foreach'
                                                ),
                                              prompt = FALSE,
                                              stop = FALSE))
        
        try(zz <- insight::check_if_installed(c("parallel"), 
                                              minimum_version =
                                                get_package_minversion(
                                                  'parallel'
                                                ),
                                              prompt = FALSE,
                                              stop = FALSE))
      }
      
      callfuns     <- TRUE
      setmarginals <- FALSE
      setpreparms  <- FALSE
      
      if(!is.null(marginals) & !is.null(preparms)) {
        stop("Please specify either marginals or preparms, not both")
      }
      
      if(!is.null(marginals)) {
        setmarginals <- TRUE
        if(method == 'custom') callfuns <- FALSE
        if(method == 'pkg')    callfuns <- FALSE
      }
      
      if(!is.null(preparms)) {
        setmarginals <- TRUE
        setpreparms  <- TRUE
        if(method == 'custom') callfuns <- FALSE
        if(method == 'pkg')    callfuns <- FALSE
      }
      
      
      if(is.null(envir)) {
        envir <- model$model_info$envir
      } else {
        envir <- envir
      }
      
      if(!is.null(transform) & !is.null(transform_draws)) {
        stop("Please specify either transform or transform_draws, not both")
      }
      
      
      # 20.03.2025
      assign_function_to_environment(transform_draws, 'transform_draws', 
                                     envir = NULL)
      model$model_info[['transform_draws']] <- transform_draws
      
      
      # IMP - as nlpar can't be defined along with dpar, dpar must be NULL
      if(!is.null(dpar)) {
        stop("Cannot use 'dpar' and 'nlpar' at the same time for", 
             " posterior_epred.brmsprep().")
      }
      
     
      
      model <- getmodel_info(model = model, 
                             dpar = dpar, 
                             resp = resp, 
                             deriv = NULL, 
                             verbose = verbose)
      
     
      
      
      if(is.null(usesavedfuns)) {
        if(!is.null(model$model_info$exefuns[[1]])) {
          usesavedfuns <- TRUE
        } else if(is.null(model$model_info$exefuns[[1]])) {
          if(expose_function) {
            model <- expose_model_functions(model, envir = envir)
            usesavedfuns <- TRUE
          } else if(!expose_function) {
            usesavedfuns <- FALSE
          }
        }
      } else { # if(!is.null(usesavedfuns)) {
        if(!usesavedfuns) {
          if(expose_function) {
            model <- expose_model_functions(model, envir = envir)
            usesavedfuns <- TRUE
          }
        } else if(usesavedfuns) {
          check_if_functions_exists(model, checks = TRUE, 
                                    usesavedfuns = usesavedfuns)
        }
      }
      
      
      draw_ids_org <- draw_ids
      draw_ids_exe <- FALSE
      ndraws_org   <- ndraws
      ndraws_exe   <- FALSE
      
      if(!is.null(draw_ids)) {
        draw_ids_exe <- TRUE
        draw_ids <- draw_ids
        draw_ids_seq <- draw_ids
      } else if(!is.null(ndraws)) {
        ndraws_exe <- TRUE
        ndraws <- ndraws
        draw_ids_seq <- seq(1, ndraws)
      } else {
        ndraws <- brms::ndraws(model)
        draw_ids_seq <- seq(1, ndraws)
      }
      
      
      
      ########################################################
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
      if(is.null(levels_id) & is.null(idvar)) {
        idvar <- model$model_info[[groupvar_]]
        if (!is.null(model$model_info[[hierarchical_]])) {
          idvar <- model$model_info[[hierarchical_]]
        }
      } else if (!is.null(levels_id)) {
        idvar <- levels_id
      } else if (!is.null(idvar)) {
        idvar <- idvar
      }
      
      xvar  <- xvar
      idvar <- idvar
      if(length(idvar) > 1) idvar <- idvar[1]

      cov_   <- paste0('cov', resp_rev_)
      cov    <- model$model_info[[cov_]]
      uvarby <- model$model_info$univariate_by$by
      
      
      
      
      if(is.null(by)) {
        by <- idvar
      } else if(is.logical(by)) {
        if(is.na(by)) {
          by <- NULL
        } else if(!by) {
          by <- idvar
        }
      } else if(is.character(by)) {
        if(by == "") {
          by <- NULL
        } else {
          by <- by
        }
      } else {
        stop("by should be either TRUE/FALSE, NULL, or an empty string")
      }
      

      if(is.null(by)) {
        if(parameter_method == 2) {
          stop("please specify 'by' argument")
        }
      }
      
      
      check_set_fun <- check_set_fun_transform(model = model, 
                                               which = 'xfuntransform2',
                                               dpar = dpar, 
                                               resp= resp, 
                                               transform = NULL,
                                               auto = FALSE, 
                                               verbose = verbose)
      
      funx_ <- check_set_fun[['setfun']]
      if(check_set_fun[['was_null']]) {
        model$model_info[[check_set_fun[['setfunname']]]] <- funx_
      }
      
      
      
      ########################################################
      ########################################################
      
      check_set_fun <- check_set_fun_transform(model = model, 
                                               which = 'ixfuntransform2',
                                               dpar = dpar, 
                                               resp= resp, 
                                               transform = itransform,
                                               auto = TRUE, 
                                               verbose = verbose)
      
      ifunx_ <- check_set_fun[['setfun']]
      if(check_set_fun[['was_null']]) {
        model$model_info[[check_set_fun[['setfunname']]]] <- ifunx_
      }
      
      ########################################################
 
      itransform_set <- get_itransform_call(itransform = itransform,
                                            model = model, 
                                            newdata = newdata,
                                            dpar = dpar, 
                                            resp = resp,
                                            auto = FALSE,
                                            verbose = verbose)
      
      if(itransform_set == "") {
        if(!isFALSE(pdrawsp)) {
          if(!is.character(pdrawsp)) pdrawsp <- "return"
          selectchoicesr <- c("return", 'add') 
          checkmate::assert_choice(pdrawsp, choices = selectchoicesr)
          if(pdrawsp == 'return' | pdrawsp == 'add') {
            ifunx_ <- function(x)x
          } 
        } # if(!isFALSE(pdrawsp)) {
      } # if(itransform_set == "") {
      ########################################################
      
      
      
      if(is.null(model_deriv)) {
        if(is.null(deriv)) {
          model_deriv <- FALSE
        } else if(deriv == 0) {
          model_deriv <- FALSE
        } else if(deriv == 1) {
          model_deriv <- TRUE
        }
      } else if(!is.null(model_deriv)) {
        if(is.null(deriv) & !model_deriv) {
          deriv <- 0
          model_deriv <- FALSE
        } else if(is.null(deriv) & model_deriv) {
          deriv <- 1
          model_deriv <- TRUE
        }
      }
      
      
      
      
      # 15 06 2025
      allowed_methods <- c('pkg', 'custom')
      if(!method %in% allowed_methods) 
        stop("Argument 'method' should be one of the following:",
             "\n ", 
             collapse_comma(allowed_methods)
        )
      
      
      
      if(method == 'custom') {
        deriv <- 1
        model_deriv <- TRUE
        if(verbose) {
          if(!setpreparms) {
            message(" For method = 'custom', deriv is set to TRUE.\n")
          }
        }
      }
      
      
      
      if (is.null(idata_method)) {
        idata_method <- 'm2'
      }
      
      if(idata_method == 'm1') {
        stop("For marginaleffects based functions, the " ,
             " \n",
             " 'idata_method' argument must be either NULL or 'm2'" )
      }
      
      
      if (is.null(eps)) eps <- 1e-6
      
      
      # Initiate non formalArgs()
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
      drawid <- NULL;
      
      draw <- NULL;
      j <- NULL;
      i <- NULL;
      
      
      peak <- NULL;
      d1 <- NULL;
      d0 <- NULL;
      rowdf <- NULL;
      row_num <- NULL;
      max_pos <- NULL;
      xtm <- NULL;
      ytm <- NULL;
      x <- NULL;
      
      xid <- NULL;
      
      
      allowed_parms <- c(
        'apgv',
        'atgv'
      )
      
      if(verbose) {
        if(setpreparms) {
          message(" For 'preparms', the argument 'parameter' is ignored.",
                  "\n All levels of parameter variable are summarised.",
                  "\n To get summary of a single variable (such as 'apgv'), you can",
                  "\n subset the parameter variable before calling the function\n")
        }
      }
      
      
      
      
      if(is.null(parameter)) {
        parameter <- c('apgv') # c('apgv', 'atgv')
      } else {
        if(length(parameter) > 1) {
          stop("please specify either 'apgv' or 'atgv', not both")
        }
      }
      
      parameter <- base::tolower(parameter)
      
      if (is.null(parameter)) {
        parm <- 'apgv' 
      } else if(length(parameter) == 1 && parameter == 'all') {
        parm <- allowed_parms 
      } else if(length(parameter) == 1) {
        parm <- parameter
      } else if(length(parameter) > 1) {
        # parameter <- base::tolower(parameter)
        for (parameteri in parameter) {
          if(!parameteri %in% allowed_parms) {
            allowed_parms_err <- c(allowed_parms, 'all')
            stop("Allowed parameter options are ", 
                 paste(paste0("'", allowed_parms_err, "'"), collapse = ", ")
            )
          }
        }
        parm <- parameter
      }
      parm <- base::tolower(parm)
      
      
      
      if(length(parm) > 1) {
        if(plot) stop("Please specify only one parameter when plot = TRUE")
      }
      
      
      conf <- conf_level
      probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
      probtitles <- probs[order(probs)] * 100
      probtitles <- paste("Q", probtitles, sep = "")
      set_names_  <- c('Estimate', probtitles)
      
      # if(!is.null(model$model_info$decomp)) {
      #   if(model$model_info$decomp == "QR") model_deriv<- FALSE
      # }
      
      expose_method_set <- model$model_info[['expose_method']]
      
      model$model_info[['expose_method']] <- 'NA' # Over ride 'R'
      
      
      setxcall_   <- match.call()
      post_processing_checks_args <- list()
      post_processing_checks_args[['model']]    <- model
      post_processing_checks_args[['xcall']]    <- setxcall_
      post_processing_checks_args[['resp']]     <- resp
      post_processing_checks_args[['envir']]    <- envir
      post_processing_checks_args[['deriv']]    <- deriv
      post_processing_checks_args[['all']]      <- FALSE
      post_processing_checks_args[['verbose']]  <- verbose
      post_processing_checks_args[['check_d0']] <- FALSE
      post_processing_checks_args[['check_d1']] <- TRUE
      post_processing_checks_args[['check_d2']] <- FALSE
      
      o    <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      
      post_processing_checks_args[['all']]      <- TRUE
      oall <- CustomDoCall(post_processing_checks, post_processing_checks_args)
      post_processing_checks_args[['all']]      <- FALSE
      
      
      
      if(!is.null(funlist)) {
        if(!is.list(funlist)) {
          stop("funlist must be a list")
        } else {
          o <- funlist
        }
      }
      
      
      test <- setupfuns(model = model, resp = resp,
                        o = o, oall = oall, 
                        usesavedfuns = usesavedfuns, 
                        deriv = deriv, envir = envir, 
                        model_deriv = model_deriv, 
                        ...)
      
      if(is.null(test)) return(invisible(NULL))
      
      call_predictions <- TRUE
      call_slopes      <- FALSE
      
      available_d1 <- o[['available_d1']]
      
      if(!available_d1) {
        deriv       <- 0
        model_deriv <- FALSE
        call_predictions <- FALSE
        call_slopes      <- TRUE
        # re-get o[[2]] as _do
        post_processing_checks_args[['deriv']]    <- 0
        o    <- CustomDoCall(post_processing_checks, 
                             post_processing_checks_args)
      }
      
      post_processing_checks_args[['deriv']]    <- deriv
      
      
      # 20.03.2025
      if(!is.null(model$model_info[['sigma_fun_mode']])) {
        sigma_fun_mode <- model$model_info[['sigma_fun_mode']]
        if(dpar == "sigma") {
          if(deriv > 0) {
            if(sigma_fun_mode == "inline") {
              check_fun    <- TRUE
              available_d1 <- FALSE
              model_deriv  <- FALSE
              call_slopes  <- TRUE
            }
          }
        }
      }
      
      
      if(!isTRUE(
        check_pkg_version_exists('brms', 
                                 minimum_version = get_package_minversion('brms'), 
                                 prompt = FALSE,
                                 stop = FALSE,
                                 verbose = FALSE))) {
        if(is.null(check_if_functions_exists(model, o, model$xcall,
                                             usesavedfuns = usesavedfuns))) {
          return(invisible(NULL))
        }
      }
      
      
      
      if(!is.null(model$xcall)) {
        if(grepl("marginal_growthparameters", model$xcall)) {
          xcall <- "marginal_growthparameters"
        }
      } else {
        rlang_trace_back <- rlang::trace_back()
        check_trace_back.bgmfit <- grepl(".bgmfit", rlang_trace_back[[1]])
        if(all(!check_trace_back.bgmfit)) {
          # nothing
        } else {
          rlang_trace_back.bgmfit_i <- min(which(check_trace_back.bgmfit == TRUE))
          rlang_trace_back.bgmfit <- rlang_trace_back[[1]][[rlang_trace_back.bgmfit_i]]
          rlang_call_name <- rlang::call_name(rlang_trace_back.bgmfit)
          xcall <- rlang_call_name
        }
      }
      
      
      check_if_package_installed(model, xcall = xcall)
      
      model$xcall <- xcall
      
      arguments <- get_args_(as.list(match.call())[-1], xcall)
      arguments$model <- model
      arguments$usesavedfuns <- usesavedfuns
      
      
      get.cores_ <- get.cores(arguments$cores)
      
      if(is.null(get.cores_[['max.cores']])) {
        if(is.null(arguments$cores)) 
          get.cores_[['max.cores']] <- future::availableCores() - 1
      }
      
      arguments$cores <- setincores <-  get.cores_[['max.cores']]
      .cores_ps <- get.cores_[['.cores_ps']]
      
      if (future) {
        if(call_function == "Stan") {
          stop("The future = TRUE is not yet suported for call_function = 'Stan'")
        }
        getfutureplan <- future::plan()
        if (future_session == 'multisession') {
          set_future_session <- future_session
        } else if (future_session == 'multicore') {
          set_future_session <- future_session
        } else if (future_session == 'sequential') {
          set_future_session <- future_session
        }
        
        setplanis <- set_future_session
        if(set_future_session == 'sequential') {
          future::plan(setplanis)
        } else {
          future::plan(setplanis, workers = setincores)
        }
        on.exit(future::plan(getfutureplan), add = TRUE)
        
        if (future_session == 'multicore') {
          multthreadplan <- getOption("future.fork.multithreading.enable")
          options(future.fork.multithreading.enable = TRUE)
          on.exit(options("future.fork.multithreading.enable" = multthreadplan), 
                  add = TRUE)
        }
      }
      
      
      
      
      future_splits_exe <- FALSE
      if(!is.null(future_splits)) {
        future_splits_exe <- TRUE
        
        if(is.logical(future_splits)) {
          if(future_splits) {
            if(ndraws_exe) {
              chunk_size_den <- setincores
              ndraws_seq <- sample.int(ndraws)
              chunk_size <- ndraws / chunk_size_den
              future_splits_at <- split(ndraws_seq, 
                                        ceiling(seq_along(ndraws_seq)/chunk_size))
              future_splits_at <- unname(future_splits_at)
            } else if(draw_ids_exe) {
              chunk_size_den <- setincores
              ndraws_seq <- draw_ids
              chunk_size <- length(draw_ids) / chunk_size_den
              future_splits_at <- split(ndraws_seq, 
                                        ceiling(seq_along(ndraws_seq)/chunk_size))
              future_splits_at <- unname(future_splits_at)
            }
          } else if(!future_splits) {
            future_splits_exe <- FALSE
            future_splits <- future_splits
          }
          
        } else if(is.list(future_splits)) {
          future_splits_at <- future_splits
          
        } else if(is.vector(future_splits)) {
          if(!is.numeric(future_splits)) {
            stop("future_splits must be a numeric vector of lenghth 2")
          } else if(length(future_splits) == 1) {
            if(draw_ids_exe) ndraws_exe <- FALSE
            if(ndraws_exe) {
              chunk_size_den <- future_splits
              ndraws_seq <- sample.int(ndraws)
              chunk_size <- ndraws / chunk_size_den
              future_splits_at <- split(ndraws_seq, 
                                        ceiling(seq_along(ndraws_seq)/chunk_size))
              future_splits_at <- unname(future_splits_at)
            } else if(draw_ids_exe) {
              chunk_size_den <- future_splits
              ndraws_seq <- draw_ids
              chunk_size <- length(draw_ids) / chunk_size_den
              future_splits_at <- split(ndraws_seq, 
                                        ceiling(seq_along(ndraws_seq)/chunk_size))
              future_splits_at <- unname(future_splits_at)
            } else if(is.null(ndraws_org)) {
              chunk_size_den <- future_splits
              ndraws_seq <- sample.int(ndraws)
              chunk_size <- ndraws / chunk_size_den
              future_splits_at <- split(ndraws_seq, 
                                        ceiling(seq_along(ndraws_seq)/chunk_size))
              future_splits_at <- unname(future_splits_at)
            } else if(!is.null(draw_ids_org)) {
              chunk_size_den <- future_splits
              ndraws_seq <- draw_ids
              chunk_size <- length(draw_ids) / chunk_size_den
              future_splits_at <- split(ndraws_seq, 
                                        ceiling(seq_along(ndraws_seq)/chunk_size))
              future_splits_at <- unname(future_splits_at)
            }
          } else if(length(future_splits) != 2) {
            stop("future_splits must be a numeric vector of lenghth 2")
          } else {
            future_splits_at <- parallel::splitIndices(future_splits[1], 
                                                       future_splits[2])
          }
        }
      }
      
      
      if(future_splits_exe) {
        if(plot) {
          stop("future_splits can not be used when plot = TRUE")
        }
        if(method == 'pkg') {
          stop("future_splits can not be used when method = 'pkg'")
        }
      }
      
      
      
      
      if(!future_splits_exe) {
        future_splits_exe_future <- FALSE
        future_splits_exe_dofuture <- FALSE
      } else if(future_splits_exe) {
        if(future_method == 'future') {
          future_splits_exe_future <- TRUE
          future_splits_exe_dofuture <- FALSE
        }
        if(future_method == 'dofuture') {
          future_splits_exe_future <- FALSE
          future_splits_exe_dofuture <- TRUE
        }
      }
      
      
      
      re_expose <- FALSE
      if (future) {
        need_future_re_expose_cpp <- FALSE
        if(any(grepl("pstream__",
                     deparse(model$model_info$exefuns[[1]])))) {
          need_future_re_expose_cpp <- TRUE
        }
        
        if(is.null(future_re_expose)) {
          if(setplanis == "multisession") {
            if(need_future_re_expose_cpp) {
              re_expose <- TRUE
              if(verbose) {
                message("For multisession plan, argument 'future_re_expose' has been set as TRUE")
              }
            } else if(!need_future_re_expose_cpp) {
              if(verbose) {
                message("To speed up the calulations, it is advised to set future_re_expose = TRUE")
              }
            }
          }
        } else if(!is.null(future_re_expose)) {
          if(future_re_expose) {
            re_expose <- TRUE
          } else if(!future_re_expose) {
            if(!need_future_re_expose_cpp) {
              # if(expose_method_set == "R") {
              if(verbose) {
                message("To speed up the calulations, it is advised ",
                        "to set 'future_re_expose = TRUE'")
              }
            } 
            if(need_future_re_expose_cpp & setplanis == "multisession") {
              # if(expose_method_set != "R") {
              stop("For plan 'multisession', the functions need to be ",
                   "\n ",
                   "re_exposed by setting 'future_re_expose = TRUE'")
            }
          }
        }
      } # if (future) {
      
      
      
      
      if (!future) {
        future_splits_at <- NULL
        future_splits_exe <- FALSE
        future_splits_exe_future <- FALSE
        future_splits_exe_dofuture <- FALSE
      }
      
      
      
      
      
      
      full.args <- evaluate_call_args(cargs = as.list(match.call())[-1], 
                                      fargs = arguments, 
                                      dargs = list(...), 
                                      verbose = verbose)
      
      full.args$model       <- model
      full.args$model_deriv <- model_deriv
      
      
      if(is.null(full.args$hypothesis) & is.null(full.args$equivalence)) {
        plot <- plot
      } else {
        plot <- FALSE
        if(verbose & plot) {
          message("Argument plot = TRUE is not allowed when either hypothesis ", 
                  "or equivalence is not NULL",
                  "\n ",
                  "Therefor, setting 'plot = FALSE'") 
        }
      }
      
      
      
      
      full.args$newdata <- newdata
      
      full.args <- 
        sanitize_CustomDoCall_args(what = "CustomDoCall", 
                                   arguments = full.args, 
                                   # check_formalArgs = marginal_growthparameters.bgmfit,
                                   check_formalArgs = modelbased_growthparameters.bgmfit,
                                   check_formalArgs_exceptions = NULL,
                                   check_trace_back = NULL,
                                   envir = parent.frame())
      
      
      if(parameter_method == 1) {
        full.args[['newdata_fixed']] <- 2
      } else if(parameter_method == 2) {
        full.args[['newdata_fixed']] <- 0
      }
      
      
      full.args$dpar    <- dpar
      
      get.newdata_args <- list()
      for (i in methods::formalArgs(get.newdata)) {
        get.newdata_args[[i]] <- full.args[[i]]
      }
      
      # For modelbased_growthparameters, dpar argument is not allowed.
      # get.newdata_args$ipts <- full.args$ipts <- ipts <- 
      #   set_for_check_ipts(ipts = ipts, nipts = 50, dpar = dpar, verbose = verbose)
      
      full.args$newdata <- newdata <- CustomDoCall(get.newdata, 
                                                   get.newdata_args)

      # Interpolation points
      if(!exists('check_fun'))    check_fun    <- FALSE
      if(!exists('available_d1')) available_d1 <- FALSE

      if(!setpreparms) {
        full.args$ipts <- ipts <- check_ipts(ipts = full.args$ipts, 
                                             nipts = NULL, 
                                             check_fun  = check_fun, 
                                             available_d1 = available_d1, 
                                             xcall = NULL, verbose = verbose)
      }
      
      
      
      
      if(!is.na(uvarby)) {
        uvarby_ind <- paste0(uvarby, resp)
        varne <- paste0(uvarby, resp)
        if(usedtplyr) {
          newdata <- newdata %>% dtplyr::lazy_dt() %>% 
            dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
        } else if(!usedtplyr) {
          newdata <- newdata %>% 
            dplyr::mutate(!! uvarby_ind := 1) %>% droplevels()
        }
      }
      full.args$newdata <- newdata
      full.args[["..."]] <- NULL
      
      
      if(!is.null(full.args[['transform_draws']])) {
        full.args[['transform']] <- full.args[['transform_draws']]
        if(verbose) message("'transform' set based on 'transform_draws'")
      }
      
      modelbased_arguments <- full.args
      
      model$xcall <- 'modelbased_growthparameters'
      
      #######################################################################
      
      allknots_ <- paste0('knots', resp_rev_)
      allknots   <- model$model_info[[allknots_]]
      knots    <- allknots[2:(length(allknots)-1)]
      bknots   <- c(allknots[1], allknots[length(allknots)])
      xknots   <- allknots
      
      create_abcd_names_vector   <- c('a', 'b', 'c', 'd')
      create_s_names_vector      <- paste0("s", 2:length(allknots)-1)
      create_abcd_s_names_vector <- c(create_abcd_names_vector, create_s_names_vector)
      
      
      fixed_ <- paste0('fixed', resp_rev_)
      check_for_fixed   <- model$model_info[[fixed_]]
      
      random_ <- paste0('random', resp_rev_)
      check_for_random   <- model$model_info[[random_]]
      
      onlyfixed <- setdiff(check_for_fixed, check_for_random)
      
      
      posterior_linpred_args <- list()
      posterior_linpred_args[['object']]                <- model
      posterior_linpred_args[['newdata']]               <- newdata
      
      posterior_linpred_args[['re_formula']]            <- re_formula
      posterior_linpred_args[['resp']]                  <- resp
      posterior_linpred_args[['dpar']]                  <- dpar
      posterior_linpred_args[['ndraws']]                <- ndraws
      posterior_linpred_args[['draw_ids']]              <- draw_ids_seq
      posterior_linpred_args[['allow_new_levels']]      <- allow_new_levels
      posterior_linpred_args[['sample_new_levels']]     <- sample_new_levels
      posterior_linpred_args[['transform']]             <- transform
      posterior_linpred_args[['incl_autocor']]          <- incl_autocor
      posterior_linpred_args[['nlpar']]                 <- NULL
      posterior_linpred_args[['incl_thres']]            <- NULL
      posterior_linpred_args[['sort']]                  <- FALSE
      posterior_linpred_args[['nug']]                   <- NULL
      posterior_linpred_args[['smooths_only']]          <- FALSE
      posterior_linpred_args[['offset']]                <- TRUE
      posterior_linpred_args[['newdata2']]              <- NULL
      posterior_linpred_args[['new_objects']]           <- NULL
      posterior_linpred_args[['point_estimate']]        <- NULL
      posterior_linpred_args[['ndraws_point_estimate']] <- 1
      posterior_linpred_args[['group_vars']]            <- NULL
      posterior_linpred_args[['req_vars']]              <- NULL
      
      
      # dpar must be NULL as we are using nlpar
      posterior_linpred_args[['dpar']]                  <- NULL
      
      
      if(parameter_method == 1) {
      ##############################################################
        
        if(add_xtm) {
          stop("For add_xtm = TRUE, please use parameter_method = 1")
        }
        
      nlpar_fixed <- array(NA, dim = c(length(eval(draw_ids_seq)),
                                       nrow(newdata),
                                       length(create_abcd_s_names_vector)))
      
      nlpar_random <- nlpar_coeffs <- nlpar_fixed
      
      
      # CustomDoCall
      
      ikl <- 0
      for (i in create_abcd_s_names_vector) {
        ikl <- ikl + 1
        if(i %in% names(model$formula$ pforms)) {
          posterior_linpred_args[['nlpar']]                 <- i
          if(is.null(posterior_linpred_args[['re_formula']])) {
            nlpar_coeffs[ , , ikl] <- do.call(brms::posterior_linpred, posterior_linpred_args)
          } else if(!is.null(posterior_linpred_args[['re_formula']])) {
            nlpar_fixed [ , , ikl] <- do.call(brms::posterior_linpred, posterior_linpred_args)
          }
        } else {
          if(is.null(posterior_linpred_args[['re_formula']])) {
            nlpar_coeffs[ , , ikl] <- 0
          } else if(!is.null(posterior_linpred_args[['re_formula']])) {
            nlpar_fixed [ , , ikl] <- 0
          }
        }
      }
      dimnames(nlpar_fixed) [[3]] <- create_abcd_s_names_vector
      dimnames(nlpar_coeffs)[[3]] <- create_abcd_s_names_vector
      
      ikl <- 0
      for (i in create_abcd_s_names_vector) {
        ikl <- ikl + 1
        if(!i %in% onlyfixed) {
          nlpar_random[ , , ikl] <- nlpar_coeffs[ , , ikl] - nlpar_fixed [ , , ikl]
        } else {
          nlpar_random[ , , ikl] <- 0
        }
      }
      dimnames(nlpar_random)[[3]] <- create_abcd_s_names_vector
      
      
      
      if(is.null(posterior_linpred_args[['re_formula']])) {
        set_frame <- nlpar_coeffs
      } else if(!is.null(posterior_linpred_args[['re_formula']])) {
        set_frame <- nlpar_fixed
      }
      
       # print(str(set_frame))
      
      splinenames         <- create_s_names_vector
      get_dims            <- dim(set_frame)
      set_draws_n         <- get_dims[1]
      get_data_cols.org   <- colnames(newdata)
      set_frame_rows      <- cbind(newdata, set_frame[1, ,])
      set_frame_rows      <- set_frame_rows %>% dplyr::mutate(row_index = dplyr::row_number())
      set_frame_rows_cols <- c(get_data_cols.org, 'row_index')
      
      if(is.null(subset_by)) {
        subset_data_by <- "acrossrows"
      } else {
        subset_data_by <- subset_by
      }
      
      if(!is.null(subset_data_by)) {
        if(subset_data_by == "acrossrows") {
          subset_data_by_names <- create_abcd_s_names_vector
        } else if(subset_data_by == "fixed") {
          subset_data_by_names <- onlyfixed
        } else if(subset_data_by == "random") {
          subset_data_by_names <- create_abcd_names_vector
        } else {
          if(!is.character(subset_by)) {
           # stop("subset_by must be a single character string")
          } else 
          subset_data_by_names <- subset_data_by
        }
        # data.table:::unique.data.table
        set_dataf_m   <- unique(data.table::as.data.table(set_frame_rows),
                                                        by = subset_data_by_names) 
      } else if(is.null(subset_data_by)) {
        set_dataf_m <- set_frame_rows
        subset_data_by_names <- NULL
      }
      
      set_nrows_n   <- nrow(set_dataf_m)
      set_dataf_m   <- set_dataf_m %>% dplyr::select(dplyr::all_of(set_frame_rows_cols))
      set_dataf_m   <- set_dataf_m %>% dplyr::mutate(rowdf = dplyr::row_number())
      set_dataf_m   <- set_dataf_m %>% data.frame()
     
      set_dataf_m_collapse <- set_dataf_m %>% 
        collapse::roworderv(c(idvar, xvar), decreasing = F) 
      
      SplineCall <- model$model_info$SplineCall
      
      if(call_function == "R") {
        GS_gps_parms_stan   <- NULL;
      }
      
      
      if(call_function == "Stan") {
        GS_gps_parms_stan_str_get_function_scode <- GS_gps_parms_stan_str_get()
        support_GS_gps_parms_stan_str_get_function_scode <- support_GS_gps_parms_stan_str_get()
        
        full_GS_gps_parms_stan_str_get_function_scode <- 
          paste0(support_GS_gps_parms_stan_str_get_function_scode,
                 "\n",
                 GS_gps_parms_stan_str_get_function_scode)
        
        GS_gps_parms_stan_str_get_function_scode <-
          paste0("functions {",
                 "\n",
                 full_GS_gps_parms_stan_str_get_function_scode,
                 "\n",
                 "} // end functions block")
        
        if(verbose) message("Prepraring Stan function...")
        rmodel <- rstan::stanc(model_code = GS_gps_parms_stan_str_get_function_scode)
        rstan::expose_stan_functions(rmodel)
        if(verbose) message("Ready Stan function...")
      } # if(call_function == "Stan") {
      
      
      
      if(call_function == "R") {
        GS_gps_parms_assign <- GS_gps_parms_R
      } else if(call_function == "Stan") {
        GS_gps_parms_assign <- GS_gps_parms_stan
      }
      
      
      drawni = 1;
      spline_precomputed_indicator = 1;
      
      degree = 3;
      set_spread = 2;
      pieces_dim = length(xknots)-1;
      degree_dim = degree + 1;
      degree_dim_set_spread = degree_dim*set_spread;
      shift_indicator = 1;
      return_indicator = 1;
      spline_subset_indicator = 1;
      parm_mat_dim = 7;
      deriv_mat_dim = 7;
      if(drawni == 0) {
        parm_mat_dim = 6;
        deriv_mat_dim = 6;
      }
      if(drawni > 0) {
        parm_mat_dim = 6+1;
        deriv_mat_dim = 6+1;
      }
      
    
      array_dim = set_nrows_n
      
      if(return_indicator == 1) {
        array_mat_dim_1 = pieces_dim;
        array_mat_dim_2 = parm_mat_dim;
      }
      if(return_indicator == 2) {
        array_mat_dim_1 = pieces_dim*degree_dim_set_spread;
        array_mat_dim_2 = deriv_mat_dim;
      }
      
      
      if(call_function == "R") {
        spline_eval_array <- array(NA, dim = c(degree_dim, pieces_dim, pieces_dim ))
        xg_array          <- array(NA, dim = c(degree_dim, pieces_dim))
        xg_curve_array    <- array(NA, dim = c(degree_dim_set_spread, pieces_dim))
      } 
      
      if(call_function == "Stan") {
        spline_eval_array <- list()
        xg_array          <- list()
        xg_curve_array    <- list()
      } 
      
      
      SplineCall[[2]]        <- quote(xg)
      
     
    
      for (i in 1:pieces_dim) {
        if(call_function == "R") {
          xg_array[, i]            <- seq_fun_R(xknots[i], xknots[i+1], degree_dim)
          xg_curve_array[, i]      <- seq_fun_R(xknots[i], xknots[i+1], degree_dim_set_spread)
          xg                       <- unlist(xg_array[,i])
          spline_eval_array[, , i] <- eval(SplineCall);
        }
        if(call_function == "Stan") {
          xg                     <- seq_fun_R(xknots[i], xknots[i+1], degree_dim)
          xg_array[[i]]          <- xg
          xg_curve_array[[i]]    <- seq_fun_R(xknots[i], xknots[i+1], degree_dim_set_spread)
          spline_eval_array[[i]] <- eval(SplineCall);
        }
      } # end for (i in 1:pieces_dim) {
      
    
      
      if (future) {
        future_globals_list = list(spline_eval_array = spline_eval_array,
                                   xg_array = xg_array,
                                   xg_curve_array = xg_curve_array,
                                   subset_data_by = subset_data_by,
                                   `%>%` = bsitar::`%>%`,
                                   subset_data_by_names = subset_data_by_names,
                                   create_abcd_names_vector = create_abcd_names_vector,
                                   create_s_names_vector = create_s_names_vector,
                                   spline_subset_indicator = spline_subset_indicator,
                                   spline_precomputed_indicator = spline_precomputed_indicator,
                                   shift_indicator = shift_indicator,
                                   set_spread = set_spread,
                                   call_function = call_function,
                                   xknots = xknots,
                                   degree = degree,
                                   my_counter = my_counter,
                                   GS_gps_parms_assign = GS_gps_parms_assign,
                                   wraper_for_drawni = wraper_for_drawni)
      }
      
      return_indicator      <- 1
      
      if(!future) {
        collect_draws_parm    <- list()
        for (drawni in 1:set_draws_n) {
          setdat_mat <- set_frame[drawni, ,]
          collect_draws_parm[[drawni]] <- wraper_for_drawni(setdat_mat = setdat_mat, 
                                                            drawni = drawni,
                                                            callvia = 'base',
                                                            return_indicator = return_indicator, 
                                                            subset_data_by = subset_data_by,
                                                            subset_data_by_names = subset_data_by_names,
                                                            create_abcd_names_vector = create_abcd_names_vector,
                                                            create_s_names_vector = create_s_names_vector,
                                                            xknots = xknots,
                                                            degree = degree,
                                                            spline_subset_indicator = spline_subset_indicator,
                                                            spline_precomputed_indicator = spline_precomputed_indicator,
                                                            shift_indicator = shift_indicator,
                                                            set_spread = set_spread,
                                                            spline_eval_array = spline_eval_array,
                                                            xg_array = xg_array,
                                                            xg_curve_array = xg_curve_array,
                                                            call_function = call_function,
                                                            GS_gps_parms_assign = GS_gps_parms_assign)
      
        } # for (drawni in 1:1) {
      } else if(future) {
        setdat_mat_list <- list()
        for (drawni in 1:set_draws_n) {
          setdat_mat <- set_frame[drawni, ,]
          setdat_mat_list[[drawni]] <- setdat_mat
        }
        
        environment(wraper_for_drawni) <- environment()
        
        wraper_for_drawni_future_mapply <- function(.x, 
                                                    drawni, 
                                                    callvia,
                                                    return_indicator #, 
                                                    # subset_data_by,
                                                    # subset_data_by_names,
                                                    # create_abcd_names_vector,
                                                    # create_s_names_vector,
                                                    # xknots,
                                                    # degree,
                                                    # spline_subset_indicator = spline_subset_indicator,
                                                    # spline_precomputed_indicator = spline_precomputed_indicator,
                                                    # shift_indicator,
                                                    # set_spread,
                                                    # spline_eval_array,
                                                    # xg_array,
                                                    # xg_curve_array,
                                                    # call_function,
                                                    # GS_gps_parms_assign
                                                    ) {
          wraper_for_drawni(setdat_mat = .x, 
                            drawni = NULL,
                            callvia = 'future',
                            return_indicator = return_indicator, 
                            subset_data_by = subset_data_by,
                            subset_data_by_names = subset_data_by_names,
                            create_abcd_names_vector = create_abcd_names_vector,
                            create_s_names_vector = create_s_names_vector,
                            xknots = xknots,
                            degree = degree,
                            spline_subset_indicator = spline_subset_indicator,
                            spline_precomputed_indicator = spline_precomputed_indicator,
                            shift_indicator = shift_indicator,
                            set_spread = set_spread,
                            spline_eval_array = spline_eval_array,
                            xg_array = xg_array,
                            xg_curve_array = xg_curve_array,
                            call_function = call_function,
                            GS_gps_parms_assign = GS_gps_parms_assign)
        }
        my_counter$reset()
        collect_draws_parm <- future.apply::future_mapply(setdat_mat_list,  
                                                          FUN = wraper_for_drawni_future_mapply, 
                                                          MoreArgs = list(
                                                            drawni = NULL,
                                                            callvia = 'future',
                                                            return_indicator = return_indicator
                                                          ) ,
                                                          future.globals = future_globals_list,
                                                          future.seed = TRUE)
      } # end else if(future) {
      
      organize_draws_parm <- array(unlist(collect_draws_parm), dim = c(pieces_dim, 
                                                                       parm_mat_dim, 
                                                                       array_dim,
                                                                       set_draws_n)) 
      
      bind_draws_parm <- apply(organize_draws_parm, 2, identity)
      
      if(future) {
        bind_draws_parm <- assign_new_sequence(mat = bind_draws_parm, col = 7)
      }
      
      
      names_parm <- c("x", "d0", "d1", "peak", "piece", "rowdf")
      names_parm <- c(names_parm, "drawid")
      
      parm_sort_keys      <- c('drawid','rowdf','piece') # could be c(7, 6) etc
      parm_is.finite_keys <- c('x','d0','d1', 'peak')
      parm_is.na_keys     <- c('x','d0','d1', 'peak')
      
      dt <- data.table::as.data.table(bind_draws_parm)
      dimnames(dt)[[2]] <- names_parm
      
      # Sort rows
      dt <- dt[order(dt[, .SD, .SDcols = parm_sort_keys], na.last = FALSE)]

      
      if(call_function == "Stan") {
        dt <- dt[, lapply(.SD, function(x) ifelse(!is.finite(x), NA, x))]
      }
        
     
      dt[, 'x']         <- model$model_info$ixfuntransform2(dt[, 'x'])
      
      peak_takeoff_data_draw     <- dt
      peak_takeoff_data_draw     <- na.omit(peak_takeoff_data_draw, cols=parm_is.na_keys, invert=FALSE)

      peak_indices    <- peak_takeoff_data_draw[, .I[peak == 1 & d1 == max(d1) ], by = c('drawid', 'rowdf', 'piece')]$V1
      peak_data_draw  <- collapse::fsubset(peak_takeoff_data_draw, peak_indices)
      
      if(nrow(peak_data_draw) > 0) {
        apgv_draw    <- peak_data_draw %>% 
          collapse::fmutate(rowdf = as.integer(rowdf)) %>% 
          collapse::fmutate(estimate = x) %>% 
          collapse::fmutate(parameter = 'apgv') 
        pgv_draw    <- peak_data_draw %>% 
          collapse::fmutate(rowdf = as.integer(rowdf)) %>% 
          collapse::fmutate(estimate = d1) %>% 
          collapse::fmutate(parameter = 'pgv') 
        spgv_draw    <- peak_data_draw %>% 
          collapse::fmutate(rowdf = as.integer(rowdf)) %>% 
          collapse::fmutate(estimate = d0) %>% 
          collapse::fmutate(parameter = 'spgv') 
        
        apgv_draw <- set_dataf_m_collapse %>% collapse::join(apgv_draw, on = 'rowdf', 
                                                             how = "right", 
                                                             verbose = FALSE)
        pgv_draw  <- set_dataf_m_collapse %>% collapse::join(pgv_draw,  on = 'rowdf', 
                                                             how = "right", 
                                                             verbose = FALSE)
        spgv_draw <- set_dataf_m_collapse %>% collapse::join(spgv_draw, on = 'rowdf', 
                                                             how = "right", 
                                                             verbose = FALSE)
        
        all_peak_data_draw <- collapse::rowbind(apgv_draw, pgv_draw, spgv_draw) 
        
       
        marginal_growthparameters_args <- modelbased_arguments
        marginal_growthparameters_args[['preparms']] <- all_peak_data_draw
        marginal_growthparameters_args[['by']] <- idvar
        
        peak_parameters <- CustomDoCall(marginal_growthparameters, 
                                        marginal_growthparameters_args)

      } # if(nrow(peak_data_draw) > 0) {
      
      
      takeoff_data_draw <- peak_takeoff_data_draw %>%
        collapse::fgroup_by(drawid, rowdf) %>%  # 1st fgroup_by
        collapse::fmutate(max_pos = which.max(d1)) %>%
        # fmutate(group_key = finteraction(drawid, rowdf)) %>%
        # fgroup_by(group_key) %>%             # 2nd fgroup_by
        collapse::fmutate(row_num = collapse::seqid(d1)) %>%
        collapse::fungroup() %>%                       # 1st fungroup
        collapse::fsubset(peak == 0 & row_num <= max_pos) %>%
        collapse::fgroup_by(drawid, rowdf) %>%  # 3rd fgroup_by
        collapse::fslice(n = 1, how = "last") %>%
        collapse::fungroup()                           # 2nd fungroup
      
      if(nrow(takeoff_data_draw) > 0) {
        atgv_draw    <- takeoff_data_draw %>% 
          collapse::fmutate(rowdf = as.integer(rowdf)) %>% 
          collapse::fmutate(estimate = x) %>% 
          collapse::fmutate(parameter = 'atgv') 
        tgv_draw    <- takeoff_data_draw %>% 
          collapse::fmutate(rowdf = as.integer(rowdf)) %>% 
          collapse::fmutate(estimate = d1) %>% 
          collapse::fmutate(parameter = 'tgv') 
        stgv_draw    <- takeoff_data_draw %>% 
          collapse::fmutate(rowdf = as.integer(rowdf)) %>% 
          collapse::fmutate(estimate = d0) %>% 
          collapse::fmutate(parameter = 'stgv') 
        
        atgv_draw <- set_dataf_m_collapse %>% collapse::join(atgv_draw, on = 'rowdf', 
                                                             how = "right", 
                                                             verbose = FALSE)
        tgv_draw  <- set_dataf_m_collapse %>% collapse::join(tgv_draw,  on = 'rowdf', 
                                                             how = "right", 
                                                             verbose = FALSE)
        stgv_draw <- set_dataf_m_collapse %>% collapse::join(stgv_draw, on = 'rowdf', 
                                                             how = "right", 
                                                             verbose = FALSE)
        
        all_takeoff_data_draw <- collapse::rowbind(atgv_draw, tgv_draw, stgv_draw)
        
        marginal_growthparameters_args <- modelbased_arguments
        marginal_growthparameters_args[['preparms']] <- all_takeoff_data_draw
        marginal_growthparameters_args[['by']] <- idvar
        
        takeoff_parameters <- CustomDoCall(marginal_growthparameters, 
                                           marginal_growthparameters_args)
        
      } # if(nrow(takeoff_data_draw) > 0) {
      
      if(parm == 'apgv') {
        out <- peak_parameters
      }
      
      if(parm == 'atgv') {
        out <- takeoff_parameters
      }
      
      peak_names.ors__2 <- base::tolower(colnames(out))
      indices_lastn     <- 2 # for est and q1 qu
      firstup_indices   <- (length(peak_names.ors__2)-indices_lastn):length(peak_names.ors__2)
      lower_case__2 <- peak_names.ors__2[1:(length(peak_names.ors__2)-indices_lastn-1)]
      upper_case__2 <- firstup(peak_names.ors__2[firstup_indices])
      
      peak_names.ors__2 <- c(lower_case__2, upper_case__2)
      
      out <- data.table::setnames(out, peak_names.ors__2)
      
      return(out)
      
      ##############################################################
    } else if(parameter_method == 2) { # if(parameter_method == 1)
      ##############################################################
      
      if(is.null(subset_by)) {
        subset_data_by <- by
      } else if(!is.null(subset_by)) {
        if(is.na(subset_by) | subset_by == "") {
          subset_data_by <- NULL
        } else {
          subset_data_by <- subset_by
        }
      }
      
      xoffset    <- model$model_info$xoffset
      
      SplineCall_d0             <- model$model_info$SplineCall
      SplineCall_d0[[2]]        <- quote(setx0)
      SplineCall_d1             <- SplineCall_d0
      SplineCall_d1[['derivs']] <- 1
      
      
      if(call_function == "R") {
        GS_gps_parms_assign <- GS_gps_parms_R
      } else if(call_function == "Stan") {
        stop("please use call_function == 'R'")
        GS_gps_parms_assign <- GS_gps_parms_stan
      }
      
     
      nlpar_fixed <- array(NA, dim = c(length(eval(draw_ids_seq)),
                                       nrow(newdata),
                                       length(create_abcd_s_names_vector)))
      
      nlpar_random <- nlpar_coeffs <- nlpar_fixed

      posterior_linpred_args_re_NA <-
        posterior_linpred_args_re_NULL <-
        posterior_linpred_args

      posterior_linpred_args_re_NA  [['re_formula']] <- NA
      posterior_linpred_args_re_NULL[['re_formula']] <-NULL

      ikl <- 0
      for (i in create_abcd_s_names_vector) {
        ikl <- ikl + 1
        if(i %in% names(model$formula$ pforms)) {
          posterior_linpred_args_re_NULL[['nlpar']] <- i
          posterior_linpred_args_re_NA  [['nlpar']] <- i
          nlpar_coeffs[ , , ikl] <- do.call(brms::posterior_linpred, posterior_linpred_args_re_NULL)
          nlpar_fixed [ , , ikl] <- do.call(brms::posterior_linpred, posterior_linpred_args_re_NA)
        } else {
          nlpar_coeffs[ , , ikl] <- 0
          nlpar_fixed [ , , ikl] <- 0
        }
      }
      dimnames(nlpar_fixed) [[3]] <- create_abcd_s_names_vector
      dimnames(nlpar_coeffs)[[3]] <- create_abcd_s_names_vector

      ikl <- 0
      for (i in create_abcd_s_names_vector) {
        ikl <- ikl + 1
        if(!i %in% onlyfixed) {
          nlpar_random[ , , ikl] <- nlpar_coeffs[ , , ikl] - nlpar_fixed [ , , ikl]
        } else {
          nlpar_random[ , , ikl] <- 0
        }
      }
      dimnames(nlpar_random)[[3]] <- create_abcd_s_names_vector

      set_frame_fixed     <- nlpar_fixed
      set_frame_random    <- nlpar_random

      set_frame           <- set_frame_random

      splinenames         <- create_s_names_vector
      get_dims            <- dim(set_frame)
      set_draws_n         <- get_dims[1]
      get_data_cols.org   <- colnames(newdata)
      set_frame_rows      <- cbind(newdata, set_frame[1, ,])
      set_frame_rows      <- set_frame_rows %>% dplyr::mutate(row_index = dplyr::row_number())
      set_frame_rows_cols <- c(get_data_cols.org, 'row_index')

      set_dataf_m          <- set_frame_rows
      subset_data_by_names <- NULL
      
      set_nrows_n   <- nrow(set_dataf_m)
      set_dataf_m   <- set_dataf_m %>% dplyr::select(dplyr::all_of(set_frame_rows_cols))
      set_dataf_m   <- set_dataf_m %>% dplyr::mutate(rowdf = dplyr::row_number())
      set_dataf_m   <- set_dataf_m %>% data.frame()
      # For later collapse::join(on = 'rowdf')
      set_dataf_m_collapse <- set_dataf_m %>% 
        collapse::roworderv(c(idvar, xvar), decreasing = F)

      
      if(is.null(re_formula)) {
        set_pdrawsp <- 'return'
        set_pdraws  <- FALSE
      } else if(!is.null(re_formula)) {
        set_pdrawsp <- FALSE
        set_pdraws  <- 'adds'
      }
      
      marginal_growthparameters_args <- modelbased_arguments
      marginal_growthparameters_args[['re_formula']] <- NA
      marginal_growthparameters_args[['pdrawsp']]    <- set_pdrawsp
      marginal_growthparameters_args[['pdraws']]     <- set_pdraws
      marginal_growthparameters_args[['parameter']]  <- parm
      marginal_growthparameters_args[['newdata']]    <- newdata
      marginal_growthparameters_args[['draw_ids']]   <- draw_ids_seq
      marginal_growthparameters_args[['by']]         <- by
      
      marginal_growthparameters_args[['newdata_fixed']] <- NULL
      
      onex0 <- CustomDoCall(marginal_growthparameters, 
                            marginal_growthparameters_args)
      

      if(!is.null(re_formula)) {
        if(add_xtm) { 
          stop("Please set 're_formula = NULL' for 'add_xtm = TRUE'")
        }
        onex0             <- onex0[['estimate']]
        peak_names.ors__2 <- base::tolower(colnames(onex0))
        colnames(onex0)   <- peak_names.ors__2
        set0_newdata                 <- newdata
        attr(set0_newdata, "list_c") <- NULL
        set0_newdata                 <- set0_newdata %>% 
                                        dplyr::left_join(., 
                                                         onex0, by = by)        
        set0_newdata[[xvar]]         <- NULL
        
        # below using newdata_fixed = 0, so apply fun here 
        
        set0_newdata[[xvar]]         <- funx_(set0_newdata[['estimate']])
        
        marginal_draws_args <- modelbased_arguments
        marginal_draws_args[['newdata']]       <- set0_newdata
        marginal_draws_args[['newdata_fixed']] <- 0
         marginal_draws_args[['by']]           <- by
        
        marginal_draws_args[['deriv']] <- 0
        get_size <- CustomDoCall(marginal_draws, 
                                 marginal_draws_args)
        
        marginal_draws_args[['deriv']] <- 1
        get_velc <- CustomDoCall(marginal_draws, 
                                 marginal_draws_args)
       
        
        colnames(get_size) <- base::tolower(colnames(get_size))
        colnames(get_velc) <- base::tolower(colnames(get_velc))
        
        get_velc <- cbind.data.frame(set0_newdata[, c("parameter", by)], 
                                     get_velc)
        get_velc[["parameter"]] <- 'PGV'
        get_velc <- get_velc %>% 
          dplyr::distinct(!! as.name(dplyr::all_of(by)), .keep_all = T)
        
        get_size <- cbind.data.frame(set0_newdata[, c("parameter", by)], 
                                     get_size)
        get_size[["parameter"]] <- 'SPGV'
        
        get_size <- get_size %>% 
          dplyr::distinct(!! as.name(dplyr::all_of(by)), .keep_all = T)
        
        
        onex0 <- dplyr::bind_rows(onex0, get_velc, get_size)
        
        marginal_draws(model, newdata = set0_newdata, 
                       newdata_fixed = 0,
                       deriv = 1) %>% nrow()
        
        
        if(!is.null(subset_data_by)) {
          group_by_indices <- "parameter"
          if(subset_data_by == "one-row") {
            group_by_indices <- group_by_indices
          } else {
            group_by_indices <- c(group_by_indices, subset_data_by)
          }
          onex0 <- onex0 %>% data.table::as.data.table()
          onex0 <- onex0[onex0[, .I[1:1], by = group_by_indices]$V1]
          onex0 <- onex0 %>% data.frame()
        }
       
        peak_names.ors__2 <- base::tolower(colnames(onex0))
        indices_lastn     <- 2 # for est and q1 qu
        firstup_indices   <- (length(peak_names.ors__2)-indices_lastn):length(peak_names.ors__2)
        lower_case__2 <- peak_names.ors__2[1:(length(peak_names.ors__2)-indices_lastn-1)]
        upper_case__2 <- firstup(peak_names.ors__2[firstup_indices])
        
        peak_names.ors__2 <- c(lower_case__2, upper_case__2)
        # change the case, after bind with xtm
        onex0 <- data.table::setnames(onex0, peak_names.ors__2)
        
        return(onex0)
      } # if(!is.null(re_formula)) {
      
      
      wraper_for_drawni_2 <- function(drawni, 
                                      nlpar_fixed, 
                                      nlpar_random,
                                      create_s_names_vector,
                                      add_xtm,
                                      callvia) {
        
        if(callvia == 'base') {
          drawniid <- drawni
        } else if(callvia == 'future') {
          drawniid <- my_counter$next_value()
          pid <- Sys.getpid()
          time_us <- as.numeric(Sys.time()) * 1e6
          drawniid <- drawniid+time_us+drawniid
        }
        
        setdat_mat_fixed  <- nlpar_fixed [drawni, ,]
        setdat_mat_random <- nlpar_random[drawni, ,]
        spmat             <- setdat_mat_fixed[, create_s_names_vector]
        setx0             <- setdat_mat_fixed[, 'Xestimate']
        
        # atgv might be Na
        if(all(is.na(setx0))) {
          x.adj <- NA_real_
          y.adj <- NA_real_
          v.adj <- NA_real_
        } else {
          setx              <- funx_(setx0)
          setx              <- setx - setdat_mat_fixed[, 'b']
          x.adj             <- setx/exp(setdat_mat_random[,"c"]) + 
            setdat_mat_random[,"b"] + 
            setdat_mat_fixed[, 'b']
          x.adj             <- ifunx_(x.adj)
          setx0             <- setx * exp(setdat_mat_fixed[, 'c'])
          y.adj             <- rowSums(eval(SplineCall_d0) * spmat) + 
            setdat_mat_random[,"a"] + 
            setdat_mat_random[,"d"] * x.adj +
            setdat_mat_fixed[, 'a']
          v.adj             <- rowSums(eval(SplineCall_d1) * spmat) * 
            exp(setdat_mat_fixed[, 'c'] + 
                  setdat_mat_random[,"c"]) + 
            setdat_mat_random[,"d"] 
        } # end else if(all(is.na(setx0))) {
          
        if(add_xtm) {
          setxx.adj_xtm <- funx_(setdat_mat_fixed[, 'xvar'])
          x.adj_xtm <- (setxx.adj_xtm - setdat_mat_random[,"b"]) * exp(setdat_mat_random[,"c"])
          x.adj_xtm <- ifunx_(x.adj_xtm)
          
          y.adj_xtm <- setdat_mat_fixed[, 'yvar'] - 
            setdat_mat_random[,"a"] - 
            setdat_mat_random[,"d"] * setxx.adj_xtm
        }
        
        
        mat.adj[, 1] <- x.adj
        mat.adj[, 2] <- y.adj
        mat.adj[, 3] <- v.adj
        if(add_xtm) {
          mat.adj[, 4] <- x.adj_xtm
          mat.adj[, 5] <- y.adj_xtm
          mat.adj[, 6] <- drawniid
          mat.adj[, 7] <- setdat_mat_fixed[, 'fomerge']
          mat.adj[, 8] <- setdat_mat_fixed[, 'xid']
        } else {
          mat.adj[, 4] <- drawniid
          mat.adj[, 5] <- setdat_mat_fixed[, 'fomerge']
          mat.adj[, 6] <- setdat_mat_fixed[, 'xid']
        }
        return(mat.adj)
      }
     
      
      onex00 <- set_dataf_m_collapse %>% collapse::join(onex0, on = by,
                                              how = "left",
                                              multiple = TRUE,
                                              verbose = FALSE)
      
      onex00 <- data.table::as.data.table(onex00)

      array_dim     <- set_nrows_n
      pieces_dim    <- 1
      parm_mat_dim  <- 6
      
      onex00 <- onex00 %>% 
        collapse::fmutate(fomerge =  
                            collapse::finteraction(onex00 %>% 
                                                     collapse::fselect(c(by)), 
                                                   factor = FALSE))
      
      xid_by_onex00 <- c("drawid", "parameter", 'fomerge')
      
      onex00$xid <- setorderv(onex00, xid_by_onex00)[, .(xid=seq_len(.N)), 
                                              by = xid_by_onex00]$xid
      
      
      which_dim                  <- 3
      nlpar_fixed_names_dim3     <- attr(nlpar_fixed, "dimnames")[[which_dim]]
      nlpar_fixed_names_dim3_add <- c('Xestimate', 'fomerge', 'xid')

      extend_array <- cbind(onex00[['estimate']], 
                            onex00[['fomerge']], 
                            onex00[['xid']])
      
      if(add_xtm) {
        xvar_ <- model$model_info$xvar
        yvar_ <- model$model_info$yvar
        extend_array <- cbind(extend_array, newdata[[xvar_]], newdata[[yvar_]])
        nlpar_fixed_names_dim3_add <- c(nlpar_fixed_names_dim3_add, "xvar", 'yvar')
        parm_mat_dim <- parm_mat_dim + 2
      }
      
      
      Sliced <- aperm(`dim<-`(t(extend_array), 
                              c(ncol(extend_array), 
                                dim(nlpar_fixed)[2], 
                                dim(nlpar_fixed)[1])), 
                      c(3, 2, 1))
      
      nlpar_fixed <- abind::abind(nlpar_fixed, Sliced, along = 3)
      
      attr(nlpar_fixed, "dimnames")[[which_dim]] <- c(nlpar_fixed_names_dim3, 
                                                      nlpar_fixed_names_dim3_add)
      
      mat.adj            <- matrix(NA_real_, nrow = array_dim, ncol = parm_mat_dim )
      # set as data.table for rowbind()
      mat.adj            <- data.table::as.data.table(mat.adj)
      if(!future) {
        collect_draws_parm <- list()
        for (drawni in 1:set_draws_n) {
          collect_draws_parm[[drawni]] <- 
            wraper_for_drawni_2(drawni = drawni, 
                                nlpar_fixed = nlpar_fixed,
                                nlpar_random = nlpar_random,
                                create_s_names_vector = create_s_names_vector,
                                add_xtm = add_xtm,
                                callvia = 'base')
        } # for (drawni in 1:set_draws_n) {
      } # if(!future) {
      
      

      if(future) {
        # setup future
        environment(wraper_for_drawni_2) <- environment()
        future_globals_list = list( mat.adj = mat.adj,
                                   `%>%` = bsitar::`%>%`,
                                   my_counter = my_counter)
        # call future
        my_counter$reset()
        collect_draws_parm <- future.apply::future_lapply(1:set_draws_n, 
                                      FUN = function(drawni, ...) 
                                        wraper_for_drawni_2(drawni = drawni, 
                                                            nlpar_fixed = nlpar_fixed,
                                                            nlpar_random = nlpar_random,
                                                            create_s_names_vector = create_s_names_vector,
                                                            add_xtm = add_xtm,
                                                            callvia = 'future'),
          future.globals = future_globals_list)
    } # end else if(future) {
      
      
      if(add_xtm) {
        names_parm      <- c("x", "d0", "d1", "xtm", "ytm", "drawid")
      } else {
        names_parm      <- c("x", "d0", "d1", "drawid")
      }
      names_parm_temp <- c(names_parm, "fomerge", "xid")
      

      bind_draws_parm <- collect_draws_parm %>% collapse::rowbind() %>%
        collapse::setrename(names_parm_temp)
      

      if(future) {
        which_cols <- 4
        if(add_xtm) {
          which_cols <- which_cols + 2
        }
        bind_draws_parm <- assign_new_sequence(mat = bind_draws_parm %>% 
                                                 as.matrix(), col = which_cols) %>% 
          data.table::as.data.table()
      }
      
     peak_data_draw <- bind_draws_parm %>% 
       collapse::join(onex00, on = c("drawid", 'fomerge', "xid"),
                      how = "right",
                      multiple = FALSE,
                      verbose = FALSE) 
     

     peak_data_draw_select <- c(names_parm, get_data_cols.org)
     peak_data_draw        <- collapse::fselect(peak_data_draw, 
                                                peak_data_draw_select)
     
     # One can subset peak_data_draw but not xtm_data_draw
     xtm_data_draw <- peak_data_draw
     
     if(!is.null(subset_data_by)) {
       group_by_indices <- c("drawid", subset_data_by) 
       peak_data_draw <- peak_data_draw[peak_data_draw[, .I[1:1], 
                                                       by = group_by_indices]$V1]
     }
     
     xid_by <- c("drawid", "parameter", "id") 
     
     if(parm == 'apgv') {
       parameter_names_vec <- c('apgv', 'pgv', 'spgv')
     }
     if(parm == 'atgv') {
       parameter_names_vec <- c('atgv', 'tgv', 'stgv')
     }
     
      if(nrow(peak_data_draw) > 0) {
        apgv_draw    <- peak_data_draw %>% 
          collapse::fmutate(estimate = x) %>% 
          collapse::fmutate(parameter = parameter_names_vec[1]) 
        pgv_draw    <- peak_data_draw %>% 
          collapse::fmutate(estimate = d1) %>% 
          collapse::fmutate(parameter = parameter_names_vec[2]) 
        spgv_draw    <- peak_data_draw %>% 
          collapse::fmutate(estimate = d0) %>% 
          collapse::fmutate(parameter = parameter_names_vec[3])
        
        all_peak_data_draw <- collapse::rowbind(apgv_draw, pgv_draw, spgv_draw)
        
        all_peak_data_draw$xid <- setorderv(all_peak_data_draw, 
                                            xid_by)[, .(xid=seq_len(.N)), 
                                                    by = xid_by]$xid
        
        marginal_growthparameters_args <- modelbased_arguments
        marginal_growthparameters_args[['preparms']] <- all_peak_data_draw
        marginal_growthparameters_args[['by']] <- c(by, 'xid')
        
        peak_parameters <- CustomDoCall(marginal_growthparameters,
                                        marginal_growthparameters_args)
      
        
        peak_names.ors__ <- colnames(peak_parameters)
        data.table::setnames(peak_parameters, tolower(names(peak_parameters)))
        
        peak_roworderv_vars <- 'parameter'
        peak_roworderv_vars <- c(peak_roworderv_vars, by, 'xid')
        peak_parameters <- collapse::roworderv(peak_parameters, peak_roworderv_vars)

        if(add_xtm) {
          xtm_draw    <- xtm_data_draw %>% 
            collapse::fmutate(estimate = xtm) %>% 
            collapse::fmutate(parameter = 'xtm') 
          ytm_draw    <- xtm_data_draw %>% 
            collapse::fmutate(estimate = ytm) %>% 
            collapse::fmutate(parameter = 'ytm') 
          all_tm_data_draw <- collapse::rowbind(xtm_draw, ytm_draw)
          
          all_tm_data_draw$xid <- setorderv(all_tm_data_draw,
                                              xid_by)[, .(xid=seq_len(.N)),
                                                      by = xid_by]$xid
          
          marginal_growthparameters_args <- modelbased_arguments
          marginal_growthparameters_args[['preparms']] <- all_tm_data_draw
          marginal_growthparameters_args[['by']] <- c(by, 'xid')
          
          tm_parameters <- CustomDoCall(marginal_growthparameters,
                                          marginal_growthparameters_args)
 
          tm_names.ors__ <- colnames(tm_parameters)
          data.table::setnames(tm_parameters, tolower(names(tm_parameters)))
          tm_roworderv_vars <- 'parameter'
          tm_roworderv_vars <- c(tm_roworderv_vars, by, 'xid')
          tm_parameters <- collapse::roworderv(tm_parameters, tm_roworderv_vars)
          tm_parameters <- data.table::setnames(tm_parameters, tm_names.ors__)
          peak_parameters <- collapse::rowbind(peak_parameters, tm_parameters)
          
        } # if(add_xtm) {
      } # if(nrow(peak_data_draw) > 0) {
      
     peak_parameters   <- peak_parameters %>% collapse::fselect(-xid)
     
     peak_names.ors__2 <- peak_names.ors__[ !grepl('xid', peak_names.ors__)]
     
     indices_lastn     <- 2 # for est and q1 qu
     firstup_indices   <- (length(peak_names.ors__2)-indices_lastn):length(peak_names.ors__2)
     lower_case__2 <- peak_names.ors__2[1:(length(peak_names.ors__2)-indices_lastn-1)]
     upper_case__2 <- firstup(peak_names.ors__2[firstup_indices])
     
     peak_names.ors__2 <- c(lower_case__2, upper_case__2)

     peak_parameters <- data.table::setnames(peak_parameters, peak_names.ors__2)
     
      return(peak_parameters) 
      ##############################################################
    } # if(parameter_method == 1) { else if(parameter_method == 2) {
    
}



#' @rdname modelbased_growthparameters
#' @export
modelbased_growthparameters <- function(model, ...) {
  UseMethod("modelbased_growthparameters")
}


