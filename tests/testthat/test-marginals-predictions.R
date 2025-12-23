

# Skip test for local R CMD Check but run on GitHub

# skip_if_not_ci()



###############################################################################
# Test marginals vs marginaleffects
###############################################################################

test_that("test-marginals-predictions", {
  skip_on_cran()
  
  
  # devtools::load_all()
  
  if(test_univariate_fit_cov) {
    fit               = readRDS(testthat::test_path("models", 
                                                    "univariate_fit_cov.rds")) 
    resp              = uvar_resp
  } else if(test_univariate_fit_cov) {
    fit               = readRDS(testthat::test_path("models", 
                                                    "multivariate_fit_cov.rds"))  
    resp              = mvar_resp
  } else {
    skip(message = 
           "Both test_univariate_fit_cov and test_univariate_fit_cov FALSE")
  }
  
  
  # Need to re-assign functions to this local test environment
  fit <- bsitar::expose_model_functions(fit, expose = F)
  
  
  
  model = fit
  draw_ids = draw_ids
  ndraws = NULL
  re_formula = NA
  xvar <- 'age'
  newdata <- bsitar:::get.newdata(model)
  
  ###############################################################################
  # marginaleffects::predictions vs marginal_draws - average = FALSE
  ###############################################################################
  
  variables = NULL
  vcov = TRUE
  conf_level = 0.95
  type = NULL
  by = FALSE
  byfun = NULL
  wts = FALSE
  transform = NULL
  hypothesis = NULL
  equivalence = NULL
  df = Inf
  numderiv = "fdforward"
  
  # marginaleffects::predictions
  marginaleffects_predictions_args <- list()
  marginaleffects_predictions_args[['model']] <- model
  marginaleffects_predictions_args[['newdata']] <- newdata
  marginaleffects_predictions_args[['variables']] <- variables
  marginaleffects_predictions_args[['vcov']] <- vcov
  marginaleffects_predictions_args[['conf_level']] <- conf_level
  marginaleffects_predictions_args[['type']] <- type
  marginaleffects_predictions_args[['by']] <- by
  marginaleffects_predictions_args[['byfun']] <- byfun
  marginaleffects_predictions_args[['wts']] <- wts
  marginaleffects_predictions_args[['transform']] <- transform
  marginaleffects_predictions_args[['hypothesis']] <- hypothesis
  marginaleffects_predictions_args[['equivalence']] <- equivalence
  marginaleffects_predictions_args[['df']] <- df
  marginaleffects_predictions_args[['numderiv']] <- numderiv
  #
  marginaleffects_predictions_args[['draw_ids']] <- draw_ids
  marginaleffects_predictions_args[['ndraws']] <- ndraws
  marginaleffects_predictions_args[['re_formula']] <- re_formula
  
  
  marginaleffects_out <- do.call(marginaleffects::predictions, 
                                 marginaleffects_predictions_args) %>%
    data.frame()
  
  
  marginal_draws_args <- marginaleffects_predictions_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- FALSE
  marginal_draws_args[['newdata_fixed']] <- 0
  
  
  marginal_out <- do.call(marginal_draws, 
                          marginal_draws_args) %>% data.frame()
  
  
  marginal_draws_args_custom                  <- marginal_draws_args
  marginal_draws_args_custom[['method']]      <- 'custom'
  marginal_draws_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  
  marginal_draws_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  
  out_1 <- round(marginaleffects_out$estimate, 2)
  out_2 <- round(marginal_out$estimate, 2)
  out_3 <- out_mdF
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  
  ###############################################################################
  # marginaleffects::avg_predictions vs marginal_draws - average = TRUE
  ###############################################################################
  # check if any new argument needed for avg_predictions than predictions
  setdiff(methods::formalArgs(marginaleffects::predictions),
          methods::formalArgs(marginaleffects::avg_predictions))
  
  
  marginaleffects_avg_predictions_args <- marginaleffects_predictions_args
  
  marginaleffects_out <- do.call(marginaleffects::predictions, 
                                 marginaleffects_avg_predictions_args) %>% 
    data.frame()
  
  
  marginal_draws_args <- marginaleffects_avg_predictions_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- TRUE
  marginal_draws_args[['newdata_fixed']] <- 0
  
  marginal_out <- do.call(marginal_draws, 
                          marginal_draws_args) %>% data.frame()
  
  
  marginal_draws_args_custom                  <- marginal_draws_args
  marginal_draws_args_custom[['method']]      <- 'custom'
  marginal_draws_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  
  marginal_draws_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  out_1 <- round(marginaleffects_out$estimate, 2)
  out_2 <- round(marginal_out$estimate, 2)
  out_3 <- out_mdF
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  ###############################################################################
  # marginaleffects::plot_predictions vs marginal_draws - average = FALSE
  # with by 
  ###############################################################################
  
  # check if any new argument needed for avg_predictions than predictions
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::predictions),
            methods::formalArgs(marginaleffects::plot_predictions))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_predictions),
            methods::formalArgs(marginaleffects::predictions))
  
  
  marginaleffects_plot_predictions_args <- marginaleffects_predictions_args
  for (i in args_remove) {
    marginaleffects_plot_predictions_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_predictions_args[['condition']] <- NULL
  marginaleffects_plot_predictions_args[['points']] <- 0
  marginaleffects_plot_predictions_args[['rug']] <- FALSE
  marginaleffects_plot_predictions_args[['gray']] <- FALSE
  marginaleffects_plot_predictions_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_predictions_args[['by']] <- xvar
  marginaleffects_plot_predictions_args[['condition']] <- NULL
  
  
  marginaleffects_plot <- do.call(marginaleffects::plot_predictions, 
                                  marginaleffects_plot_predictions_args)
  
  marginaleffects_plot_predictions_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_predictions, 
                                 marginaleffects_plot_predictions_args) %>% 
    data.frame()
  
  
  marginal_draws_args <- marginaleffects_plot_predictions_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- TRUE
  marginal_draws_args[['newdata_fixed']] <- 0
  marginal_draws_args[['plot']]          <- TRUE
  
  marginal_plot <- do.call(marginal_draws, 
                           marginal_draws_args)
  
  
  marginal_draws_args[['plot']]          <- FALSE
  marginal_out <- do.call(marginal_draws, 
                          marginal_draws_args) %>% data.frame()
  
  
  marginal_draws_args_custom                  <- marginal_draws_args
  marginal_draws_args_custom[['method']]      <- 'custom'
  marginal_draws_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  
  marginal_draws_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  
  # Note that for plot, length differs between marginaleffects and marginals
  # But rounded mean should be same
  out_1 <- round(mean(marginaleffects_out$estimate), 1)
  out_2 <- round(mean(marginal_out$estimate), 1)
  out_3 <- round(mean(marginal_out_custom_mdF$estimate), 1)
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  # # for testthat, add tolerance
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  
  ###############################################################################
  # marginaleffects::plot_predictions vs marginal_draws - average = FALSE
  # with condition 
  ###############################################################################
  
  # check if any new argument needed for avg_predictions than predictions
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::predictions),
            methods::formalArgs(marginaleffects::plot_predictions))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_predictions),
            methods::formalArgs(marginaleffects::predictions))
  
  
  marginaleffects_plot_predictions_args <- marginaleffects_predictions_args
  for (i in args_remove) {
    marginaleffects_plot_predictions_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_predictions_args[['condition']] <- NULL
  marginaleffects_plot_predictions_args[['points']] <- 0
  marginaleffects_plot_predictions_args[['rug']] <- FALSE
  marginaleffects_plot_predictions_args[['gray']] <- FALSE
  marginaleffects_plot_predictions_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_predictions_args[['by']] <- NULL
  marginaleffects_plot_predictions_args[['condition']] <- xvar
  
  marginaleffects_plot <- do.call(marginaleffects::plot_predictions, 
                                  marginaleffects_plot_predictions_args)
  
  marginaleffects_plot_predictions_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_predictions, 
                                 marginaleffects_plot_predictions_args) %>% 
    data.frame()
  
  marginal_draws_args <- marginaleffects_plot_predictions_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- TRUE
  marginal_draws_args[['newdata_fixed']] <- 0
  marginal_draws_args[['plot']]          <- TRUE
  
  
  marginal_plot <- do.call(marginal_draws, 
                           marginal_draws_args)
  
  marginal_draws_args[['plot']]          <- FALSE
  
  marginal_out <- do.call(marginal_draws, 
                          marginal_draws_args) %>% data.frame()
  
  
  marginal_draws_args_custom                  <- marginal_draws_args
  marginal_draws_args_custom[['method']]      <- 'custom'
  marginal_draws_args_custom[['model_deriv']] <- T
  marginal_out_custom_mdT <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  
  marginal_draws_args_custom[['model_deriv']] <- F
  marginal_out_custom_mdF <- do.call(marginal_draws, 
                                     marginal_draws_args_custom) %>% data.frame()
  
  
  # Note that for plot, length differs between marginaleffects and marginals
  # But rounded mean should be same
  out_1 <- round(mean(marginaleffects_out$estimate), 1)
  out_2 <- round(mean(marginal_out$estimate), 1)
  out_3 <- round(mean(marginal_out_custom_mdF$estimate), 1)
  
  expect_equal(out_1,   out_2,   tolerance = 0.01)
  expect_equal(out_mdT, out_mdF, tolerance = 0.01)
  expect_equal(out_2,   out_3,   tolerance = 0.01)
  
  # # for testthat, add tolerance
  # if(!identical(out_1, out_2)) {
  #   stop()
  # }
  # 
  # if(!identical(out_mdT, out_mdF)) {
  #   stop()
  # }
  # 
  # if(!identical(out_2, out_3)) {
  #   stop()
  # }
  
  
  
  
  
})
