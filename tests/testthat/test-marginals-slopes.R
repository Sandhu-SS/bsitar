


# Skip test for local R CMD Check but run on GitHub

# skip_if_not_ci()



###############################################################################
# Test marginals vs marginaleffects
###############################################################################

test_that("test-marginals-slopes", {
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
  
  marginal_deriv_args <- 1
  
  ###############################################################################
  # marginaleffects::slopes vs marginal_draws - average = FALSE
  ###############################################################################
  
  variables = xvar # note variables set as xvar for slope
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
  
  # marginaleffects::slopes
  marginaleffects_slopes_args <- list()
  marginaleffects_slopes_args[['model']] <- model
  marginaleffects_slopes_args[['newdata']] <- newdata
  marginaleffects_slopes_args[['variables']] <- variables
  marginaleffects_slopes_args[['vcov']] <- vcov
  marginaleffects_slopes_args[['conf_level']] <- conf_level
  marginaleffects_slopes_args[['type']] <- type
  marginaleffects_slopes_args[['by']] <- by
  marginaleffects_slopes_args[['byfun']] <- byfun
  marginaleffects_slopes_args[['wts']] <- wts
  marginaleffects_slopes_args[['transform']] <- transform
  marginaleffects_slopes_args[['hypothesis']] <- hypothesis
  marginaleffects_slopes_args[['equivalence']] <- equivalence
  marginaleffects_slopes_args[['df']] <- df
  marginaleffects_slopes_args[['numderiv']] <- numderiv
  #
  marginaleffects_slopes_args[['draw_ids']] <- draw_ids
  marginaleffects_slopes_args[['ndraws']] <- ndraws
  marginaleffects_slopes_args[['re_formula']] <- re_formula
  
  
  marginaleffects_out <- do.call(marginaleffects::slopes, 
                                 marginaleffects_slopes_args) %>%
    data.frame()
  
  
  marginal_draws_args <- marginaleffects_slopes_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- FALSE
  marginal_draws_args[['newdata_fixed']] <- 0
  marginal_draws_args[['deriv']]         <- marginal_deriv_args
  
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
  # marginaleffects::avg_slopes vs marginal_draws - average = TRUE
  ###############################################################################
  # check if any new argument needed for avg_slopes than slopes
  setdiff(methods::formalArgs(marginaleffects::slopes),
          methods::formalArgs(marginaleffects::avg_slopes))
  
  
  marginaleffects_avg_slopes_args <- marginaleffects_slopes_args
  
  marginaleffects_out <- do.call(marginaleffects::slopes, 
                                 marginaleffects_avg_slopes_args) %>% 
    data.frame()
  
  
  marginal_draws_args <- marginaleffects_avg_slopes_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- TRUE
  marginal_draws_args[['newdata_fixed']] <- 0
  marginal_draws_args[['deriv']]         <- marginal_deriv_args
  
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
  # marginaleffects::plot_slopes vs marginal_draws - average = FALSE
  # with by 
  ###############################################################################
  
  # check if any new argument needed for avg_slopes than slopes
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::slopes),
            methods::formalArgs(marginaleffects::plot_slopes))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_slopes),
            methods::formalArgs(marginaleffects::slopes))
  
  
  marginaleffects_plot_slopes_args <- marginaleffects_slopes_args
  for (i in args_remove) {
    marginaleffects_plot_slopes_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_slopes_args[['condition']] <- NULL
  marginaleffects_plot_slopes_args[['points']] <- NULL # 0
  marginaleffects_plot_slopes_args[['rug']] <- FALSE
  marginaleffects_plot_slopes_args[['gray']] <- FALSE
  marginaleffects_plot_slopes_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_slopes_args[['by']] <- xvar
  marginaleffects_plot_slopes_args[['condition']] <- NULL
  
  
  marginaleffects_plot <- do.call(marginaleffects::plot_slopes, 
                                  marginaleffects_plot_slopes_args)
  
  marginaleffects_plot_slopes_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_slopes, 
                                 marginaleffects_plot_slopes_args) %>% 
    data.frame()
  
  
  marginal_draws_args <- marginaleffects_plot_slopes_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- TRUE
  marginal_draws_args[['newdata_fixed']] <- 0
  marginal_draws_args[['plot']]          <- TRUE
  marginal_draws_args[['deriv']]         <- marginal_deriv_args
  
  
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
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  # out_1 <- round(marginaleffects_out$estimate, 2)
  # out_2 <- round(marginal_out$estimate, 2)
  # out_3 <- out_mdT
  
  out_1 <- round(mean(marginaleffects_out$estimate), 1)
  out_2 <- round(mean(marginal_out$estimate), 1)
  out_3 <- round(mean(marginal_out_custom_mdF$estimate), 1)
  
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
  # marginaleffects::plot_slopes vs marginal_draws - average = FALSE
  # with condition 
  ###############################################################################
  
  # check if any new argument needed for avg_slopes than slopes
  args_remove <-
    setdiff(methods::formalArgs(marginaleffects::slopes),
            methods::formalArgs(marginaleffects::plot_slopes))
  
  args_add <-
    setdiff(methods::formalArgs(marginaleffects::plot_slopes),
            methods::formalArgs(marginaleffects::slopes))
  
  
  marginaleffects_plot_slopes_args <- marginaleffects_slopes_args
  for (i in args_remove) {
    marginaleffects_plot_slopes_args[[i]] <- NULL
  }
  
  # args_add
  marginaleffects_plot_slopes_args[['condition']] <- NULL
  marginaleffects_plot_slopes_args[['points']] <- NULL # 0
  marginaleffects_plot_slopes_args[['rug']] <- FALSE
  marginaleffects_plot_slopes_args[['gray']] <- FALSE
  marginaleffects_plot_slopes_args[['draw']] <- TRUE
  # need to change by
  marginaleffects_plot_slopes_args[['by']] <- NULL
  marginaleffects_plot_slopes_args[['condition']] <- xvar
  
  marginaleffects_plot <- do.call(marginaleffects::plot_slopes, 
                                  marginaleffects_plot_slopes_args)
  
  marginaleffects_plot_slopes_args[['draw']] <- FALSE
  marginaleffects_out <- do.call(marginaleffects::plot_slopes, 
                                 marginaleffects_plot_slopes_args) %>% 
    data.frame()
  
  marginal_draws_args <- marginaleffects_plot_slopes_args
  
  marginal_draws_args[['method']]        <- 'pkg'
  marginal_draws_args[['reformat']]      <- FALSE
  marginal_draws_args[['average']]       <- TRUE
  marginal_draws_args[['newdata_fixed']] <- 0
  marginal_draws_args[['plot']]          <- TRUE
  marginal_draws_args[['deriv']]         <- marginal_deriv_args
  
  
  
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
  
  marginal_draws_args_custom_plot <- marginal_draws_args_custom
  marginal_draws_args_custom_plot[['plot']] <- TRUE
  do.call(marginal_draws, 
          marginal_draws_args_custom_plot)
  
  out_mdT <- round(marginal_out_custom_mdT$estimate, 2)
  out_mdF <- round(marginal_out_custom_mdF$estimate, 2)
  
  # out_1 <- round(marginaleffects_out$estimate, 2)
  # out_2 <- round(marginal_out$estimate, 2)
  # out_3 <- out_mdT
  
  out_1 <- round(mean(marginaleffects_out$estimate), 1)
  out_2 <- round(mean(marginal_out$estimate), 1)
  out_3 <- round(mean(marginal_out_custom_mdF$estimate), 1)
  
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
  
  
  
  
  
  
  
  
  
  
  
})
