
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


###############################################################################
# Test bsitar with rcs settings
###############################################################################

test_that("bsitar works fully with rcs settings", {
  skip_on_cran()
 
  test_scode <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = TRUE,
                       get_standata = FALSE, 
                       chains = 1, cores = 1, iter = 10, 
                       backend = "rstan",  
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refres = 0, silent = 2,
                       seed = 123)
  
  expect_type(test_scode, "character")
  
  
  test_sdata <- bsitar(x = age, y = height, id = id, 
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = TRUE, 
                       chains = 1, cores = 1, iter = 10, 
                       backend = "rstan",  
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refres = 0, silent = 2,
                       seed = 123)
  
  expect_type(test_sdata, "list")
  
  
  suppressWarnings(suppressMessages({
    test_fit <- bsitar(x = age, y = height, id = id,
                       data = test_data_male,  df = 4,
                       stype = list(type = 'rcs', normalize = TRUE),
                       get_stancode = FALSE,
                       get_standata = FALSE,
                       chains = 1, cores = 1, iter = 10,
                       backend = "rstan",
                       sample_prior = "no",
                       threads = threading(NULL),
                       # init = '0',
                       init = NULL, # Don't use default random with init_r = 0.5
                       vcov_init_0 = TRUE,
                       refres = 0, silent = 2,
                       seed = 123)
  }))
  
  

  true_sbetas <- c(128.05, 0.00, 0.00, 5.02, 4.28, -12.27, -13.26)
  
  test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2)
  
  expect_equal(true_sbetas, test_sbetas, tolerance = 0.01)
  
  test_gparms <- get_growthparameters(test_fit, re_formula = NA)
  
  expect_equal(round(test_gparms$Estimate[1], 2), 12.86, tolerance = 0.01)
  expect_equal(round(test_gparms$Estimate[2], 2), 6.38,  tolerance = 0.01)
  # 6.45 -> 6.38 to68fix
  
  
  # Also test get_growthparameters vs growthparameters
  
  test_gparms_gp <- growthparameters(test_fit, re_formula = NA)
  test_gparms_gp[['Est.Error']] <- NULL
  expect_equal(test_gparms_gp, test_gparms, 
               tolerance = 0.01)
  
  
  ##############################################################################
  # get_growthparameters
  ##############################################################################
  
  test_gparms_re <- get_growthparameters(test_fit, by = 'id', re_formula = NULL)
  test_gparms_gp_re <- growthparameters(test_fit, re_formula = NULL)
  # This will test mean of APGV and PGV
  expect_equal(mean(test_gparms_re$Estimate), 
               mean(test_gparms_gp_re$Estimate), 
               tolerance = 0.01)
  
  ##############################################################################
  # modelbased_growthparameters_call.bgmfit
  ##############################################################################
  
  mgc <- modelbased_growthparameters_call.bgmfit(test_fit, method = 'custom')
  expect_equal(mean(mgc$distance$Estimate), 158.2923, tolerance = 0.01)
  expect_equal(mean(mgc$velocity$Estimate), 6.4, tolerance = 0.01)
  
  
  ##############################################################################
  # optimize_model
  ##############################################################################
  omc <- optimize_model(test_fit, 
                        newdata = NULL,
                        optimize_df = 3,
                        optimize_x = list(NULL, log, sqrt),
                        optimize_y = list(NULL, log, sqrt),
                        transform_prior_class = c('beta', 'sd', 'rsd', 'sigma', 'dpar'),
                        transform_beta_coef = c('b', 'c'),
                        transform_sd_coef = c('b', 'c'),
                        exclude_default = TRUE,
                        add_fit_criteria = c("loo", "waic", "bayes_R2"),
                        byresp = FALSE,
                        model_name = NULL,
                        overwrite = FALSE,
                        file = NULL,
                        force_save = FALSE,
                        save_each = FALSE,
                        digits = 2,
                        cores = 1,
                        verbose = FALSE,
                        expose_function = FALSE,
                        usesavedfuns = FALSE,
                        clearenvfuns = NULL,
                        envir = NULL)
  
  expect_true(is.data.frame(omc$optimize_summary))
  
  
  
  ##############################################################################
  # hypothesis_test
  ##############################################################################
  
  
  ##############################################################################
  # summary_table
  ##############################################################################
  # stab <- summary_table(test_fit)
  
  
})



