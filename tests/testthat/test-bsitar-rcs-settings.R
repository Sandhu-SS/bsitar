

# Skip test for local R CMD Check but run on GitHub

# skip_if_not_ci()


###############################################################################
# Test bsitar with rcs settings
###############################################################################

test_that("bsitar works fully with nsk settings", {
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
                       vcov_init_0 = TRUE,
                       refres = 0, silent = 2,
                       seed = 123)
  }))
  
  
  # test_fit <- test_fit_rcs
  
  true_sbetas <- c(128.05, 0.00, 0.00, 5.02, 4.28, -12.27, -13.26)
  
  test_sbetas <- round(unname(brms::fixef(test_fit)[,1]), 2)
  
  expect_equal(true_sbetas, test_sbetas, tolerance = 0.01)
  
  test_gparms <- marginal_growthparameters(test_fit, re_formula = NA)
  
  expect_equal(round(test_gparms$Estimate[1], 2), 12.86, tolerance = 0.01)
  expect_equal(round(test_gparms$Estimate[2], 2), 6.470,  tolerance = 0.01)
  
  # marginal_draws(test_fit, re_formula = NA, deriv = 0, by = 'age', plot = T)
  # marginal_draws(test_fit, re_formula = NA, deriv = 1, by = 'age', plot = T)
  
  # expect_error(bsitar(x=xx, y=y, id=id, data = dat, backend = "rstan",
  #                    get_stancode = TRUE, sample_prior = "only"), 
  #              "variable xx not in the dataframe")
  
})



