
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  #  skip_local_run_ci()
}

if(set_skip_run_ci) {
  # skip_run_ci()
}


###############################################################################
# Test marginals vs marginaleffects
###############################################################################

test_that("test-hypothesis_test", {
  skip_on_cran()

  
  testthat::test_that("prior_sensitivity validates plot argument", {
    fake_fit <- structure(list(), class = "bgmfit")
    
    testthat::skip_if_not_installed("priorsense")
    testthat::skip_if_not_installed("brms")
    
    # expect_error(
    #   prior_sensitivity(fake_fit, plot = "bad-plot-type"),
    #   NA
    # )
  })
  
  
  
  
  testthat::test_that("prior_sensitivity_conflict returns expected names", {
    x <- list(
      sensitivity = data.frame(
        variable = c("a", "a", "b", "b"),
        component = c("prior", "likelihood", "prior", "likelihood"),
        prior = c(0.10, 0.10, 0.10, 0.01),
        likelihood = c(0.10, 0.10, 0.10, 0.01),
        sensitivity = c(0.10, 0.10, 0.10, 0.01)
      )
    )
    class(x) <- c("prior_sensitivity", "bgmfit")
    
    out <- prior_sensitivity_conflict(x)
    
    testthat::expect_named(
      out,
      c(
        "prior_flagged",
        "likelihood_flagged",
        "prior_likelihood_flagged"
      )
    )
    testthat::expect_true("a" %in% out$prior_flagged)
    testthat::expect_true("b" %in% out$likelihood_flagged)
  })
  
  
  
  
  
  
  
}) # test_that("test-marginals-slopes-bycov-byvariable", {






