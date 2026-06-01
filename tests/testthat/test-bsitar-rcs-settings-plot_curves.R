
# Skip test for local R CMD Check but run on GitHub

if(skip_test_local_rcmd_check) {
  skip_local_run_ci()
}


###############################################################################
# Test plot_curves with rcs settings
###############################################################################

test_that("test plot_curves", {
  skip_on_cran()
 
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
  
  
  expected <- list(
    d = 161.3626,
    v = 4.366173,
    D = 160.1998,
    V = 4.479366,
    a = 160.5832
  )
  
  # plots_c <- c("d", "v", "D", "V", "a", "u", "O", "Dd", "Ddv", "DVdv", "dvDVauO")
  
  plots_c <- c("d", "v", "D", "V", "a", "u", "O", "dvDVauO")
  plots_list <- list()
  plots_est_list <- list()
  for (plots_ci in plots_c) {
    plots_list[[plots_ci]] <- 
      plot_curves(test_fit, plots_ci, bands = '', print = FALSE)
    if(!inherits(plots_list[[plots_ci]], "patchwork")) {
      if(plots_ci == "d" | plots_ci == "v" |
         plots_ci == "D" | plots_ci == "V" |
         plots_ci == "a") {
        plots_est_list[[plots_ci]] <- mean(plots_list[[plots_ci]]$data$Estimate,
                                           na.rm = TRUE)
      }
    }
  }
  
  for (i in seq_along(expected)) {
    expect_equal(plots_est_list[[i]], expected[[i]], tolerance = 0.01)
  }
  
  
})



