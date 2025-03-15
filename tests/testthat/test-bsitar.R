test_that("bsitar works fully with default backend", {
  skip_on_cran()
  dat <- data.frame(x = rep(1:5, 2),
                    y = rep(seq(100, 150, length.out = 5), 2)+rnorm(10, 0, 0.5),
                    id = factor(rep(c("a", "b"),times=c(5,5))))
  
  expect_type(bsitar(x=x, y=y, id=id, data = dat, backend = "rstan",
                     threads = threading(NULL),
                   get_stancode = TRUE, sample_prior = "only"),
               "character")
  
  expect_type(bsitar(x=x, y=y, id=id, data = dat, backend = "rstan",
                     threads = threading(NULL),
                     get_standata = TRUE, sample_prior = "only"),
              "list")
  
  # expect_error(bsitar(x=xx, y=y, id=id, data = dat, backend = "rstan",
  #                    get_stancode = TRUE, sample_prior = "only"), 
  #              "variable xx not in the dataframe")
  
})
