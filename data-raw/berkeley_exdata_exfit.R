


# library(bsitar)

# devtools::load_all()

# berkeley_exdata

# sitar::heights

berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                         df = 3, chains = 2, cores = 2, iter = 1000, thin = 15,
                         a_prior_beta = normal(ymean, ysd, autoscale = TRUE),
                         b_prior_beta = normal(0, 2, autoscale = FALSE),
                         c_prior_beta = normal(0, 1, autoscale = FALSE),
                         d_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         b_prior_sd = normal(0, 2, autoscale = FALSE),
                         c_prior_sd = normal(0, 1, autoscale = FALSE),
                         d_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         s_prior_beta = normal(lm, lm, autoscale = FALSE),
                         control = list(adapt_delta = 0.9, max_treedepth = 15),
                         save_pars = save_pars(all = TRUE),
                         # sample_prior = "only",
                         # init = 0,
                         # init_r = 0,
                         sample_prior = "yes",
                         normalize = TRUE, seed = 123)


 berkeley_exfit$test_mode <- TRUE
 
 # priorsense::powerscale_sensitivity(berkeley_exfit)
 
# object.size(berkeley_exfit)
# 1563920 bytes

usethis::use_data(berkeley_exdata, overwrite = TRUE, compress = 'xz')
usethis::use_data(berkeley_exfit,  overwrite = TRUE, compress = 'xz')

