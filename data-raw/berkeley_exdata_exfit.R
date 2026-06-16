


# library(bsitar)

# devtools::load_all()

# berkeley_exdata

# sitar::heights

berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                         df = 3, chains = 2, cores = 2, iter = 2000, 
                         thin = 10, warmup = 1000,
                         a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                         b_prior_beta = normal(0, 2, autoscale = FALSE),
                         c_prior_beta = normal(0, 1, autoscale = FALSE),
                         d_prior_beta = normal(0, 1, autoscale = FALSE),
                         a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         b_prior_sd = normal(0, 2, autoscale = FALSE),
                         c_prior_sd = normal(0, 1, autoscale = FALSE),
                         d_prior_sd = normal(0, 1, autoscale = FALSE),
                         rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         s_prior_beta = normal(lm, lm, autoscale = 1),
                         control = list(adapt_delta = 0.95, max_treedepth = 12),
                         save_pars = save_pars(all = TRUE),
                         sample_prior = "yes",
                         normalize = TRUE, 
                         
                         init = NULL,
                         
                         seed = 123)




 berkeley_exfit$test_mode <- TRUE
 
 # priorsense::powerscale_sensitivity(berkeley_exfit)
 
# object.size(berkeley_exfit)
# 1563920 bytes

usethis::use_data(berkeley_exdata, overwrite = TRUE, compress = 'xz')
usethis::use_data(berkeley_exfit,  overwrite = TRUE, compress = 'xz')




# 
# devtools::load_all()
# 
# fit_3.3 <- bsitar(x = age, y = height, id = id, data = berkeley_exdata, 
#        df = 3, stype = "rcs", 
#        a_prior_beta = normal(ymean, ysd, autoscale = 2.5), 
#        b_prior_beta = normal(0, 2, autoscale = FALSE), 
#        c_prior_beta = normal(0, 1, autoscale = FALSE), 
#        s_prior_beta = normal(lm, lm, autoscale = FALSE), 
#        a_prior_sd = normal(0, ysd, autoscale = FALSE), 
#        b_prior_sd = normal(0, 2, autoscale = FALSE), 
#        c_prior_sd = normal(0, 1, autoscale = FALSE), 
#        rsd_prior_sigma = normal(0, ysd, autoscale = FALSE), 
#        init = NULL,  init_r = 0.0, 
#        chains = 2, iter = 1000, cores = 2, thin = 12,
#        sample_prior = "yes", save_pars = save_pars(all = TRUE),
#        parameterization = "ncp", 
#       #  threads = 2,
#       # vcov_init_0 = FALSE,
#        threads = list(threads = 1, grainsize = NULL, static = FALSE, force = FALSE),
#        seed = 123)
# 
# 
# usethis::use_data(fit_3.3,  overwrite = TRUE, compress = 'xz')
# 
# 

