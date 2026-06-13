


# library(bsitar)

# devtools::load_all()

# berkeley_exdata

# sitar::heights

berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                         df = 3, chains = 2, cores = 2, iter = 1000, thin = 15,
                         
                         
                         
                         a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                         b_prior_beta = normal(0, 2, autoscale = FALSE),
                         c_prior_beta = normal(0, 1, autoscale = FALSE),
                         d_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         s_prior_beta = normal(lm, lm, autoscale = FALSE),
                         a_cov_prior_beta = normal(0, 5.0, autoscale = FALSE),
                         b_cov_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                         d_cov_prior_beta = normal(0, 1.0, autoscale = FALSE),
                         s_cov_prior_beta = normal(lm, lm, autoscale = FALSE),
                         a_prior_sd = normal(0, ysd, autoscale = FALSE),
                         b_prior_sd = normal(0, 2, autoscale = FALSE),
                         c_prior_sd = normal(0, 1, autoscale = FALSE),
                         d_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         a_cov_prior_sd = normal(0, 5.0, autoscale = FALSE),
                         b_cov_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         c_cov_prior_sd = normal(0, 0.1, autoscale = FALSE),
                         d_cov_prior_sd = normal(0, 1.0, autoscale = FALSE),
                         a_prior_sd_str = NULL,
                         b_prior_sd_str = NULL,
                         c_prior_sd_str = NULL,
                         d_prior_sd_str = NULL,
                         a_cov_prior_sd_str = NULL,
                         b_cov_prior_sd_str = NULL,
                         c_cov_prior_sd_str = NULL,
                         d_cov_prior_sd_str = NULL,
                         sigma_prior_beta = normal(0, 1, autoscale = FALSE),
                         sigma_cov_prior_beta = normal(0, 0.5, autoscale = FALSE),
                         sigma_prior_sd = normal(0, 0.25, autoscale = FALSE),
                         sigma_cov_prior_sd = normal(0, 0.15, autoscale = FALSE),
                         sigma_prior_sd_str = NULL,
                         sigma_cov_prior_sd_str = NULL,
                         rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         dpar_prior_sigma = normal(0, ysd, autoscale = FALSE),
                         dpar_cov_prior_sigma = normal(0, 1, autoscale = FALSE),
                         autocor_prior_acor = uniform(-1, 1, autoscale = FALSE),
                         autocor_prior_unstr_acor = lkj(1),
                         gr_prior_cor = lkj(1),
                         gr_prior_cor_str = lkj(1),
                         sigma_prior_cor = lkj(1),
                         sigma_prior_cor_str = lkj(1),
                         mvr_prior_rescor = lkj(1),
                         init = NULL,
                         init_r = 0.5,
                         a_init_beta = 'lm',
                         b_init_beta = 0,
                         c_init_beta = 0,
                         d_init_beta = 0,
                         s_init_beta = 'lm',
                         a_cov_init_beta = 0,
                         b_cov_init_beta = 0,
                         c_cov_init_beta = 0,
                         d_cov_init_beta = 0,
                         s_cov_init_beta = 0,
                         a_init_sd = 'random',
                         b_init_sd = 'random',
                         c_init_sd = 'random',
                         d_init_sd = 'random',
                         
                         control = list(adapt_delta = 0.95, max_treedepth = 15),
                         
                         
                         
                         
                         
                         save_pars = save_pars(all = TRUE),
                         # sample_prior = "only",
                         # init = random,
                         # init_r = 0.1,
                         sample_prior = "yes",
                         normalize = TRUE, seed = 123)


 berkeley_exfit$test_mode <- TRUE
 
 # priorsense::powerscale_sensitivity(berkeley_exfit)
 
# object.size(berkeley_exfit)
# 1563920 bytes

usethis::use_data(berkeley_exdata, overwrite = TRUE, compress = 'xz')
usethis::use_data(berkeley_exfit,  overwrite = TRUE, compress = 'xz')





devtools::load_all()

fit_3.3 <- bsitar(x = age, y = height, id = id, data = berkeley_exdata, 
       df = 3, stype = "rcs", 
       a_prior_beta = normal(ymean, ysd, autoscale = 2.5), 
       b_prior_beta = normal(0, 2, autoscale = FALSE), 
       c_prior_beta = normal(0, 1, autoscale = FALSE), 
       s_prior_beta = normal(lm, lm, autoscale = FALSE), 
       a_prior_sd = normal(0, ysd, autoscale = FALSE), 
       b_prior_sd = normal(0, 2, autoscale = FALSE), 
       c_prior_sd = normal(0, 1, autoscale = FALSE), 
       rsd_prior_sigma = normal(0, ysd, autoscale = FALSE), 
       init = NULL,  init_r = 0.0, 
       chains = 2, iter = 1000, cores = 2, thin = 12,
       sample_prior = "yes", save_pars = save_pars(all = TRUE),
       parameterization = "ncp", 
      #  threads = 2,
      # vcov_init_0 = FALSE,
       threads = list(threads = 1, grainsize = NULL, static = FALSE, force = FALSE),
       seed = 123)


usethis::use_data(fit_3.3,  overwrite = TRUE, compress = 'xz')



