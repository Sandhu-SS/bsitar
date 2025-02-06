## code to prepare `berkeley_fit` dataset goes here

library(bsitar)


berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                        df = 3, 
                        stype = list('nsp', F),
                        chains = 4, cores = 4, iter = 2000, thin = 6,
                        a_prior_beta = normal(lm, ysd, autoscale = 2.5),
                        b_prior_beta = normal(0, 1.5),
                        c_prior_beta = normal(0, 0.5),
                        s_prior_beta = normal(lm, lm, autoscale = 2.5),
                        # a_prior_beta = flat,
                        # b_prior_beta = flat,
                        # c_prior_beta = flat,
                        # s_prior_beta = flat,
                        a_prior_sd = normal(0, ysd, autoscale = 1),
                        b_prior_sd = normal(0, 1.0),
                        c_prior_sd = normal(0, 0.25),
                        rsd_prior_sigma = normal(0, ysd, autoscale = 1),
                        a_init_beta = lm,
                        b_init_beta = 0,
                        c_init_beta = 0,
                        s_init_beta = lm,
                        backend = "rstan", 
                        normalize = FALSE, 
                        threads = brms::threading(NULL),
                        seed = 123)


usethis::use_data(berkeley_exdata, overwrite = TRUE, compress = 'xz')
usethis::use_data(berkeley_exfit,  overwrite = TRUE, compress = 'xz')

