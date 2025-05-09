## code to prepare `berkeley_fit` dataset goes here

#library(bsitar)

devtools::load_all()

berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                        df = 3, 
                        stype = list('nsp', F),
                        chains = 2, cores = 2, iter = 2000, thin = 5,
                        # a_prior_beta = normal(lm, ysd, autoscale = 2.5),
                        # b_prior_beta = normal(0, 1.5),
                        # c_prior_beta = normal(0, 0.5),
                        # s_prior_beta = normal(lm, lm, autoscale = 2.5),
                        # # a_prior_beta = flat,
                        # # b_prior_beta = flat,
                        # # c_prior_beta = flat,
                        # # s_prior_beta = flat,
                        # a_prior_sd = normal(0, ysd, autoscale = 1),
                        # b_prior_sd = normal(0, 1.0),
                        # c_prior_sd = normal(0, 0.25),
                        # rsd_prior_sigma = normal(0, ysd, autoscale = 1),
                        # a_init_beta = lm,
                        # b_init_beta = 0,
                        # c_init_beta = 0,
                        # s_init_beta = lm,
                        # backend = "rstan", 
                        sample_prior = "no",
                        normalize = FALSE, 
                        genquant_xyadj = FALSE,
                        # threads = brms::threading(NULL),
                        seed = 123)


# load("C:/Users/drsat/OneDrive/Documents/GitHub/bsitar/data/berkeley_exfit.rda")

# save_file_exdata      <- "berkeley_exdata_temp.rds"
# save_file_exfit       <- "berkeley_exfit_temp.rds"
# 
# saveRDS(berkeley_exdata, file = save_file_exdata, compress = 'xz')
# saveRDS(berkeley_exfit,  file = save_file_exfit,  compress = 'xz')
# 
# 
# rm(list=setdiff(ls(), c('save_file_exdata', 'save_file_exfit')))
# 
# 
# berkeley_exdata <- readRDS(file = save_file_exdata)
# berkeley_exfit  <- readRDS(file = save_file_exfit)
# 
# file.remove(save_file_exdata)
# file.remove(save_file_exfit)


# Moving from data 'rda' to sysdata internal

# usethis::use_data(berkeley_exdata, overwrite = TRUE)
# usethis::use_data(berkeley_exfit, overwrite = TRUE)



# for (i in names(berkeley_exfit$model_info)) {
#   print(i)
#   object.size(berkeley_exfit$model_info[[i]]) %>% print()
# }



usethis::use_data(berkeley_exdata, overwrite = TRUE, compress = 'xz')
usethis::use_data(berkeley_exfit,  overwrite = TRUE, compress = 'xz')

