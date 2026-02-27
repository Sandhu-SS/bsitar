## code to prepare `berkeley_fit` dataset goes here

#library(bsitar)

# devtools::load_all()

berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                        df = 3, 
                       #  stype = list('nsk', normalize = T),
                        
                       # stype = list('nsp', normalize = FALSE),
                        
                        chains = 2, cores = 2, iter = 1000, thin = 6,
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
                        # backend = "cmdstanr", 
                       
                       ######################################
                       # MCMC Setting - See SECTION 1: SETUP
                       ######################################
                       
                       a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                       b_prior_beta = normal(0, 2, autoscale = FALSE),
                       c_prior_beta = normal(0, 1, autoscale = FALSE),
                       
                       s_prior_beta = normal(lm, lm, autoscale = FALSE),
                       
                       # a_cov_prior_beta = normal(0, .10, autoscale = FALSE),
                       # b_cov_prior_beta = normal(0, .2, autoscale = FALSE),
                       # c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                       
                       a_prior_sd = normal(0, ysd, autoscale = FALSE),
                       b_prior_sd = normal(0, 2, autoscale = FALSE),
                       c_prior_sd = normal(0, 1, autoscale = FALSE),
                       
                       rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                       
                       # cp -> multi_normal_cholesky multi_normal, default cp is multi_normal_cholesky
                       parameterization = "ncp",        # Suppress verbose output,
                       stype = "rcs",        # Suppress verbose output,
                       init = 0.5,
                       
                       ######################################
                       
                       
                        sample_prior = "no",
                        normalize = FALSE, 
                        genquant_xyadj = FALSE,
                        # threads = brms::threading(NULL),
                        seed = 123)


 berkeley_exfit$test_mode <- TRUE
 
#  object.size(berkeley_exfit)
# 1563920 bytes

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

