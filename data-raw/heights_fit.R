## code to prepare `heights_fit` dataset goes here


library(magrittr)
data(heights, package = "sitar")
data <- heights
rm(heights)


data <- data %>% 
  dplyr::select(id, age, height, men) %>% 
  tidyr::drop_na(height) %>% 
  dplyr::rename(y = height)

# data %>% dplyr::glimpse()

# sfit <- sitar::sitar(x = age, y = y, id = id, data = data, df = 4, random = 'a+b+c') 

heights_fit <- bsitar::bgm(x = age, y = y, id = id, data = data, df = 4,
                           chains = 2, iter = 2000,
                           select_model = sitar,
                           backend = 'cmdstanr',
                           thin = 2)

# heights_fit <- bsitar::bgm(x = age, y = y, id = id, data = data, 
#                            chains = 2, iter = 1000,
#                            select_model = sitar,
#                            a_prior_beta = normal(ymean, ysd, autoscale = 2.5),
#                            b_prior_beta = normal(0, 2),
#                            c_prior_beta = normal(0, 0.25),
#                            s_prior_beta = normal(0, lm, autoscale = 2.5),
#                            a_prior_sd = normal(0, ysd, autoscale = 2.5),
#                            b_prior_sd = normal(0, 2),
#                            c_prior_sd = normal(0, 0.15),
#                            d_prior_sd = normal(0, 1),
#                            rsd_prior_sigma = normal(0, ysd, autoscale = 2.5),
#                            
#                            threads = brms::threading(44),
#                            control = list(max_treedepth = 15),
#                            backend = 'cmdstanr',
#                            
#                            sample_prior = 'no')


usethis::use_data(heights_fit, overwrite = TRUE)


