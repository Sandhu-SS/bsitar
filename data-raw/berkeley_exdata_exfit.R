## code to prepare `berkeley_fit` dataset goes here

library(bsitar)

# data(berkeley, package = "sitar")
# berkeley_data <- berkeley
# rm(berkeley)
# 
# berkeley_exdata <- berkeley_data %>%
#   dplyr::select(id, age, height, sex) %>%
#   dplyr::filter(age %in% c(8:18) ) %>%
#   tidyr::drop_na(height) %>%
#   dplyr::mutate(sex =
#                   dplyr::recode_factor(sex, "1" = "Male", "2" = "Female")) %>%
#   dplyr::select(id, age, sex, height) %>%
#   tidyr::drop_na() %>%
#   dplyr::filter(sex == "Female") %>%
#   dplyr::select(-sex) %>%
#   droplevels() %>% data.frame()




# sample_n_of_groups <- function(data, size, ...) {
#   dots <- rlang::quos(...)
#   group_ids <- data %>%
#     dplyr::group_by(!!! dots) %>%
#     dplyr::group_indices()
#   sampled_groups <- sample(unique(group_ids), size)
#   data %>%
#     dplyr::filter(group_ids %in% sampled_groups) %>%
#     droplevels()
# }
# 
# set.seed(1234)
# berkeley_exdata <- berkeley_exdata %>% sample_n_of_groups(size = 20, id)


# sitar_fit <- sitar::sitar(x = age, y = height, id = id, df = 5,
#                           data = berkeley_exdata)


berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                        df = 3,
                        stan_model_args = list(include_paths = "C:/Users/drsat/OneDrive/Documents/GitHub/bsitar/data-raw"),
                        chains = 2, cores = 2, iter = 1000, thin = 5,
                        a_prior_beta = normal(lm, ysd, autoscale = 2.5),
                        b_prior_beta = normal(0, 1.0),
                        c_prior_beta = normal(0, 0.25),
                        s_prior_beta = normal(lm, lm, autoscale = 2.5),
                        a_prior_sd = normal(0, ysd, autoscale = 2.5),
                        b_prior_sd = normal(0, 1.0),
                        c_prior_sd = normal(0, 0.25),
                        a_init_beta = lm, 
                        b_init_beta = 0,
                        c_init_beta = 0,
                        s_init_beta = lm,
                        backend = "rstan", 
                        normalize = FALSE, 
                        threads = brms::threading(NULL),
                        seed = 123)


# plot(berkeley_exfit, nvariables = 1, combo = c('dens', 'trace'))

save_file_exdata      <- "berkeley_exdata_temp.rds"
save_file_exfit       <- "berkeley_exfit_temp.rds"

saveRDS(berkeley_exdata, file = save_file_exdata, compress = 'xz')
saveRDS(berkeley_exfit,  file = save_file_exfit,  compress = 'xz')

rm(list=setdiff(ls(), c('save_file_exdata', 'save_file_exfit')))


berkeley_exdata <- readRDS(file = save_file_exdata)
berkeley_exfit  <- readRDS(file = save_file_exfit)

file.remove(save_file_exdata)
file.remove(save_file_exfit)


# Moving from data 'rda' to sysdata internal

usethis::use_data(berkeley_exdata, overwrite = TRUE)
usethis::use_data(berkeley_exfit, overwrite = TRUE)


# usethis::use_data(berkeley_exdata, berkeley_exfit, 
#                   overwrite = TRUE, internal = TRUE, compress = 'xz')




