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
                        stype = 'nsp',
                        sample_prior = 'no', 
                        chains = 2, cores = 2, iter = 2000, thin = 10,
                        a_prior_beta = normal(ymean, ysd, autoscale = 1),
                        b_prior_beta = normal(0, 1.0, autoscale = FALSE),
                        c_prior_beta = normal(0, 0.5, autoscale = FALSE),
                        s_prior_beta = normal(lm, lm, autoscale = 1),
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
                        init_r = 0.25,
                        backend = "rstan", 
                        normalize = FALSE, 
                        threads = threading(NULL),
                        seed = 123)

# berkeley_exfitx <- berkeley_exfit

# berkeley_exfitxx <- berkeley_exfit

# berkeley_exfitxxx <- berkeley_exfit

# plot(berkeley_exfit, nvariables = 1, combo = c('dens', 'trace'))

save_file_exdata      <- "berkeley_exdata_temp.rds"
save_file_exfit       <- "berkeley_exfit_temp.rds"

saveRDS(berkeley_exdata, file = save_file_exdata, compress = 'xz')

saveRDS(berkeley_exfit,  file = save_file_exfit,  compress = 'xz')

# saveRDS(reduce_size(berkeley_exfit), file = save_file_exfit,  compress = 'xz')


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



# https://stackoverflow.com/questions/75924112/r-brms-saving-models-to-file-inside-a-function-call-saves-the-entire-local-env

# hack_size <- function(x, ...) {
#   UseMethod("hack_size")
# }
# hack_size.stanfit <- function(x) {
#   x@stanmodel <- structure(numeric(0), class="stanmodel")
#   x@.MISC <- new.env()
#   return(x)
# }
# 
# hack_size.brmsfit <- function(x) {
#   x$fit <- hack_size(x$fit)
#   return(x)
# }
# 
# hack_size.stanreg <- function(x) {
#   x$stanfit <- hack_size(x$stanfit)
#   return(x)
# }
# 
# 
# hack_size.data.frame <- function(x) {
#   environment(attr(x, "terms")) <- new.env(parent = baseenv())
#   return(x)
# }

# hack_size.brmsfit <- function(x, formual  = FALSE, data = TRUE, fit = TRUE) {
#   
#   hack_size.data.frame <- function(x) {
#     environment(attr(x, "terms")) <- new.env(parent = baseenv())
#     return(x)
#   }
#   
#   hack_size.brmsfit <- function(x) {
#     x$fit <- hack_size(x$fit)
#     return(x)
#   }
#   
#   if(formual) x$formula <- hack_size(x$formula)
#   if(data)    x$data    <- hack_size(x$data)
#   if(fit)     x$fit     <- hack_size(x$fit)
#   return(x)
# }
# 
# hack_size <- function(x, ...) {
#   UseMethod("hack_size")
# }


