## code to prepare `berkeley_fit` dataset goes here

library(bsitar)

data(berkeley, package = "sitar")
berkeley_data <- berkeley
rm(berkeley)

berkeley_exdata <- berkeley_data %>%
  dplyr::select(id, age, height, sex) %>%
  dplyr::filter(age %in% c(8:18) ) %>%
  tidyr::drop_na(height) %>%
  dplyr::mutate(sex =
                  dplyr::recode_factor(sex, "1" = "Male", "2" = "Female")) %>%
  dplyr::select(id, age, sex, height) %>%
  tidyr::drop_na() %>%
  dplyr::filter(sex == "Female") %>%
  dplyr::select(-sex) %>%
  droplevels() %>% data.frame()




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
                        df = 5, 
                        threads = brms::threading(NULL), 
                        
                        # backend = 'cmdstanr',
                        
                        a_prior_beta = normal(ymean, ysd, autoscale = TRUE),
                        b_prior_beta = normal(0, 3.5),
                        c_prior_beta = normal(0, 1.5),
                        
                        a_prior_sd = normal(0, ysd),
                        b_prior_sd = normal(0, 2.0),
                        c_prior_sd = normal(0, 1.0),
                        
                        rsd_prior_sigma = normal(0, ysd),
                        
                        sample_prior = 'no',
                        expose_function = FALSE,
                         # chains = 2, cores = 2, iter = 2000, thin = 15,
                        chains = 2, cores = 2, iter = 1000, thin = 7,
                        seed = 123)



# https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/

# breakItDown = function(mod) {
#   sapply(mod, FUN=function(x){length(serialize(x, NULL))}, simplify=T)
# }
# cleanModel1 = function(cm) {
#   cm$model_info$bgmfit.data = c()
#   cm$bmodel  = c()
#   cm$model_info$emodel  = c()
#   cm$model_info$exefuns  = c()
#   cm$model_info$bgmfit.data   = c()
#   cm
# }
# 
# berkeley_exfit_checks = breakItDown(berkeley_exfit)
# berkeley_exfit_reduce = cleanModel1(berkeley_exfit)
# 
# berkeley_exfit <- berkeley_exfit_reduce




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




