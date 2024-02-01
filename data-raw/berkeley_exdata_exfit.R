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




sample_n_of_groups <- function(data, size, ...) {
  dots <- rlang::quos(...)
  group_ids <- data %>%
    dplyr::group_by(!!! dots) %>%
    dplyr::group_indices()
  sampled_groups <- sample(unique(group_ids), size)
  data %>%
    dplyr::filter(group_ids %in% sampled_groups) %>%
    droplevels()
}

set.seed(1234)
#berkeley_exdata <- berkeley_exdata %>% sample_n_of_groups(size = 20, id)


# sitar_fit <- sitar::sitar(x = age, y = height, id = id, df = 5,
#                           data = berkeley_exdata)


berkeley_exfit <- bsitar(x = age, y = height, id = id, data = berkeley_exdata,
                        df = 5, 
                        threads = brms::threading(NULL), 
                        # backend = 'cmdstanr',
                        # b_prior_beta = student_t(3, 0, 2.5),
                        # c_prior_beta = student_t(3, 0, 1.0),
                        sample_prior = 'no',
                        expose_function = FALSE,
                        chains = 2, cores = 2, iter = 2000, thin = 10)


usethis::use_data(berkeley_exdata, overwrite = TRUE)

save_file       <- "berkeley_exfittemp.rds"
saveRDS(berkeley_exfit, file = save_file, compress = 'xz')

rm(list = ls())

save_file      <- "berkeley_exfittemp.rds"
berkeley_exfit <- readRDS(file = save_file)

file.remove(save_file)

usethis::use_data(berkeley_exfit, overwrite = TRUE)


