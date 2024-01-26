## code to prepare `berkeley_fit` dataset goes here

library(bsitar)

data(berkeley, package = "sitar")
data <- berkeley
rm(berkeley)

berkeley_mdata <- data %>%
  dplyr::select(id, age, height, sex) %>%
  dplyr::filter(age %in% c(6:20) ) %>%
  tidyr::drop_na(height) %>%
  dplyr::mutate(sex =
                  dplyr::recode_factor(sex, "1" = "Male", "2" = "Female")) %>%
  dplyr::select(id, age, sex, height) %>%
  tidyr::drop_na() %>%
  dplyr::filter(sex == "Male") %>%
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
berkeley_mdata <- berkeley_mdata %>% sample_n_of_groups(size = 20, id)


sitar_fit <- sitar::sitar(x = age, y = height, id = id, df = 4,
                          data = berkeley_mdata)


berkeley_mfit <- bsitar(x = age, y = height, id = id, data = berkeley_mdata,
                        df = 4, xoffset = mean, fixed = a+b+c, random = a+b+c,
                        a_formula = ~1, b_formula = ~1, c_formula = ~1,
                        threads = brms::threading(NULL),
                        chains = 2, cores = 2, iter = 5000, thin = 15)





usethis::use_data(berkeley_mdata, overwrite = TRUE)

usethis::use_data(berkeley_mfit, overwrite = TRUE)


