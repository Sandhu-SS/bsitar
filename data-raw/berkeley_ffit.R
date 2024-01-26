## code to prepare `berkeley_fit` dataset goes here

library(bsitar)

berkeley_fdata <- na.omit(berkeley[berkeley$sex == "Female" &
                           berkeley$age >= 8 & berkeley$age <= 18,
                         c('id', 'age', 'height')])


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

berkeley_fdata <- berkeley_fdata %>% sample_n_of_groups(size = 20, id)


sitar_ffit <- sitar::sitar(x = age, y = height, id = id, data = berkeley_fdata, df = 5)


berkeley_ffit <- bsitar(x = age, y = height, id = id, data = berkeley_fdata,
                     df = 5, xoffset = mean, fixed = a+b+c, random = a+b+c,
                     a_formula = ~1, b_formula = ~1, c_formula = ~1,
                     threads = brms::threading(NULL),
                     sample_prior = "no",
                     chains = 2, cores = 2, iter = 5000, thin = 15)

usethis::use_data(berkeley_fdata, overwrite = TRUE)

usethis::use_data(berkeley_ffit, overwrite = TRUE)


