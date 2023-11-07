## code to prepare `berkeley_fit` dataset goes here

# library(magrittr)

library(bsitar)

data(berkeley, package = "sitar")
data <- berkeley
rm(berkeley)

data <- data %>% 
  dplyr::select(id, age, height, sex) %>% 
  dplyr::filter(age %in% c(6:20) ) %>% 
  tidyr::drop_na(height) %>% 
  dplyr::mutate(sex = 
                  dplyr::recode_factor(sex, "1" = "Male", "2" = "Female")) %>% 
  dplyr::select(id, age, sex, height) %>% 
  tidyr::drop_na() %>% 
  dplyr::filter(sex == "Male")

data %>% dplyr::glimpse()

sitar_fit <- sitar::sitar(x = age, y = height, id = id, data = data, df = 4)

berkeley_fit <- bsitar::bgm(x = age, y = height, id = id, data = data, df = 4,
                            sample_prior = 'only',
                            chains = 2, iter = 1000, thin = 8)


# usethis::use_data(berkeley_fit, overwrite = TRUE)
