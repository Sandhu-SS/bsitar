## code to prepare `berkeley_fit` dataset goes here

# library(magrittr)

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
  droplevels()


sitar_fit <- sitar::sitar(x = age, y = height, id = id, df = 4,
                          data = berkeley_mdata)

berkeley_mfit <- bsitar::bgm(x = age, y = height, id = id, df = 4,
                            data = berkeley_mdata,
                            sample_prior = 'only',
                            chains = 2, iter = 1000, thin = 10)


usethis::use_data(berkeley_mfit, overwrite = TRUE)
