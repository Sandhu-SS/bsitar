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

berkeley_mdata %>% dplyr::glimpse()

usethis::use_data(berkeley_mdata, overwrite = TRUE)
