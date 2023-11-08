## code to prepare `berkeley_fit` dataset goes here

# library(magrittr)

library(bsitar)

data(berkeley, package = "sitar")
data <- berkeley
rm(berkeley)

berkeley <- data %>% 
  dplyr::mutate(sex = 
                  dplyr::recode_factor(sex, "1" = "Male", "2" = "Female")) %>% 
  tidyr::drop_na() %>% 
  droplevels()

berkeley %>% dplyr::glimpse()

usethis::use_data(berkeley, overwrite = TRUE)
