## code to prepare `heights_fit` dataset goes here


library(magrittr)
data(heights, package = "sitar")
data <- heights
rm(heights)


data <- data %>% 
  dplyr::select(id, age, height, men) %>% 
  tidyr::drop_na(height) %>% 
  dplyr::rename(y = height)

# data %>% dplyr::glimpse()

heights_fit <- bsitar::bgm(x = age, y = y, id = id, data = data, 
                               chains = 2, iter = 100,
                               silent = 0,
                               sample_prior = 'only')


usethis::use_data(heights_fit, overwrite = TRUE)


