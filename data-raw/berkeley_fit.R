## code to prepare `berkeley_fit` dataset goes here

library(magrittr)
data(berkeley, package = "sitar")
data <- berkeley
rm(berkeley)

data <- data %>% 
  dplyr::select(id, age, height, sex) %>% 
  dplyr::filter(age %in% c(6:20) ) %>% 
  tidyr::drop_na(height)

data$sex <- dplyr::recode_factor(data$sex, "1" = "Male", "2" = "Female")

data <- data %>% 
  sjmisc::to_dummy(sex, suffix = "label") %>% 
  dplyr::bind_cols(data,.) %>% 
  dplyr::mutate(Male = height, Female = height) %>% 
  dplyr::select(-height) %>% 
  dplyr::arrange(id, age, sex)


data <- data %>% dplyr::mutate(y = Male, id = as.factor(id) ) 
data <- data %>% dplyr::select(id, age, sex, y)

data <- data %>% dplyr::filter(sex == "Male")

# data %>% dplyr::glimpse()

berkeley_fit <- bsitar::bsitar(age, y, id, data=data, 
                       chains = 2, iter = 100,
                       sample_prior = 'only')


usethis::use_data(berkeley_fit, overwrite = TRUE)
