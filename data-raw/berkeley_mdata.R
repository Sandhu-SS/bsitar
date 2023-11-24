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

berkeley_mdata %>% dplyr::glimpse()

usethis::use_data(berkeley_mdata, overwrite = TRUE)
