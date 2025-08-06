

###############################################################################
# Create test data
# he data will allow testing all three forms of models
#   univariate model
#   univariate_by model
#   multivariate model
###############################################################################

# set.seed(123)
# test_data <- data.frame(x = rep(1:5, 4),
#                   y = rep(seq(100, 150, length.out = 5), 2)+rnorm(20, 0, 0.5),
#                   y2 = rep(seq(100, 150, length.out = 5), 2)+rnorm(20, 0, 0.25),
#                   id = factor(rep(c("1", "2", "3", "4"),times=c(5,5,5,5))),
#                   gender = factor(rep(c("Male", "Female"),times=c(10,10))))


  test_data <- bsitar::berkeley %>% 
    dplyr::select(id, age, height, weight, sex) %>% 
    dplyr::relocate(height, .before = weight) %>% 
    dplyr::relocate(sex, .before = id) %>% 
    dplyr::filter(age <= 20)
  
  if(nlevels(test_data$sex) != 2) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(nlevels(test_data$id) != 136) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(round(mean(test_data$age), 2) != 13.36) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(round(mean(test_data$height), 2) != 157.56) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  if(round(mean(test_data$weight), 2) != 50.16) {
    stop("data 'bsitar::berkeley' used as test data has changed")
  }
  
  test_data_male   <- test_data %>% dplyr::filter(sex == 'Male') %>% 
    droplevels()
  test_data_female <- test_data %>% dplyr::filter(sex == 'Female') %>% 
    droplevels()


