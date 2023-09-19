


#' An internal function to evaluate and see priors specified
#'
#' @param model An object of class \code{bsitar}
#' @keywords internal
#' @return A data frame object. 
#' @noRd
#'
priors_to_textdata <- function(model,
                               spriors = NULL, 
                               sdata = NULL,
                               prior_name_asit = F, 
                               gsub_coef = NULL,
                               gsub_group = NULL,
                               sort_response = NULL,
                               # sort_parameter = NULL,
                               # sort_coefficient = NULL,
                               # sort_class = NULL,
                               sort_group = NULL,
                               sort_parameter = c(letters[1:26], "sigma"),
                               sort_coefficient = c("Intercept"),
                               sort_class = c("b", "sd", "cor")
                               ) {
  arguments <- as.list(match.call())[-1]
  
  if(missing(model)) {
    model <- NULL
  }
  # stop()
  
  if(is.null(model) & is.null(spriors) & is.null(sdata)) {
    stop("Supply either model or spriors and sdata arguments")
  } else if(!is.null(model) & !is.null(spriors) & !is.null(sdata)) {
    stop("Supply only model or spriors and sdata arguments")
  } else if(!is.null(model)) {
    spriors <- prior_summary(model)
    sdata <- standata(model)
  } else if(is.null(model)) {
    if(is.null(spriors) & is.null(sdata)) {
      stop("Supply spriors and sdata arguments")
    }
    if(is.null(spriors) & is.null(sdata)) {
      stop("Supply spriors and sdata arguments")
    }
  }
  
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  spriors <- spriors %>% filter(source == 'user')
  
  env_ <- environment()
  list2env(sdata, envir =  env_)
  # a_cov_b_scale %>% print() 
  
  for(i in 1:nrow(spriors)) {
    getxit <- spriors[i,]$prior
    prior_name <- strsplit(getxit, "\\(")[[1]][1]
    
    if(!prior_name_asit) {
      if(prior_name == 'lkj') {
        prior_name_case <- toupper(prior_name)
      } else if(prior_name == 'lkj_corr_cholesky') {
        prior_name_case <- 'LKJ'
      } else {
        prior_name_case <- firstup(prior_name)
      }
    }
    
    if(prior_name_asit) {
      prior_name_case <- prior_name
    }
    
    getxit_2 <- regmatches(getxit, gregexpr("(?<=\\().*?(?=\\))", getxit, perl=T))[[1]]
    getxit_3 <- strsplit(getxit_2, ",")[[1]]
    getxit_4 <- sapply(getxit_3, function(x) eval(parse(text=x)))
    getxit_5 <- paste(getxit_4, collapse = ", ")
    getxit_6 <- paste0("(", getxit_5, ")")
    getxit_7 <- paste0(prior_name_case, getxit_6)
    spriors[i,]$prior <- getxit_7
  }
  
  # for(i in 1:nrow(spriors)) {
  #   spriors[i,] %>% print()
  # }
  
  
  spriors <- spriors %>% data.frame() %>% dplyr::select(-c('lb', 'ub', 'source'))
  spriors <- spriors %>% `rownames<-`( NULL )
  # spriors <- spriors %>%  dplyr::mutate(class =  dplyr::if_else(class == 'b', 'Beta', class))
  # spriors <- spriors %>%  dplyr::mutate(class =  dplyr::if_else(class == 'sd', 'Std.dev', class))
  # spriors <- spriors %>%  dplyr::mutate(class =  dplyr::if_else(class == 'cor', 'Corr', class))
  spriors <- spriors %>%  dplyr::mutate(class =  dplyr::if_else(class == 'L', 'cor', class))
 
  
  if(!is.null(gsub_coef)) {
    for (gsub_coefi in gsub_coef) {
      spriors <- spriors %>%  dplyr::mutate(coef = gsub(gsub_coefi, "" , coef))
    }
  }
  
  if(!is.null(gsub_group)) {
    for (gsub_groupi in gsub_group) {
      spriors <- spriors %>%  dplyr::mutate(.data$group = gsub(gsub_groupi, "" , .data$group))
    }
  }
  
  
  spriors <- spriors %>% relocate(.data$nlpar, .data$coef, 
                                  .data$class, .data$prior, 
                                  .data$group, .data$resp, 
                                  .data$dpar)
  #stop()
  # for sigma betas
  spriors <- spriors %>%  dplyr::mutate(coef =  dplyr::if_else(coef == '' &
                                                 class == 'Intercept',
                                               class, coef))
  
  spriors <- spriors %>%  dplyr::mutate(class =  dplyr::if_else(class == 'Intercept' &
                                                  dpar == 'sigma' &
                                                  class == 'Intercept',
                                                'b', class))
  
  
  
  spriors <- spriors %>%  dplyr::mutate(nlpar =  dplyr::if_else(nlpar == '' & dpar != '', 
                                                dpar, nlpar)) %>% 
    dplyr::select(-'dpar')
  
  
  
  spriors <- spriors %>% rename(.data$Parameter = .data$nlpar,
                                .data$Coefficient = .data$coef,
                                .data$Class = .data$class, 
                                .data$Prior = .data$prior,
                                .data$Group = .data$group,
                                .data$Response = .data$resp)
  
  
  spriors <- spriors %>% 
    dplyr::arrange(match(Response, sort_response  )) %>% 
    dplyr::arrange(match(Coefficient, sort_coefficient  )) %>% 
    dplyr::arrange(match(Parameter, sort_parameter  )) %>% 
    dplyr::arrange(match(Group, sort_group)) %>% 
    dplyr::arrange(match(Class, sort_class)) 
  
  if(!is.null(model)) {
    if(is.na(model$model_info$univariate_by) |
       !model$model_info$multivariate) {
      spriors <- spriors %>%  dplyr::select(-'Response')
    }
  }
  
  # spriors <- spriors %>% 
  #   gt::gt()  %>% 
  #   gt::cols_align(
  #     align = "left",
  #     columns = everything()
  #   )
  
  spriors
}

