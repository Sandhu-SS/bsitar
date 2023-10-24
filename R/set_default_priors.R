


#' An internal function to set default prior for model specific parameters
#'
#' @param select_model A character string specifying the model fitted. Please 
#' see \code{bgm} function for details. 
#' 
#' @param prior A character string specifying the parameter name. Please see
#' \code{bgm} function for details. 
#' 
#' @param class A character string specifying the parameter name. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param parameter A character string specifying the parameter name. Options
#'  are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'}, \code{'e'}, \code{'f'},
#'  \code{'g'}, \code{'h'}, and \code{'i'}. Default \code{NULL} indicates that
#'  parameter name in infered automatically. 
#'
#' @return A character string.
#' @keywords internal
#' @noRd
set_default_priors <- function(select_model,
                               prior,
                               class = NULL,
                               parameter = NULL) {
  if (is.null(parameter)) {
    parameter <- strsplit(deparse(substitute(prior)), "_")[[1]][1]
  }
  if (is.null(class)) {
    get_suffix <- strsplit(deparse(substitute(prior)), "_")[[1]]
    get_suffix <- get_suffix[length(get_suffix)]
    if (grepl("^beta", get_suffix))
      class <- 'b'
    if (grepl("^sd", get_suffix))
      class <- 'sd'
  }
  
  ##############################################################
  # class b
  ##############################################################
  
  # parameter a class b
  if (parameter == 'a' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter b class b
  if (parameter == 'b' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymaxs, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymaxs, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 0.75, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter c class b
  if (parameter == 'c' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0.1, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0.1, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter d class b
  if (parameter == 'd' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter e class b
  if (parameter == 'e' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(7, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter f class b
  if (parameter == 'f' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(13, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class b
  if (parameter == 'g' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class b
  if (parameter == 'h' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter i class b
  if (parameter == 'i' & class == 'b') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  
  
  
  ##############################################################
  # class sd
  ##############################################################
  
  # parameter a class sd
  if (parameter == 'a' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(ymax, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(ymaxs, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(ymean, ysd, autoscale = 2.5)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter b class sd
  if (parameter == 'b' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(ymaxs, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 0.75, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter c class sd
  if (parameter == 'c' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(0.1, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter d class sd
  if (parameter == 'd' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter e class sd
  if (parameter == 'e' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter f class sd
  if (parameter == 'f' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'g' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'h' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  # parameter i class sd
  if (parameter == 'i' & class == 'sd') {
    if (prior == 'NA' | prior == '') {
      if (grepl('^sitar', select_model)) {
        prior_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        prior_out <- "normal(ymean, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        prior_out <- "normal(1.2, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          prior_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          prior_out <- "normal(1.2, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      prior_out <- gsub_space(prior_out)
    } else {
      prior_out <- prior
    }
  }
  
  
  
  
  # print(parameter)
  # print(prior_out)
  return(prior_out)
}

