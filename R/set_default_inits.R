


#' An internal function to set default initials for model specific parameters
#'
#' @param select_model A character string specifying the model fitted. Please 
#' see \code{bgm} function for details. 
#' 
#' @param init A character string specifying the initials. Please see
#' \code{bgm} function for details. 
#' 
#' @param class A character string specifying the parameter class. Options
#' are \code{'b'}, \code{'sd'} and \code{'cor'}. Default \code{NULL} indicates 
#' that class name in infered automatically. 
#' 
#' @param parameter A character string specifying the parameter name. Options
#'  are \code{'a'}, \code{'b'}, \code{'c'}, \code{'d'}, \code{'e'}, \code{'f'},
#'  \code{'g'}, \code{'h'}, and \code{'i'}. Default \code{NULL} indicates that
#'  parameter name in infered automatically. 
#'
#' @return A character string.
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' @keywords internal
#' @noRd
set_default_inits <- function(select_model,
                               init,
                               class = NULL,
                               parameter = NULL) {
  if (is.null(parameter)) {
    parameter <- strsplit(deparse(substitute(init)), "_")[[1]][1]
  }
  if (is.null(class)) {
    get_suffix <- strsplit(deparse(substitute(init)), "_")[[1]]
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
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "ymean"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "ymean"
      } else if (grepl('^pb', select_model)) {
        init_out <- "ymax"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "ymax"
        }
        if (select_model == 'logistic2') {
          init_out <- "ymax"
        }
        if (select_model == 'logistic3') {
          init_out <- "ymin"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter b class b
  if (parameter == 'b' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- "ymaxs"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 0.1
        }
        if (select_model == 'logistic2') {
          init_out <- "ymaxs"
        }
        if (select_model == 'logistic3') {
          init_out <- 1.5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter c class b
  if (parameter == 'c' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- 0.1
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 5.0
        }
        if (select_model == 'logistic2') {
          init_out <- 0.1
        }
        if (select_model == 'logistic3') {
          init_out <- 0.1
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter d class b
  if (parameter == 'd' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- 1.0
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- 1.2
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- 1.0
        }
        if (select_model == 'logistic2') {
          init_out <- 1.2
        }
        if (select_model == 'logistic3') {
          init_out <- "ymeanxmid"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter e class b
  if (parameter == 'e' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        if(select_model == 'pb1') init_out <- 13
        if(select_model == 'pb2') init_out <- 13
        if(select_model == 'pb3') init_out <- 13
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- 7
        }
        if (select_model == 'logistic3') {
          init_out <- 0.15
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter f class b
  if (parameter == 'f' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        if(select_model == 'pb1') init_out <- NULL
        if(select_model == 'pb2') init_out <- 2
        if(select_model == 'pb3') init_out <- 1
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- 13
        }
        if (select_model == 'logistic3') {
          init_out <- 5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter g class b
  if (parameter == 'g' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 
            "ymeanxmidxmaxdiff"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter h class b
  if (parameter == 'h' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 1.5
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter i class b
  if (parameter == 'i' & class == 'b') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- NULL
      } else if (grepl('^rcs', select_model)) {
        init_out <- NULL
      } else if (grepl('^pb', select_model)) {
        init_out <- NULL
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- NULL
        }
        if (select_model == 'logistic2') {
          init_out <- NULL
        }
        if (select_model == 'logistic3') {
          init_out <- 13
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  
  
  
  ##############################################################
  # class sd
  ##############################################################
  
  # parameter a class sd
  if (parameter == 'a' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, ysd, autoscale = 2.5)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmin, autoscale = 2.5)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter b class sd
  if (parameter == 'b' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 2.5)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, ysd, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter c class sd
  if (parameter == 'c' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.1, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 3, autoscale = 1)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.1, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 1, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter d class sd
  if (parameter == 'd' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmid, autoscale = 2.5)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter e class sd
  if (parameter == 'e' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 0.15, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter f class sd
  if (parameter == 'f' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter g class sd
  if (parameter == 'g' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, ysdxmidxmaxdiff, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter h class sd
  if (parameter == 'h' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  # parameter i class sd
  if (parameter == 'i' & class == 'sd') {
    if (init == 'NA' | init == '') {
      if (grepl('^sitar', select_model)) {
        init_out <- "normal(0, 2, autoscale = 1)"
      } else if (grepl('^rcs', select_model)) {
        init_out <- "normal(0, ysd, autoscale = 1)"
      } else if (grepl('^pb', select_model)) {
        init_out <- "normal(0, 0.25, autoscale = 1)"
      } else if (grepl('^logistic', select_model)) {
        if (select_model == 'logistic1') {
          init_out <- "normal(0, 2, autoscale = 2.5)"
        }
        if (select_model == 'logistic2') {
          init_out <- "normal(0, 0.25, autoscale = 1)"
        }
        if (select_model == 'logistic3') {
          init_out <- "normal(0, 2, autoscale = 1)"
        }
      } else {
        
      }
      init_out <- gsub_space(init_out)
    } else {
      init_out <- init
    }
  }
  
  
  
  
  # print(parameter)
  # print(init_out)
  return(init_out)
}

